# Snakemake Coding Guidelines

> Snakemake v8.30.0 | Snap-HiC-Diff project conventions

**⚠️ Snakemake v8 Executor Plugin Architecture (Breaking Change):**

In Snakemake v8, executors (`cluster-generic`, `slurm`, `drmaa`, etc.) are **separate pip packages**. The base `snakemake` package no longer includes them. This project uses `cluster-generic` which requires:

```yaml
# envs/SNAP.yaml
- pip:
  - snakemake-executor-plugin-cluster-generic
```

Without this plugin, `--executor cluster-generic` raises:
```
snakemake: error: argument --executor/-e: invalid choice: 'cluster-generic'
```

**Login node OOM with cluster-generic:** Even `--executor cluster-generic -n` (dry-run) loads all registered executor plugins on the login node, causing OOM on resource-constrained nodes. Test cluster submission via `dsub -s test.sh` on a compute node instead.

## Table of Contents

1. [Snakefile Structure](#1-snakefile-structure)
2. [Rule Anatomy](#2-rule-anatomy)
3. [Wildcards](#3-wildcards)
4. [Input / Output](#4-input--output)
5. [params](#5-params)
6. [threads & resources](#6-threads--resources)
7. [conda](#7-conda)
8. [Execution methods: shell / script / wrapper](#8-execution-methods-shell--script--wrapper)
9. [expand() — Wildcard Expansion](#9-expand--wildcard-expansion)
10. [Input functions](#10-input-functions)
11. [Semantic helpers](#11-semantic-helpers)
12. [Config & Configuration](#12-config--configuration)
13. [Modules & Include](#13-modules--include)
14. [Common.smk: Shared helpers](#14-commonsmk-shared-helpers)
15. [Best practices from upstream docs](#15-best-practices-from-upstream-docs)
16. [DAG & execution](#16-dag--execution)
17. [Project-specific patterns](#17-project-specific-patterns)

---

## 1. Snakefile Structure

### Standard layout (per upstream recommendation)

```
project/
├── Snakefile               # Entry point (discovers workflow/Snakefile)
├── workflow/
│   ├── Snakefile           # Main snakefile (recommended)
│   ├── rules/
│   │   ├── common.smk      # Shared Python vars, helper functions
│   │   ├── rule1.smk
│   │   └── rule2.smk
│   ├── scripts/
│   │   ├── rule1.py
│   │   └── rule2.py
│   └── envs/
│       ├── rule1.yaml
│       └── rule2.yaml
├── config/
│   ├── config.yaml         # Main config
│   ├── threads.yaml        # Per-rule thread counts
│   └── resources.yaml      # Per-rule memory specs
└── results/
```

### Entry point Snakefile pattern (used in this project)

```python
import os, yaml
from snakemake.utils import validate

configfile: os.path.join(os.path.dirname(workflow.snakefile), "config", "config.yaml")

# Resolve outdir to absolute path
_outdir_raw = config.get("outdir", "results")
config["outdir"] = os.path.abspath(...)
```

**Key:** `configfile:` path is relative to the **working directory**, not the Snakefile location.
Use `os.path.dirname(workflow.snakefile)` as base for relative paths within included rules.

**⚠ cluster-generic Working Directory Trap:**

When using `--executor cluster-generic`, Snakemake generates a Python jobscript that is executed by dsub on remote nodes. The remote node's working directory is **not guaranteed** to be the project root — even though `--cwd "$PROJECT_ROOT"` is set in dsub_submit.sh, it controls the shell wrapper's cwd, not the Python jobscript execution context.

**Consequence:** Relative paths in `--configfile`, `configfile:`, and `include:` statements resolve from the remote node's actual cwd, causing `FileNotFoundError`.

**Required pattern (always):**

```bash
# Always pass absolute paths and --directory for cluster-generic
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_CONFIG_ABS="$PROJECT_ROOT/config/config.yaml"

snakemake \
    --executor cluster-generic \
    --directory "$PROJECT_ROOT" \      # ← locks Snakemake's cwd on remote nodes
    -s Snakefile \
    --configfile "$_CONFIG_ABS" \    # ← absolute path
    ...
```

**Local dry-run cannot detect this issue** — `snakemake -n` always runs from the local shell's cwd. Test against the cluster executor with a single job before full submission.

### cluster-generic-submit-cmd placeholders (critical)

The `--cluster-generic-submit-cmd` template receives the **jobscript path as a positional argument** — the Snakemake cluster-generic plugin always auto-appends it as the 4th argument (`$3` in zero-indexed bash). The `{jobscript}` **placeholder does not exist** in Snakemake's formatting system and will cause `NameError`.

**Correct pattern — positional argument only (no `{jobscript}` placeholder):**

```bash
--cluster-generic-submit-cmd 'bash dsub_submit.sh {threads} {rulename}'
```

> The plugin internally appends the jobscript path, so `$3` in the wrapper script receives it automatically.

**Invalid — `{jobscript}` is NOT a valid placeholder in Snakemake v8:**

```bash
# WRONG: {jobscript} does not exist in Snakemake's .format() context
--cluster-generic-submit-cmd 'bash dsub_submit.sh {threads} {rulename} {jobscript}'
# This causes: NameError: The name 'jobscript' is unknown in this context
```

**Valid placeholders for `--cluster-generic-submit-cmd`:**

| Placeholder | Description |
|-------------|-------------|
| `{threads}` | Number of threads requested |
| `{rule}`, `{rulename}` | Rule name |
| `{jobid}` | Internal Snakemake job ID |
| `{profile}` | Profile name |
| `{resources.*}` | Resource values (e.g. `{resources.mem_mb}`) |
| `{attempt}` | Job attempt number |

**Why `{jobscript}` fails**: Snakemake's cluster-generic executor plugin (snakemake-executor-plugin-cluster-generic v1.0.9) uses Python's `str.format()` internally with a fixed set of variables. The `{jobscript}` variable is not in that set. The plugin additionally writes the jobscript to a temp file and passes its path directly as a quoted positional argument to the shell command, bypassing the template entirely.

The wrapper script must capture it with positional parameters:

```bash
_threads="${1:?missing threads}"
_rulename="${2:?missing rulename}"
_jobscript="${3:?missing jobscript}"
```

**Snakemake profile YAML pattern (preferred over CLI flags):**

```yaml
# workflow/profile/dsub/config.yaml
executor: cluster-generic
jobs: 200
keep-going: true
rerun-incomplete: true
notemp: true
cluster-generic-submit-cmd: "bash {DEPLOYMENT_SOURCE}/dsub_submit.sh {threads} {rulename}"
cluster-generic-cancel-cmd: "dkill"
```

> **Note**: Do NOT append `{jobscript}` to the `cluster-generic-submit-cmd` template. The plugin auto-appends it as a quoted positional argument.

> **CI/Audit**: Before each deployment, run `grep -rn '{jobscript}' *.sh` to detect this pattern. The only valid placeholders are: `{threads}`, `{rulename}`, `{rule}`, `{jobid}`, `{profile}`, `{resources.*}`, `{attempt}`.

**⚠ Flag name typo trap:** The flag is `--cluster-generic-submit-cmd` (with `-cmd` suffix). Writing `--cluster-generic-submit` (without `-cmd`) is silently ignored by Snakemake.

---

## 2. Rule Anatomy

Every rule must declare all of the following:

```python
rule rule_name:
    """
    One-line description.

    Input:
        - description

    Output:
        - description
    """
    input:
        file_or_func,
    output:
        file_or_func,
    log:
        os.path.join(LOGS_DIR, "rule_name", "{wildcard}", "rule_name.log"),
    benchmark:
        os.path.join(BENCHMARK_DIR, "rule_name", "{wildcard}", "rule_name.tsv"),
    params:
        param1=config.get("section", {}).get("param1", default),
        threads=config.get("threads", {}).get("rule_name", 1),
    threads: lambda wildcards, resources, params: params.threads,
    resources:
        mem_mb=config.get("resources", {}).get("rule_name", {}).get("mem_mb", 4096),
    conda:
        "envs/rule_name.yaml"
    script:
        "workflow/scripts/rule_name.py"
```

**Order of directives does not matter**, but grouping as above improves readability.

### Required directives for every rule in this project

| Directive | Required | Notes |
|-----------|----------|-------|
| `input` | Yes | Can be a string, list, dict, lambda, or function |
| `output` | Yes | Supports named outputs: `cool=..., done=...` |
| `log` | Yes | Pattern must match output wildcards |
| `benchmark` | Yes | TSV benchmark output |
| `params` | Yes | Load from `config`, include `threads` here |
| `threads` | Yes | Lambda form: `lambda wildcards, resources, params: params.threads` |
| `conda` | Yes | Relative to the `.smk` file location |
| `script` / `shell` / `wrapper` | Yes | Execution method |

### Output with side-effect marker

Use `touch()` to create a sentinel `.done` file alongside real output:

```python
output:
    cool=os.path.join(OUTDIR, "chr_mcool_kr", "{sample_id}", "{chromosome}.cool"),
    done=touch(os.path.join(OUTDIR, "chr_mcool_kr", "{sample_id}", "{chromosome}.kr.done")),
```

This enables DAG tracking for operations that don't produce a primary output.

---

## 3. Wildcards

Wildcards are `{name}` placeholders in input/output paths. Snakemake automatically resolves them by matching output patterns against requested files.

### Basic example

```python
rule extract_chr:
    input:
        mcool=lambda wildcards: os.path.join(RAWDATA_DIR, f"{wildcards.sample_id}.mcool"),
    output:
        cool=os.path.join(OUTDIR, "chr_mcool", "{sample_id}", "{chromosome}.cool"),
```

### Wildcard constraints

Constrain wildcards with regex to avoid ambiguity:

```python
# Inline constraint
output: r"{sample_id,[A-Za-z0-9_-]+}/{chromosome}.cool"

# Or global constraint (in common.smk)
wildcard_constraints:
    sample_id="[A-Za-z0-9_-]+"
```

### Wildcard names must match between input and output

The same `{wildcard}` name in input and output must refer to the same value. A wildcard can appear only in output (not input) but not vice versa without matching.

---

## 4. Input / Output

### Named inputs/outputs

Prefer named inputs for clarity, especially for rules with multiple files:

```python
input:
    cool=os.path.join(OUTDIR, "chr_mcool", "{sample_id}", "{chromosome}.cool"),
    manifest=os.path.join(OUTDIR, "input_check", "manifest.tsv"),
output:
    cool=os.path.join(OUTDIR, "chr_mcool_kr", "{sample_id}", "{chromosome}.cool"),
```

Access in script: `snakemake.input.cool`, `snakemake.output.cool`

### Positional access

```python
input: "file1.txt", "file2.txt"
output: "out.txt"
# Access: snakemake.input[0], snakemake.output[0]
```

### Lists and dicts

```python
input:
    files=["a.txt", "b.txt"],       # space-separated in shell: "a.txt b.txt"
    meta={"config": "cfg.yaml"},    # named access only
output: multiext("results/plot", ".pdf", ".svg", ".png")
```

### Formatting in shell commands

```python
shell:
    "command {input} {output}"          # space-separated list
    "command {input:q} {output:q}"       # quote each element
    "command {wildcards.sample_id}"       # direct wildcard access
    "command {params.threshold}"          # param access
    "{{input}}"                          # escaped — literal "{input}", not replaced
```

---

## 5. params

`params` is for non-file parameters passed to the execution script/shell.

### Loading from config

```python
params:
    resolution=RESOLUTION,                          # global constant
    threads=config.get("threads", {}).get("rule_name", 20),
    gradient_threshold=config.get("gap", {}).get("gradient_threshold", 10.0),
```

### Callable params (evaluated at job execution)

```python
params:
    threads=config.get("threads", {}).get("rule_name", 8),
threads: lambda wildcards, resources, params: params.threads
```

### Params functions signature

Functions assigned to `params` (or `resources`, `threads`) can take:
- `wildcards` — required first argument
- `input` — optional second
- `output` — optional
- `threads` — optional
- `resources` — optional
- `attempt` — optional (for retry-based resource scaling)

---

## 6. threads & resources

### Pattern: define in params, use in threads directive

```python
params:
    threads=config.get("threads", {}).get("rule_name", 8),
threads: lambda wildcards, resources, params: params.threads
```

### Dynamic threads based on available cores

```python
threads: workflow.cores * 0.75   # reserve 75% of available cores
```

### Standard resources

| Resource | Meaning | Unit |
|----------|---------|------|
| `mem_mb` | Memory | MB |
| `disk_mb` | Disk I/O | MB |
| `runtime` | Wall clock time | minutes or string (e.g. `"2h30m"`) |
| `tmpdir` | Temp directory | path (sets `$TMPDIR`) |
| `gpu` | Number of GPUs | int |

### Default resources

Snakemake auto-computes: `mem_mb=max(2*input.size_mb, 1000)` unless overridden.

### Callable resources (dynamic sizing)

```python
resources:
    mem_mb=lambda wildcards, input: max(2.5 * input.size_mb, 300)
    # For retry: attempt parameter grows memory on each retry
    mem_mb=lambda wildcards, attempt: attempt * 100
```

### Per-rule resource override via config

```python
resources:
    mem_mb=config.get("resources", {}).get("rule_name", {}).get("mem_mb", 4096),
```

### Memory follows threads — no manual mem setting needed

> **所有 Snakemake rule 的内存需求通过 `threads` 参数间接申请，不单独设置内存参数。**

本项目通过 `dsub_submit.sh` 提交作业，只传入 `{threads}`，服务器按队列比例自动分配内存：

```
dsub -R "cpu=${_THREADS}"  →  arm 队列：${_THREADS} × 5 GB
```

**核心原则**：
- **每个 rule 的 `threads` 值同时决定了 CPU 核数和内存上限**（arm 每核 5 GB）
- rule 作者应估算内存需求，合理设置 `threads`（不是随意设大）
- 申请 10 核 → 自动 50 GB；申请 40 核 → 自动 200 GB
- 内存需求 > 500 GB 时改用 `fat_arm` 或 `fat_x86` 队列（修改 `dsub_submit.sh` 的 `-q` 参数）

配置方式：

```python
# threads.yaml — 决定 CPU 核数，同时决定内存上限
threads:
  extract_chr_mcool: 20   # → arm 队列 100 GB
  kr_normalize: 8         # → arm 队列 40 GB
  detect_gaps: 4          # → arm 队列 20 GB
```

详细内存规则见 [dsub-guidelines.md](./dsub-guidelines.md#7-内存限制)

---

## 7. conda

Path is relative to the `.smk` file location. In this project, use relative paths from `workflow/rules/`:

```python
conda:
    "../../envs/SNAP.yaml"
```

### Environment definition pattern

```yaml
channels:
  - conda-forge
  - bioconda
  - nodefaults
dependencies:
  - python=3.12
  - cooler=0.10
```

### Offline conda environment management (无网络集群方案)

> **问题场景**：计算节点无网络，无法在运行时创建 conda 环境。
>
> **解决方案**：在有网节点预创建环境，计算节点自动复用。

#### 架构概览

```
┌─────────────────────────────────────────────────────────────────┐
│  共享存储 (NFS)                                                  │
│  project/                                                       │
│     ├── envs/SNAP.yaml         ← 环境定义文件 (git 可管理)       │
│     ├── .snakemake/conda/      ← 预创建的 conda 环境             │
│     │      └── snakemake-snaphic/                                │
│     ├── run_snakemake.sh        ← 主运行脚本                      │
│     ├── setup_conda_envs.sh     ← 环境预创建脚本 (新增)           │
│     └── workflow/profile/                                                       │
│         └── dsub_submit.sh      ← 作业提交脚本 (已更新)           │
└─────────────────────────────────────────────────────────────────┘
         ↑ 有网节点创建                    ↓ 无网计算节点使用
         环境                               自动复用
```

#### 文件清单

| 文件 | 用途 | 可移植性 |
|------|------|---------|
| `envs/SNAP.yaml` | 环境定义 | ✓ 随项目迁移 |
| `.snakemake/conda/` | 预创建的环境 | ✗ 需重新创建 |
| `setup_conda_envs.sh` | 环境创建脚本 | ✓ 随项目迁移 |
| `dsub_submit.sh` | 作业提交脚本 | ✓ 已配置相对路径 |

#### 使用流程

```bash
# Step 1: 在有网节点预创建环境（一次性）
cd project/
bash setup_conda_envs.sh

# Step 2: 提交作业到无网计算节点
bash run_snakemake.sh

# Step 3: 恢复失败的作业
bash run_resume.sh
```

#### 关键配置点

**1. run_snakemake.sh 必须包含 --use-conda:**

```bash
snakemake \
    --use-conda \
    --executor cluster-generic \
    ...
```

**2. dsub_submit.sh 使用 ${cwd} 实现可移植路径:**

```bash
# 使用 dsub 的 ${cwd} 变量（在 wrapper 脚本中展开）
_project_root="\${cwd}"
_conda_envs_path="\${_project_root}/.snakemake/conda"

export CONDA_ENVS_PATH="${_conda_envs_path}"
conda activate "${_conda_envs_path}/snakemake-snaphic"
```

**3. setup_conda_envs.sh 预创建所有环境:**

```bash
snakemake \
    --use-conda \
    --conda-create-envs-only \
    -s Snakefile \
    --configfile config/config.yaml
```

#### 扩展多环境

当需要多个 conda 环境时：

```yaml
# envs/env_A.yaml
name: env-a
channels:
  - conda-forge
dependencies:
  - tool-a=1.0
```

```python
# workflow/rules/rule_a.smk
rule rule_a:
    conda:
        "../../envs/env_A.yaml"   # 引用新环境
```

运行 `bash setup_conda_envs.sh` 会自动创建所有规则引用的环境。

#### 可移植性验证

当项目迁移到新位置时：

```bash
# 方案 1: 重新创建环境
cp -r /old/project /new/project
cd /new/project
bash setup_conda_envs.sh  # 重新创建环境

# 方案 2: 整体迁移（如果目标机器 conda 路径相同）
cp -r /old/project /new/project
# 直接运行，无需重新创建环境
bash run_snakemake.sh
```

### conda + container combo

```python
container: "docker://continuumio/miniconda3:4.4.10"
conda: "envs/tool.yaml"   # executed inside container
# Run with: snakemake --sdm conda apptainer
```

---

## 8. Execution methods: shell / script / wrapper

### script (preferred for complex logic)

```python
script: "workflow/scripts/rule_name.py"
```

Inside the Python script, access Snakemake object:

```python
from snakemake.script import snakemake

cool_path = snakemake.input.cool
out_path  = snakemake.output.cool_kr
threads   = snakemake.threads
config    = snakemake.config
wildcards = snakemake.wildcards
log       = snakemake.log
```

### shell (for simple commands)

```python
shell:
    "cooler cat {input.cool} | cooler load {params} - {output.cool}"
```

### wrapper (for curated tools)

```python
wrapper: "0.0.8/bio/samtools/sort"
# Supports: git tags, commit IDs, local file:// URLs
```

---

## 9. expand() — Wildcard Expansion

`expand()` generates all combinations of wildcard values (default: `itertools.product`).

```python
expand(
    os.path.join(OUTDIR, "chr_mcool", "{sid}", "{chr}.cool"),
    sid=get_all_sample_ids(),
    chr=get_all_chromosomes()
)
# → all combinations of sample × chromosome
```

### zip mode (pairwise)

```python
expand("{dataset}/a.{ext}", zip, dataset=DATASETS, ext=FORMATS)
# → pairs dataset[i] with ext[i]
```

### Masked wildcards (keep as placeholder)

```python
expand("{{dataset}}/a.{ext}", ext=FORMATS)
# → "dataset/a.txt", "dataset/a.csv" (literal {dataset})
```

### ⚠️ CRITICAL: `.format()` in expand() strings — Python str.format() runs BEFORE Snakemake wildcard resolution

When a path string uses `.format()` for OUTDIR/RESOLUTION substitution, the `.format()` call happens during Python module load — **before** Snakemake's `expand()` evaluates wildcards.

**Decision table:**

| Path string has `.format()`? | Wildcard in expand() args? | Use |
|-------------------------------|---------------------------|-----|
| No (literal string) | Yes | `{sample}` |
| Yes `"{OUTDIR}/out/{sample}"` | Yes | `{{sample}}` (double brace so `.format()` produces literal `{sample}`) |
| Yes `"{OUTDIR}/out/{chrom}"` | Yes | `{chrom}` (only if `chrom` IS in expand kwargs — then double-brace NOT needed) |

**Examples from this project:**

```python
# GOOD: literal string, wildcards in expand kwargs → single brace
expand(OUTDIR + "/balanced/{sample}__{chrom}_iced.cool",
       sample=config["samples"], chrom=get_chroms())

# GOOD: .format() string, wildcard in expand kwargs → double brace for .format() targets
expand("{}/smooth/{{sample}}_smooth.cool".format(OUTDIR),
       sample=config["samples"])

# GOOD: .format() string, wildcard NOT in expand kwargs → single brace
expand("{}/windows/windows_{}.parquet".format(OUTDIR, RESOLUTION))

# BAD: .format() string with single brace — Python KeyError at module load
expand("{}/smooth/{sample}_smooth.cool".format(OUTDIR), sample=config["samples"])
#         ↑ Python tries to format with no 'sample' kwarg → KeyError: 'sample'
```

**Why `{{sample}}` works in `.format()` strings:**
`"{{sample}}".format()` → `"{sample}"` (double brace `{{` → literal `{`, so `.format()` produces the Snakemake wildcard placeholder).

**Why single brace `{sample}` works without `.format()`:**
No `.format()` to intercept it — `expand()` sees `{sample}` and fills it from kwargs.

**Rule of thumb:** If the path string has `.format(...)`, double-brace every `{}` that represents a Snakemake wildcard but is NOT a `.format()` kwarg. Use `os.path.join()` instead of `.format()` to avoid this whole class of bugs.

### Functions as expand arguments

```python
expand(
    os.path.join(OUTDIR, "normalized_windows", "{sid}", "{chr}", "{wid}.norm.done"),
    sid=[p[0] for p in get_comparison_pairs()],   # list comprehension
    chr=get_all_chromosomes(),
    wid=get_window_ids_for_chr(sample_id, chromosome)  # evaluated later
)
```

---

## 10. Input functions

Input can be a callable returning file(s) based on wildcards:

```python
input:
    mcool=lambda wildcards: os.path.join(RAWDATA_DIR, f"{wildcards.sample_id}.mcool"),
```

### Input functions that depend on checkpoints

When a rule's output wildcards need to be dynamically determined by a `checkpoint`, use `checkpoints` (plural, the module-level object) inside an **input function**:

```python
def get_metric_chunks(wildcards):
    """Return all chunk metric files for a given pair_id (reads checkpoint output)."""
    manifest_path = checkpoints.chunk_windows.get(pair_id=wildcards.pair_id).output.manifest
    df = pd.read_csv(manifest_path, sep="\t")
    return [f"{OUTDIR}/metrics/{wildcards.pair_id}/chunk_{cid:03d}_metrics.tsv"
            for cid in df["chunk_id"].tolist()]

rule aggregate_metrics:
    input:
        get_metric_chunks  # No parentheses — the function object itself
```

**Critical rules for `checkpoints` usage:**

1. **Always inside an input function** — `checkpoints.chunk_windows.get(...)` must be inside `lambda wildcards:` or a named function. Never at module level.
2. **Never wrap in `ancient()`** — `ancient()` evaluates at module level; `checkpoint()` is not available there.
3. **Use `.output.<name>`** — access checkpoint output attributes by the name defined in the checkpoint's `output:` directive (e.g., `.output.manifest`, not `.output[0]`).
4. **No parentheses in `input:`** — pass the function object (`get_metric_chunks`), not its result (`get_metric_chunks(wildcards)`).

**Wrong patterns:**

```python
# ❌ Module-level checkpoint call — NameError
input:
    ancient(checkpoint(OUTDIR + "/chunks/{pair_id}/chunk_manifest.tsv")),

# ❌ Wrapping checkpoint in ancient() — checkpoint not available at module level
input:
    ancient(checkpoints.chunk_windows.get(pair_id="{pair_id}").output.manifest),

# ❌ Calling the function in input:
input:
    get_metric_chunks(wildcards)  # wildcards not defined here!
```

### Returns a dict (use unpack())

```python
def myfunc(wildcards):
    return {"foo": f"{wildcards.token}.txt"}

rule:
    input:
        unpack(myfunc),
    output:
        "out.txt"
```

### Returns a list (auto-flattened, no unpack needed)

```python
def myfunc(wildcards):
    return ["a.txt", "b.txt"]

input:
    *myfunc()      # inline unpacking
```

---

## 11. Semantic helpers

### multiext — Multiple extensions

```python
output: multiext("results/plot", ".pdf", ".svg", ".png")
# Equivalent to: expand("results/plot{ext}", ext=[".pdf", ".svg", ".png"])
```

### touch() — Sentinel files

```python
output:
    cool="out.cool",
    done=touch("out.done"),   # creates empty marker file on success
```

### branch — Conditional inputs

```python
input:
    branch(
        config.get("use_kr", False),
        then="input_kr.txt",
        otherwise="input_ice.txt"
    )
```

### lookup — Config/dataframe lookup by wildcard

```python
input:
    lookup(
        query="sample == '{sample}'",
        within=samples_df,
        cols=["bam_path"]
    )
```

### exists — Conditional file existence

```python
input:
    "out.txt" if exists("test.txt") else []
```

### subpath — Extract path components

```python
params:
    outdir=subpath(output.cool, parent=True),
    basename=subpath(output.cool, basename=True),
```

---

## 12. Config & Configuration

### YAML config loading

```python
configfile: "config/config.yaml"

# Access in rules
params:
    resolution=config.get("resolution", 4000)
    section=config.get("section", {})
```

### ⚠️ Chromosome naming: git HEAD vs working tree

> **CRITICAL**: `dsub` jobs run from git HEAD, NOT the working tree. Working tree changes to `config/chroms.yaml` are invisible to dsub.

```
git checkout chroms.yaml       → changes to Ghir_* (committed version)
# working tree: HC04_*         → invisible to dsub!
git add chroms.yaml
git commit                     → now dsub sees HC04_*
```

**Failure mode**: All 36 `extract_and_balance_sample` jobs fail with:

```
[ERROR] Chromosome Ghir_A01 not found in /path/to/sample.mcool::/resolutions/4000
```

The mcool files use `HC04_*`, `chroms.yaml` in git HEAD has `Ghir_*`, `hictk dump` output doesn't match `awk '$1 == chr'`.

**Rule of thumb**: Any config change that affects all samples must be committed before `dsub -s run_snakemake.sh`. Always run `git status` before submission.

**Prevention**: After any `config/chroms.yaml` change, immediately commit before running dsub. Add chromosome naming column to manifest TSV for automatic inference.

### Tabular config (samples.tsv, comparisons.tsv)

```python
import pandas as pd

samples_df = pd.read_csv(config["samples_tsv"], sep="\t", comment="#")
pairs = list(zip(df["sample_a"], df["sample_b"]))
```

### Validation with JSON schema

```python
from snakemake.utils import validate

validate(config, "config/schema.yaml")
validate(samples_df, "config/samples.schema.yaml")
```

### Overwrite via CLI

```bash
snakemake --config outdir=/custom/path resolution=8000
snakemake --configfile config/local.yaml
```

### Thread / resource config files (this project pattern)

```python
# threads.yaml
threads:
  extract_chr_mcool: 20
  kr_normalize: 8
  detect_gaps: 4

# resources.yaml
resources:
  extract_chr_mcool:
    mem_mb: 16384
  kr_normalize:
    mem_mb: 8192
```

Load in Snakefile before includes:

```python
_threads_cfg = yaml.safe_load(open(_threads_cfg_path))
config["threads"] = _threads_cfg["threads"]
```

---

## 13. Modules & Include

### include (flat merging)

```python
include: "workflow/rules/common.smk"
include: "workflow/rules/extract_chr_mcool.smk"
```

All rules from included files share the same namespace. No rule name conflicts allowed.

### module + use (Snakemake 6.0+, named isolation)

```python
module bwa_workflow:
    snakefile: github("owner/repo", path="workflow/Snakefile", tag="v1.0.0")
    config: config["bwa"]

use rule * from bwa_workflow as bwa_*
use rule align from bwa_workflow with:
    input:
        get_input
```

Key: imported module evaluated in separate namespace. Module cannot be named `workflow`.

### module prefix (automatic path isolation)

```python
module rna_seq:
    snakefile: "rna_seq/Snakefile"
    prefix: "rna-seq"      # all relative outputs go under rna-seq/

use rule * from rna_seq as rna_seq_*
```

---

## 14. common.smk: Shared helpers

`common.smk` is included by all rules and must contain **only Python code** — no Snakemake directives (no `rule`, `output:`, etc.).

### Required content

```python
# === Project root and paths ===
PROJECT_ROOT = os.path.abspath(os.path.dirname(workflow.snakefile))
OUTDIR      = config.get("outdir", os.path.join(PROJECT_ROOT, "results"))
LOGS_DIR    = config.get("logs_dir", os.path.join(PROJECT_ROOT, "logs"))
BENCHMARK_DIR = config.get("benchmark_dir", os.path.join(PROJECT_ROOT, "benchmark"))

# === Config defaults (fallbacks) ===
RESOLUTION   = config.get("resolution", 4000)
WINDOW_SIZE  = config.get("window_size", 300000)

# === TSV path constants ===
SAMPLES_TSV      = config.get("samples_tsv", os.path.join(PROJECT_ROOT, "config", "samples.tsv"))
COMPARISONS_TSV  = config.get("comparisons_tsv", ...)
```

### Helper function patterns

**Read TSV without pandas** (for compatibility in module mode):

```python
def _read_tsv_no_pandas(path, sep="\t"):
    result = []
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                cols = line.split(sep)
                result.append(cols)
        return result
    except Exception:
        return []

def get_all_sample_ids():
    if pd is not None:
        try:
            df = pd.read_csv(SAMPLES_TSV, sep="\t", comment="#")
            return df["sample_id"].tolist()
        except Exception:
            pass
    return _get_all_sample_ids_raw()
```

**Wildcard expansion helpers** (used in `rule all`):

```python
def get_all_sample_ids():         # used in expand() for all samples
def get_comparison_pairs():       # returns [(sample_a, sample_b), ...]
def get_all_chromosomes():        # respects test_region config
def get_window_ids_for_chr(sample_id, chromosome):  # reads windows.tsv
```

---

## 15. Best practices from upstream docs

### From Snakemake best practices page

1. **Linting** — Run `snakemake --lint` before publishing
2. **Formatting** — Use [Snakefmt](https://github.com/snakemake/snakefmt) to format workflows
3. **Separation** — Keep Python helper functions in `workflow/rules/common.smk`, not inside rules
4. **Named parameters** — Prefer named `input:`/`output:` keys over positional; avoid lambdas inside rules
5. **Conda/containers** — Annotate every rule with versioned environments for reproducibility
6. **Config** — Use config files for metadata; use CLI args (`--set-threads`, `--set-resources`) for runtime
7. **Reports** — Use Snakemake's [interactive reports](https://snakemake.readthedocs.io/en/v8.30.0/snakefiles/reporting.html) for results
8. **Avoid duplication** — Use [Snakemake wrappers](https://snakemake-wrappers.readthedocs.io/) when available

### Readable rule docstrings

Every rule docstring should follow this template:

```python
rule rule_name:
    """
    Brief one-liner.

    Input:
        - description of input

    Output:
        - description of output

    Params:
        - param: description
    """
```

---

## 16. DAG & execution

### Key execution options

| Flag | Purpose |
|------|---------|
| `-n`, `--dry-run` | Show plan without executing |
| `-p` | Print shell commands |
| `--dag` | Print DAG in dot format |
| `--touch` | Touch output files (mark as complete) |
| `--forcerun -r rule` | Force re-run of specific rule |
| `--until target` | Build up to target only |
| `--jobs N`, `-j N` | Number of parallel jobs |
| `--retries N` | Retry failed jobs N times |

### Profiles

Use profiles for cluster/HPC execution (replaces deprecated cluster config):

```
profiles/
└── slurm/
    └── config.yaml    # sets: executor, jobs, resources, etc.
```

Run: `snakemake --profile slurm`

---

## 17. Project-specific patterns

### Rule template (`workflow/rules/_template.smk`)

Copy this template when creating a new rule. Every new rule must:
- Define all required directives (input, output, log, benchmark, params, threads, conda, script/shell)
- Load thread count from `config.get("threads", {}).get("rule_name", N)`
- Load memory from `config.get("resources", {}).get("rule_name", {}).get("mem_mb", 4096)`
- Use named outputs for multi-file rules
- Use `touch()` sentinel files for side-effect-only outputs
- Document wildcards in the docstring

### File naming conventions

| Pattern | Example |
|---------|---------|
| Rule file | `workflow/rules/{rule_name}.smk` |
| Python script | `workflow/scripts/{rule_name}.py` |
| Conda env | `workflow/envs/{rule_name}.yaml` |
| Log dir | `logs/{rule_name}/{wildcard}/` |
| Benchmark dir | `benchmark/{rule_name}/{wildcard}/` |

### Absolute output paths

Always resolve OUTDIR to absolute path in Snakefile entry point:

```python
config["outdir"] = os.path.abspath(os.path.join(PROJECT_ROOT, _outdir_raw))
```

### Conditional rules via commented includes

Commented-out rule includes serve as a feature gate for future phases:

```python
# include: "workflow/rules/replicate_baseline.smk"
# include: "workflow/rules/background_pairs.smk"
```

### touch() sentinel in multi-output rules

When a rule produces multiple output files but DAG needs a single dependency marker:

```python
output:
    result="output/data.tsv",
    done=touch("output/data.done"),
```

### Python-script-first convention

Complex rule logic lives in `workflow/scripts/*.py` (accessed via `script:` directive).
Shell commands are reserved for one-liners or simple tool invocations.

---

## 17b. Checkpoint Output: Real Files + Two-Stage Architecture

### The problem: checkpoint cannot mix wildcard sets

All outputs of a checkpoint must share the same wildcard set. The following is **invalid**:

```python
# ❌ WRONG: manifest={sample} + gaps expand({chrom}) → wildcard conflict
checkpoint compute_gap_qc:
    output:
        manifest = OUTDIR + "/logs/gap_qc/{sample}.done",
        gaps = expand(
            OUTDIR + "/windows/{{sample}}/{{sample}}__{chrom}_gaps.tsv",
            chrom=CHROMS
        ),
```

Snakemake raises:
```
Not all output, log and benchmark files of rule compute_gap_qc contain the same wildcards.
```

### Correct pattern: two-stage architecture

**Stage 1 — checkpoint: execution coordinator + manifest only**

```python
checkpoint compute_gap_qc:
    input:
        ice_done = lambda w: checkpoints.extract_and_balance_sample.get(
            sample=w.sample
        ).output.manifest,
        script = SCRIPT_DIR + "/scripts/compute_gap_qc.py"
    output:
        touch = OUTDIR + "/logs/gap_qc/{sample}.done"
    shell:
        r"""
        # ... computation ...
        touch "{output.touch}"
        """
```

**Stage 2 — ordinary rule: declares real output files**

```python
# Declares H5 files produced by sliding_window checkpoint.
# Downstream rules depend on this rule's output to get mtime tracking.
rule sliding_window_h5:
    input:
        lambda w: checkpoints.sliding_window.get(
            sample=w.sample
        ).output.manifest
    output:
        h5_files = expand(
            OUTDIR + "/windows/{{sample}}/{{sample}}__{{chrom}}_{{stage}}.h5",
            chrom=CHROMS,
            stage=["ln_oe", "gap_filled", "zscore", "gaussian"]
        )
```

**Stage 3 — downstream: depends on the ordinary rule**

```python
rule compare_gaussian_smoothed_windows:
    input:
        h5_a = lambda w: rules.sliding_window_h5.output.h5_files.get(
            sample=get_pair_samples(w.comparison)[0]
        ),
```

### Why this works

1. `checkpoint` only has `{sample}` in its output → wildcard set = `{sample}` → valid
2. Ordinary `rule sliding_window_h5` depends on checkpoint's `.done` via lambda input
3. H5 files declared in `rule sliding_window_h5.output.h5_files` → Snakemake tracks their mtime
4. Downstream rule `compare_gaussian_smoothed_windows` depends on `h5_a` → Snakemake knows if H5 files change

### Checklist before modifying any checkpoint output

- [ ] Does the checkpoint's output contain `expand()` that introduces extra wildcards (e.g., `{chrom}`)?
  - **Yes** → Must use two-stage architecture; cannot add directly to checkpoint output
  - **No** → Can add real files directly to checkpoint output
- [ ] Does the Python script internally construct paths without CLI args?
  - **Yes** → Downstream rule's `input:` declaration is for mtime tracking only; script keeps constructing paths internally
  - **Yes + want content tracking** → Must also modify the Python script to accept explicit `--xxx-files` args

---

## 17c. Real Output Files + Python Script Parameter Passing

### Goal

Let Snakemake track every file that a Python script actually produces, so manual edits or script changes trigger downstream reruns.

### Pattern A: input declaration for tracking only (script constructs paths internally)

Use when: script logic already constructs paths internally; you only want Snakemake to track file existence.

```python
# rule: list files in input, pass as shell arg
input:
    qc_a = lambda w: expand(
        OUTDIR + "/windows/{sample}__{chrom}_window_qc.tsv",
        sample=get_pair_samples(w.comparison)[0], chrom=CHROMS
    ),

shell:
    r"""
    QC_A=""
    for f in {input.qc_a}; do
        QC_A="$QC_A $f"
    done
    python {input.script} --qc-a-files $QC_A ...
    """
```

```python
# script.py: accept explicit file list via argparse
parser.add_argument(
    "--qc-a-files", nargs="+", default=None,
    help="Explicit QC file list (from Snakemake input)"
)

def compute_covariates(..., qc_a_files=None):
    if qc_a_files:
        for f in qc_a_files:
            if os.path.exists(f):
                df = pd.concat([df, pd.read_csv(f, sep="\t")])
    else:
        # fallback: auto-construct paths from outdir
        for chrom in chroms:
            path = f"{outdir}/windows/{sample}__{chrom}_window_qc.tsv"
```

### Pattern B: fully pass via CLI args (recommended for new rules)

Use when: writing a new rule / willing to modify the script.

```python
shell:
    r"""
    python {input.script} \
        --qc-a-files {input.qc_a} \
        --qc-b-files {input.qc_b}
    """
```

```python
# script.py: argparse with nargs="+"
parser.add_argument("--qc-a-files", nargs="+", required=True)
parser.add_argument("--qc-b-files", nargs="+", required=True)
```

### Checklist when changing output declarations

- [ ] Do downstream Python scripts need new `--xxx-files` args?
- [ ] Does the script prefer explicit file lists over auto-constructed paths?
- [ ] After adding a file to `output:`, are all downstream `input:` declarations updated?

---

## 17d. Pre-Modification Checklist for All Rules

Before committing any change to `rules/*.smk`, check:

**output declaration**
- [ ] Does this checkpoint's output contain `expand()` introducing extra wildcards?
  - **Yes** → Use two-stage architecture
- [ ] What real files does this rule produce? Are they constructed internally by the script or declared in `output:`?
  - **Internal construction** → Downstream rules need explicit `input:` for mtime tracking
- [ ] Do all real output files need mtime tracking for downstream reruns?
  - **Yes** → Must be in downstream rule's `input:`

**`expand()` + `.format()` mixing**
- [ ] Does the output/input path string contain `.format()`?
  - **Yes** → Does the wildcard exist in `.format()` kwargs or in `expand()` kwargs?
    - In `.format()` kwargs only → single brace `{}`
    - In `expand()` kwargs only → double brace `{{}}`
    - In both → single brace `{}` (expand resolves first)

**Python script parameters**
- [ ] Are new input files passed via `--` args or constructed internally by the script?
- [ ] Does the script need new `--xxx-files` arguments?

**Cross-rule consistency**
- [ ] After changing this rule's `output:`, are all downstream `input:` declarations consistent?
- [ ] Do Python scripts need CLI arg updates?

---

## 18. Common Bug Patterns

### Bug: NameError in standalone mode (shell directive)

**Problem:** A Python script uses `config.get()` or `snakemake.config.get()` assuming Snakemake will inject these objects. However, when the rule uses `shell:` directive (instead of `script:`), Snakemake does NOT inject `snakemake` or `config` objects.

**Symptom:**
```
NameError: name 'config' is not defined
```

**Root Cause:** The script's `main()` function calls `config.get("rawdata_dir", "")` directly, but `config` is only defined inside the `if "snakemake" in dir():` block. When running via `shell:` directive (standalone mode), this code path is taken and `config` is undefined.

**Affected Code:**
```python
def main(args=None):
    args = sys.argv[1:]
    ap = argparse.ArgumentParser()
    # ...
    # BUG: config is not defined here!
    mcool_path = _read_manifest_mcool_path(
        parsed.manifest, parsed.sample_id,
        config.get("rawdata_dir", "")  # NameError!
    )
```

**Fix:** Define `config = {}` in `main()` before any `config.get()` calls:

```python
def main(args=None):
    if args is None: args = sys.argv[1:]
    # In shell mode, config is not injected by Snakemake.
    # Define an empty dict so that config.get() calls don't raise NameError.
    config = {}
    ap = argparse.ArgumentParser()
    # ...
```

**Prevention:**
- Use `script:` directive instead of `shell:` for Python scripts that need Snakemake context
- If using `shell:`, ensure the script can run standalone without Snakemake injection
- Add a lint check: `grep -n "config\." workflow/scripts/*.py` to find potential issues

### Snakemake injection modes

| Directive | `snakemake` object | `config` object |
|-----------|-------------------|-----------------|
| `script:` | ✅ Injected | ✅ Injected |
| `shell:` | ❌ Not injected | ❌ Not injected |
| `wrapper:` | ✅ Injected | ✅ Injected |

**Rule:** If a script uses `snakemake` or `config`, use `script:` directive.

