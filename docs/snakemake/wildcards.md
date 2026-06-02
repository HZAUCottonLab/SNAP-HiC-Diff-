# Snakemake 通配符与 Shell 占位符规范

> 本文档是 `rules/*.smk` 中 Snakemake wildcard、params、bash 变量、转义变量的权威规范。
>
> **核心原则：params 只用来存 Python 常量，wildcard 的值统一在 shell 开头提取。**

---

## 1. 总览

Snakemake 的 `shell:` block 涉及三类变量的交织：

| 层次 | 示例 | 展开时机 | 展开者 |
|------|------|----------|--------|
| Snakemake wildcard | `{wildcards.sample}` | 执行前 | Snakemake |
| Snakemake params/input/output | `{params.X}`, `{input.X}`, `{output.X}`, `{threads}` | 执行前 | Snakemake |
| Bash 变量 | `$SAMPLE`, `${{CHR}}` | 执行时 | Bash |

---

## 2. 写法分类表

| # | 类别 | 推荐写法 | 展开者 | 注意事项 |
|---|------|----------|--------|----------|
| 1 | 路径 wildcard pattern | `{sample}` / `{chrom}` | Snakemake | 双大括号 `{{sample}}` 在 `.format()` 中转义 |
| 2 | shell 中访问 wildcard | `{wildcards.sample}` | Snakemake | **唯一推荐方式**，shell block 开头提取 |
| 3 | shell 中访问 params | `{params.X}` | Snakemake | 仅存 Python 常量 |
| 4 | shell 中访问 input/output | `{input.X}`, `{output.X}` | Snakemake | named input/output |
| 5 | shell 中访问 threads | `{threads}` | Snakemake | Python 脚本并行参数必须来自此处 |
| 6 | bash 变量 | `$VAR` / `"$VAR"` | Bash | 大写命名，加引号 |
| 7 | 转义 bash 变量 | `${{VAR}}` | Snakemake→Bash | 防止变量边界歧义 |
| 8 | lambda wildcards | `lambda w: w.sample` | Python | 动态路径构造 |
| 9 | 全局常量 | `{OUTDIR}`, `{SCRIPT_DIR}` | Snakemake | Snakefile 顶部定义的常量 |

---

## 3. 每类写法的适用场景

### 3.1 路径 wildcard pattern（类别 1）

用于 `input:`、`output:`、`log:` 的路径声明。

```python
# 单 wildcard
output:
    done = OUTDIR + "/logs/{sample}.done"

# .format() + 双大括号转义
output:
    manifest = "{}/logs/expected/{{sample}}_cooltools.done".format(OUTDIR)
```

### 3.2 shell 中访问 wildcard（类别 2）

**在 shell block 开头统一提取，不要裸写 `{sample}`。**

```bash
shell:
    r"""
    set -euo pipefail

    SAMPLE="{wildcards.sample}"
    CHR="{wildcards.chrom}"

    # 后续使用 $SAMPLE, $CHR
    """
```

### 3.3 params 存 Python 常量（类别 3）

**params 只存 Python 常量，不存 wildcard 值。**

```python
params:
    resolution = RESOLUTION,           # Python 常量
    chroms = " ".join(CHROMS),        # 派生常量
    window_size = WINDOW_SIZE,         # Python 常量
    sigma = config["smoothing"]["sigma"]  # config 值

# ❌ 错误：params.sample = lambda w: w.sample
```

### 3.4 threads 用于并行参数（类别 5）

**Python 脚本的并行参数必须从 `{threads}` 获取，禁止硬编码。**

```python
threads: config["parallel"]["threads_per_job"]

shell:
    r"""
    python {SCRIPT_DIR}/scripts/sliding_window.py \
        --n-jobs {threads}    # ✓ 正确
        # --n-jobs 80         # ❌ 硬编码，禁止
    """
```

### 3.5 转义 bash 变量 `${{VAR}}`（类别 7）

Snakemake 会把 `${{VAR}}` 转成 `${VAR}` 再交给 bash。

| 变量后接字符 | 写法 | 原因 |
|-------------|------|------|
| 空格或引号 | `$SAMPLE` | 边界清晰 |
| `_` | `${{SAMPLE}}` | 防止 `$SAMPLE_` 被识别为变量 |
| 字母/数字 | `${{CHR}}` | 防止 `$CHRSomething` 被误识别 |
| `.` | `$BASE` | `.` 不是变量名合法字符 |

```bash
# ✅ 正确
BASE="${{SAMPLE}}__${{CHR}}"
RAW="{params.tmpdir}/${{BASE}}.raw.cool"

# ❌ 错误：$SAMPLE__$CHR 会导致 bash 把 $SAMPLE__ 当作变量名
```

### 3.6 checkpoint output 必须用二级架构

**checkpoint 的 output 不能包含 `expand()` 引入的额外 wildcard（如 `{chrom}`）。**

正确做法：

```python
# Step 1 — checkpoint 只保留简单 wildcard 的 manifest
checkpoint compute_gap_qc:
    output:
        touch = OUTDIR + "/logs/gap_qc/{sample}.done"

# Step 2 — 普通 rule 声明真实文件
rule sliding_window_h5:
    input:
        lambda w: checkpoints.compute_gap_qc.get(sample=w.sample).output.touch
    output:
        h5_files = expand(
            OUTDIR + "/windows/{{sample}}/{{sample}}__{{chrom}}_{{stage}}.h5",
            chrom=CHROMS,
            stage=["ln_oe", "gap_filled", "zscore", "gaussian"]
        )
```

---

## 4. 不推荐写法

| 写法 | 原因 | 示例 |
|------|------|------|
| `{sample}` 在 shell block 中裸写 | 歧义：可能是 Snakemake 路径 pattern，也可能是 bash 变量 | ❌ `mkdir $SAMPLE` |
| `params.sample = lambda w: w.sample` | params 应该存常量，wildcard 在 shell 开头提取 | ❌ `params.sample = lambda w: w.sample` |
| 硬编码 CPU 核数 | threads 应该从 `{threads}` 获取 | ❌ `--n-jobs 80` |
| f-string 在 Snakemake rule 中 | 与 Snakemake `{}` 语法冲突 | ❌ `f"{OUTDIR}/logs/{sample}.done"` |
| `$SAMPLE__$CHR` | bash 会把 `$SAMPLE__` 当作变量名 | ❌ `file_$SAMPLE__$CHR.txt` |

---

## 5. Agent 编写 rule 时的检查清单

在提交任何 Snakemake rule 修改前，请逐项核对：

- [ ] **config 访问**：使用 `config["section"]["key"]` 直接访问，不使用 `.get()`
- [ ] **wildcard 提取**：shell block 开头有 `SAMPLE="{wildcards.sample}"` 或类似提取
- [ ] **params 用途**：params 只存 Python 常量，没有 `params.sample = lambda w: w.sample`
- [ ] **threads 来源**：Python 脚本的并行参数来自 `{threads}`，不是硬编码
- [ ] **bash 变量转义**：变量后接 `_` 或字母数字时使用 `${{VAR}}`
- [ ] **不使用 f-string**：使用 `.format()` 或字符串拼接
- [ ] **checkpoint 格式**：`output:` 声明 `touch`，shell 中 `touch "{output.touch}"`
- [ ] **checkpoint 依赖**：必须显式声明所有上游依赖到 `input:`，否则并行执行时竞态

---

## 附录：Snakemake vs Bash 变量对照

| Snakemake | Bash | 示例 |
|-----------|------|------|
| `{wildcards.sample}` | `$SAMPLE` | `SAMPLE="{wildcards.sample}"` |
| `{params.resolution}` | `$RESOLUTION` | `{params.resolution}` 直接传递 |
| `{input.mcool}` | `$INPUT` | `URI="{input.mcool}::/resolutions/..."` |
| `{output.manifest}` | `$OUTPUT` | `touch "{output.manifest}"` |
| `{threads}` | `$NJOBS` | `--n-jobs {threads}` |
| `{OUTDIR}` | `$OUTDIR` | `mkdir -p "{OUTDIR}/logs"` |
