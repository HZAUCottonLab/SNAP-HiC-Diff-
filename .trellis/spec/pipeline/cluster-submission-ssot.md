# Cluster Submission SSOT — Snakemake + dsub 提交锁定规范

> **Purpose**: 锁定 Snakemake 集群提交方式，禁止通过硬编码或临时命令行绕过 SSOT 配置。
>
> **Scope**: `run_snakemake.sh`、`config/hpc/dsub_submit.sh`、所有 Snakemake 规则文件。
>
> **Status**: Active — 修改前必须审阅 REVIEW_CHECKLIST_CLUSTER.md

---

## 1. 核心原则

### 1.1 唯一真相来源 (SSOT)

```
config/hpc/dsub_submit.sh  ← 集群提交的单一 wrapper
run_snakemake.sh           ← 动态构建 snakemake 命令，通过 PROFILE_DIR 变量管理 wrapper 路径
```

**所有 snakemake 集群提交参数** 必须通过读取 `cluster_config.yaml` 注入到 `run_snakemake.sh`，禁止在 shell 脚本中硬编码 `--executor`、`--cluster-generic-submit-cmd` 等标志。

### 1.1.1 YAML 格式约束

`cluster_config.yaml` 中所有配置值**禁止使用 inline `#` 注释**：

```yaml
# ✅ 正确：无注释
latency_wait: 300
rerun_incomplete: true

# ❌ 错误：grep/sed 解析会吞掉 inline 注释
latency_wait: 300  # seconds to wait for NFS filesystem consistency
rerun_incomplete: true     # re-run incomplete jobs on resume
```

run_snakemake.sh 使用 `grep | sed` 解析值，`;s/\s*#.*//` 会正确剥离注释。但测试/调试脚本若缺少此步骤会导致：
```
snakemake: error: argument --latency-wait/-w: invalid int value: '300          # seconds to wait...'
```

所有描述性注释必须放在值上方的行内（如 SSOT 文件顶部），不能放在行尾。

### 1.2 Executor 写死

```yaml
# cluster_config.yaml — 禁止修改
executor: cluster-generic
```

`run_snakemake.sh` 启动时**强制验证**：
- 如果 `executor` 不是 `cluster-generic`，立即报错退出
- 如果 `submit_cmd` 不包含 `dsub_submit.sh`，立即报错退出

### 1.3 提交命令模板锁定

```bash
# run_snakemake.sh — 通过 PROFILE_DIR 变量管理
PROFILE_DIR="${PROFILE_DIR:-config/hpc}"
submit_cmd: "bash ${PROFILE_DIR}/dsub_submit.sh {threads} {rulename}"
```

有效 placeholder（Snakemake v8 cluster-generic）：
| Placeholder | 说明 |
|-------------|------|
| `{threads}` | CPU 核数（来自 rule 的 `threads:` 指令） |
| `{rule}`, `{rulename}` | 规则名称 |
| `{jobid}` | Snakemake 内部作业 ID |
| `{resources.*}` | 资源值（如 `{resources.mem_mb}`） |
| `{attempt}` | 作业尝试次数 |

**⚠️ `{jobscript}` 不是有效 placeholder**，Snakemake cluster-generic 插件自动追加 jobscript 路径作为第3个位置参数。

---

## 2. 文件清单与约束

### 2.1 SSOT 文件（只读，修改需审阅）

| 文件 | 约束 |
|------|------|
| `config/hpc/dsub_submit.sh` | **只可读取** — dsub 包装脚本，包含 conda 环境激活、嵌套并行抑制逻辑 |
| `run_snakemake.sh` | **只可读取** — 通过 PROFILE_DIR 变量管理 wrapper 路径 |
| `REVIEW_CHECKLIST_CLUSTER.md` | **审阅清单** — 修改 SSOT 文件前必须填写的检查项 |

### 2.2 执行文件

| 文件 | 说明 |
|------|------|
| `run_snakemake.sh` | 从 `cluster_config.yaml` 读取配置，动态构建 snakemake 命令 |
| `Snakefile` | 入口，通过 `--directory` 和绝对路径确保远程节点正确解析 |

### 2.3 规则文件

所有 `workflow/rules/*.smk` 必须：
- 使用绝对路径（`OUTDIR`、`LOGS_DIR`、`BENCHMARK_DIR`）
- 从 `config.get("threads", ...)` 加载线程数
- 定义完整的 `input/output/log/benchmark/params/threads/conda` 指令
- **禁止使用 `localrules:`**，除非该 rule 明确只需要本地执行

---

## 3. run_snakemake.sh 架构

```
┌─────────────────────────────────────────────────────────┐
│  run_snakemake.sh                                        │
│                                                          │
│  1. 读取 cluster_config.yaml (grep/sed 解析)             │
│  2. 验证 executor == "cluster-generic"                   │
│  3. 验证 submit_cmd 包含 "dsub_submit.sh"               │
│  4. 动态构建 snakemake 命令                              │
│     --executor "$_executor"                              │
│     --cluster-generic-submit-cmd "$_submit_cmd"          │
│     --cluster-generic-cancel-cmd "$_cancel_cmd"          │
│     --latency-wait "$_latency_wait"                      │
│     $([ "$_rerun_incomplete" = "true" ] && echo ...)    │
│  5. 执行 snakemake                                       │
└─────────────────────────────────────────────────────────┘
         ↓
┌─────────────────────────────────────────────────────────┐
│  config/hpc/dsub_submit.sh
│                                                          │
│  参数: $1=threads, $2=rulename, $3=jobscript           │
│  1. 设置 OMP_NUM_THREADS=1 等嵌套并行抑制               │
│  2. 设置 CONDA_ENVS_PATH 指向预创建环境                  │
│  3. 激活 conda 环境 snakemake-snaphic                   │
│  4. 包装 jobscript → 自删除 wrapper                      │
│  5. dsub 提交: -q arm, -R "cpu=${threads}"              │
└─────────────────────────────────────────────────────────┘
```

---

## 4. 审阅流程

### 4.1 修改前必须完成

1. 阅读 `REVIEW_CHECKLIST_CLUSTER.md`
2. 确认修改的文件在 SSOT 范围内
3. 确认 executor 仍为 `cluster-generic`
4. 确认 submit_cmd 仍调用 `dsub_submit.sh`
5. 用 `bash run_snakemake.sh --help` 验证配置加载

### 4.2 紧急回滚

```bash
# Step 1: 恢复 SSOT 文件
git restore config/hpc/dsub_submit.sh
git restore run_snakemake.sh

# Step 2: 验证恢复成功
bash run_snakemake.sh 2>&1 | head -20
```

---

## 5. Snakemake 规则约束

### 5.1 禁止 `localrules:` 滥用

只有以下情况允许 `localrules:`：
- Phase 聚合规则（如 `rule phase2`, `rule all`）
- 纯粹的数据汇总/报告规则（不调用外部脚本）
- DAG 目标定义规则

**禁止**：将实际计算 rule 设为 localrule，导致无法提交到集群。

### 5.2 线程数必须与实际并行匹配

```python
# 正确：线程数通过 config 声明
params:
    threads=config.get("threads", {}).get("extract_chr_mcool", 20),
threads: lambda wildcards, resources, params: params.threads

# 错误：硬编码线程数
threads: 20
```

线程数同时决定 dsub 的 CPU 核数和内存上限（arm 队列每核 5 GB）。

### 5.3 路径必须为绝对路径

```python
# 正确
OUTDIR = config.get("outdir", os.path.join(PROJECT_ROOT, "results"))

# 错误（相对路径在远程节点会失败）
output: "results/sample_1/output.txt"
```

---

## 6. 关联文档

| 文档 | 说明 |
|------|------|
| [snakemake-guidelines.md](./snakemake-guidelines.md) | Snakemake 语法、规则编写、资源管理 |
| [dsub-guidelines.md](./dsub-guidelines.md) | dsub CLI、队列限制、MPI/NPU 作业、Snakemake 集成 |
| `REVIEW_CHECKLIST_CLUSTER.md` | 修改 SSOT 前的审阅清单（项目根目录） |
| `config/hpc/dsub_submit.sh` | dsub 包装脚本（实际提交逻辑） |

---

## 7. CI/Audit 检查项

每次修改后自动检查：

```bash
# 检查是否有硬编码的 cluster-generic 配置（不应有）
grep -n 'cluster-generic-submit-cmd' run_snakemake.sh

# 检查 cluster_config.yaml 是否被正确读取
grep -n 'cluster_config.yaml' run_snakemake.sh

# 检查 executor 是否被验证
grep -n 'cluster-generic' run_snakemake.sh

# 检查 {jobscript} 误用（应不存在）
grep -rn '{jobscript}' run_snakemake.sh config/hpc/
```
