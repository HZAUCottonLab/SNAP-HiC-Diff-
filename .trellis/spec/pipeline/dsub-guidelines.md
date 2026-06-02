# 多瑙调度器 (dsub) — 华中农业大学超算中心

> dsub CLI v1 | HPC 集群调度器 | Snakemake 集成参考

## Table of Contents

1. [认证与初始化](#1-认证与初始化)
2. [作业提交 (dsub)](#2-作业提交-dsub)
3. [脚本提交 (推荐)](#3-脚本提交-推荐)
4. [DSUB 指令速查](#4-dsub-指令速查)
5. [作业队列](#5-作业队列)
6. [资源申请规则](#6-资源申请规则)
7. [内存限制](#7-内存限制)
8. [作业依赖 (-D)](#8-作业依赖--d)
9. [交互作业](#9-交互作业)
10. [MPI 并行作业](#10-mpi-并行作业)
11. [NPU 作业](#11-npu-作业)
12. [批量作业提交](#12-批量作业提交)
13. [查看作业 (djob)](#13-查看作业-djob)
14. [作业控制](#14-作业控制)
15. [日志查看](#15-日志查看)
16. [节点查看 (dnode)](#16-节点查看-dnode)
17. [Snakemake 集成](#17-snakemake-集成)
18. [常见故障](#18-常见故障)

---

## 1. 认证与初始化

> **首次使用前必须执行**

```bash
dlogin
# 输入用户名密码获取 token
```

---

## 2. 作业提交 (dsub)

### 命令行直接提交

```bash
dsub -n jobname -q arm -aa -R 'cpu=4' -o jobname_%J.out -e jobname_%J.err sleep 100
```

### 常用选项

| 选项 | 说明 | 示例 |
|------|------|------|
| `-n` | 作业名称 | `-n align_bwa` |
| `-q` | 作业队列 | `-q arm`, `-q x86` |
| `-aa` | 任意 CPU 架构 | 与 `-q` 必须同时使用 |
| `-R 'cpu=N'` | 申请 CPU 核数 | `-R 'cpu=40'` |
| `-R 'npu=N'` | 申请 NPU 资源 | `-R 'npu=1'` |
| `-o FILE` | stdout 日志（追加模式） | `-o job_%J.out` |
| `-oo FILE` | stdout 日志（覆盖模式） | `-oo job.out` |
| `-e FILE` | stderr 日志路径 | `-e job_%J.err` |
| `-N N` | MPI 总进程数 | `-N 256` |
| `--mpi TYPE` | MPI 类型 | `--mpi hmpi` |
| `-nl 'NODE'` | 指定运行节点 | `-nl 'agent-ARM-01'` |
| `-nl 'RANGE'` | 节点范围 | `-nl 'agent-ARM-0[1-9]'` |
| `-I` | 交互模式 | `-I bash` |
| `-w` | 阻塞提交（等待结束） | `-w` |

### 日志文件名变量

| 变量 | 含义 |
|------|------|
| `%U` | 用户名 (userName) |
| `%J` | 作业 ID (JobID) |
| `%G` | TaskGroupName |
| `%I` | Index |
| `%A` | UserName_JobID_TaskGroupName_Index |

---

## 3. 脚本提交 (推荐)

### 脚本格式

```bash
#!/bin/bash
#DSUB -n jobname
#DSUB -R 'cpu=2'
#DSUB -o jobname_%J.out
#DSUB -aa
#DSUB -q arm

date
echo "this is script job"
sleep 10
date
```

### 提交方式

```bash
chmod +x submit.sh
dsub -s submit.sh
```

> **注意**：脚本必须有可执行权限 (`chmod +x`)。

---

## 4. DSUB 指令速查

| 指令 | 说明 | 示例 |
|------|------|------|
| `#DSUB -n NAME` | 作业名 | `#DSUB -n bwa_align` |
| `#DSUB -R 'cpu=N'` | CPU 核数 | `#DSUB -R 'cpu=40'` |
| `#DSUB -R 'npu=N'` | NPU 数量 | `#DSUB -R 'npu=1'` |
| `#DSUB -o FILE` | stdout 日志（追加） | `#DSUB -o align_%J.out` |
| `#DSUB -oo FILE` | stdout 日志（覆盖） | `#DSUB -oo align.out` |
| `#DSUB -e FILE` | stderr 日志 | `#DSUB -e align_%J.err` |
| `#DSUB -aa` | 任意架构 | `#DSUB -aa` |
| `#DSUB -q QUEUE` | 作业队列 | `#DSUB -q arm` |
| `#DSUB -N N` | MPI 总进程数 | `#DSUB -N 256` |
| `#DSUB -nn N` | 每节点进程数 | `#DSUB -nn 2` |
| `#DSUB --mpi TYPE` | MPI 类型 | `#DSUB --mpi hmpi` |
| `#DSUB -nl NODE` | 指定节点 | `#DSUB -nl 'agent-ARM-01'` |

---

## 5. 作业队列

### 查看队列

```bash
dqueue
```

### 队列与节点对应关系

| 队列 | 节点 | CPU 核数 | 每核内存限制 | 适用场景 |
|------|------|---------|------------|---------|
| `arm` | agent-ARM-{01-43} | 128 | 5 GB | 普通串行/多线程作业（推荐） |
| `x86` | agent-X86-{01-65} | 64 | 10 GB | 普通串行/多线程作业 |
| `fat_arm` | fat-agent-ARM-{01-03} | 128 | 16 GB | 大内存作业 (>500 GB) |
| `fat_x86` | fat-agent-X86-{01-03} | 64 | 32 GB | 大内存作业 (>500 GB) |
| `AI` | AI-01 | 192 | 10 GB | AI/NPU 作业（需申请权限） |
| `interactive` | — | — | — | 交互调试 |
| `default` | 混合 | — | — | 默认队列（不推荐） |

> **重要**：`#DSUB -q` 和 `#DSUB -aa` 必须同时使用，否则作业无法调度到计算节点。

---

## 6. 资源申请规则

### 程序类型

- **串行程序**：只能在单个节点运行
  - 单线程程序：只用 1 个 CPU 核
  - 多线程程序：可使用多个核（需用专用选项设置线程数）
- **并行程序 (MPI)**：可跨节点运行（如 vasp、gromacs、lammps）

### 核心原则

> **申请的 CPU 核数应与程序使用的线程数相等**，避免资源浪费或节点压力过大。

### 一般规则

- 串行作业申请的 CPU 核数 **不能超过节点核数**（如 arm 队列 ≤ 128 核）
- 内存需求 **< 500 GB**：通过申请更多 CPU 核来满足（arm 每核 5 GB）
- 内存需求 **≥ 500 GB**：使用 `fat_x86` 或 `fat_arm` 队列
- 内存需求 **< 500 GB** 的作业**不要**使用 fat 队列，以免排队时间过长

---

## 7. 内存限制

系统按核数强制限制每个作业的最大内存。超过限制会被系统杀掉，日志中提示：

```
EXIT_MESSAGE: Job was terminated due to reaching the mem limit, errCode: 143
```

### 各队列内存限制

| 队列 | 每核内存限制 | 计算方式 |
|------|------------|---------|
| arm | 5 GB | 需要 50 GB → 申请 `cpu=10` |
| x86 | 10 GB | 需要 50 GB → 申请 `cpu=5` |
| AI | 10 GB | — |
| fat_x86 | 32 GB | — |
| fat_arm | 16 GB | — |

---

## 8. 作业依赖 (-D)

用于指定作业间的执行顺序依赖。

### 依赖条件格式

```bash
dsub -D "key=value" [command]
```

- **key**：`jobid` 或 `JobName`
- **value**：状态 — `STARTED` / `RUNNING` / `SUCCEEDED` / `FAILED` / `ENDED`
  - `STARTED` = `RUNNING` + `STOPPED` + `SSTOPPED` + `SUCCEEDED` + `FAILED`
  - `ENDED` = `SUCCEEDED` + `FAILED`

### 作业 ID 约束

- 1~12 位正整数
- 支持同一 jobid 指定多种状态

### JobName 约束

- 特殊字符 `%=()|` 需用 `\` 转义
- 支持 `*` 通配符（开头=后缀匹配，中间=关键字，结尾=前缀匹配）
- 只能依赖**本用户**提交的作业

### 逻辑组合

支持 `&&`（与）和 `||`（或），或用 `;` 表示与；`&&、||、;` 不能混用。

### 纯数字名称

- 不加单引号 → 按 jobid 处理
- 加单引号 → 按 JobName 处理

### 示例

```bash
# 按 jobid 指定单个依赖
dsub -D "1=RUNNING" "echo 'hello'"

# 按 JobName 指定
dsub -D "job1=RUNNING" "echo 'hello'"

# 多作业依赖（与/或）
dsub -D "(1=RUNNING||2=SUCCEEDED)&&3=ENDED" "echo 'hello'"

# 通配符匹配
dsub -D "job_*=SUCCEEDED" "echo 'hello'"

# 特殊字符转义
dsub -D "job\(202306\)=SUCCEEDED" "echo 'hello'"

# 纯数字 JobName（加单引号）
dsub -D "'12'=RUNNING" "echo 'hello'"
```

---

## 9. 交互作业

```bash
# 普通交互作业
dsub -q interactive -I bash

# 进入 AI 节点（需向管理员申请权限）
dsub -aa -q AI -I bash
```

---

## 10. MPI 并行作业

### 必需指令

```bash
#DSUB -N 256       # 总进程数
#DSUB -nn 2        # 每节点进程数
#DSUB --mpi hmpi   # MPI 类型：openmpi / intelmpi / hmpi / mpich
```

### 环境变量

MPI 作业中 `$CCS_MPI_OPTIONS` 包含 `mpirun` 所需的 hostfile 参数。

```bash
# 打印 hostfile（用于 debug）
echo $CCS_MPI_OPTIONS
hostfile=$(echo $CCS_MPI_OPTIONS | sed 's/-hostfile //g')
cat $hostfile

# 运行
mpirun $CCS_MPI_OPTIONS -x OMP_NUM_THREADS=1 lmp_mpi -in in.lj
```

---

## 11. NPU 作业

### 申请 NPU 资源

```bash
#!/bin/bash
#DSUB -n jobname
#DSUB -R 'npu=1'
#DSUB -o jobname_%J.out
#DSUB -aa
#DSUB -q AI

python npu.py
```

### 查看 NPU 使用情况

```bash
dnode --npu
```

---

## 12. 批量作业提交

用于同一程序处理多个样本或不同参数，提高效率。

```bash
for sample in /share/home/username/work/lsf_bwait/raw_data/*_R1.fastq.gz; do
    index=$(basename $sample | sed 's/_R1.fastq.gz//')
    prefix=$(dirname $sample)
    dsub -n "bwa_${prefix}" -aa -q arm -R 'cpu=40' \
         -o "bwa_${prefix}_%J.out" \
         "module load arm/bwa/0.7.18 arm/samtools/1.21; \
          bwa mem -t 40 -R '@RG\tID:${prefix}\tPL:illumina\tLB:library\tSM:${prefix}' \
          hg38.fa ${prefix}_R1.fastq.gz ${prefix}_R2.fastq.gz | \
          samtools sort -@ 40 -o ${prefix}_srt.bam"
done
```

---

## 13. 查看作业 (djob)

### 基本使用

```bash
djob           # 查看所有作业
djob jobid     # 查看指定作业
djob -l jobid  # 查看调度详情（排队原因）
djob -ll jobid # 查看详细运行信息（CPU 利用率、内存消耗）
```

### 字段说明

| 字段 | 含义 |
|------|------|
| ID | 作业 ID |
| NAME | 作业名称 |
| STATE | 状态 |
| USER | 用户名 |
| QUEUE | 所在队列 |
| START_TIME / END_TIME | 开始/结束时间 |
| EXEC_NODES | 执行节点 |

### 作业状态

| 状态 | 含义 |
|------|------|
| `WAITING` | 调度器已接收，等待调度 |
| `PENDING` | 排队等待资源 |
| `RUNNING` | 正在执行 |
| `STOPPED` / `SSTOPPED` | 已停止 |
| `SUCCEEDED` | 执行成功 |
| `FAILED` | 执行失败 |

### CPU 利用率（djob -ll）

`cpuUtil` 格式为 `avg/recent`：

- `avg`：从运行到采样的平均 CPU 利用率
- `recent`：最近一次上报的 CPU 利用率

> 建议关注实际 CPU 使用率，以准确申请核数，**避免浪费核时产生不必要费用**。

### 定制输出 alias

```bash
# 写入 ~/.bashrc
alias dbs="djob --output 'jobId:8 name:8 user:8 state:8 queue:8 startTime:20 execNodes:15 totalMaxMem:15 totalUtime:15 totalStime:15'"
```

---

## 14. 作业控制

| 命令 | 说明 |
|------|------|
| `dkill jobid` | 终止作业 |
| `dkill --force jobid` | 强制终止 |
| `djctl stop jobid` | 挂起作业 |
| `djctl resume jobid` | 恢复作业 |
| `djctl requeue jobid` | 重启作业 |
| `djctl resubmit jobid` | 重新提交作业 |

---

## 15. 日志查看

### 实时查看运行作业输出

```bash
dpeek
```

### 日志文件关键信息

日志末尾包含作业执行信息：

| 字段 | 含义 |
|------|------|
| `EXEC_NODE` | 计算节点 |
| `EXIT_CODE` | 退出码（0=成功） |
| `EXIT_MESSAGE` | 退出信息 |
| `MRUN_TIME` | 运行时间（秒） |
| `MEM_MAX` | 最大内存 |
| `USER_CPU_TIME` | 用户态 CPU 时间 |
| `SYSTEM_CPU_TIME` | 系统态 CPU 时间 |
| `CPU_UTIL_AVG(%)` | 平均 CPU 利用率 |

### 内存不足退出

```
EXIT_MESSAGE: Job was terminated due to reaching the mem limit, errCode: 130
```

---

## 16. 节点查看 (dnode)

```bash
dnode        # 查看所有节点及资源
dqueue       # 查看可用队列
```

### 节点状态

| 状态 | 说明 |
|------|------|
| `OK` | 正常，允许关闭 |
| `CLOSE_ADMIN` | 管理员关闭 |
| `CLOSE_LOCK` | 作业批量失败被隔离 |
| `CLOSE_BUSY` | 资源使用达阈值，停止下发新作业 |
| `CLOSE_FAULT` | 局部故障（通信/网卡/磁盘），需管理员处理 |
| `UNAVAILABLE` | 无法连接或严重故障 |
| `REGISTERING` | 节点注册中 |
| `UNLICENSED` | 无可用 License |
| `SUSPENDED_ADMIN` | 管理员手动休眠 |
| `SUSPENDED_AUTO` | 自动休眠 |
| `CLOSE_EXIT` | 即将删除（运行中作业不受影响） |
| `CLOSE_FULL` | CPU/GPU 资源全部分配完 |

---

## 17. Snakemake 集成

### Snakemake v8 之前的版本

```bash
snakemake --cluster 'dsub -R "cpu={threads}"' -s Snakefile
```

### Snakemake v8+ (cluster-generic executor)

本项目使用 `dsub_submit.sh` 包装脚本，通过 Snakemake 的 `cluster-generic` executor 调用：

```bash
snakemake --executor cluster-generic \
          --cluster-generic-submit-cmd 'bash dsub_submit.sh {threads} {rulename}' \
          -s Snakefile
```

### dsub_submit.sh 核心逻辑

1. 接收 Snakemake 传来的 `{threads}`, `{rulename}` 参数（插件自动追加 `$3` = jobscript 路径）
2. 设置环境变量抑制嵌套并行（OMP_NUM_THREADS=1 等）
3. 将原 jobscript 包装进一个自删除脚本（避免共享文件系统残留）
4. 解析 jobscript 为绝对路径
5. 调用 `dsub` 提交作业：

```bash
dsub \
    --cwd "$_PROJECT_ROOT" \
    -n "smk_${_RULENAME}" \
    -aa \
    -q arm \
    -e "${_LOG_ROOT}/${_RULENAME}_%J.err" \
    -o "${_LOG_ROOT}/${_RULENAME}_%J.out" \
    -R "cpu=${_THREADS}" \
    -- "$_WRAPPER"
```

> **⚠ Important distinction:** `--cwd` sets the shell wrapper's working directory, but the Snakemake jobscript (Python) may resolve relative paths from a different cwd depending on the cluster's job script execution context. **Always pair `--cwd` in dsub_submit.sh with `--directory` in the Snakemake invocation to guarantee all relative paths resolve correctly on remote nodes.**


### Snakemake profile 配置

```yaml
# config/profiles/dsub/config.yaml
executor: cluster-generic
jobs: 200                    # 最大并发 dsub 作业数
keep-going: true
rerun-incomplete: true
cluster-generic-submit-cmd: "bash {DEPLOYMENT_SOURCE}/dsub_submit.sh {threads} {rulename}"
```

---

## 18. 常见故障

### 作业提交后立即结束，无日志输出

**原因**：节点根目录写满

**排查**：

```bash
djob -ll jobid
```

查找关键报错：

```
Job execution environment build failed, err: container save failed failureHost: agent-ARM-03
```

**解决**：联系管理员清理节点根目录

### 内存超限被杀掉

日志中：

```
EXIT_MESSAGE: Job was terminated due to reaching the mem limit, errCode: 143
```

**解决**：增加申请核数（arm 每核 5 GB）或换用 fat 队列

---

## 19. 测试作业约束 (强制)

> **所有测试作业必须提交到服务器执行，禁止在本机（登录节点）直接运行测试。**

### 触发条件

当需要运行任何测试（单元测试、集成测试、端到端测试）时，**必须**使用 dsub 提交脚本方式执行，不得在本地终端直接运行。

### CPU 规格分级

根据测试规模选择对应 CPU 核数，**默认提交到 `arm` 队列**：

| 规模 | CPU 核数 | `#DSUB -R` 参数 | 适用场景 |
|------|---------|-----------------|---------|
| 小测试 | 4 核 | `-R 'cpu=4'` | 单个模块的单元测试、快速回归 |
| 中测试 | 8 核 | `-R 'cpu=8'` | 多模块集成测试、较大数据集 |
| 大测试 | 16 核 | `-R 'cpu=16'` | 全量回归测试、端到端测试 |

> **不需要手动设置内存参数**，服务器按队列比例自动分配：
> - `arm` 队列：每核 5 GB（小测试 20 GB，中测试 40 GB，大测试 80 GB）
> - `x86` 队列：每核 10 GB

### 提交流程（必须按顺序执行）

```bash
# Step 1: 编写测试脚本（如 tests/run_all.sh）
# Step 2: 赋予可执行权限
chmod 777 tests/run_all.sh
# Step 3: 提交到服务器
dsub -s tests/run_all.sh
```

### 项目级作业提交流程

```bash
# Step 1: 赋予可执行权限
chmod 777 run_snakemake.sh
# Step 2: 通过 dsub 提交到服务器（作业名 smk_snakemake_run）
dsub -s run_snakemake.sh
```

> **注意**：项目级运行脚本和测试脚本都必须通过 `dsub -s <script>` 提交，**禁止在本机直接运行 `bash run_snakemake.sh`**。

### 测试脚本模板

```bash
#!/bin/bash
#DSUB -n unittest_<项目名>
#DSUB -R 'cpu=8'         # 根据规模选择：4(小)/8(中)/16(大)
#DSUB -o logs/test_%J.out
#DSUB -e logs/test_%J.err
#DSUB -aa
#DSUB -q arm             # 默认 arm 队列

# 设置环境（如需要）
source activate snakemake

# 运行测试（pytest/pytest-xdist 等）
pytest tests/ -v
```

### 常见错误

| 错误 | 原因 | 正确做法 |
|------|------|----------|
| 在本机直接运行 `pytest` | 违反约束 | 使用 `dsub -s` 提交 |
| 脚本无执行权限 | 未执行 `chmod 777` | 提交前必须赋予权限 |
| 手动指定 `--mem` | 不必要且可能浪费资源 | 依赖服务器自动分配 |
| 小测试申请 16 核 | 资源浪费 | 按规模分级：小=4核，中=8核 |
| 忘记 `-q arm` | 默认队列可能调度慢 | 必须加 `-q arm` |

### 为什么禁止在本机运行测试

1. **资源限制**：登录节点资源有限，长时间测试会影响其他用户
2. **环境一致性**：服务器节点与本地环境一致，测试结果更可靠
3. **并行效率**：服务器可申请多核并行，测试速度更快
