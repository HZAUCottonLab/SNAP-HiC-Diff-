# CWorld Rscript 与 dsub wrapper 修复实施计划

**目标：** 让每个 `compute_IS` dsub 子任务从其独立 CWorld Conda 前缀获得 `Rscript`，并让 dsub 准确报告 jobscript 失败；验证流程能继续生成供 Meltron 使用的 bedGraph。

**约束：** 计算节点离线；CPU 数决定内存配额；不改变当前 `compute_IS` 的 4 CPU 配置；不触碰未提交的用户文件。

## 1. 写出失败回归测试

修改 `test/cworld_is/test_offline_conda_execution.py` 与 `test/cworld_is/test_compute_is_oom_script.py`，断言：

- `env/cworld.yaml` 需要 `r-base`；
- bootstrap 验证 CWorld 前缀中的 `Rscript`；
- dsub wrapper 在 jobscript 前启用 `set -euo pipefail`；
- 冒烟测试使用不继承环境 PATH 的受控 PATH，并检查 `Rscript`。

运行聚焦 pytest，确认这些新断言在修复前失败。

## 2. 补全并验证 CWorld 环境

在 `env/cworld.yaml` 添加最小依赖 `r-base`；在 `scripts/conda_envs.sh` 中将 Rscript 加入 CWorld prefix 完整性和可执行性检查。

运行聚焦 pytest、`bash -n scripts/conda_envs.sh`，确认绿色。

## 3. 修复 dsub 退出码传播

在 `config/dsub_submit.sh` 生成 wrapper 的 shebang 后写入 `set -euo pipefail`，确保 jobscript 的非零退出码不会被尾部清理命令覆盖。

再次运行静态测试及 `bash -n config/dsub_submit.sh`。

## 4. 让计算节点冒烟测试自包含

在 `test/compute_is_oom_20260609/run_test.sh` 中将 PATH 固定为 CWorld prefix 加系统基础命令目录，显式检查 prefix 内 `Rscript` 后再调用验证脚本。

运行该脚本的静态测试与 shell 语法检查。

## 5. 环境与集群验收

1. 登录节点执行 `bash scripts/conda_envs.sh bootstrap`，仅此一次重建因 YAML 哈希变化而产生的 CWorld 前缀。
2. 执行 `CONDA_OFFLINE=true bash scripts/conda_envs.sh verify` 和 `CONDA_OFFLINE=true snakemake --dry-run --use-conda`，确认不触发联网或重建。
3. 以 dsub 提交 `test/compute_is_oom_20260609/run_test.sh`，检查 CWorld 产出非空 bedGraph/TSV 且任务成功。
4. 仅在步骤 3 成功后提交 `dsub -s run_snakemake.sh`，观察至少十分钟并检查 controller、子任务状态与首个 `compute_IS` 日志。

**验收：** 真实 CWorld 计算产生 bedGraph；dsub 不再把失败任务报告为成功；主流程无 Rscript/环境/wrapper 类错误。
