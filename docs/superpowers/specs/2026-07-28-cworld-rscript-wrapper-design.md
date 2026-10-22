# CWorld Rscript 与 dsub 失败传播修复设计

## 背景与证据

作业 `38534914` 的首个 `compute_IS` 子任务在 CWorld 已完成 cool-to-matrix 转换后失败。`results/logs/compute_IS/H1397_compute_IS.log` 显示 `matrix2insulation.pl` 的依赖检查找不到 `Rscript`；当前 `.snakemake/conda/.../bin/` 确无该可执行文件。

这是环境不完整，不是 bedtools、输入矩阵、Meltron 或内存不足。该任务申请 4 CPU（约 20 GB 内存），矩阵加载峰值约 18.3 GB；Rscript 检查发生在真正 insulation 计算前。

此外，`config/dsub_submit.sh` 生成的 wrapper 在 Snakemake jobscript 失败后仍执行 `rm`，其退出码会掩盖 jobscript 的非零状态。因此失败的工作节点会被 dsub 误报为 `SUCCEEDED`。

## 已批准的最小修复

1. 在 `env/cworld.yaml` 显式加入 `r-base`，使 CWorld 前缀自包含 `Rscript`。
2. 在 `scripts/conda_envs.sh` 的 CWorld 前缀检查中验证可执行的 `Rscript` 并运行 `Rscript --version`。
3. 在 dsub 生成 wrapper 的 shebang 后加入 `set -euo pipefail`，保留 Snakemake jobscript 的失败退出码。
4. 将 CWorld 冒烟测试改为 `CWORLD_PREFIX/bin:/usr/bin:/bin` 的受控 PATH，并显式检查该前缀中的 `Rscript`。
5. 增加针对以上四项的回归测试；重建登录节点环境后，以离线模式验证，并在计算节点执行真实 CWorld 冒烟任务，随后再提交主流程观察至少十分钟。

## 明确不在本次范围内

- 不为当前问题提高 `compute_IS` 的 CPU 数。Rscript 修复通过后，才以真实 CWorld 运行的峰值内存决定是否另行调整。
- 不加入参考项目中未被当前 IS 产物需要的 `r-cairo` 或 ImageMagick。当前 vendored CWorld 缺少绘图 R 脚本；其绘图命令即使失败也不影响 pipeline 所需的 `.bedGraph`、`.insulation` 和 TSV。`r-base` 是满足 CWorld 强制依赖的最小项。
- 不修改 Meltron 算法、样本并行度或输入矩阵格式。

## 验收标准

- 新 CWorld Conda 前缀中 `bin/Rscript --version` 成功。
- 任何 jobscript 失败使其 dsub 工作节点状态为失败，而非 `SUCCEEDED`。
- 受控 PATH 的 CWorld 冒烟测试能生成非空 bedGraph/TSV。
- `CONDA_OFFLINE=true` 下 bootstrap、verify 与 Snakemake dry-run 不触发环境创建或联网。
- 真实计算节点冒烟任务成功后，主流程提交至少十分钟内无新的环境、Rscript 或 wrapper 失败。
