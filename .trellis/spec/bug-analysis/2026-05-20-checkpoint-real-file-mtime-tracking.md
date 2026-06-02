# Bug Analysis: checkpoint Output Real File Mtime Not Tracked (2026-05-20)

## 1. Root Cause Category

- **Category**: C — Change Propagation Failure（变更传播失败）
- **Specific Cause**: checkpoint 的 `output:` 策略性地只用 `.done` marker，而 Python 脚本在内部独立构造输出文件路径（不走 Snakemake `input:` 传递），导致 Snakemake DAG 无法感知下游 Python 脚本实际产出的文件的 mtime 变化。

## 2. Why Fixes Failed

1. **First attempt — 直接在 checkpoint output 加 expand() 文件**：`output: manifest + gaps = expand(...)` 导致 Snakemake 报错 `Not all output ... contain the same wildcards` — checkpoint 要求所有 output 共享同一 wildcard 集合，`{sample}` 和 `expand({chrom})` 不兼容
2. **Second attempt — `}{stage}` typo**：`.format()` 字符串中 `{{chrom}}_{{stage}}` 误写成 `{{chrom}}_{stage}`，导致单 `}` 无匹配
3. **Third attempt — 移除 checkpoint 的 gaps，只保留 manifest**：Snakemake 报错消失，但真实文件追踪问题未解决
4. **Final fix — 二级架构**：checkpoint 只保留 `manifest`（简单 wildcard），新增普通 rule `sliding_window_h5` 专门声明 H5 文件，checkpoint 触发后该 rule 立即感知

## 3. Prevention Mechanisms

| Priority | Mechanism | Specific Action | Status |
|----------|-----------|-----------------|--------|
| P0 | Rule Template | 新增 Section 17b：checkpoint 二级架构模式及检查清单 | DONE |
| P0 | Rule Template | 新增 Section 17c：真实文件 + Python 脚本参数传递模式 | DONE |
| P0 | Rule Template | 新增 Section 17d：预修改检查清单 | DONE |
| P1 | Bug Analysis Doc | 创建本分析文档 | DONE |
| P2 | Trellis Sync | 同步模板到 `src/templates/markdown/spec/` | TODO |

## 4. Systematic Expansion

- **Similar Issues**: `extract_expected.smk`、`compute_oe_matrix.smk` 等其他 checkpoint 是否也有同样问题？Python 脚本内部是否也独立构造了路径字符串？
- **Design Improvement**: checkpoint 只做"执行协调"，真实文件声明全部在普通 rule 中——符合 Snakemake 的"每个真实文件都有对应 output"哲学
- **Process Improvement**: 写 rule 前先画数据流图（Python 脚本读什么，写什么），再写 `input:` / `output:`

## 5. Knowledge Capture

### Files Changed

| 文件 | 改动 |
|------|------|
| `rules/compute_gap_qc.smk` | `output:` 新增 `gaps` + `qc` (expand) |
| `rules/sliding_window.smk` | checkpoint 只保留 `manifest`；新增 `rule sliding_window_h5` |
| `rules/window_compare.smk` | `input:` 新增真实文件引用，移除 `checkpoint_done` lambda |
| `scripts/compute_window_pair_covariates.py` | 新增 `--qc-a-files` + `--qc-b-files` 参数 |
| `scripts/compare_gaussian_smoothed_windows.py` | 新增 `--h5-a-files` + `--h5-b-files` 参数 |

### Bug Analysis Document Created

- [x] Created `.trellis/spec/bug-analysis/2026-05-20-checkpoint-real-file-mtime-tracking.md` (this file)

### Spec Updates

- [x] Updated `.trellis/spec/pipeline/snakemake-guidelines.md` 新增 Section 17b, 17c, 17d
- [ ] 检查其他 checkpoint rule（`extract_expected.smk`、`compute_oe_matrix.smk`）是否有相同问题
- [ ] Trellis sync to `src/templates/markdown/spec/`
