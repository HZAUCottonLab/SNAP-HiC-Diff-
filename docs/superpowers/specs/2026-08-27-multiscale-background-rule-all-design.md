# 多尺度 KNN 背景接入与完整 `rule all` 设计

日期：2026-08-27
状态：已由用户批准

## 目标

把 `background_pool` 已生成并验收通过的五个 4 kb KNN 背景表，以绝对路径配置到 SNAP-HiC-Diff。一次 `dsub` 提交构建一个 Snakemake DAG，在同一轮运行中展开五个窗口尺度，完成尺度独立的窗口比较、KNN 和 Meltron 分析，再汇总候选并调用现有 IS-change 流程。最终必须完成原有完整 `rule all`，包括 TAD 等其他默认目标。

## 已确认决策

- 活动矩阵分辨率为 4,000 bp/bin。
- 窗口为 `[25, 50, 75, 100, 125]` bins，即 100、200、300、400 和 500 kb。
- 步长为 `ceil(window_bins × 0.2)` bins，即 20、40、60、80 和 100 kb。
- 旧的 `window_size: 300000` 和 `step_size: 32000` 不再是活动配置，也不能满足新 DAG 的任何尺度目标。
- 背景文件保留在 `background_pool`；本项目只配置绝对路径，不复制数据、不建立软链接，也不引入背景生成规则。
- 五个尺度由一次 `dsub -s run_snakemake.sh` 提交的同一个 DAG 展开；禁止按尺度在外层循环运行或分别提交。
- 可以复用窗口尺度无关且仍有效的共享产物。所有窗口相关产物必须使用新尺度路径，不能复用旧单窗口结果。
- 任一窗口尺度通过现有 KNN 与 Meltron 联合阈值即可进入跨尺度候选并集，不要求至少两个窗口尺度支持。
- 有效空候选和空 IS-change 结果可以作为成功结果；不得为了非空而放宽统计阈值。

## 窗口换算

令 `R` 为分辨率、`W` 为窗口 bin 数、`F` 为步长比例：

```text
window_bp = R × W
step_bins = ceil(W × F)
step_bp   = R × step_bins
```

活动参数映射为：

| 路径标签 | 窗口 bins | 窗口 bp | 步长 bins | 步长 bp |
|---|---:|---:|---:|---:|
| `w25` | 25 | 100,000 | 5 | 20,000 |
| `w50` | 50 | 200,000 | 10 | 40,000 |
| `w75` | 75 | 300,000 | 15 | 60,000 |
| `w100` | 100 | 400,000 | 20 | 80,000 |
| `w125` | 125 | 500,000 | 25 | 100,000 |

换算只在配置边界完成。现有计算脚本继续接收明确的 bp 参数，不各自重新解释比例。

## 配置契约

`config/config.yaml` 使用一个窗口配置源和五个显式背景 preset：

```yaml
resolution: 4000

window_analysis:
  window_bins: [25, 50, 75, 100, 125]
  step_fraction: 0.2

knn_empirical_pvalue:
  background_presets:
    - resolution_bp: 4000
      window_bins: 25
      step_bins: 5
      file: /share/home/xhhuang24/project/PanTE/background_pool/results/multiscale_4k/window_scales/100kb/window_compare/HiC_window_background.tsv
    - resolution_bp: 4000
      window_bins: 50
      step_bins: 10
      file: /share/home/xhhuang24/project/PanTE/background_pool/results/multiscale_4k/window_scales/200kb/window_compare/HiC_window_background.tsv
    - resolution_bp: 4000
      window_bins: 75
      step_bins: 15
      file: /share/home/xhhuang24/project/PanTE/background_pool/results/multiscale_4k/window_scales/300kb/window_compare/HiC_window_background.tsv
    - resolution_bp: 4000
      window_bins: 100
      step_bins: 20
      file: /share/home/xhhuang24/project/PanTE/background_pool/results/multiscale_4k/window_scales/400kb/window_compare/HiC_window_background.tsv
    - resolution_bp: 4000
      window_bins: 125
      step_bins: 25
      file: /share/home/xhhuang24/project/PanTE/background_pool/results/multiscale_4k/window_scales/500kb/window_compare/HiC_window_background.tsv
```

配置解析必须拒绝非正整数分辨率或窗口、空或重复窗口、非有限或非正步长比例、仍存在的旧窗口键，以及活动分辨率下缺失、重复或未使用的 preset。每个推导出的 `(resolution_bp, window_bins, step_bins)` 必须恰好匹配一个 preset。

## 背景消费边界与校验

`background_pool` 负责背景队列和表生成；SNAP-HiC-Diff 负责 preset 选择、消费前校验、目标 comparison 的 KNN 校准和下游汇总。本项目不调用 `background_pool` 工作流。

每个尺度在 KNN 前生成一份小型校验报告。校验要求：

- 文件存在、可读且非空，并作为显式 Snakemake 输入参与依赖跟踪；
- provenance/坐标列 `comparison`、`chrom`、`window_start`、`window_end` 存在；
- KNN 所需列 `EC`、`IB`、`SP`、`MSE`、`SSIM`、`Spearman`、`tw_pvalue` 和 `NCC` 存在；
- 数值清理后至少保留 `knn_k` 行；
- 完整窗口长度等于对应 `window_bp`；染色体末端窗口只按生产端现有的至少半窗、染色体截断规则接受；
- `window_start` 位于对应 `step_bp` 网格，允许上游 QC 删除中间窗口；
- 不存在重复的 `(comparison, chrom, window_start, window_end)` 行；不同背景 comparison 复用相同坐标是合法的。

校验报告记录绝对路径、预期尺度元组、行数和观察到的坐标属性。缺失映射在 DAG 构建前失败；表内容错误在 KNN 分配大量资源前失败，并指出具体尺度和路径。2026-08-27 的只读验收已确认五个现有背景均通过 `background_pool/scripts/validate_multiscale_backgrounds.py`。

## 工作流架构与调度

Snakemake 使用一个尺度通配符从配置展开五个分支：

```text
manifest / comparisons
        |
        +--> convert / balanced / expected / O-E（共享）
        |
        +--> 十个 CWorld IS 距离（每样本共享）
        |
        +--> w25  --> QC --> windows --> compare --> KNN --> Meltron
        +--> w50  --> QC --> windows --> compare --> KNN --> Meltron
        +--> w75  --> QC --> windows --> compare --> KNN --> Meltron
        +--> w100 --> QC --> windows --> compare --> KNN --> Meltron
        +--> w125 --> QC --> windows --> compare --> KNN --> Meltron
        |
        +--> 五尺度候选拼接 --> 区间并集 --> 现有 IS-change caller
        |
        +--> 原有 TAD 和其他 rule-all 目标
```

只提交一次 `dsub -s run_snakemake.sh`。cluster-generic executor 把所有依赖已满足的作业交给 Donau，因此不同尺度的作业可以交错并发；实际同时运行的节点数由队列容量决定。设计不引入“完成一个尺度后再提交下一个尺度”的串行控制。

共享且依赖未变化的 convert、balanced、expected、O-E、CWorld IS 和 TAD 产物可以由 Snakemake 复用。尺度相关 QC、滑窗、比较、KNN、热图、Meltron、候选和日志必须使用新路径，旧的未限定单窗口路径不属于新 DAG。

## 路径与身份

尺度相关根目录为：

```text
results/window_scales/r4000/w25/
results/window_scales/r4000/w50/
results/window_scales/r4000/w75/
results/window_scales/r4000/w100/
results/window_scales/r4000/w125/
```

每个根目录保留现有逻辑阶段，并包含尺度独立的日志和 benchmark。稳定的跨阶段窗口身份为：

```text
(resolution_bp, window_bins, chrom, window_start, window_end)
```

`window_id` 只在单尺度内使用，不得单独用于跨尺度连接。KNN 与 Meltron 必须以完整尺度身份连接，并继续执行现有起点一致性和终点差异审计。

## 每尺度筛选

每个 `comparison × window_bins` 分支保持现有统计定义：

1. 从原始指标派生 `MSE`、`d_SSIM`、`d_Spearman`、`tw_score` 和 `d_NCC`。
2. 用精确匹配的尺度背景计算五个原始 KNN 经验 P 值。
3. 任一原始经验 P 值 `< 0.05` 时 KNN 通过。
4. `meltron_best_score >= 5` 时 Meltron 通过。
5. 两者同时通过时该窗口进入本尺度候选。

不改变原始阈值，不新增跨窗口尺度的多重检验校正。跨尺度并集属于候选筛选，不能描述为新的全局显著性结论。

## 跨尺度汇总和 IS-change

先无损拼接五个尺度的通过窗口，形成逐窗口证据表；再按染色体排序并合并重叠或 BED 相邻区间。单尺度支持即可保留。每个合并区域记录：

- `support_window_bins`：排序后的唯一支持窗口 bins；
- `support_scale_count`：支持窗口尺度数；
- `source_window_count`：来源窗口数；
- 尺度限定的来源窗口与 Meltron 标识；
- `best_meltron_score`；
- 五项 KNN 经验 P 值的区域最小值。

聚合输出为：

```text
results/multiscale/r4000/<comparison>/all_scale_candidates.tsv
results/multiscale/r4000/<comparison>/candidate_union.bed
results/multiscale/r4000/<comparison>/candidate_union.tsv
results/multiscale/r4000/<comparison>/is_change/<comparison>_is_change_loci.tsv
results/multiscale/r4000/<comparison>/is_change/<comparison>_is_change_loci_preliminary.tsv
results/multiscale/r4000/<comparison>/is_change/<comparison>_is_change_masked.bed
results/multiscale/r4000/<comparison>/is_change/figures/manifest.tsv
```

逐窗口表是证据账本，区域表是 IS-change 输入。现有 IS-change 评分、transition/focal-pulse/mixed/broad-shift 分类、局部极值恢复、十个 CWorld IS 距离和三个 IS 尺度带保持不变。primary call 仍要求至少两个 IS 尺度带支持；`support_scale_count` 仅表示窗口尺度支持，不能替代该条件。

## 空结果与失败处理

- 非法窗口算术、旧配置键或 preset 映射错误在 DAG 构建前失败。
- 背景文件缺失或内容不兼容在对应 KNN 前失败。
- 尺度身份重复、KNN/Meltron 坐标不一致或多对多连接直接失败。
- 五个配置尺度缺少任一输出时拒绝聚合，不发布部分结果。
- 单尺度零候选生成带表头的空表；全部尺度零候选时，候选并集和 IS-change 仍生成结构正确的空结果。
- 尺度相关最终文件原子发布，失败的临时文件不能满足 DAG。
- 旧单窗口结果不读取、不覆盖、不删除。

## 验证

最小但充分的验证包括：

1. 窗口和步长算术及向上取整；
2. 非法、重复和旧配置入口拒绝；
3. 五个精确背景 preset 的匹配，以及缺失、重复、空、错表头、错跨度、错网格、重复逻辑行和有效末端窗口测试；
4. 尺度限定路径和完整窗口身份无碰撞；
5. KNN 任一 P 值 `< 0.05` 与 Meltron `>= 5` 的交集保持不变；
6. 五尺度候选拼接、重叠/相邻合并、证据摘要、单尺度支持和全空结果；
7. 小型合成 Snakemake fixture 在一次 DAG 中展开五个尺度；
8. 现有 KNN、Meltron、IS-change、绘图和 TAD 回归测试；
9. 真实配置 dry-run，证明旧单窗口路径不能满足新目标；
10. 一次 dsub 真实运行并持续监控到完整 `rule all` 成功。

## 运行验收

任务只有同时满足以下条件才完成：

- 只进行一次顶层 dsub 提交；
- 同一 Snakemake DAG 包含五个尺度，并向 Donau 提交来自多个尺度的就绪作业；
- 五份背景校验报告全部通过；
- 五套尺度结果全部完成，且没有 incomplete 文件；
- `data/window_compare_pairs.yaml` 中所有 comparison 均生成结构正确的聚合和 IS-change 结果，允许空结果；
- 原有完整 `rule all` 的 TAD 等其他默认目标也全部满足；
- Snakemake 最终正常退出。

## 非目标

- 不在本项目生成或复制 KNN 背景。
- 不修改 QC、Gaussian、Tracy-Widom、NCC、KNN、Meltron 或 IS-change 算法。
- 不要求两个以上窗口尺度支持候选。
- 不删除或迁移旧结果。
- 不增加第二套调度器、五份重复规则或按尺度循环提交脚本。
