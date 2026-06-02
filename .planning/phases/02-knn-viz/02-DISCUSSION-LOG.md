# Phase 02 Discussion Log

**Phase:** 02-knn-viz
**Date:** 2026-05-30

## Areas Discussed

### 1. 可视化类型组合
- **Question:** 哪些图是核心必需的？
- **Options presented:** 极简核心（2~3 张）/ 标准分析套件（5~7 张）/ 全面探索（10+ 张）
- **Decision:** 全面探索（10+ 张）— 包含所有候选图类型

### 2. 交互性深度
- **Question:** 静态图、轻微交互还是丰富交互？
- **Options presented:** 纯静态 / 轻度交互（hover tooltip）/ 丰富交互（滑块/联动）
- **Decision:** 纯静态图（matplotlib/seaborn）

### 3. Notebook 结构
- **Question:** 如何组织 Notebook 代码结构？
- **Options presented:** 线性流程 / 模块化（Section 分组）/ 工具类封装
- **Decision:** 工具类封装 — Plotter 类封装所有绘图方法，Notebook 只调用方法

### 4. 多 comparison 支持
- **Question:** 如何处理多 comparison 场景？
- **Options presented:** 硬编码单一 / 参数化 / 自动发现
- **Decision:** 自动发现 — glob 扫描所有 `results/window_compare/` 下的 comparison

### 5. 输出形式
- **Question:** 图的输出方式？
- **Options presented:** 纯 inline / inline + 保存 PNG / inline + SVG + PNG + 路径记录
- **Decision:** inline + SVG + PNG + 路径记录 — 保存到 `results/figures/{comparison}/`

---

*Discussion completed: 2026-05-30*
