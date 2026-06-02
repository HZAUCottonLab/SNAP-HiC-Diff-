# Bug Analysis: hictk Table Parsing + Expected Value Computation (May 2026)

## 1. Root Cause Category

- **Category**: B — Cross-Layer Contract + C — Change Propagation Failure

两个独立 bug 属于不同类别：

| Bug | Category | 具体原因 |
|-----|----------|---------|
| `lines[1:]` 跳过表头 | B — Cross-Layer Contract | hictk CLI 输出格式（无表头）与代码假设（有表头）不一致 |
| `--matrix-type expected` 对 mcool 无效 | C — Change Propagation Failure | 代码从 hic 文档复制到 mcool 场景，`.hic` 和 `.mcool` 的能力差异未更新 |

---

## 2. Bug 1: hictk table 输出误判表头

### 2.1 Symptom

测试报告 `L128_rep1.mcool` 缺少 4000 分辨率，实际文件包含 `[4000, 8000, 16000, 32000, 64000, 128000, 256000]`。

### 2.2 Root Cause

`hictk dump --table resolutions` 和 `hictk dump --table chroms` **均无表头**，但代码用 `lines[1:]` 跳过了第一行：

```python
# scripts/validate_inputs.py
for line in lines[1:]:  # BUG: 跳过了第一个分辨率 4000
    parts = line.split()
    if len(parts) >= 1:
        resolutions.append(int(parts[0]))
```

```bash
$ hictk dump L128_rep1.mcool --table resolutions
4000    # ← 被当成 header 跳掉了
8000
16000
...

$ hictk dump L128_rep1.mcool::/resolutions/32000 --table chroms
HC04_A01 119507322   # ← 同样被跳过
```

### 2.3 Why Fixes Failed

1. **表面修复**：用户最初认为分辨率确实不存在，让我验证文件
2. **直接验证**：用 `hictk dump --table resolutions` 直接确认了所有分辨率都存在
3. **Python 重现**：用 Python subprocess 重现了相同的解析逻辑，精确定位到 `lines[1:]` 是罪魁祸首
4. **最终修复**：改为 `lines[:]`（或直接 `lines`），不对无表头 table 使用 header-skip 逻辑

### 2.4 Prevention Mechanisms

| Priority | Mechanism | Specific Action | Status |
|----------|-----------|-----------------|--------|
| P0 | Documentation | 在 `hictk.md` 的 table 输出格式中标注每个 table 是否有表头 | DONE |
| P1 | Code review | 审查所有 hictk dump 解析代码，确认是否正确处理表头 | TODO |
| P2 | Test coverage | 单元测试验证 `hictk dump` 解析结果与 `hictk dump --table resolutions` 输出完全匹配 | TODO |

---

## 3. Bug 2: `--matrix-type expected` 不支持 mcool

### 3.1 Symptom

测试脚本运行 `hictk dump L128_rep1.mcool --matrix-type expected` 失败（返回非零退出码）。

### 3.2 Root Cause

`--matrix-type expected` 仅对 `.hic` 文件有意义，hictk 文档中针对 `.hic` 的示例被复制到了 mcool 场景：

```bash
# hictk.md 中示例：
# Expected values
hictk dump sample.hic --resolution 1000 --matrix-type expected --join
```

`.mcool` 文件没有 `expected` 数据集，`--matrix-type expected` 对 mcool 无效。

### 3.3 Why Fixes Failed

用户主动确认不需要 hictk 计算 expected——`extract_expected.smk` 已用 Python/cooler 正确实现。测试脚本应删除该 hictk 命令。

### 3.4 Prevention Mechanisms

| Priority | Mechanism | Specific Action | Status |
|----------|-----------|-----------------|--------|
| P0 | Documentation | 在 `hictk.md` 明确标注 `--matrix-type expected` 适用范围 | DONE |
| P1 | Architecture | expected value 计算统一使用 `extract_expected.smk`（Python/cooler），不在测试中混用 hictk | DONE |
| P2 | Test coverage | 测试脚本的 expected 验证使用 cooler 分块流式读取，不调用 hictk | DONE |

---

## 4. Systematic Expansion

### 4.1 Similar Issues

- **Bug A 的同类问题**：其他所有 hictk `--table` 参数的输出格式（有无表头）需要逐一验证

```bash
# 已有确认：
--table resolutions     # 无表头 ✓
--table chroms          # 无表头 ✓
--table bins            # 有表头（bin_id, chrom, start, end, weight）
--table pixels          # 有表头（bin1_id, bin2_id, count）
```

- **Bug B 的同类问题**：hictk 文档中的示例混合了 `.hic` 和 `.cool` 命令，应在每个示例前标注文件类型

### 4.2 Design Improvement

- **统一的 hictk 解析工具函数**：所有解析 `hictk dump` 输出的代码应使用同一个 helper function，统一处理表头判断

```python
def hictk_dump_table(uri, table_name):
    """Parse hictk dump table output.
    Tables WITHOUT header: resolutions, chroms
    Tables WITH header: bins, pixels, normalizations, cells
    """
    tables_without_header = {"resolutions", "chroms"}
    cmd = ["hictk", "dump", uri, "--table", table_name]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    lines = result.stdout.strip().split("\n")
    if table_name not in tables_without_header:
        lines = lines[1:]  # skip header
    return [line.split() for line in lines if line.strip()]
```

### 4.3 Process Improvement

- **dsub 命令语法精确性**：dsub 命令格式为 `dsub ... "command arg1 arg2"`（无 `--` 分隔符），测试脚本注释中记录的语法需与实际一致
- **当前节点 vs 计算节点**：需要区分可在当前节点执行的轻量测试（如语法检查）和必须通过 dsub 提交的重量测试（如 expected value 计算）

---

## 5. Knowledge Capture

- [x] **更新 `hictk.md` 表头格式说明**：在 `--table` 参数表格中标注每个 table 是否有表头
- [x] **更新 `3d-genomics-guide.md` expected value 计算方式**：明确使用 Python/cooler，不使用 hictk
- [x] **创建 `scripts/test_pipeline.py`**：覆盖所有 pipeline 阶段的集成测试，dsub 提交方式正确
- [x] **更新 `config/config.yaml`**：默认分辨率改为 4kb
- [ ] **为 `validate_inputs.py` 添加单元测试**：用 mock 验证 `lines[1:]` vs `lines[:]` 的差异
- [ ] **CI/Audit 检查**：添加 `grep -n 'lines\[1:\]' scripts/validate_inputs.py` 检查是否有其他类似 bug

---

## 6. Additional Fixes (Same Session)

本次对话中同时修复了以下问题：

### 6.1 Bug C: config import 在循环内

```python
# BEFORE (scripts/validate_inputs.py)
for idx, row in manifest.iterrows():
    from config import config  # 在循环内 import，scripts/ 目录找不到
    target_res = config.get("resolution", 32000)

# AFTER
import yaml
config_path = os.path.join(os.path.dirname(__file__), "..", "config", "config.yaml")
with open(config_path) as f:
    cfg = yaml.safe_load(f)
target_res = cfg.get("resolution", 32000)  # 循环外加载
```

### 6.2 Bug D: 测试脚本 OOM

4kb 分辨率下 119Mb 染色体 = 29877 bins，一次性 `clr.matrix().fetch()` 会展开 29877² 的稀疏矩阵导致 OOM。

```python
# BEFORE
mat = clr.matrix(balance=True, sparse=True).fetch(chrom).tocsr()  # OOM

# AFTER: 分块流式读取
for start in range(0, n_bins, chunk_size):
    end = min(start + chunk_size, n_bins)
    mat_chunk = clr.matrix(balance=True, sparse=True).fetch(chrom)[start:end, :].tocoo()
    # 累加 dist_sums, dist_counts
```

### 6.3 Bug E: dsub 命令语法

```bash
# BEFORE
dsub ... -- python3 scripts/test_pipeline.py

# AFTER (正确语法)
dsub ... "python3 scripts/test_pipeline.py"
```

---

## 7. Summary Table

| Bug | Category | File | Fix |
|-----|----------|------|-----|
| `lines[1:]` 误跳首行 | B — Cross-Layer Contract | `scripts/validate_inputs.py` | `lines[1:]` → `lines[:]` |
| `--matrix-type expected` 对 mcool 无效 | C — Change Propagation Failure | 测试脚本设计 | 删除该命令，用 Python/cooler |
| config import 在循环内 | D — Test Coverage Gap | `scripts/validate_inputs.py` | 循环外 yaml.load |
| 4kb 矩阵 OOM | E — Implicit Assumption | `scripts/test_pipeline.py` | 分块流式读取 |
| dsub 命令语法错误 | E — Implicit Assumption | `scripts/test_pipeline.py` | `--` → `"..."` |
