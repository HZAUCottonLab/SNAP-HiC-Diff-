# Bug Analysis: `rule all` WildcardError — Python `.format()` 与 Snakemake `{{}}` 转义冲突 (2026-05-09)

## 1. Root Cause Category

- **Category**: B — Cross-Layer Contract (Python str.format() ↔ Snakemake expand())
- **Specific Cause**: Snakemake DAG 解析时，`rule all` 的 `input:` 中的 `expand()` 路径字符串会经过 Python 的 `str.format()` 预处理。`.format()` 先于 Snakemake 的 wildcard 解析执行，因此 `{sample}` / `{pair_id}` / `{chrom}` 在 `expand()` 中被当作 Python 格式化占位符处理，导致 `KeyError`。

## 2. Why Fixes Failed

1. **First attempt**: 改 `{chrom}` → `{{chrom}}` —— 表面正确但遗漏了其他 wildcard，且未理解 `{{}}` 在 `.format()` 字符串中的真实行为
2. **Second attempt**: 修复 `rule all` 第 127 行 `{chrom}` → `{{chrom}}` —— 对了一半，但遗漏了在 `.format()` 字符串中的 `{sample}`、`{pair_id}`
3. **Third attempt**: 改回 `{sample}`、`{pair_id}`（去掉 `{{}}`）—— 导致 Python `KeyError`，因为这些在 `"{}/smooth/{sample}_smooth.cool".format(OUTDIR)` 中被 `.format()` 解析
4. **Final fix**: 区分两种情况：
   - 有 `.format()` 的字符串 → 用 `{{sample}}`（让 `.format()` 产出字面量 `{sample}`）
   - 无 `.format()` 的字符串 → 用 `{sample}`（直接由 `expand()` 填充）

## 3. Prevention Mechanisms

| Priority | Mechanism | Specific Action | Status |
|----------|-----------|-----------------|--------|
| P0 | Documentation | 更新 `snakemake-guidelines.md` Section 9 `expand()` 规则，明确 `.format()` 字符串中的转义行为 | TODO |
| P0 | Documentation | 添加本 bug 分析文档 | DONE |
| P1 | Code Review Checklist | 审查 `rule all` 和 `report` 的 `expand()` 调用时，检查路径字符串是否含 `.format()` 以决定是否需要 `{{}}` | TODO |
| P2 | Lint Rule | 添加 Snakemake lint 检查：检测 `expand()` 中同时使用 `{{}}` 和 `.format()` 的模式 | TODO |

## 4. Systematic Expansion

- **Similar Issues**: 所有在 Snakefile 中使用 `"{path}/{wildcard}".format(...)` + `expand()` 组合的地方，都面临同样的 `.format()` vs Snakemake wildcard 冲突风险。本项目的 `rule all`、`report`、`validate_inputs` 等入口规则都可能受影响。
- **Design Improvement**: 统一使用 `os.path.join()` 而非字符串插值 `.format()` 来构建路径，可以彻底避免此冲突。
- **Process Improvement**: 在 `snakemake-guidelines.md` 的 `rule all` 模式中添加明确的路径构建规则说明。

## 5. Knowledge Capture

### Files Changed

- `Snakefile` line 127, 130, 133, 137, 149: `{{sample}}` → `{sample}` (`.format()` 字符串中), `{chrom}` 保留

### Bug Analysis Document Created

- [x] Created `.trellis/spec/bug-analysis/2026-05-09-expand-format-wildcard-conflict.md` (this file)

### Spec Updates

- [ ] 更新 `snakemake-guidelines.md` Section 9 `expand()` 规则
- [ ] 创建 issue for Snakemake linting rule
- [ ] 创建 issue for path-building convention standardization
