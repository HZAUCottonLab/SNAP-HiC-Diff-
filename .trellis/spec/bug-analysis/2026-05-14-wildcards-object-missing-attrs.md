# Bug: Wildcards 对象缺失未出现在 output 路径中的属性

> 日期: 2026-05-14
> 标签: snakemake, wildcards, params
> 状态: 已解决

## 问题描述

在实现 `rules/window_compare.smk` 时，`compute_window_pair_covariates` rule 需要 `sample_a` 和 `sample_b` wildcards，但这两个值没有出现在 output 路径中（只有 `{comparison}` 和 `{chrom}`）。

## 报错信息

```
AttributeError: 'Wildcards' object has no attribute 'sample_a'
Wildcards:
  comparison=L128_rep1_vs_L128_rep2
  chrom=HC04_A01
```

## 根因分析

Snakemake 的 wildcard 推断机制：
- **只有 output 路径中的 `{wildcard}` 才会被 Snakemake 解析并加入 Wildcards 对象**
- `input` 路径中的 `{wildcard}` 不会被自动解析（除非是 `lambda w:` 函数动态构造）

当在 `input` lambda 函数中访问 `w.sample_a` 时，由于 `sample_a` 不在 output 路径中，Snakemake 不知道这个值，导致 `AttributeError`。

## 解决方案

通过 `params` + lambda function + 外部 YAML 配置推导 `sample_a`/`sample_b`：

```python
def get_pair_samples(comparison):
    """从外部 YAML 配置中查找比较对的样本名称。"""
    pairs_file = config.get("window_compare", {}).get("pairs_file", "data/window_compare_pairs.yaml")
    if os.path.exists(pairs_file):
        with open(pairs_file) as f:
            data = yaml.safe_load(f)
            for p in data.get("comparisons", []):
                if p["name"] == comparison:
                    return p["sample_a"], p["sample_b"]
    return None, None

rule compute_window_pair_covariates:
    params:
        sample_a = lambda w: get_pair_samples(w.comparison)[0],
        sample_b = lambda w: get_pair_samples(w.comparison)[1],
        qc_a = lambda w: "{}/windows/{}/{}__{}_window_qc.tsv".format(
            OUTDIR, get_pair_samples(w.comparison)[0], get_pair_samples(w.comparison)[0], w.chrom),
    shell:
        r"""
        SAMPLE_A="{params.sample_a}"   # ← 从 params 获取，而非 wildcards
        """
```

## 关键点

1. **params lambda 在执行时求值**：当 Snakemake 执行 rule 时，`params` 中的 lambda 函数会被调用，此时 `w.comparison` 已经被解析
2. **外部配置作为桥梁**：通过 YAML 配置文件中存储的 `comparison → sample_a/sample_b` 映射，实现间接访问
3. **output 路径决定 wildcards**：这是 Snakemake 的核心约束，必须在设计 rule 时就考虑

## 相关文件

- `rules/window_compare.smk`
- `data/window_compare_pairs.yaml`
