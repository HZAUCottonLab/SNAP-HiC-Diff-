# 3D Genomics Thinking Guide

> Concepts for understanding chromatin architecture analysis with Hi-C data

## Table of Contents

1. [Core Concepts](#1-core-concepts)
2. [Analysis Levels](#2-analysis-levels)
3. [Normalization Strategies](#3-normalization-strategies)
4. [Key Algorithms](#4-key-algorithms)
5. [Common Pitfalls](#5-common-pitfalls)
6. [Project-Specific Considerations](#6-project-specific-considerations)

---

## 1. Core Concepts

### 1.1 Chromatin Contact Matrix

```
           A   B   C   D   E   F   G   H
        ┌───────────────────────────────┐
    A   │   •   •   •   •   •   •   •   •  │
    B   │   •   •   •   •   •   •   •   •  │
    C   │   •   •   •   •   •   •   •   •  │
    D   │   •   •   •   •   •   •   •   •  │
    E   │   •   •   •   •   •   •   •   •  │
    F   │   •   •   •   •   •   •   •   •  │
    G   │   •   •   •   •   •   •   •   •  │
    H   │   •   •   •   •   •   •   •   •  │
        └───────────────────────────────┘

• = contact frequency between genomic bins
Matrix is symmetric (A-B contact = B-A contact)
Diagonal = self-contacts (always highest)
```

### 1.2 Distance Decay

**Key insight**: Contact frequency decreases with genomic distance.

```
Contact Frequency ∝ 1 / distance^s

where s ≈ 0.5-1.0 (varies by species/resolution)
```

This is **NOT noise** — it's fundamental biology. Must account for it in normalization.

### 1.3 Resolution Considerations

| Resolution | Bin Size | Best For |
|------------|----------|----------|
| Low | 100kb+ |Compartments, large-scale structure |
| Medium | 10-50kb | TADs, large loops |
| High | 1-5kb | Fine loops, CTCF-CTCF contacts |
| Ultra-high | <1kb | Nucleosome-level, single-molecule |

**Rule**: Resolution should match biological feature size:
- TADs (~100kb-1Mb) → 10-50kb bins
- Loops (~10-100kb) → 5-10kb bins

---

## 2. Analysis Levels

### 2.1 Hierarchical Organization

```
┌─────────────────────────────────────────────────────────────┐
│  Chromosome Territories                                      │
│  ┌─────────┐ ┌─────────┐ ┌─────────┐                       │
│  │  Chr 1  │ │  Chr 2  │ │  Chr 3  │  (inter-chromosomal)   │
│  └─────────┘ └─────────┘ └─────────┘                        │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│  A/B Compartments                                           │
│  ┌──┬────────┬──┬────────┬──┬────────┬──┬────────┬──┐      │
│  │A │   B    │A │   B    │A │   B    │A │   B    │A │     │
│  └──┴────────┴──┴────────┴──┴────────┴──┴────────┴──┘      │
│  Active (A) vs Repressive (B) chromatin                    │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│  TADs (Topologically Associating Domains)                   │
│  ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐       │
│  │  TAD 1   │ │  TAD 2   │ │  TAD 3   │ │  TAD 4   │       │
│  │  ┌────┐  │ │  ┌────┐  │ │  ┌────┐  │ │  ┌────┐  │       │
│  │  │    │  │ │  │    │  │ │  │    │  │ │  │    │  │       │
│  │  └────┘  │ │  └────┘  │ │  └────┘  │ │  └────┘  │       │
│  └──────────┘ └──────────┘ └──────────┘ └──────────┘       │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│  Loops / CTCF Peaks                                        │
│         ●─────────────────●                                │
│         │                 │                                 │
│         │     ○ loop      │                                 │
│         │                 │                                 │
│         ●─────────────────●                                │
└─────────────────────────────────────────────────────────────┘
```

### 2.2 Compartment Analysis (A/B)

**What**: Segregate genome into active (A) and repressive (B) compartments.

**How**:
1. Compute Pearson correlation matrix
2. Eigenvector decomposition
3. PC1 separates A (positive) from B (negative)
4. Flip sign based on gene density or histone marks

**Interpretation**:
- A compartments: gene-dense, H3K36me3, early replication
- B compartments: gene-sparse, H3K9me3, late replication

### 2.3 TAD Analysis

**What**: Identify self-interacting genomic domains.

**Metrics**:
- **Insulation Score**: Low at TAD boundaries
- **TAD Separation Score**: Maximum at boundaries
- **Boundary Strength**: Delta between flanking TADs

**hicFindTADs Algorithm**:
```
1. Compute z-score matrix
2. For each bin, calculate TAD-separation score:
   - Left region: bins i to i-depth
   - Right region: bins i to i+depth
3. Local minima = TAD boundaries
4. Statistical test (Wilcoxon) + delta filter
```

### 2.4 Loop Detection

**What**: Identify specific pairwise contacts (chromatin loops).

**Key insight**: Loops appear as off-diagonal peaks above background.

**hicDetectLoops Algorithm**:
```
1. Compute O/E matrix (remove distance decay)
2. Fit negative binomial per distance
3. Preselect candidates with p < threshold
4. For each candidate neighborhood:
   - Compare peak region vs background
   - Wilcoxon rank-sum test
5. Keep strongest per neighborhood
```

---

## 3. Normalization Strategies

### 3.1 Why Normalize?

```
Raw contacts contain:
├── Coverage bias (GC content, mappability)
├── Distance dependence (power law decay)
├── Technical artifacts (PCR duplicates, bubbles)
└── Biological signal of interest

Goal: Remove biases, preserve biological signal
```

### 3.2 ICE (Iterative Correction and Eigenvector)

```python
def ice_normalize(matrix):
    weights = ones(n_bins)
    
    for iteration in range(max_iter):
        # 1. Compute row sums
        row_sums = (matrix * weights).sum(axis=1)
        
        # 2. Target (mean of non-zero sums)
        target = mean(row_sums[row_sums > 0])
        
        # 3. Update weights
        weights = weights * sqrt(target / row_sums)
        
        # 4. Apply to matrix
        matrix = matrix / outer(weights, weights)
        
        # 5. Check convergence
        if max_change(weights) < tolerance:
            break
    
    return weights
```

**Properties**:
- Assumes symmetric coverage bias
- Assumes uniform interaction probability
- Sensitive to zero/missing bins

### 3.3 Knight-Ruiz (KR) Balancing

```python
def kr_balance(matrix):
    # Fast eigenvector-based balancing
    # 1. Initialize scaling vector
    x = ones(n)
    
    # 2. Iterate
    for i in range(max_iter):
        # Normalize to unit diagonal
        d = sqrt(diag(matrix @ diag(x)))
        y = x / d
        
        # Update
        x = y * sqrt(target / (y.T @ matrix @ y))
        
        if converged(x, y):
            break
    
    return x  # balance weights
```

**Properties**:
- Faster than ICE
- Memory efficient
- Good for high-resolution matrices

### 3.4 O/E (Observed/Expected)

```python
def compute_oe(matrix):
    n = len(matrix)
    expected = zeros(n)
    
    for d in range(n):
        # All bins at same diagonal distance
        diag_vals = diag(matrix, d)
        expected[d] = median(diag_vals[diag_vals > 0])
    
    # O/E for each entry
    oe = matrix / expected[abs(row - col)]
    
    return oe
```

**Properties**:
- Removes distance dependence
- Preserves relative differences
- Good for comparison across conditions

### 3.5 When to Use Which

| Scenario | Recommended |
|----------|-------------|
| General QC + visualization | KR or ICE |
| A/B compartment calling | O/E then Pearson |
| TAD boundary detection | ICE (preserves boundaries) |
| Loop calling | O/E (removes distance bias) |
| Differential analysis | Same normalization for all samples |

---

## 4. Key Algorithms

### 4.1 TAD Separation Score

```
For bin i at depth d:
  
  TAD_sep(i, d) = mean(z-score[left region]) - mean(z-score[right region])
  
Where:
  left region = bins [i-d, i)
  right region = bins (i, i+d]
  
Boundary = local minimum of TAD_sep across multiple depths
```

### 4.2 Insulation Score

```
For bin i at window size w:

  IS(i) = log2( sum(contacts crossing boundary at i) )
  
Low IS = strong boundary
High IS = inside TAD
```

### 4.3 Loop Detection (HiCCUPS-style)

```
For each pixel (i, j) at distance d:

  1. Compute O/E(i, j)
  2. Get p-value from NB distribution at distance d
  3. If p < threshold_preselect:
     - Compare peak region vs background
     - Wilcoxon test
     - If significant: call loop

Peak region: (2*peakWidth)² around (i, j)
Background: neighborhood minus peak region
```

### 4.4 Tracy-Widom for Differential

```
For differential signal Δ:

  1. Normalize Δ to standard normal
  2. Compute TW statistic from largest eigenvalue
  3. p-value from TW distribution
  
Intuition: Large deviations from random noise → significant
```

---

## 5. Common Pitfalls

### 5.1 Resolution Mismatch

**Wrong**: Using 100kb bins to detect 10kb loops
**Right**: Match bin size to feature size

### 5.2 Ignoring Coverage Bias

**Wrong**: Comparing raw counts
**Right**: Normalize with ICE/KR first

### 5.3 Distance Decay Confusion

**Wrong**: Expecting uniform contact frequency
**Right**: Account for 1/distance^s relationship

### 5.4 Zero Bin Problems

**Problem**: ICE amplifies low-coverage bins

**Solutions**:
- Filter low-coverage bins before normalization
- Use KR instead (more robust)
- Set `--filterThreshold` in hicCorrectMatrix

### 5.5 Asymmetric Normalization

**Problem**: Different samples normalized differently

**Solution**: Apply same parameters to all samples

### 5.6 Diagonal Contamination

**Problem**: Self-contacts dominate signal

**Solution**: `--ignoreDiags 2` or more

---

## 6. Project-Specific Considerations

### 6.1 Expected Value Computation

**工具选择**：`extract_expected.smk`（Python/cooler）计算 expected 值，**不使用 hictk**。

原因：`hictk dump --matrix-type expected` 只对 `.hic` 文件有效，`.mcool` 文件无 expected 数据集。

```python
# extract_expected.smk 实现：分块流式读取，避免 OOM
import cooler
import numpy as np

uri = "sample_iced.mcool::/resolutions/4000"
clr = cooler.Cooler(uri)

chunk_size = 5000
max_dist = 10000
dist_sums = np.zeros(max_dist)
dist_counts = np.zeros(max_dist)

n_bins = clr.bins().fetch(chrom).shape[0]
for start in range(0, n_bins, chunk_size):
    end = min(start + chunk_size, n_bins)
    mat = clr.matrix(balance=True, sparse=True).fetch(chrom)[start:end, :].tocoo()
    for i, j, val in zip(mat.row, mat.col, mat.data):
        row_idx = start + i
        if row_idx < j:  # 上三角，避免重复计算
            d = j - row_idx
            if d < max_dist and np.isfinite(val) and val != 0:
                dist_sums[d] += val
                dist_counts[d] += 1

expected = np.zeros(max_dist)
for d in range(max_dist):
    if dist_counts[d] > 0:
        expected[d] = dist_sums[d] / dist_counts[d]
```

### 6.2 hictk dump 解析注意事项

不同 `--table` 输出格式不同（有/无表头）：

| `--table` | 有表头 | 解析方式 |
|-----------|--------|---------|
| `resolutions` | 无 | `lines[:]` |
| `chroms` | 无 | `lines[:]` |
| `bins` | 有 | `lines[1:]` |
| `pixels` | 有 | `lines[1:]` |

**常见错误**：

```python
# ❌ 对 resolutions/chroms 用 lines[1:]，会丢失第一个值
for line in lines[1:]:  # 4000 被当成 header 跳过

# ✅ 使用 lines[:]
for line in lines:
    resolutions.append(int(line.split()[0]))
```

### 6.3 SNAP-HiC-Diff Pipeline Comparison

| Feature | HiCExplorer | SNAP-HiC-Diff |
|---------|-------------|---------------|
| Input | BAM → matrix | .mcool (pre-processed) |
| Normalization | ICE/KR | ICE + O/E + z-score |
| Differential | TAD-based (Wilcoxon) | Window-based (Tracy-Widom) |
| Gap Handling | Mask bins | Coverage gradient detection |

### 6.2 Key Design Decisions

**Why ICE + O/E + z-score?**
- ICE: Standard Hi-C normalization
- O/E: Removes distance dependence
- Z-score: Enables cross-condition comparison

**Why Tracy-Widom?**
- Robust to non-normal distributions
- Computationally efficient (GPU)
- Well-suited for matrix comparisons

### 6.3 Gap Detection Strategy

```
Coverage gradient approach:
1. Bin coverage by balanced matrix marginals
2. Sort by coverage (descending)
3. Find breakpoint with largest gradient
4. Classify low-coverage bins as gaps

This is different from:
- Assembly gaps (N regions)
- Blacklisted regions (ENCODE)
```

### 6.4 Window-Based vs TAD-Based

| Aspect | Window-Based (SNAP) | TAD-Based (HiCExplorer) |
|--------|---------------------|-------------------------|
| Granularity | Uniform windows | Variable TAD sizes |
| Sensitivity | Good for focal changes | Good for domain changes |
| Resolution | Fixed by window size | Flexible |
| Statistics | TW for global significance | Wilcoxon per TAD region |

---

## Quick Reference: Decision Tree

```
                    ┌─────────────────────┐
                    │ Start with .mcool   │
                    └──────────┬──────────┘
                               │
                    ┌──────────▼──────────┐
                    │ Normalize (ICE/KR)  │
                    └──────────┬──────────┘
                               │
              ┌────────────────┼────────────────┐
              │                │                │
       ┌──────▼──────┐  ┌──────▼──────┐  ┌──────▼──────┐
       │ Compartments │  │    TADs    │  │    Loops    │
       │ (PCA/Eigen)  │  │hicFindTADs │  │hicDetectLoops│
       └──────────────┘  └─────────────┘  └─────────────┘
              │                │                │
       ┌──────▼──────┐  ┌──────▼──────┐  ┌──────▼──────┐
       │ Compare A/B │  │ Boundary    │  │ Anchor      │
       │ between     │  │ strength    │  │ motif       │
       │ conditions  │  │ changes    │  │ analysis    │
       └──────────────┘  └─────────────┘  └─────────────┘
```

---

## Further Reading

- Lieberman-Aiden et al. (2009). "Comprehensive mapping of long-range interactions reveals folding principles of the human genome." Science.
- Dixon et al. (2012). "Topological domains in mammalian genomes identified by analysis of chromatin interactions." Nature.
- Rao et al. (2014). "A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping." Cell.
- Imakaev et al. (2012). "Iterative correction of Hi-C data reveals hallmarks of chromosome organization." Nature Methods.
