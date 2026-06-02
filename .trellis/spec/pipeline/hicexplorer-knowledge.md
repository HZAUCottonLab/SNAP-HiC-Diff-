# HiCExplorer Knowledge Spec

> Distilled from [HiCExplorer documentation](https://hicexplorer.readthedocs.io/en/latest/) | deepTools team, MPI-IE & University of Freiburg

## Table of Contents

1. [Overview](#1-overview)
2. [File Formats](#2-file-formats)
3. [Preprocessing Tools](#3-preprocessing-tools)
4. [Normalization Methods](#4-normalization-methods)
5. [TAD Detection](#5-tad-detection)
6. [Loop Detection](#6-loop-detection)
7. [Differential Analysis](#7-differential-analysis)
8. [Visualization](#8-visualization)
9. [Quality Control](#9-quality-control)
10. [Comparison with Project Pipeline](#10-comparison-with-project-pipeline)

---

## 1. Overview

**HiCExplorer** is a comprehensive tool suite for processing, normalizing, analyzing, and visualizing Hi-C and cHi-C data.

### Citation

```
Wolff et al. (2020). "Galaxy HiCExplorer 3: a web server for reproducible Hi-C,
capture Hi-C and single-cell Hi-C data analysis, quality control and visualization."
Nucleic Acids Research, 48(W1), W177-W184.
```

### Installation

```bash
# Via conda (recommended)
conda install -c bioconda hicexplorer

# Or via pip
pip install hicexplorer
```

### Tool Categories

| Category | Tools |
|----------|-------|
| **Preprocessing** | `hicBuildMatrix`, `hicFindRestSite` |
| **QC** | `hicQuickQC`, `hicQC` |
| **Normalization** | `hicCorrectMatrix`, `hicNormalize`, `hicTransform` |
| **TAD Analysis** | `hicFindTADs`, `hicMergeDomains`, `hicMergeTADbins` |
| **Loop Detection** | `hicDetectLoops`, `hicMergeLoops`, `hicValidateLocations` |
| **Differential** | `hicDifferentialTAD`, `hicCompareMatrices`, `hicCorrelate` |
| **Visualization** | `hicPlotMatrix`, `hicPlotTADs`, `hicPlotDistVsCounts`, `hicAggregateContacts` |
| **PCA/Compartments** | `hicPCA`, `hicCompartmentalization` |

---

## 2. File Formats

### 2.1 Cooler Format (.cool / .mcool)

```bash
# List available resolutions in .mcool
cooler list sample.mcool
```

Output files use HDF5-based cooler format:
- `.cool` — single resolution
- `.mcool` — multi-resolution (multiple bin sizes)
- `.scool` — single-cell Hi-C

### 2.2 hic Format (.h5)

Legacy format from HiCExplorer. Can be converted:

```bash
hicConvertFormat -m input.cool --outFileName output.h5 --inputFormat cool --outputFormat h5
```

### 2.3 BED/Interval Formats

| File | Format | Description |
|------|--------|-------------|
| `*_domains.bed` | BED6 | TAD regions |
| `*_boundaries.bed` | BED | TAD boundaries |
| `*_boundaries.gff` | GFF3 | Boundaries with p-values |
| `*_score.bedgraph` | bedGraph | TAD separation scores |
| `loops.bedgraph` | Custom | Loop coordinates (6+ columns) |

---

## 3. Preprocessing Tools

### 3.1 hicFindRestSite

Find restriction enzyme cut sites in the genome.

```bash
hicFindRestSite \
    --genomeAssembly dm6 \
    --searchPattern GATC \
    --outputFileName dm6_DpnII.bed
```

**Parameters:**
- `--genomeAssembly` — Genome name (for UCSC download) or FASTA file
- `--searchPattern` — Restriction enzyme motif (e.g., `GATC` for DpnII, `AAGCTT` for HindIII)
- `--outputFileName` — Output BED file

### 3.2 hicBuildMatrix

Build Hi-C matrix from aligned BAM files.

```bash
hicBuildMatrix \
    --samFiles forward.bam reverse.bam \
    --outFileName hic_matrix.cool \
    --QCfolder QC \
    --restrictionCutFile dm6_DpnII.bed \
    --restrictionSequence GATC \
    --danglingSequence GATC \
    --binSize 10000 \
    --threads 16
```

**Key Parameters:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--binSize` | Bin size in bp. Multiple values create .mcool | - |
| `--minDistance` | Min distance between restriction sites to merge | 300 bp |
| `--maxLibraryInsertSize` | Max insert size for valid pairs | 1000 bp |
| `--minMappingQuality` | Minimum MAPQ | 15 |
| `--region` | Chromosome region (chr:start-end) | all |
| `--doTestRun` | Test with first 1M reads | False |

**Multi-resolution output:**
```bash
hicBuildMatrix \
    --binSize 10000 20000 50000 100000 \
    --outFileName multi_res.mcool
```

---

## 4. Normalization Methods

### 4.1 hicCorrectMatrix

Matrix balancing via Knight-Ruiz (KR) or Iterative Correction (ICE).

#### Step 1: Diagnostic Plot

```bash
hicCorrectMatrix diagnostic_plot \
    --matrix hic_matrix.cool \
    --outFileName plot.png
```

#### Step 2: Correct Matrix

**KR (recommended for high resolution):**
```bash
hicCorrectMatrix correct \
    --matrix hic_matrix.cool \
    --correctionMethod KR \
    --outFileName corrected_KR.cool
```

**ICE:**
```bash
hicCorrectMatrix correct \
    --matrix hic_matrix.cool \
    --correctionMethod ICE \
    --filterThreshold -1.5 5.0 \
    --iterNum 500 \
    --outFileName corrected_ICE.cool
```

**Parameters:**

| Parameter | Description | ICE only |
|-----------|-------------|----------|
| `--filterThreshold` | Remove low/high coverage bins (z-score) | Yes |
| `--iterNum` | Max iterations | Yes |
| `--inflationCutoff` | Max bin scaling factor | Yes |
| `--transCutoff` | Clip top trans contacts (0.05) | Yes |
| `--skipDiagonal` | Exclude diagonal from correction | Yes |

### 4.2 hicTransform

Convert matrix to Obs/Exp, Pearson, or covariance.

```bash
# Observed/Expected
hicTransform \
    --matrix hic_matrix.cool \
    --method obs_exp \
    --outFileName obs_exp.cool

# Pearson correlation
hicTransform \
    --matrix hic_matrix.cool \
    --method pearson \
    --outFileName pearson.cool

# Lieberman method (A/B compartments)
hicTransform \
    --matrix hic_matrix.cool \
    --method obs_exp_lieberman \
    --outFileName lieberman.cool
```

**Methods:**

| Method | Formula | Use Case |
|--------|---------|----------|
| `obs_exp` | Sum(diagonal) / count | General normalization |
| `obs_exp_lieberman` | Sum / (chr_len - distance) | A/B compartments |
| `obs_exp_non_zero` | Sum / non-zero count | Sparse matrices |
| `pearson` | C_ij / sqrt(C_ii × C_jj) | Correlation analysis |
| `covariance` | E[M_i×M_j] - μ_i×μ_j | Covariance analysis |

---

## 5. TAD Detection

### 5.1 hicFindTADs

Detect Topologically Associating Domains using TAD-separation score.

```bash
hicFindTADs \
    --matrix corrected.cool \
    --outPrefix TADs_min3k_max100k \
    --minDepth 30000 \
    --maxDepth 100000 \
    --step 5000 \
    --thresholdComparisons 0.02 \
    --delta 0.01 \
    --correctForMultipleTesting fdr \
    --numberOfProcessors 32
```

**Parameters:**

| Parameter | Description | Recommendation |
|-----------|-------------|----------------|
| `--minDepth` | Min window size (3× bin size) | 3× bin size |
| `--maxDepth` | Max window size (6-10× bin size) | 6-10× bin size |
| `--step` | Step between window sizes | 5000-15000 |
| `--thresholdComparisons` | p-value or q-value threshold | 0.01-0.05 |
| `--delta` | Min delta for boundary separation | 0.01-0.05 |
| `--correctForMultipleTesting` | bonferroni / fdr / None | fdr |
| `--TAD_sep_score_prefix` | Reuse precomputed z-score matrix | For parameter testing |

**Output files:**
- `TADs_domains.bed` — TAD regions
- `TADs_boundaries.bed` — TAD boundaries
- `TADs_score.bedgraph` — TAD separation score track
- `TADs_zscore_matrix.h5` — Z-score matrix (reusable)

**Fast parameter testing:**
```bash
# Reuse existing z-score matrix to test different thresholds
hicFindTADs \
    --matrix corrected.cool \
    --outPrefix new_threshold \
    --TAD_sep_score_prefix old_TADs_zscore_matrix \
    --thresholdComparisons 0.01 \
    --delta 0.02
```

### 5.2 hicMergeDomains

Merge TADs from different resolutions.

```bash
hicMergeDomains \
    --domainFiles 10kb_domains.bed 50kb_domains.bed \
    --proteinFile ctcf_peaks.broadPeak \
    --minimumNumberOfPeaks 2 \
    --value 5000 \
    --percent 0.5 \
    --outputMergedList merged.bed \
    --outputRelationList relations.txt
```

---

## 6. Loop Detection

### 6.1 hicDetectLoops

Detect chromatin loops using negative binomial + Wilcoxon test.

```bash
hicDetectLoops \
    --matrix corrected.cool \
    --outFileName loops.bedgraph \
    --windowSize 10 \
    --peakWidth 6 \
    --pValuePreselection 0.05 \
    --pValue 0.05 \
    --maxLoopDistance 2000000 \
    --threads 16
```

**Parameters:**

| Parameter | Description | Recommendation |
|-----------|-------------|----------------|
| `--windowSize` | Neighborhood size in bins | 10 for 10kb |
| `--peakWidth` | Peak region size (≤ windowSize) | 6 for 10kb |
| `--pValuePreselection` | NB distribution p-value filter | 0.05-0.1 |
| `--peakInteractionsThreshold` | Min peak contacts | 10 |
| `--obsExpThreshold` | Min O/E for peak | 1.5 |
| `--maxLoopDistance` | Max genomic distance | 2MB |
| `--expected` | Expected computation method | mean_nonzero |

**Hyperparameter optimization:**
```bash
# Find optimal parameters
hicHyperoptDetectLoops \
    --matrix corrected.cool \
    --outFileName optimal_loops \
    --maxLoopDistance 2000000
```

---

## 7. Differential Analysis

### 7.1 hicDifferentialTAD

Compare TADs between two conditions.

```bash
hicDifferentialTAD \
    --targetMatrix treatment.cool \
    --controlMatrix control.cool \
    --tadDomains treatment_domains.bed \
    --outFileNamePrefix differential \
    --pValue 0.01 \
    --mode all \
    --modeReject one \
    --threads 8
```

**Parameters:**

| Parameter | Description | Options |
|-----------|-------------|---------|
| `--mode` | Regions to compare | intra-TAD, left-inter-TAD, right-inter-TAD, all |
| `--modeReject` | Rejection criterion | all (all must fail), one (any fails) |

**Output:**
- `differential_accepted.diff_tad` — TADs with similar patterns
- `differential_rejected.diff_tad` — Differential TADs

### 7.2 hicCompareMatrices

Direct matrix comparison (difference, ratio, log2ratio).

```bash
hicCompareMatrices \
    --matrices treatment.cool control.cool \
    --operation log2ratio \
    --outFileName comparison.cool
```

**Operations:**
- `diff` — Simple difference
- `ratio` — Division
- `log2ratio` — log2(treatment/control)

### 7.3 hicCorrelate

Correlation analysis between multiple samples.

```bash
hicCorrelate \
    --matrices wt1.cool wt2.cool kd1.cool kd2.cool \
    --method pearson \
    --log1p \
    --labels WT1 WT2 KD1 KD2 \
    --range 5000:200000 \
    --outFileNameHeatmap heatmap.png \
    --outFileNameScatter scatter.png \
    --plotNumbers
```

---

## 8. Visualization

### 8.1 hicPlotMatrix

Plot Hi-C matrix as heatmap.

```bash
hicPlotMatrix \
    --matrix corrected.cool \
    --outFileName heatmap.png \
    --title "Hi-C Contact Map" \
    --log1p \
    --region chrX:5000000-10000000 \
    --colorMap RdYlBu_r \
    --vMin 1 \
    --vMax 80 \
    --loops loops.bedgraph \
    --tads domains.bed
```

**Parameters:**

| Parameter | Description |
|-----------|-------------|
| `--region` | Genomic region (chr:start-end) |
| `--region2` | Asymmetric region for off-diagonal |
| `--log1p` | Log transform |
| `--clearMaskedBins` | Remove masked bins |
| `--perChromosome` | Plot each chromosome separately |
| `--bigwig` | Add track below matrix |
| `--loops` | Overlay loop calls |
| `--tads` | Overlay TAD boundaries |

### 8.2 hicPlotTADs

Plot TAD tracks with multiple data tracks.

```bash
hicPlotTADs \
    --tracks tracks.ini \
    --region chrX:5000000-8000000 \
    --outFileName tad_tracks.png
```

**tracks.ini example:**
```ini
[x-axis]
fontsize=10

[hic]
file = hic_matrix.cool
title = Hi-C
colormap = Spectral_r
depth = 400000
transform = log1p
file_type = hic_matrix

[tads]
file = TADs_domains.bed
file_type = domains
border_color = black
color = none

[spacer]
height = 0.5

[bigwig]
file = chipseq.bw
title = H3K27ac
color = red
height = 2
file_type = bigwig
```

### 8.3 hicAggregateContacts

Aggregate contacts around specific regions.

```bash
hicAggregateContacts \
    --matrix corrected.cool \
    --BED HiCpeaks.bed \
    --outFileName aggregated.png \
    --plotType pixels
```

---

## 9. Quality Control

### 9.1 hicQuickQC

Quick QC from BAM files (first 1M reads).

```bash
hicQuickQC \
    --samFiles forward.bam reverse.bam \
    --QCfolder QC_test \
    --restrictionCutFile DpnII.bed \
    --restrictionSequence GATC \
    --danglingSequence GATC
```

### 9.2 hicQC

Generate QC report from multiple samples.

```bash
hicQC \
    --logfiles sample1/QC.log sample2/QC.log \
    --labels "Sample 1" "Sample 2" \
    --outputFolder QC_report
```

**QC Metrics:**
- Sequenced reads
- Mappable reads
- Unique reads
- Hi-C contacts
- Dangling ends (re-ligation issues)
- Duplicated pairs (PCR artifacts)
- Same fragment (invalid ligation)
- Self-circles (size < 25kb)
- Inter- vs intra-chromosomal ratio

---

## 10. Comparison with Project Pipeline

### 10.1 Functional Equivalents

| HiCExplorer | SNAP-HiC-Diff | Notes |
|-------------|---------------|-------|
| `hicBuildMatrix` | Input .mcool | Already processed |
| `hicCorrectMatrix` KR/ICE | `extract_chr_mcool.py` ICE | Same balancing |
| `hicTransform --method obs_exp` | `normalize_window.py` O/E | Distance normalization |
| `hicPCA` | - | A/B compartment analysis |
| `hicFindTADs` | - | TAD boundary detection |
| `hicDetectLoops` | - | Loop detection |
| `hicDifferentialTAD` | `compute_diff_metrics.py` + `tracy_widom.py` | Differential analysis |

### 10.2 Key Differences

| Aspect | HiCExplorer | SNAP-HiC-Diff |
|--------|-------------|---------------|
| **Input** | BAM files | Pre-processed .mcool |
| **Resolution** | User-defined | Fixed from .mcool |
| **Differential** | TAD-based | Window-based with TW |
| **Statistics** | Wilcoxon rank-sum | Tracy-Widom |
| **Output** | BED domains | Differential windows |
| **Visualization** | Built-in plots | External tools |

### 10.3 Integration Opportunities

1. **Preprocessing**: Use `hicBuildMatrix` if starting from BAM
2. **QC**: Use `hicQC` for sample quality comparison
3. **TAD Detection**: Use `hicFindTADs` for boundary analysis
4. **Visualization**: Use `hicPlotMatrix` for contact maps
5. **Compartments**: Use `hicPCA` for A/B compartment tracking

---

## Quick Reference: Common Workflows

### Full Preprocessing Pipeline

```bash
# 1. Find restriction sites
hicFindRestSite -g genome.fa -r GATC -o rest_sites.bed

# 2. Build matrix
hicBuildMatrix -s R1.bam R2.bam -o matrix.cool \
    --binSize 10000 --QCfolder QC

# 3. QC check
hicQC -l sample1/QC.log sample2/QC.log -o report/

# 4. Diagnostic
hicCorrectMatrix diagnostic_plot -m matrix.cool -o diag.png

# 5. Correct
hicCorrectMatrix correct -m matrix.cool -o corrected.cool \
    --correctionMethod KR
```

### TAD Analysis Pipeline

```bash
# 1. Find TADs
hicFindTADs -m corrected.cool -o TADs \
    --minDepth 30000 --maxDepth 100000 --step 5000

# 2. Visualize
hicPlotTADs --tracks tracks.ini -r chrX:1M-5M -o tad_view.png

# 3. Differential
hicDifferentialTAD -tm treated.cool -cm control.cool \
    -td TADs_domains.bed -o diff_tads
```

---

## References

- HiCExplorer GitHub: https://github.com/deeptools/HiCExplorer
- Documentation: https://hicexplorer.readthedocs.io/
- Galaxy Instance: http://hicexplorer.usegalaxy.eu
- Citation: Wolff et al. (2020), NAR, 48(W1), W177-W184
