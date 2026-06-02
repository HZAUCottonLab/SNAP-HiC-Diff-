# hictk CLI Specification

Blazing fast toolkit to work with .hic and .cool files.

**Documentation**: https://hictk.readthedocs.io/en/stable/

## Overview

hictk is a C++ tool for Hi-C data analysis supporting:
- `.cool` / `.mcool` / `.scool` - Cooler format
- `.hic` - HiC format (including hic v9)

## Subcommands

| Command | Purpose |
|---------|---------|
| `balance` | Balance Hi-C matrices (ICE, SCALE, VC) |
| `convert` | Convert between Hi-C formats |
| `dump` | Read interactions and data to stdout |
| `fix-mcool` | Fix corrupted .mcool files |
| `load` | Build .cool/.hic from text formats |
| `merge` | Merge multiple files |
| `metadata` | Print file metadata |
| `rename-chromosomes` | Rename chromosomes in Cooler files |
| `validate` | Validate file integrity |
| `zoomify` | Create multi-resolution files |

---

## hictk validate

**Purpose**: Verify file integrity, detect corruption.

### Signatures

```bash
hictk validate [OPTIONS] uri
hictk validate FILE [--validate-index] [--validate-pixels]
```

### Key Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--validate-index` | flag | false | Check Cooler index (slow) |
| `--validate-pixels` | flag | false | Validate pixel data (slow) |
| `--output-format` | enum | json | Output format: json/toml/yaml |
| `--exhaustive` | flag | false | Don't stop at first error |
| `--quiet` | flag | false | Suppress stdout output |

### Common Issues Detected

1. **Cooler index corruption** - Older cooler versions had bugs causing duplicate pixels
2. **File truncation** - Failed downloads
3. **Unsorted pixels** - Data integrity issues
4. **Invalid counts** - Zero or negative values

### Examples

```bash
# Quick validation
hictk validate sample.cool

# Check index (detect duplicate pixels)
hictk validate sample.mcool::/resolutions/1000 --validate-index

# Full validation with pixel check
hictk validate sample.cool --validate-pixels --exhaustive

# JSON output
hictk validate sample.hic --output-format json
```

### Validation Matrix

| File Type | Quick Check | --validate-index | --validate-pixels |
|-----------|-------------|------------------|-------------------|
| .hic | ✓ | N/A | N/A |
| .cool | ✓ | ✓ | ✓ |
| .mcool | ✓ (per resolution) | ✓ | ✓ |
| .scool | ✓ | ✓ | ✓ |

---

## hictk dump

**Purpose**: Read interactions and data tables.

### Signatures

```bash
hictk dump [OPTIONS] uri
```

### Key Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--resolution` | uint | auto | Matrix resolution (for .hic) |
| `--matrix-type` | enum | observed | Type: observed/oe/expected | **只对 `.hic` 有效**，`.cool/.mcool/.scool` 无 expected 数据集 |
| `--matrix-unit` | enum | BP | Unit: BP/FRAG |
| `-t, --table` | enum | pixels | Table: chroms/bins/pixels/normalizations/resolutions/cells/weights |
| `-r, --range` | text | all | Genomic region (chr:start-end) |
| `--range2` | text | - | Second region for trans queries |
| `--query-file` | file | - | BEDPE file for batch queries |
| `--cis-only` | flag | false | Intra-chromosomal only |
| `--trans-only` | flag | false | Inter-chromosomal only |
| `-b, --balance` | text | NONE | Normalization method |
| `--join` | flag | false | Output BG2 format |
| `--sorted` | flag | false | Sort output |

### Examples

```bash
# Basic dump (COO format: bin1 bin2 count)
hictk dump sample.cool

# BedGraph2 format with coordinates
hictk dump sample.cool --join

# Dump specific region
hictk dump sample.hic --resolution 1000 --range chr2L:0-100000 --join

# Dump with chromosome names
hictk dump sample.mcool --resolution 1000 --join | head

# Balanced interactions
hictk dump sample.hic --resolution 1000 --balance SCALE --join

# Expected values (只对 .hic 文件有效)
hictk dump sample.hic --resolution 1000 --matrix-type expected --join

# Trans interactions only
hictk dump sample.hic --resolution 1000 --trans-only --join

# Dump chromosome list
hictk dump sample.hic --table chroms

# Dump available normalizations
hictk dump sample.hic --table normalizations

# Batch query from BEDPE file
hictk dump sample.cool --query-file regions.bedpe
```

### Output Formats

**COO format** (default):
```
bin1 bin2 count
7 7 1745
7 12 2844
```

**BG2 format** (`--join`):
```
chr1 start1 end1 chr2 start2 end2 count
chr2L 7000 8000 chr2L 7000 8000 1745
```

---

## hictk dump `--table` 输出格式：表头的有无

> ⚠️ **重要**：不同 `--table` 选项的输出格式（有无表头）不同。解析代码必须区分处理，否则会丢失第一行数据。

| `--table` | 是否有表头 | 输出格式 | 示例 |
|-----------|-----------|---------|------|
| `resolutions` | **无** | 每行一个整数（分辨率值） | `4000\n8000\n16000` |
| `chroms` | **无** | `name\tsize`（两列，无 header 行） | `HC04_A01\t119507322` |
| `bins` | 有 | `bin_id\tchrom\tstart\tend\t...` | `0\tHC04_A01\t0\t4000\t1.23` |
| `pixels` | 有 | `bin1_id\tbin2_id\tcount` | `7\t7\t1745` |
| `normalizations` | 有 | `name\tvalue` | `ICE\t1.0` |
| `cells` | 有 | `cell_id\t...` | - |

**解析示例（Python）**：

```python
TABLES_WITHOUT_HEADER = {"resolutions", "chroms"}

def parse_hictk_table(stdout_text, table_name):
    lines = stdout_text.strip().split("\n")
    if table_name not in TABLES_WITHOUT_HEADER:
        lines = lines[1:]  # 跳过表头
    return [line.split() for line in lines if line.strip()]
```

**常见错误**：

```python
# ❌ 错误：对 resolutions/chroms 使用 lines[1:]，会丢失第一个分辨率
for line in lines[1:]:
    resolutions.append(int(line.split()[0]))

# ✅ 正确：对 resolutions/chroms 使用 lines[:]
for line in lines:
    resolutions.append(int(line.split()[0]))
```

---

## hictk balance

**Purpose**: Normalize Hi-C matrices using ICE, SCALE, or VC.

### Subcommands

| Subcommand | Method | Description |
|------------|--------|-------------|
| `hictk balance ice` | ICE | Iterative correction |
| `hictk balance scale` | SCALE | Single-step correction |
| `hictk balance vc` | VC | Vanilla coverage |

### hictk balance ice

```bash
hictk balance ice [OPTIONS] input
```

**Key Options**:

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--mode` | enum | gw | gw/trans/cis |
| `--tmpdir` | path | - | Temporary directory |
| `--ignore-diags` | uint | 2 | Diagonals to mask |
| `--mad-max` | float | 5 | MAD-max threshold |
| `--min-nnz` | uint | 10 | Min non-zero entries |
| `--min-count` | uint | 0 | Min interaction count |
| `--tolerance` | float | 1e-05 | Convergence threshold |
| `--max-iters` | uint | 500 | Max iterations |
| `--rescale-weights` | flag | false | Rescale to sum=2 |
| `--name` | text | auto | Weight dataset name |
| `--in-memory` | flag | false | Faster but memory-heavy |
| `--stdout` | flag | false | Write to stdout |
| `-t, --threads` | uint | 1 | Thread count |
| `-f, --force` | flag | false | Overwrite |

### Examples

```bash
# Basic ICE balancing
hictk balance ice sample.cool

# Multi-threaded with temp dir
hictk balance ice sample.mcool --threads 4 --tmpdir /tmp

# In-memory mode (faster)
hictk balance ice sample.cool --in-memory

# Cis-only mode
hictk balance ice sample.cool --mode cis

# Custom parameters
hictk balance ice sample.cool --max-iters 1000 --tolerance 1e-06

# Output weights to stdout
hictk balance ice sample.cool --stdout > weights.tsv
```

---

## hictk convert

**Purpose**: Convert between Hi-C file formats.

### Signatures

```bash
hictk convert [OPTIONS] input output
```

### Key Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-r, --resolutions` | uint[] | all | Target resolutions |
| `--normalization-methods` | text[] | ALL | Normalization vectors |
| `-g, --genome` | text | auto | Genome assembly |
| `--chunk-size` | uint | 10M | Batch size |
| `-t, --threads` | uint | 2 | Thread count |
| `-f, --force` | flag | false | Overwrite |
| `--count-type` | enum | auto | int/float/auto |

### Supported Conversions

| Input | Output | Notes |
|-------|--------|-------|
| .hic | .cool | Single resolution |
| .hic | .mcool | All resolutions |
| .cool | .hic | Requires bin size |
| .cool | .mcool | Multi-resolution via zoomify |
| .mcool | .hic | Extract resolution |

### Examples

```bash
# hic to cool (single resolution)
hictk convert sample.hic sample.cool

# hic to mcool (all resolutions)
hictk convert sample.hic sample.mcool

# Extract specific resolution
hictk convert sample.mcool extracted.cool -r 10000

# Specify genome
hictk convert sample.hic sample.cool -g hg38

# Copy specific normalizations only
hictk convert sample.hic sample.cool --normalization-methods SCALE
```

---

## hictk load

**Purpose**: Build .cool/.hic from text interaction files.

### Signatures

```bash
hictk load [OPTIONS] interactions output-path
```

### Input Formats

| Format | Option | Description |
|--------|--------|-------------|
| 4DN pairs | `--format 4dn` | Standard 4DN format |
| Valid pairs | `--format validpairs` | UHGG format |
| BG2 | `--format bg2` | BedGraph2 |
| COO | `--format coo` | Sparse matrix |

### Key Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-c, --chrom-sizes` | file | required* | Chromosome sizes |
| `-b, --bin-size` | uint | required* | Bin size in bp |
| `--bin-table` | file | required* | Custom BED file |
| `-f, --format` | enum | **required** | Input format |
| `--assembly` | text | unknown | Assembly name |
| `--one-based` | flag | - | Coordinates are 1-based |
| `--count-as-float` | flag | false | Float counts |
| `-t, --threads` | uint | 2 | Thread count |

*Either `--bin-size` or `--bin-table` required.

### Examples

```bash
# Load from 4DN format
hictk load pairs.tsv output.cool -c hg38.chrom.sizes -b 10000 --format 4dn

# Load from stdin
cat pairs.tsv | hictk load - output.cool -c hg38.sizes -b 5000 --format 4dn

# With custom bins
hictk load interactions.tsv output.cool --bin-table custom_bins.bed --format bg2

# Float counts
hictk load pairs.tsv output.cool -c hg38.sizes -b 10000 --format 4dn --count-as-float

# Compressed input
hictk load pairs.tsv.gz output.cool -c hg38.sizes -b 10000 --format 4dn
```

---

## hictk merge

**Purpose**: Merge multiple Hi-C files.

### Signatures

```bash
hictk merge [OPTIONS] input-files... -o output
```

### Key Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-o, --output-file` | path | **required** | Output file |
| `--output-fmt` | enum | auto | cool/hic |
| `--resolution` | uint | - | For .hic input |
| `--count-type` | enum | int | int/float |
| `--skip-all-vs-all` | flag | false | Skip symmetric matrix |
| `-t, --threads` | uint | 1 | Thread count |

### Examples

```bash
# Merge two coolers
hictk merge sample1.cool sample2.cool -o merged.cool

# Merge hic files
hictk merge rep1.hic rep2.hic -o merged.hic -r 10000

# Float counts
hictk merge sample1.cool sample2.cool -o merged.cool --count-type float
```

---

## hictk metadata

**Purpose**: Print file metadata in structured format.

### Signatures

```bash
hictk metadata [OPTIONS] uri
```

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-f, --output-format` | enum | json | json/toml/yaml |
| `--include-file-path` | flag | false | Include URI |
| `--recursive` | flag | false | All resolutions |

### Examples

```bash
# JSON metadata
hictk metadata sample.hic

# YAML format
hictk metadata sample.mcool --output-format yaml

# All resolutions (for mcool)
hictk metadata sample.mcool --recursive

# Include file path
hictk metadata sample.cool --include-file-path
```

---

## hictk zoomify

**Purpose**: Create multi-resolution files by coarsening.

### Signatures

```bash
hictk zoomify [OPTIONS] input output
```

### Key Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--resolutions` | uint[] | auto | Target resolutions |
| `--nice-steps` | flag | true | Nice (vs pow2) progression |
| `--copy-base-resolution` | flag | false | Include base resolution |
| `-t, --threads` | uint | 1 | Thread count |

### Resolution Generation

| Mode | Example (base=1000) |
|------|---------------------|
| `--nice-steps` | 1000, 2000, 5000, 10000... |
| `--pow2-steps` | 1000, 2000, 4000, 8000... |

### Examples

```bash
# Basic zoomify (nice steps)
hictk zoomify sample.cool sample.mcool

# Custom resolutions
hictk zoomify sample.cool sample.mcool --resolutions 1000 5000 10000 25000

# Power of 2 steps
hictk zoomify sample.cool sample.mcool --pow2-steps

# Include base resolution
hictk zoomify sample.cool sample.mcool --copy-base-resolution
```

---

## hictk fix-mcool

**Purpose**: Repair corrupted .mcool files (index corruption).

### Signatures

```bash
hictk fix-mcool [OPTIONS] input output
```

### Key Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--skip-balancing` | flag | false | Skip balancing step |
| `--check-base-resolution` | flag | false | Verify base resolution |
| `--in-memory` | flag | false | Faster balancing |
| `-t, --threads` | uint | 1 | Thread count |

### Examples

```bash
# Fix corrupted mcool
hictk fix-mcool corrupted.mcool fixed.mcool

# Skip balancing (faster, if already balanced)
hictk fix-mcool corrupted.mcool fixed.mcool --skip-balancing

# Check base resolution integrity
hictk fix-mcool corrupted.mcool fixed.mcool --check-base-resolution
```

---

## hictk rename-chromosomes

**Purpose**: Change chromosome naming conventions.

### Signatures

```bash
hictk rename-chromosomes [OPTIONS] uri
```

### Options

| Option | Description |
|--------|-------------|
| `--add-chr-prefix` | Add "chr" prefix |
| `--remove-chr-prefix` | Remove "chr" prefix |
| `--name-mappings FILE` | Custom TSV mappings |

### Examples

```bash
# Add chr prefix
hictk rename-chromosomes sample.cool --add-chr-prefix

# Remove chr prefix
hictk rename-chromosomes sample.cool --remove-chr-prefix

# Custom mappings (TSV: old_name new_name)
hictk rename-chromosomes sample.cool --name-mappings rename.tsv
```

---

## File Formats

### Supported Formats

| Extension | Type | Description |
|-----------|------|-------------|
| `.cool` | Single resolution | Cooler format |
| `.mcool` | Multi-resolution | Multiple resolutions |
| `.scool` | Single-cell | Per-cell matrices |
| `.hic` | Multi-resolution | Juicer format (incl. v9) |

### Cooler URI Syntax

Access specific resolution in mcool:
```
file.mcool::/resolutions/10000
```

---

## Common Workflows

### 1. Validate and Inspect

```bash
# Validate file
hictk validate sample.cool

# Check metadata
hictk metadata sample.hic

# View chromosomes
hictk dump sample.hic --table chroms

# Check available resolutions
hictk dump sample.hic --table resolutions
```

### 2. Extract Region

```bash
# Extract region from hic
hictk dump sample.hic --resolution 10000 --range chr1:0-1000000 --join

# Extract trans interactions
hictk dump sample.hic --resolution 10000 --range chr1 --range2 chr2 --join
```

### 3. Balance and Export

```bash
# Balance matrix
hictk balance ice sample.cool --in-memory

# Export balanced data
hictk dump sample.cool --balance ICE --join > balanced.tsv
```

### 4. Format Conversion

```bash
# hic to mcool
hictk convert sample.hic sample.mcool

# Extract resolution
hictk convert sample.mcool extracted.cool -r 10000

# hic to single cool
hictk convert sample.hic sample.cool -r 10000
```

### 5. Create from Pairs

```bash
# Load 4DN pairs
hictk load pairs.tsv output.cool -c hg38.sizes -b 10000 --format 4dn

# Balance
hictk balance ice output.cool

# Zoomify to multi-resolution
hictk zoomify output.cool output.mcool
```

---

## Error Handling

| Error | Cause | Solution |
|-------|-------|----------|
| `Invalid resolution` | Resolution not in file | Check available: `hictk dump file --table resolutions` |
| `Chromosome not found` | Wrong naming | Check: `hictk dump file --table chroms` |
| `Duplicate pixels` | Corrupted index | Use `hictk validate --validate-index` |
| `Bins not sorted` | Data corruption | May need `hictk fix-mcool` |

---

## Performance Tips

1. **Use `--in-memory`** for `balance` when RAM allows
2. **Increase `--threads`** for convert/zoomify/balance
3. **Use `--tmpdir`** on fast storage for large files
4. **Use `--chunk-size`** tuning for memory/speed tradeoff
5. **Pre-validate** with `hictk validate` before processing
