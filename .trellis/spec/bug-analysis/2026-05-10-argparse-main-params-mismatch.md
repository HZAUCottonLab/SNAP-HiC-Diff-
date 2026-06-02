# Bug: argparse 参数名与 main() 函数签名不匹配

**Date**: 2026-05-10
**Severity**: HIGH (blocks execution)
**Status**: FIXED (but RECURRING)

## Problem

运行 `compute_gap_qc` 和 `sliding_window` 失败：

```
TypeError: main() got an unexpected keyword argument 'log'
TypeError: main() got an unexpected keyword argument 'threads'
```

## Root Cause

Python 脚本中 argparse 参数名与 `main()` 函数参数名不一致：

```python
# argparse 定义（问题：argparse 名与 main() 名不同）
parser.add_argument("--log", required=True)      # → vars(args) 产生 key: "log"
parser.add_argument("--threads", type=int)      # → vars(args) 产生 key: "threads"

# main() 函数签名（期望不同的名字）
def main(..., log_file, n_jobs):  # 期望 key: "log_file", "n_jobs"
```

当 `main(**vars(args))` 时，argparse 转换后的 key 与 `main()` 期望的参数名不匹配。

## Affected Scripts

- `scripts/compute_gap_qc.py`
- `scripts/sliding_window.py`

## Solution

统一参数名：argparse 参数名必须与 `main()` 签名完全一致

1. `--log` → `--log-file` (main 期望 `log_file`)
2. `--threads` → `--n-jobs` (main 期望 `n_jobs`)

## Changes

1. `scripts/compute_gap_qc.py`: `--log` → `--log-file`, `--threads` → `--n-jobs`
2. `scripts/sliding_window.py`: `--log` → `--log-file`, `--threads` → `--n-jobs`
3. `rules/compute_gap_qc.smk`: `--log` → `--log-file`, `--threads` → `--n-jobs`
4. `rules/sliding_window.smk`: `--log` → `--log-file`, `--threads` → `--n-jobs`

## Prevention

**Project Rule**: argparse 参数名必须与 `main()` 函数签名完全一致

- 推荐做法 1：直接使用与 `main()` 一致的参数名（推荐）
- 推荐做法 2：如果 argparse 需要不同名字，使用 `dest` 参数或手动映射

```python
# 推荐做法 1（简单一致）
parser.add_argument("--log-file", required=True)  # main() 也用 log_file
parser.add_argument("--n-jobs", type=int)          # main() 也用 n_jobs

# 推荐做法 2（手动映射）
args = parser.parse_args()
sys.exit(0 if main(log_file=args.log, n_jobs=args.threads, **vars(args)) else 1)
```

## Related Files

- `rules/compute_gap_qc.smk`
- `rules/sliding_window.smk`
- `scripts/compute_gap_qc.py`
- `scripts/sliding_window.py`
