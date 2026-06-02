# dsub Direct Command Template

> 方式一：命令行直接提交（适合简单任务或一次性作业）

## 模板

```bash
#!/bin/bash
#DSUB -q arm
#DSUB -aa
#DSUB -n <job-name>
#DSUB -e %J.err
#DSUB -o %J.out
#DSUB -R "cpu=10"

set -euo pipefail

# Disable nested parallelism
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# Your commands
python3 -m pytest tests/ -v
```

## 提交流程

```bash
chmod +x submit.sh
dsub -s submit.sh
```

## CPU 选择参考

| 规模 | CPU | 内存（arm 队列） |
|------|-----|-----------------|
| 小 | 4 核 | 20 GB |
| 中 | 8 核 | 40 GB |
| 大 | 16 核 | 80 GB |

> arm 队列：每核 5 GB。内存需求 > 500 GB 时改用 `fat_arm` 或 `fat_x86` 队列。
