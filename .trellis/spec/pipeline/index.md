# Pipeline Layer Specs

> Specs for the computational pipeline layer (Snakemake workflows, data processing).

---

## Snakemake / dsub

| Document | Purpose |
|----------|---------|
| [snakemake-guidelines.md](./snakemake-guidelines.md) | Full syntax reference: rules, wildcards, expand(), config, modules, conda, best practices — distilled from official docs v8.30.0 |
| [dsub-guidelines.md](./dsub-guidelines.md) | 多瑙调度器 (dsub) CLI reference for HZAU HPC: queues, resource limits, job dependencies, MPI/NPU, Snakemake v8 integration via cluster-generic, **以及测试作业提交约束** |
| [cluster-submission-ssot.md](./cluster-submission-ssot.md) | 锁定 Snakemake 集群提交方式，禁止硬编码绕过 SSOT 配置 |

## Hi-C Tools

| Document | Purpose |
|----------|---------|
| [hictk.md](./hictk.md) | hictk CLI 速查手册：validate、dump、balance、convert、load、merge、zoomify |
| [hicexplorer-knowledge.md](./hicexplorer-knowledge.md) | HiCExplorer 知识库：normalize、findTADs、detectLoops、hicCorrelate、hicPlotMatrix |
