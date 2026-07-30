from pathlib import Path
import re
import subprocess


REPO = Path(__file__).resolve().parents[2]
RULE = REPO / "rules" / "compute_IS.smk"
MELTRON_RULE = REPO / "rules" / "meltron.smk"
SAMPLES = [
    "H1394",
    "H1395",
    "H1396",
    "H1397",
    "H1398",
    "H1406",
    "H1409",
    "H1410",
    "H1436",
    "L021",
]
TARGETS = [f"results/IS/{sample}/{sample}_IS.tsv" for sample in SAMPLES]
MELTRON_COMPARISON = "H1398_vs_L134"
MELTRON_TARGET = (
    f"results/meltron/{MELTRON_COMPARISON}/"
    f"{MELTRON_COMPARISON}_meltron_merged.tsv"
)


def test_compute_is_schedules_one_cworld_job_per_sample_chromosome():
    result = subprocess.run(
        [
            "snakemake",
            "-n",
            "-p",
            "--configfile",
            "config/config.yaml",
            "--",
            *TARGETS,
        ],
        cwd=REPO,
        text=True,
        capture_output=True,
        check=False,
    )
    plan = result.stdout + result.stderr

    assert result.returncode == 0, plan
    assert re.search(r"^\s*compute_IS_chrom\s+260\s*$", plan, re.MULTILINE), plan
    assert re.search(r"^\s*merge_IS\s+10\s*$", plan, re.MULTILINE), plan
    assert "checkpoint compute_IS" not in RULE.read_text()
    assert 'config["meltron"]["is_distances"]' not in RULE.read_text()
    assert plan.count("matrix2insulation.pl") == 260, plan


def test_compute_is_rule_uses_chromosome_outputs_and_configured_cpus():
    rule = RULE.read_text()

    assert "rule compute_IS_chrom:" in rule
    assert "rule merge_IS:" in rule
    assert '"--expected-bin-size"' in rule
    assert "tmp/cworld_matrix" in rule
    assert "CWORLD_IS.chrom_job_threads" in rule
    assert "CWORLD_IS.merge_job_threads" in rule
    assert 'MANIFEST_TMP="$MANIFEST.tmp.$$"' in rule
    assert "find" in rule and "-newer" in rule


def test_meltron_schedules_per_sample_is_merges_without_compute_is_checkpoint():
    rule = MELTRON_RULE.read_text()

    assert "checkpoints.compute_IS" not in rule

    result = subprocess.run(
        [
            "snakemake",
            "-n",
            "-p",
            "--forceall",
            "--configfile",
            "config/config.yaml",
            "--",
            MELTRON_TARGET,
        ],
        cwd=REPO,
        text=True,
        capture_output=True,
        check=False,
    )
    plan = result.stdout + result.stderr

    assert result.returncode == 0, plan
    assert re.search(r"^\s*compute_IS_chrom\s+52\s*$", plan, re.MULTILINE), plan
    assert re.search(r"^\s*merge_IS\s+2\s*$", plan, re.MULTILINE), plan
    assert re.search(r"^\s*meltron_per_comparison\s+1\s*$", plan, re.MULTILINE), plan
    assert "results/IS/H1398/H1398_IS.tsv" in plan
    assert "results/IS/L134/L134_IS.tsv" in plan


def test_meltron_cores_are_params_and_precede_logging():
    rule = MELTRON_RULE.read_text()

    params_start = rule.index("    params:\n")
    conda_start = rule.index("    conda:\n")
    core_param = rule.index("        meltron_cores = 20,\n")
    assert params_start < core_param < conda_start

    for direction in ("decondensation", "condensation"):
        run_start = rule.index(f"--comparison-direction {direction}")
        log_start = rule.index('--log-file "{log}"', run_start)
        assert run_start < rule.index("--cutoff-na", run_start, log_start)
        assert run_start < rule.index("--cores {params.meltron_cores}", run_start, log_start)


def test_rule_all_uses_final_is_tsv_without_removed_done_markers():
    snakefile = (REPO / "Snakefile").read_text()

    assert 'expand(OUTDIR + "/logs/compute_IS/{sample}.done", sample=SAMPLES),' not in snakefile
    assert 'expand(OUTDIR + "/IS/{sample}/{sample}_IS.tsv", sample=SAMPLES),' in snakefile
