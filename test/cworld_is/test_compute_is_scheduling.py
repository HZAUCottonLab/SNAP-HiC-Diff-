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
MELTRON_WINDOW_BINS = 25
MELTRON_TARGET = (
    f"window_scales/r4000/w{MELTRON_WINDOW_BINS}/meltron/{MELTRON_COMPARISON}/"
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


def test_merge_is_passes_configured_source_resolution_to_converter():
    rule = (REPO / "rules" / "compute_IS.smk").read_text()
    params_start = rule.index("rule merge_IS:")
    shell_start = rule.index("    shell:", params_start)
    assert "expected_bin_size=CWORLD_IS.source_resolution_bp" in rule[params_start:shell_start]
    assert '--expected-bin-size "{params.expected_bin_size}"' in rule[shell_start:]


def write_config_with_temporary_backgrounds(tmp_path):
    output_root = tmp_path / "results"
    scratch_root = tmp_path / "tmp"
    config_text = (REPO / "config" / "config.yaml").read_text()
    assert 'outdir: "results"' in config_text
    config_text = config_text.replace('outdir: "results"', f'outdir: "{output_root}"')
    assert 'tmpdir: "tmp"' in config_text
    config_text = config_text.replace('tmpdir: "tmp"', f'tmpdir: "{scratch_root}"')
    for window_bins in (25, 50, 75, 100, 125):
        background = tmp_path / f"background_r4000_w{window_bins}.tsv"
        background.write_text("comparison\tchrom\twindow_start\twindow_end\n")
        configured = (
            "/share/home/xhhuang24/project/PanTE/background_pool/results/"
            f"window_scales/r4000/w{window_bins}/window_compare/HiC_window_background.tsv"
        )
        assert configured in config_text
        config_text = config_text.replace(configured, str(background))
    config_path = tmp_path / "config.yaml"
    config_path.write_text(config_text)
    return config_path, output_root

def test_meltron_schedules_per_sample_is_merges_without_compute_is_checkpoint(tmp_path):
    rule = MELTRON_RULE.read_text()
    config_path, output_root = write_config_with_temporary_backgrounds(tmp_path)
    knn_input = output_root / (
        f"window_scales/r4000/w{MELTRON_WINDOW_BINS}/knn_empirical_pvalue/"
        f"{MELTRON_COMPARISON}/{MELTRON_COMPARISON}_knn_enriched.tsv"
    )
    knn_input.parent.mkdir(parents=True)
    knn_input.write_text("chrom\twindow_id\n")
    target = str(output_root / MELTRON_TARGET)

    assert "checkpoints.compute_IS" not in rule

    result = subprocess.run(
        [
            "snakemake",
            "-n",
            "-p",
            "--configfile",
            str(config_path),
            "--",
            target,
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
    assert f"{output_root}/IS/H1398/H1398_IS.tsv" in plan
    assert f"{output_root}/IS/L134/L134_IS.tsv" in plan


def test_meltron_cpu_request_and_cli_cores_share_one_config_value():
    rule = MELTRON_RULE.read_text()
    assert 'knn_enriched = SCALE_ROOT + "/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv"' in rule
    assert 'merged = SCALE_ROOT + "/meltron/{comparison}/{comparison}_meltron_merged.tsv"' in rule
    assert 'SCALE_ROOT + "/logs/meltron/{comparison}.log"' in rule
    assert 'SCALE_ROOT + "/benchmark/meltron/{comparison}.tsv"' in rule
    assert "--raw-pvalue-threshold 0.05" in rule
    assert "--fdr-threshold" not in rule
    config_text = (REPO / "config" / "config.yaml").read_text()

    threads_start = rule.index("    threads:\n")
    params_start = rule.index("    params:\n")
    conda_start = rule.index("    conda:\n")
    core_expression = 'config["meltron"]["cores"]'
    assert rule.index(core_expression, threads_start, params_start) > threads_start
    assert f"        meltron_cores = {core_expression},\n" in rule
    assert "    chrom_job_threads: 2" in config_text
    assert "    merge_job_threads: 2" in config_text

    core_param = rule.index(f"        meltron_cores = {core_expression},\n")
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
