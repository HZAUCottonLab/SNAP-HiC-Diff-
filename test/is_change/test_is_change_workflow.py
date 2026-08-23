from pathlib import Path


REPO = Path(__file__).resolve().parents[2]


def test_is_change_rules_aggregate_all_scales_under_the_multiscale_root():
    snakefile = (REPO / "Snakefile").read_text()
    rules = (REPO / "rules" / "is_change.smk").read_text()

    assert "MULTISCALE_ROOT" in snakefile
    assert "WINDOW_BINS" in snakefile
    assert 'include: "rules/meltron.smk"' in snakefile
    assert 'include: "rules/is_change.smk"' in snakefile

    assert "rule build_is_change_candidates:" in rules
    assert "window_bins=WINDOW_BINS" in rules
    assert 'MULTISCALE_ROOT + "/{comparison}/all_scale_candidates.tsv"' in rules
    assert 'MULTISCALE_ROOT + "/{comparison}/candidate_union.tsv"' in rules
    assert 'MULTISCALE_ROOT + "/{comparison}/candidate_union.bed"' in rules
    assert "--expected-window-bins {params.expected_window_bins}" in rules
    assert "--out-all" in rules
    assert "--out-regions" in rules
    assert "--out-bed" in rules

    assert "rule call_is_change_loci:" in rules
    assert 'MULTISCALE_ROOT + "/{comparison}/candidate_union.tsv"' in rules
    assert 'MULTISCALE_ROOT + "/{comparison}/is_change/{comparison}_is_change_loci.tsv"' in rules
    assert "get_comparison_is_file(w.comparison, 0)" in rules
    assert "get_comparison_is_file(w.comparison, 1)" in rules
    assert 'OUTDIR + "/balanced/{sample}/{sample}__{chrom}_iced.cool"' in rules

    assert "rule plot_is_change_loci:" in rules
    assert 'MULTISCALE_ROOT + "/{comparison}/is_change/figures/manifest.tsv"' in rules
    assert rules.count('conda:\n        "../env/SNAP-HiC-Diff.yaml"') == 3
    assert "is_change/logs" in rules
