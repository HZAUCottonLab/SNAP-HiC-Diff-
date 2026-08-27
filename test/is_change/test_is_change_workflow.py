from pathlib import Path


REPO = Path(__file__).resolve().parents[2]


def test_is_change_rules_aggregate_ordered_scale_inputs_to_shared_targets():
    snakefile = (REPO / "Snakefile").read_text()
    rules_path = REPO / "rules" / "is_change.smk"

    assert rules_path.exists()
    rules = rules_path.read_text()

    assert "def get_comparison_is_file(comparison, member_index):" in snakefile
    assert 'include: "rules/is_change.smk"' in snakefile
    assert snakefile.index('include: "rules/meltron.smk"') < snakefile.index(
        'include: "rules/is_change.smk"'
    )
    assert "MULTISCALE_ALL_OUTPUTS" in snakefile

    assert "rule build_is_change_candidates:" in rules
    assert "MULTISCALE_ROOT" in rules
    assert "SCALE_LABELS" in rules
    assert "all_scale_candidates.tsv" in rules
    assert "candidate_union.tsv" in rules
    assert "candidate_union.bed" in rules
    assert "knn_enriched=lambda w: expand(" in rules
    assert "meltron_merged=lambda w: expand(" in rules
    assert "build_is_change_candidates.py" in rules
    assert "--knn-enriched {input.knn_enriched}" in rules
    assert "--meltron-merged {input.meltron_merged}" in rules
    assert "--out-candidates \"{output.candidates}\"" in rules
    assert "--out-union-tsv \"{output.union_tsv}\"" in rules
    assert "--out-union-bed \"{output.union_bed}\"" in rules

    assert "rule call_is_change_loci:" in rules
    assert "candidate_union.tsv" in rules
    assert "--candidate-regions \"{input.regions}\"" in rules
    assert "get_comparison_is_file(w.comparison, 0)" in rules
    assert "get_comparison_is_file(w.comparison, 1)" in rules
    assert "_iced.cool" in rules
    assert "call_is_change_loci.py" in rules

    assert "rule plot_is_change_loci:" in rules
    assert "_is_change_loci.tsv" in rules
    assert "plot_is_change_loci.py" in rules
    assert "manifest.tsv" in rules
    assert rules.count('cool_a=lambda w: expand') == 2
    assert rules.count("--cool-a {input.cool_a}") == 2
    assert rules.count("--cool-b {input.cool_b}") == 2

    assert rules.count('conda:\n        "../env/SNAP-HiC-Diff.yaml"') == 3
    assert rules.count("--config config/config.yaml") == 3
    assert "logs/is_change" in rules
    assert "_candidate_regions.tsv" not in rules
