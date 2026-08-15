from pathlib import Path


REPO = Path(__file__).resolve().parents[2]


def test_is_change_rules_wire_existing_handoffs_to_final_targets():
    snakefile = (REPO / "Snakefile").read_text()
    rules_path = REPO / "rules" / "is_change.smk"

    assert rules_path.exists()
    rules = rules_path.read_text()

    assert "def get_comparison_is_file(comparison, member_index):" in snakefile
    assert 'include: "rules/is_change.smk"' in snakefile
    assert snakefile.index('include: "rules/meltron.smk"') < snakefile.index(
        'include: "rules/is_change.smk"'
    )
    assert "_is_change_loci.tsv" in snakefile
    assert "figures/is_change" in snakefile

    assert "rule build_is_change_candidates:" in rules
    assert "_knn_enriched.tsv" in rules
    assert "_meltron_merged.tsv" in rules
    assert "build_is_change_candidates.py" in rules

    assert "rule call_is_change_loci:" in rules
    assert "_candidate_regions.tsv" in rules
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
