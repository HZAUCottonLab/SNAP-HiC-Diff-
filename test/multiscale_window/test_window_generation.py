from pathlib import Path
import re
import sys


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from build_windows import build_windows


def test_build_windows_records_the_active_scale_on_every_window():
    windows = build_windows(
        {"chr1": 140_000},
        window_size=100_000,
        step_size=20_000,
        resolution=4_000,
        window_bins=25,
        step_bins=5,
    )

    assert not windows.empty
    assert windows[
        [
            "resolution_bp",
            "window_bins",
            "window_size_bp",
            "step_bins",
            "step_size_bp",
        ]
    ].drop_duplicates().to_dict("records") == [
        {
            "resolution_bp": 4_000,
            "window_bins": 25,
            "window_size_bp": 100_000,
            "step_bins": 5,
            "step_size_bp": 20_000,
        }
    ]


def test_window_preprocessing_paths_are_scale_qualified_without_moving_shared_inputs():
    snakefile = (REPO / "Snakefile").read_text()
    assert 'SCALE_ROOT = "{}/window_scales/r{}/w{{window_bins}}".format(OUTDIR, RESOLUTION)' in snakefile

    build_windows_rule = (REPO / "rules" / "build_windows.smk").read_text()
    gap_qc_rule = (REPO / "rules" / "compute_gap_qc.smk").read_text()
    sliding_rule = (REPO / "rules" / "sliding_window.smk").read_text()

    for rule in (build_windows_rule, gap_qc_rule, sliding_rule):
        assert "SCALE_ROOT" in rule
        assert "window_bins" in rule

    assert gap_qc_rule.count('window_bins="{{window_bins}}"') == 2
    assert sliding_rule.count('window_bins="{{window_bins}}"') == 2
    assert "shared_outdir = OUTDIR" in gap_qc_rule
    assert "{params.shared_outdir}/balanced/" in gap_qc_rule
    assert "shared_outdir = OUTDIR" in sliding_rule
    assert "{params.shared_outdir}/oe/" in sliding_rule
    assert "window_bins=w.window_bins" in sliding_rule


def test_all_window_rule_declarations_are_scale_qualified():
    expected_paths = {
        "build_windows.smk": [
            ("output", "/windows/windows.parquet"),
            ("log", "/logs/build_windows.log"),
            ("benchmark", "/benchmark/build_windows/build_windows.tsv"),
        ],
        "compute_gap_qc.smk": [
            ("output", "/logs/gap_qc/{sample}.done"),
            ("output", "/windows/{{sample}}/{{sample}}__{chrom}_gaps.tsv"),
            ("output", "/windows/{{sample}}/{{sample}}__{chrom}_window_qc.tsv"),
            ("log", "/logs/gap_qc/{sample}.log"),
            ("benchmark", "/benchmark/gap_qc/{sample}.tsv"),
        ],
        "sliding_window.smk": [
            ("output", "/logs/windows/{sample}.done"),
            ("log", "/logs/windows/{sample}.log"),
            ("benchmark", "/benchmark/windows/{sample}.tsv"),
            ("output", "/windows/{{sample}}/{{sample}}__{{chrom}}_{{stage}}.h5"),
        ],
    }

    declaration_blocks = re.compile(
        r"(?ms)^    (output|log|benchmark):\n(.*?)(?=^    [a-z_]+:|\Z)"
    )
    path_literals = re.compile(r'"([^"]*/[^"]*)"')

    for filename, expected in expected_paths.items():
        rule_text = (REPO / "rules" / filename).read_text()
        actual = []
        for declaration, body in declaration_blocks.findall(rule_text):
            for path_match in path_literals.finditer(body):
                prefix = body[:path_match.start()].rstrip()
                assert prefix.endswith("SCALE_ROOT +"), (
                    f"{filename} {declaration} path is not scale-qualified: "
                    f"{path_match.group(1)}"
                )
                actual.append((declaration, path_match.group(1)))

        assert actual == expected
