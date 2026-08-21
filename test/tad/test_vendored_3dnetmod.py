from pathlib import Path
import runpy


REPO = Path(__file__).resolve().parents[2]
CODE_DIR = REPO / "scripts" / "python3_3dnetmod"


def test_vendored_entrypoints_recursive_closure_and_seed_exist():
    required = {
        "3DNetMod_preprocess.py",
        "3DNetMod_GPS_MMCP.py",
        "3DNetMod_HSVM.py",
        "outdir_helper.py",
        "read_settings.py",
        "read_seeds.py",
        "gen_louvain_100_seeds.txt",
        "find_bad_regions_mean_var_simplified_per_chrom.py",
        "load_bad_regions_simplified_per_chrom.py",
        "find_max_gamma_simplified_per_chrom.py",
        "ThreeDNetMod_GPS_MMCP_per_chr_simplified_per_chrom.py",
    }
    assert required <= {path.name for path in CODE_DIR.iterdir()}


def test_vendored_code_has_no_runtime_reference_to_source_checkout():
    source = "\n".join(
        path.read_text(errors="replace")
        for path in CODE_DIR.glob("*.py")
    )
    assert "/share/home/xhhuang24/project/HiC_analysis_pipeline" not in source
    assert "PanTE/compare_mtx" not in source


def test_consensus_filename_is_deterministic_and_filesystem_safe():
    helper = runpy.run_path(str(CODE_DIR / "outdir_helper.py"))
    name = helper["consensus_filtered_filename"](
        "sample", "long_parameter_tag_" * 20, 0.4050000000000002
    )
    assert name == helper["consensus_filtered_filename"](
        "sample", "long_parameter_tag_" * 20, 0.4050000000000002
    )
    assert len(name.encode()) <= 255
