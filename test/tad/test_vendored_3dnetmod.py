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


def test_vendored_preprocessor_does_not_hard_filter_chrM():
    source = (CODE_DIR / "slice_genomic_coordinates_genomewide_majregv.py").read_text()
    assert "if chr == 'chrM':" not in source


def test_vendored_preprocessor_defaults_to_no_chromosome_exclusions(monkeypatch):
    monkeypatch.syspath_prepend(str(CODE_DIR))
    preprocess = runpy.run_path(str(CODE_DIR / "3DNetMod_preprocess.py"))
    configured = preprocess["configured_excluded_chromosomes"]

    assert configured({}) == set()
    assert configured({"excluded_chromosomes": "chrM, chrY"}) == {"chrM", "chrY"}


def test_vendored_stages_use_settings_processor_limits_without_taskset_affinity():
    source_by_name = {
        path.name: path.read_text(errors="replace")
        for path in CODE_DIR.glob("*.py")
    }
    source = "\n".join(source_by_name.values())

    assert "taskset" not in source
    assert "multiprocessing.Pool(4)" not in source
    assert "Pool(processes=5)" not in source
    assert "processors = 23" not in source
    assert "def worker_count(settings, work_items):" in source_by_name[
        "processor_limits.py"
    ]
    for name in (
        "3DNetMod_preprocess.py",
        "3DNetMod_GPS_MMCP.py",
        "3DNetMod_GPS_MMCP_perchr.py",
        "find_max_gamma.py",
        "write_consistent_communities.py",
    ):
        assert "worker_count(settings," in source_by_name[name]


def test_consensus_filename_is_deterministic_and_filesystem_safe():
    helper = runpy.run_path(str(CODE_DIR / "outdir_helper.py"))
    name = helper["consensus_filtered_filename"](
        "sample", "long_parameter_tag_" * 20, 0.4050000000000002
    )
    assert name == helper["consensus_filtered_filename"](
        "sample", "long_parameter_tag_" * 20, 0.4050000000000002
    )
    assert len(name.encode()) <= 255


def test_variance_thresholds_handle_empty_communities():
    module = runpy.run_path(str(CODE_DIR / "get_variance_distribution_and_threshold.py"))
    threshold = module["get_variance_distribution_and_threshold"]

    distributions, thresholds, selected = threshold(
        {200: [], 400: [{"left_var": "0.1", "right_var": "0.2"}]},
        [200, 400],
        [45, 45],
    )

    assert distributions[200] == []
    assert thresholds[200] == 0
    assert selected[200] == []
    assert thresholds[400] == 0.2
    assert selected[400] == [{"left_var": "0.1", "right_var": "0.2"}]


def test_gps_skips_regions_without_counts_or_bed_files():
    source = (CODE_DIR / "heatmap_removal.py").read_text()
    guard = "if not os.path.isfile(counts_file) or not os.path.isfile(bed_file):"
    assert source.count(guard) == 3
    assert "skipping region without counts or bed file:" in source


def test_hsvm_writes_empty_merged_output_without_bin_size():
    source = (CODE_DIR / "hierarchical_spatial_variance_minimization.py").read_text()
    assert "if total_thresholded_communities:" in source
    assert "else:\n                total_merged_communities = []" in source
