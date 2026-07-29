# Meltron 4 kb Per-Chromosome CWorld IS Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (- [ ]) syntax for tracking.

**Goal:** Produce Meltron-ready insulation-score tables exclusively from 4 kb CWorld matrices, with one 4-CPU dsub/Snakemake task per sample × chromosome and one validated merge task per sample.

**Architecture:** Parse one nested meltron.cworld_is configuration object when Snakemake loads. Replace the current per-sample compute_IS checkpoint with a chromosome-level rule that creates one manifest after ten sequential CWorld calls; a separate merge_IS rule consumes every manifest and retains the existing results/IS/<sample>/<sample>_IS.tsv contract.

**Tech Stack:** Snakemake, Bash, Python 3.10+, cooler, pandas, CWorld-Dekker Perl scripts, pytest, dsub/Donau.

## Global Constraints

- The source CWorld input resolution is exactly 4000 bp; CWorld never runs on 20 kb or 100 kb matrices.
- IS distances are exactly 32000, 60000, 92000, 120000, 152000, 180000, 212000, 240000, 272000, 300000 bp.
- Every distance is positive, unique, divisible by 4000, and uses --ids = --is / 2.
- Every compute_IS_chrom job requests 4 CPUs by default. Do not set memory separately because Donau derives memory from CPUs.
- One task computes one sample × one chromosome × all ten scales sequentially. One merge task combines all chromosome manifests for a sample.
- Matrix files belong in tmp/cworld_matrix/<sample>/; bedGraphs and manifests belong in results/cworld_insulation/<sample>/.
- The final Meltron table remains results/IS/<sample>/<sample>_IS.tsv with columns chrom, viewpoint, IS_distance, IS.
- Cleanup accepts explicit samples, defaults to dry-run, and removes only the four approved output classes after --execute.
- Preserve the existing offline Conda and dsub-wrapper behavior.

---

## File Structure

| File | Action | Responsibility |
| --- | --- | --- |
| scripts/cworld_is_config.py | Create | Parse and validate fixed CWorld IS settings. |
| config/config.yaml | Modify | Store the resolution, ten scales, IDs divisor, and CPU counts. |
| Snakefile | Modify | Construct validated CWORLD_IS once before rule inclusion. |
| scripts/convert_cool_to_cworld.py | Modify | Reject a mismatched bin size when invoked by the production rule. |
| scripts/convert_cworld_IS_to_meltron_IS.py | Modify | Merge only supplied manifests and require full scale coverage. |
| rules/compute_IS.smk | Replace | Implement compute_IS_chrom and merge_IS. |
| scripts/cleanup_cworld_is.py | Create | Preview and execute explicit-sample cleanup safely. |
| test/cworld_is/test_cworld_is_config.py | Create | Validate scale and CPU configuration. |
| test/cworld_is/test_convert_cool_to_cworld.py | Create | Test bin-size guard and CWorld header. |
| test/cworld_is/test_convert_cworld_is.py | Modify | Test strict manifest merge behavior. |
| test/cworld_is/test_compute_is_scheduling.py | Modify | Test 26 chromosome jobs and one merge per sample. |
| test/cworld_is/test_cleanup_cworld_is.py | Create | Test dry-run and approved deletion scope. |
| README.md | Modify | Explain constraints, metrics, outputs, cleanup, and restart. |
| docs/CODEBASE_OVERVIEW.md | Modify | Keep architecture summary current. |
| docs/SNAKEMAKE_PIPELINE_MAP.md | Modify | Keep scheduling topology current. |

## Task 1: Establish the single 4 kb configuration contract

**Files:**

- Create: scripts/cworld_is_config.py
- Create: test/cworld_is/test_cworld_is_config.py
- Modify: config/config.yaml
- Modify: Snakefile

**Interfaces:**

- Consumes: the top-level Snakemake config mapping.
- Produces: CworldISSettings(source_resolution_bp, is_distances_bp, ids_divisor, chrom_job_threads, merge_job_threads).
- Produces: CWORLD_IS, a parsed Snakefile global visible to rules/compute_IS.smk.

- [ ] **Step 1: Write the failing parser tests**

Create test/cworld_is/test_cworld_is_config.py:

~~~python
from pathlib import Path
import sys

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
from cworld_is_config import ConfigError, parse_cworld_is_settings

DISTANCES = (
    32000, 60000, 92000, 120000, 152000,
    180000, 212000, 240000, 272000, 300000,
)
IDS = (
    16000, 30000, 46000, 60000, 76000,
    90000, 106000, 120000, 136000, 150000,
)

def valid_config():
    return {
        "meltron": {
            "cworld_is": {
                "source_resolution_bp": 4000,
                "is_distances_bp": list(DISTANCES),
                "ids_divisor": 2,
                "chrom_job_threads": 4,
                "merge_job_threads": 4,
            }
        }
    }

def test_parses_exact_scale_series_and_ids():
    settings = parse_cworld_is_settings(valid_config())

    assert settings.source_resolution_bp == 4000
    assert settings.is_distances_bp == DISTANCES
    assert settings.ids_distances_bp == IDS
    assert settings.chrom_job_threads == settings.merge_job_threads == 4

@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda c: c["meltron"]["cworld_is"].update(source_resolution_bp=20000), "4000"),
        (lambda c: c["meltron"]["cworld_is"].update(is_distances_bp=list(DISTANCES[:-1])), "exactly 10"),
        (lambda c: c["meltron"]["cworld_is"].update(is_distances_bp=[*DISTANCES[:-1], 30001]), "divisible"),
        (lambda c: c["meltron"]["cworld_is"].update(is_distances_bp=[*DISTANCES[:-1], 272000]), "unique"),
        (lambda c: c["meltron"]["cworld_is"].update(ids_divisor=3), "ids_divisor"),
        (lambda c: c["meltron"]["cworld_is"].update(chrom_job_threads=0), "positive"),
    ],
)
def test_rejects_invalid_constraints(mutate, message):
    config = valid_config()
    mutate(config)
    with pytest.raises(ConfigError, match=message):
        parse_cworld_is_settings(config)
~~~

- [ ] **Step 2: Confirm the test initially fails**

Run:

~~~bash
pytest -q test/cworld_is/test_cworld_is_config.py
~~~

Expected: collection fails because cworld_is_config does not exist.

- [ ] **Step 3: Implement the parser and configuration**

Create scripts/cworld_is_config.py:

~~~python
from dataclasses import dataclass
from typing import Any, Mapping

EXPECTED_RESOLUTION_BP = 4000
EXPECTED_SCALE_COUNT = 10

class ConfigError(ValueError):
    pass

@dataclass(frozen=True)
class CworldISSettings:
    source_resolution_bp: int
    is_distances_bp: tuple[int, ...]
    ids_divisor: int
    chrom_job_threads: int
    merge_job_threads: int

    @property
    def ids_distances_bp(self) -> tuple[int, ...]:
        return tuple(distance // self.ids_divisor for distance in self.is_distances_bp)

def parse_cworld_is_settings(config: Mapping[str, Any]) -> CworldISSettings:
    try:
        raw = config["meltron"]["cworld_is"]
    except (KeyError, TypeError) as exc:
        raise ConfigError("meltron.cworld_is is required") from exc
    if not isinstance(raw, Mapping):
        raise ConfigError("meltron.cworld_is must be a mapping")

    keys = (
        "source_resolution_bp", "is_distances_bp", "ids_divisor",
        "chrom_job_threads", "merge_job_threads",
    )
    missing = [key for key in keys if key not in raw]
    if missing:
        raise ConfigError(f"meltron.cworld_is missing keys: {', '.join(missing)}")

    resolution = raw["source_resolution_bp"]
    distances = raw["is_distances_bp"]
    divisor = raw["ids_divisor"]
    chrom_threads = raw["chrom_job_threads"]
    merge_threads = raw["merge_job_threads"]
    scalar_values = (resolution, divisor, chrom_threads, merge_threads)
    if any(isinstance(value, bool) or not isinstance(value, int) for value in scalar_values):
        raise ConfigError("resolution, ids_divisor, and thread counts must be integers")
    if resolution != EXPECTED_RESOLUTION_BP:
        raise ConfigError("source_resolution_bp must be 4000")
    if divisor != 2:
        raise ConfigError("ids_divisor must be 2 so --ids equals --is / 2")
    if chrom_threads <= 0 or merge_threads <= 0:
        raise ConfigError("thread counts must be positive")
    if not isinstance(distances, list) or any(
        isinstance(value, bool) or not isinstance(value, int) for value in distances
    ):
        raise ConfigError("is_distances_bp must be a list of integers")
    if len(distances) != EXPECTED_SCALE_COUNT:
        raise ConfigError("is_distances_bp must contain exactly 10 values")
    if len(set(distances)) != len(distances):
        raise ConfigError("is_distances_bp values must be unique")
    if any(value <= 0 or value % resolution for value in distances):
        raise ConfigError("every IS distance must be positive and divisible by 4000")
    return CworldISSettings(
        resolution, tuple(distances), divisor, chrom_threads, merge_threads
    )
~~~

Replace the old compute_is block and meltron.is_distances with:

~~~yaml
meltron:
  cworld_is:
    source_resolution_bp: 4000
    is_distances_bp: [32000, 60000, 92000, 120000, 152000,
                      180000, 212000, 240000, 272000, 300000]
    ids_divisor: 2
    chrom_job_threads: 4
    merge_job_threads: 4
  comparison_directions: [decondensation, condensation]
  cutoff_na: 0.5
  score_threshold: 5.0
  min_significant: 1
~~~

After SCRIPT_DIR is assigned in Snakefile, add:

~~~python
import sys

sys.path.insert(0, os.path.join(SCRIPT_DIR, "scripts"))
from cworld_is_config import parse_cworld_is_settings
CWORLD_IS = parse_cworld_is_settings(config)
~~~

- [ ] **Step 4: Run parser tests**

~~~bash
pytest -q test/cworld_is/test_cworld_is_config.py
~~~

Expected: 7 passed.

- [ ] **Step 5: Commit**

~~~bash
git add scripts/cworld_is_config.py test/cworld_is/test_cworld_is_config.py config/config.yaml Snakefile
git commit -m "feat: configure fixed 4kb cworld IS scales"
~~~

## Task 2: Reject a non-4 kb cooler before CWorld work begins

**Files:**

- Create: test/cworld_is/test_convert_cool_to_cworld.py
- Modify: scripts/convert_cool_to_cworld.py

**Interfaces:**

- Consumes: --cool, --chrom, --out, --log, and optional --expected-bin-size.
- Produces: a gzip CWorld matrix only after matching the requested bin size.
- Used by: compute_IS_chrom with --expected-bin-size 4000.

- [ ] **Step 1: Write failing bin-size and matrix-format tests**

Create test/cworld_is/test_convert_cool_to_cworld.py:

~~~python
from pathlib import Path
import gzip
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
import convert_cool_to_cworld as converter

def test_rejects_wrong_bin_size_before_loading_matrix(monkeypatch, tmp_path):
    fake = SimpleNamespace(chromnames=["HC04_A01"], info={"bin-size": 20000})
    monkeypatch.setattr(converter.cooler, "Cooler", lambda _: fake)

    with pytest.raises(ValueError, match="expected 4000 bp, found 20000 bp"):
        converter.convert_cool_to_cworld(
            "wrong.cool", "HC04_A01", tmp_path / "wrong.matrix.gz",
            expected_bin_size=4000,
        )

def test_writes_cworld_header_and_labels(monkeypatch, tmp_path):
    bins = pd.DataFrame(
        {"chrom": ["HC04_A01", "HC04_A01"], "start": [0, 4000], "end": [4000, 8000]}
    )

    class FakeMatrix:
        def fetch(self, chrom):
            assert chrom == "HC04_A01"
            return np.array([[1.0, np.nan], [2.0, 3.0]])

    fake = SimpleNamespace(
        chromnames=["HC04_A01"], info={"bin-size": 4000},
        bins=lambda: bins, matrix=lambda balance: FakeMatrix(),
    )
    monkeypatch.setattr(converter.cooler, "Cooler", lambda _: fake)
    output = tmp_path / "matrix.gz"

    converter.convert_cool_to_cworld(
        "ok.cool", "HC04_A01", output, expected_bin_size=4000
    )

    with gzip.open(output, "rt") as handle:
        lines = handle.read().splitlines()
    assert lines[0] == "\tHiC_1|aaa|HC04_A01:0-4000\tHiC_2|aaa|HC04_A01:4000-8000"
    assert lines[1] == "HiC_1|aaa|HC04_A01:0-4000\t1\t0"
    assert lines[2] == "HiC_2|aaa|HC04_A01:4000-8000\t2\t3"
~~~

- [ ] **Step 2: Confirm the new API fails first**

~~~bash
pytest -q test/cworld_is/test_convert_cool_to_cworld.py
~~~

Expected: failure because convert_cool_to_cworld does not accept expected_bin_size.

- [ ] **Step 3: Implement the guard**

Change the function signature and put the guard immediately after Cooler construction:

~~~python
def convert_cool_to_cworld(cool_path, chrom, output_path, log_file=None, expected_bin_size=None):
    log(f"[INFO] Opening cooler: {cool_path}", log_file)
    clr = cooler.Cooler(cool_path)
    resolution = int(clr.info["bin-size"])
    if expected_bin_size is not None and resolution != expected_bin_size:
        raise ValueError(
            f"Cooler {cool_path} has bin size {resolution} bp; "
            f"expected {expected_bin_size} bp"
        )
    if chrom not in clr.chromnames:
        raise ValueError(f"Chromosome {chrom} not found in {cool_path}")
~~~

Use that resolution variable for later logging, create output_path.parent before gzip.open, and add:

~~~python
parser.add_argument(
    "--expected-bin-size", type=int,
    help="Require this cooler bin size in bp before conversion",
)
try:
    convert_cool_to_cworld(
        args.cool, args.chrom, args.out, log_file=args.log,
        expected_bin_size=args.expected_bin_size,
    )
except ValueError as exc:
    parser.error(str(exc))
~~~

The historical validate_cworld_is smoke harness keeps working because it omits this optional argument. The production rule supplies 4000.

- [ ] **Step 4: Run focused tests**

~~~bash
pytest -q test/cworld_is/test_convert_cool_to_cworld.py test/cworld_is/test_compute_is_oom_script.py
~~~

Expected: all pass.

- [ ] **Step 5: Commit**

~~~bash
git add scripts/convert_cool_to_cworld.py test/cworld_is/test_convert_cool_to_cworld.py
git commit -m "feat: guard cworld conversion at 4kb"
~~~

## Task 3: Make IS conversion strictly manifest-driven

**Files:**

- Modify: scripts/convert_cworld_IS_to_meltron_IS.py
- Modify: test/cworld_is/test_convert_cworld_is.py

**Interfaces:**

- Consumes: --sample SAMPLE --chroms CHR... --manifests PATH... --outdir RESULTS --is-distances CSV.
- Produces: RESULTS/IS/SAMPLE/SAMPLE_IS.tsv only after all requested entries validate.
- Used by: merge_IS. The existing validate_cworld_is harness already creates one manifest and will pass it explicitly.

- [ ] **Step 1: Write failing manifest-only tests**

Change the existing write_bedgraph helper to return path after writing it, then replace the runner helper in test_convert_cworld_is.py:

~~~python
def write_bedgraph(root, sample, chrom, distance):
    path = (
        root / "cworld_insulation" / sample /
        f"{sample}__{chrom}--is{distance + 1}--nt0--ids{distance // 2 + 1}--ss1--immean.insulation.bedGraph"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"track type=bedGraph name=IS\n{chrom}\t0\t16000\t0.25\n"
    )
    return path

def write_manifest(root, sample, chrom, entries):
    path = root / "cworld_insulation" / sample / f"{sample}__{chrom}.is_manifest.tsv"
    path.write_text("".join(f"{distance}\t{filename}\n" for distance, filename in entries))
    return path

def run_converter(root, chroms, manifests, distances="32000"):
    return subprocess.run(
        [
            sys.executable, str(SCRIPT), "--sample", "S1", "--chroms", *chroms,
            "--manifests", *(str(path) for path in manifests),
            "--outdir", str(root), "--is-distances", distances,
            "--log", str(root / "converter.log"),
        ],
        text=True, capture_output=True, check=False,
    )
~~~

Keep the successful 32000 assertion, but use a manifest containing the generated is32001 filename. Add:

~~~python
def test_converter_rejects_duplicate_distance_in_manifest(tmp_path):
    bedgraph = write_bedgraph(tmp_path, "S1", "HC04_A01", 32000)
    manifest = write_manifest(
        tmp_path, "S1", "HC04_A01",
        [(32000, bedgraph.name), (32000, bedgraph.name)],
    )

    result = run_converter(tmp_path, ["HC04_A01"], [manifest])

    assert result.returncode != 0
    assert not (tmp_path / "IS" / "S1" / "S1_IS.tsv").exists()

def test_converter_rejects_missing_manifest_target(tmp_path):
    manifest = write_manifest(
        tmp_path, "S1", "HC04_A01", [(32000, "missing.bedGraph")]
    )

    result = run_converter(tmp_path, ["HC04_A01"], [manifest])

    assert result.returncode != 0
    assert "manifest target missing" in result.stderr
    assert not (tmp_path / "IS" / "S1" / "S1_IS.tsv").exists()
~~~

- [ ] **Step 2: Run tests to prove the CLI lacks manifests**

~~~bash
pytest -q test/cworld_is/test_convert_cworld_is.py
~~~

Expected: failure because --manifests is unrecognized.

- [ ] **Step 3: Implement strict manifest selection**

Remove glob imports and legacy filename matching. Add:

~~~python
def manifest_bedgraph(manifest_path, distance):
    matches = []
    for line in manifest_path.read_text().splitlines():
        requested, separator, filename = line.partition("\t")
        if separator and requested == str(distance):
            matches.append(filename)
    if len(matches) != 1:
        raise ValueError(
            f"{manifest_path}: distance {distance} requires exactly one entry; "
            f"found {len(matches)}"
        )
    root = manifest_path.parent.resolve()
    bedgraph = (root / matches[0]).resolve()
    if root not in bedgraph.parents or not bedgraph.is_file():
        raise ValueError(f"{manifest_path}: manifest target missing: {bedgraph}")
    return bedgraph

def convert_manifest(manifest_path, distances, log_file):
    frames = []
    for distance in distances:
        bedgraph = manifest_bedgraph(manifest_path, distance)
        log(f"[INFO] Converting {bedgraph} (requested dist={distance})", log_file)
        frames.append(convert_bedgraph_to_meltron(bedgraph, distance, log_file))
    return frames
~~~

Change both CLI options to required multi-value arguments: parser.add_argument("--chroms", nargs="+", required=True) and parser.add_argument("--manifests", nargs="+", required=True). Parse chroms and manifests as lists, reject unequal lengths, missing manifests, and empty frames. Iterate with zip(chrom_list, manifest_paths, strict=True), collect each convert_manifest result, sort by chromosome, write a temporary TSV, then use os.replace only after every chromosome and scale succeeds. IS_distance always comes from the requested manifest key; never parse it from CWorld output filenames.

Change scripts/validate_cworld_is.py so its converter command adds:

~~~python
"--manifests", str(manifest_path),
~~~

- [ ] **Step 4: Run converter and smoke assertion tests**

~~~bash
pytest -q test/cworld_is/test_convert_cworld_is.py test/cworld_is/test_compute_is_oom_script.py
~~~

Expected: all pass.

- [ ] **Step 5: Commit**

~~~bash
git add scripts/convert_cworld_IS_to_meltron_IS.py scripts/validate_cworld_is.py test/cworld_is/test_convert_cworld_is.py
git commit -m "feat: merge meltron IS from cworld manifests"
~~~

## Task 4: Schedule one CWorld job per sample × chromosome

**Files:**

- Modify: rules/compute_IS.smk
- Modify: test/cworld_is/test_compute_is_scheduling.py

**Interfaces:**

- Consumes: results/balanced/<sample>/<sample>__<chrom>_iced.cool, the extract_and_balance_sample checkpoint, and CWORLD_IS.
- Produces: one manifest per chromosome, ten matching bedGraphs, per-chromosome logs/benchmarks, and one per-sample final IS TSV.
- Preserves: rules/meltron.smk input paths and the four-column table schema.

- [ ] **Step 1: Replace the old topology test**

Replace test_compute_is_schedules_one_cworld_job_per_requested_sample with:

~~~python
def test_compute_is_schedules_one_cworld_job_per_sample_chromosome():
    result = subprocess.run(
        ["snakemake", "-n", "-p", "--configfile", "config/config.yaml", "--", *TARGETS],
        cwd=REPO, text=True, capture_output=True, check=False,
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
~~~

- [ ] **Step 2: Run test to demonstrate the old one-per-sample topology fails**

~~~bash
pytest -q test/cworld_is/test_compute_is_scheduling.py
~~~

Expected: failure shows compute_IS count 10 and no merge_IS rule.

- [ ] **Step 3: Replace rules/compute_IS.smk**

At module scope derive:

~~~python
CWORLD_IS_DISTANCES = tuple(str(distance) for distance in CWORLD_IS.is_distances_bp)
CWORLD_IS_DISTANCES_CSV = ",".join(CWORLD_IS_DISTANCES)
~~~

Define the resource-facing compute rule:

~~~python
rule compute_IS_chrom:
    input:
        upstream=lambda w: checkpoints.extract_and_balance_sample.get(
            sample=w.sample
        ).output.manifest,
        cool=OUTDIR + "/balanced/{sample}/{sample}__{chrom}_iced.cool"
    output:
        manifest=OUTDIR + "/cworld_insulation/{sample}/{sample}__{chrom}.is_manifest.tsv"
    log:
        "{}/logs/compute_IS/{{sample}}__{{chrom}}_compute_IS.log".format(OUTDIR)
    benchmark:
        "{}/benchmark/compute_IS/{{sample}}__{{chrom}}_compute_IS.tsv".format(OUTDIR)
    threads:
        CWORLD_IS.chrom_job_threads
    conda:
        "../env/cworld.yaml"
~~~

Its params contain OUTDIR, TMP, CWORLD_IS.source_resolution_bp, the exact ten distances, CWORLD_IS.ids_divisor, the two conversion scripts, and env/cworld-dekker. Its shell begins:

~~~bash
set -euo pipefail
SAMPLE="{wildcards.sample}"
CHR="{wildcards.chrom}"
OUTDIR="{params.outdir}"
TMP="{params.tmpdir}"
CWORLD_SAMPLE_DIR="$OUTDIR/cworld_insulation/$SAMPLE"
MATRIX_DIR="$TMP/cworld_matrix/$SAMPLE"
MATRIX="${MATRIX_DIR}/${SAMPLE}__${CHR}.matrix.gz"
PREFIX="${CWORLD_SAMPLE_DIR}/${SAMPLE}__${CHR}"
MANIFEST="$PREFIX.is_manifest.tsv"
MANIFEST_TMP="$MANIFEST.tmp.$$"
mkdir -p "$CWORLD_SAMPLE_DIR" "$MATRIX_DIR" "$OUTDIR/logs/compute_IS" "$OUTDIR/benchmark/compute_IS"
trap 'rm -f "$MANIFEST_TMP"' EXIT
: > "$MANIFEST_TMP"
python "{params.script_convert}" --cool "{input.cool}" --chrom "$CHR" --out "$MATRIX" \
  --expected-bin-size "{params.source_resolution_bp}" --log "{log}"
~~~

For each configured distance, create a unique mktemp marker in CWORLD_SAMPLE_DIR; call matrix2insulation.pl exactly once with --is "$DIST" and --ids "$((DIST / {params.ids_divisor}))"; find only files matching the exact sample/chromosome prefix and newer than the marker; fail unless that list has one bedGraph; append requested distance, a tab, and the bedGraph basename to MANIFEST_TMP; remove the marker on both success and error. After ten calls, require ten manifest lines and atomically move MANIFEST_TMP to MANIFEST. Retain MATRIX below tmp/cworld_matrix for diagnostics.

Define the merge rule:

~~~python
rule merge_IS:
    input:
        manifests=expand(
            OUTDIR + "/cworld_insulation/{{sample}}/{{sample}}__{chrom}.is_manifest.tsv",
            chrom=CHROMS,
        )
    output:
        is_file=OUTDIR + "/IS/{sample}/{sample}_IS.tsv"
    log:
        "{}/logs/compute_IS/{{sample}}_merge_IS.log".format(OUTDIR)
    benchmark:
        "{}/benchmark/compute_IS/{{sample}}_merge_IS.tsv".format(OUTDIR)
    threads:
        CWORLD_IS.merge_job_threads
    conda:
        "../env/cworld.yaml"
~~~

Its shell calls convert_cworld_IS_to_meltron_IS.py once with sample, all CHROMS, all input manifests, output directory, the CSV of ten scales, and log; then requires a nonempty output TSV. It never invokes matrix2insulation.pl.

- [ ] **Step 4: Run targeted implementation checks**

~~~bash
pytest -q test/cworld_is/test_cworld_is_config.py test/cworld_is/test_convert_cool_to_cworld.py test/cworld_is/test_convert_cworld_is.py test/cworld_is/test_compute_is_scheduling.py
CONDA_OFFLINE=true snakemake -n -p --use-conda --conda-prefix "$PWD/.snakemake/conda" --configfile config/config.yaml -- results/IS/L021/L021_IS.tsv
~~~

Expected: targeted tests pass; the direct dry run shows 26 compute_IS_chrom jobs, 1 merge_IS job, no compute_IS checkpoint, and no 20 kb or 100 kb CWorld calculation.

- [ ] **Step 5: Commit**

~~~bash
git add rules/compute_IS.smk test/cworld_is/test_compute_is_scheduling.py
git commit -m "feat: schedule cworld IS per chromosome"
~~~

## Task 5: Add safe explicit-sample cleanup

**Files:**

- Create: scripts/cleanup_cworld_is.py
- Create: test/cworld_is/test_cleanup_cworld_is.py

**Interfaces:**

- Consumes: --outdir RESULTS --samples SAMPLE [SAMPLE ...] with optional --execute.
- Produces: a printed, deterministic target list; filesystem changes only with --execute.
- Deletes only: results/cworld_insulation/<sample>/, results/IS/<sample>/, results/logs/compute_IS/<sample>*, results/benchmark/compute_IS/<sample>*.

- [ ] **Step 1: Write failing cleanup tests**

Create test/cworld_is/test_cleanup_cworld_is.py:

~~~python
from pathlib import Path
import subprocess
import sys

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "cleanup_cworld_is.py"

def populate(results, sample):
    (results / "cworld_insulation" / sample).mkdir(parents=True)
    (results / "IS" / sample).mkdir(parents=True)
    (results / "logs" / "compute_IS").mkdir(parents=True, exist_ok=True)
    (results / "benchmark" / "compute_IS").mkdir(parents=True, exist_ok=True)
    (results / "cworld_insulation" / sample / "a.bedGraph").write_text("x")
    (results / "IS" / sample / f"{sample}_IS.tsv").write_text("x")
    (results / "logs" / "compute_IS" / f"{sample}__A01.log").write_text("x")
    (results / "benchmark" / "compute_IS" / f"{sample}__A01.tsv").write_text("x")

def run_cleanup(results, *extra):
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--outdir", str(results), "--samples", "S1", *extra],
        text=True, capture_output=True, check=False,
    )

def test_cleanup_defaults_to_dry_run(tmp_path):
    results = tmp_path / "results"
    populate(results, "S1")

    result = run_cleanup(results)

    assert result.returncode == 0, result.stderr
    assert "[DRY-RUN]" in result.stdout
    assert (results / "cworld_insulation" / "S1").exists()
    assert (results / "IS" / "S1" / "S1_IS.tsv").exists()

def test_execute_removes_only_selected_scope(tmp_path):
    results = tmp_path / "results"
    populate(results, "S1")
    populate(results, "S2")
    protected = results / "balanced" / "S1" / "S1__HC04_A01_iced.cool"
    protected.parent.mkdir(parents=True)
    protected.write_text("must remain")

    result = run_cleanup(results, "--execute")

    assert result.returncode == 0, result.stderr
    assert not (results / "cworld_insulation" / "S1").exists()
    assert not (results / "IS" / "S1").exists()
    assert not list((results / "logs" / "compute_IS").glob("S1*"))
    assert not list((results / "benchmark" / "compute_IS").glob("S1*"))
    assert (results / "cworld_insulation" / "S2").exists()
    assert (results / "IS" / "S2").exists()
    assert protected.read_text() == "must remain"
~~~

- [ ] **Step 2: Confirm the helper is absent**

~~~bash
pytest -q test/cworld_is/test_cleanup_cworld_is.py
~~~

Expected: failure because scripts/cleanup_cworld_is.py is absent.

- [ ] **Step 3: Implement path-validated cleanup**

Create the script with argparse, pathlib.Path, glob, and shutil. Reject blank samples and names containing slash, backslash, dot, or dot-dot. Resolve outdir once. For each sample, construct only:

~~~python
targets = [
    root / "cworld_insulation" / sample,
    root / "IS" / sample,
    *sorted((root / "logs" / "compute_IS").glob(f"{sample}*")),
    *sorted((root / "benchmark" / "compute_IS").glob(f"{sample}*")),
]
~~~

Before deletion require target.resolve().is_relative_to(root) and require target is not root. Without --execute print "[DRY-RUN] would remove: <path>". With --execute, delete directories with shutil.rmtree and files/symlinks with unlink, printing "[DELETE] removed: <path>". Print "[SKIP] missing: <path>" for missing paths. No other paths are constructed or removed.

- [ ] **Step 4: Run cleanup test**

~~~bash
pytest -q test/cworld_is/test_cleanup_cworld_is.py
~~~

Expected: 2 passed; all fixtures reside under pytest temporary directories.

- [ ] **Step 5: Commit**

~~~bash
git add scripts/cleanup_cworld_is.py test/cworld_is/test_cleanup_cworld_is.py
git commit -m "feat: add safe cworld IS cleanup"
~~~

## Task 6: Synchronize documentation and complete local acceptance

**Files:**

- Modify: README.md
- Modify: docs/CODEBASE_OVERVIEW.md
- Modify: docs/SNAKEMAKE_PIPELINE_MAP.md

**Interfaces:**

- Consumes: the configuration and output contracts from Tasks 1-5.
- Produces: user-facing instructions that match runtime behavior exactly.

- [ ] **Step 1: Add README metrics and constraints**

Add a CWorld insulation scores for Meltron section after Run on HPC. It states all of the following verbatim facts:

- CWorld IS is solely a Meltron input and always uses a 4000 bp balanced cool file.
- The ten windows are 32000, 60000, 92000, 120000, 152000, 180000, 212000, 240000, 272000, and 300000 bp.
- --ids equals --is / 2; an output suffix such as is32001 does not change recorded distance 32000.
- A sample has 26 four-CPU compute_IS_chrom jobs; each computes ten scales sequentially; one four-CPU merge_IS job emits the final table.
- CPU count controls Donau memory allocation; change meltron.cworld_is.chrom_job_threads instead of setting a memory field.
- The matrix, bedGraph/manifest, and final TSV paths are tmp/cworld_matrix/<sample>/, results/cworld_insulation/<sample>/, and results/IS/<sample>/<sample>_IS.tsv.
- The final table columns are chrom, viewpoint, IS_distance, IS; one failed chromosome blocks only its own sample merge and downstream comparison.
- The cleanup command must first be shown without --execute, then repeated with --execute before dsub -s run_snakemake.sh.

Include these exact commands:

~~~bash
python scripts/cleanup_cworld_is.py --outdir results --samples H1394 H1395 H1396 H1397 H1398 H1406 H1409 H1410 H1436 L021
python scripts/cleanup_cworld_is.py --outdir results --samples H1394 H1395 H1396 H1397 H1398 H1406 H1409 H1410 H1436 L021 --execute
dsub -s run_snakemake.sh
~~~

- [ ] **Step 2: Update architecture maps**

In docs/CODEBASE_OVERVIEW.md, describe rules/compute_IS.smk as two rules: compute_IS_chrom creates ten bedGraphs and one manifest for one chromosome; merge_IS creates the per-sample table. In docs/SNAKEMAKE_PIPELINE_MAP.md, replace the old worker loop over all chromosomes with 26 chromosome jobs followed by one merge. Both entries state fixed 4 kb resolution and the ten-scale series.

- [ ] **Step 3: Run complete local validation**

~~~bash
python -m compileall -q scripts/cworld_is_config.py scripts/convert_cool_to_cworld.py scripts/convert_cworld_IS_to_meltron_IS.py scripts/cleanup_cworld_is.py
pytest -q test/cworld_is
CONDA_OFFLINE=true snakemake -n -p --use-conda --conda-prefix "$PWD/.snakemake/conda" --configfile config/config.yaml -- results/IS/L021/L021_IS.tsv
git diff --check
~~~

Expected: compile succeeds, all CWorld IS tests pass, the dry run contains exactly 26 compute_IS_chrom and 1 merge_IS task for L021, and git diff --check is silent.

- [ ] **Step 4: Controlled production restart after local gates pass**

First preview the exact sample scope:

~~~bash
python scripts/cleanup_cworld_is.py --outdir results --samples H1394 H1395 H1396 H1397 H1398 H1406 H1409 H1410 H1436 L021
~~~

Only after its output shows the four approved classes, execute:

~~~bash
python scripts/cleanup_cworld_is.py --outdir results --samples H1394 H1395 H1396 H1397 H1398 H1406 H1409 H1410 H1436 L021 --execute
dsub -s run_snakemake.sh
~~~

Expected: dsub exposes one four-CPU compute_IS_chrom task per sample × chromosome and starts merge_IS only after that sample has all 26 manifests. This step never removes Conda environments, balanced cool inputs, or unrelated result directories.

- [ ] **Step 5: Commit**

~~~bash
git add README.md docs/CODEBASE_OVERVIEW.md docs/SNAKEMAKE_PIPELINE_MAP.md
git commit -m "docs: describe 4kb per-chromosome meltron IS workflow"
~~~

## Plan Self-Review

| Approved requirement | Plan task |
| --- | --- |
| Fixed 4 kb source and exact ten historical scales | Task 1 |
| --ids = --is / 2 | Tasks 1 and 4 |
| Reject an actual non-4 kb input | Tasks 2 and 4 |
| One task per sample × chromosome; ten sequential scales | Task 4 |
| Separate merge with unchanged Meltron table path | Tasks 3 and 4 |
| Filename +1 never changes IS_distance | Task 3 |
| Matrix in tmp, bedGraphs/manifests in results | Task 4 and Task 6 |
| Four CPUs and CPU-driven memory | Tasks 1 and 6 |
| Dry-run-first cleanup of exactly four classes | Task 5 and Task 6 |
| README metrics and constraints | Task 6 |
| Offline Conda/dsub and Meltron interface retained | Tasks 2, 4, and 6 |

The names CworldISSettings, CWORLD_IS, compute_IS_chrom, merge_IS, and the final per-sample TSV path are used consistently. Every code task includes a failing test, implementation detail, verification command, expected result, and isolated commit.
