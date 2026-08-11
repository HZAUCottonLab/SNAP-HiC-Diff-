# Compute insulation scores with cworld-dekker and convert them to Meltron IS.
# CWorld runs once per sample and chromosome; the per-sample merge consumes
# exactly one manifest from each chromosome. Matrices remain under tmp/cworld_matrix/<sample>.

CWORLD_IS_DISTANCES = tuple(str(distance) for distance in CWORLD_IS.is_distances_bp)
CWORLD_IS_DISTANCES_CSV = ",".join(CWORLD_IS_DISTANCES)


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
    params:
        outdir=OUTDIR,
        tmpdir=TMP,
        source_resolution_bp=CWORLD_IS.source_resolution_bp,
        distances=" ".join(CWORLD_IS_DISTANCES),
        ids_divisor=CWORLD_IS.ids_divisor,
        script_convert=SCRIPT_DIR + "/scripts/convert_cool_to_cworld.py",
        script_merge=SCRIPT_DIR + "/scripts/convert_cworld_IS_to_meltron_IS.py",
        cworld_dir=SCRIPT_DIR + "/env/cworld-dekker",
    conda:
        "../env/cworld.yaml"
    shell:
        r"""
        set -euo pipefail

        SAMPLE="{wildcards.sample}"
        CHR="{wildcards.chrom}"
        OUTDIR="{params.outdir}"
        TMP="{params.tmpdir}"
        CWORLD_SAMPLE_DIR="$OUTDIR/cworld_insulation/$SAMPLE"
        MATRIX_DIR="$TMP/cworld_matrix/$SAMPLE"
        MATRIX="$MATRIX_DIR/$SAMPLE"__"$CHR.matrix.gz"
        PREFIX="$CWORLD_SAMPLE_DIR/$SAMPLE"__"$CHR"
        MANIFEST="$PREFIX.is_manifest.tsv"
        MANIFEST_TMP="$MANIFEST.tmp.$$"
        MARKER=""
        mkdir -p "$CWORLD_SAMPLE_DIR" "$MATRIX_DIR" "$OUTDIR/logs/compute_IS" "$OUTDIR/benchmark/compute_IS"
        trap 'rm -f "$MANIFEST_TMP"; if [ -n "$MARKER" ]; then rm -f "$MARKER"; fi' EXIT
        : > "$MANIFEST_TMP"

        python "{params.script_convert}" --cool "{input.cool}" --chrom "$CHR" --out "$MATRIX" \
          "--expected-bin-size" "{params.source_resolution_bp}" --log "{log}"

        for DIST in {params.distances}; do
            MARKER="$(mktemp "$CWORLD_SAMPLE_DIR/$SAMPLE"__"$CHR.is"$DIST".XXXXXX")"
            if ! PERL5LIB="{params.cworld_dir}/lib" perl -I "{params.cworld_dir}/lib" \
                "{params.cworld_dir}/scripts/perl/matrix2insulation.pl" \
                --is "$DIST" --ids "$((DIST / {params.ids_divisor}))" \
                -i "$MATRIX" -o "$PREFIX" >> "{log}" 2>&1; then
                rm -f "$MARKER"
                MARKER=""
                exit 1
            fi
            FILENAME_PREFIX="$SAMPLE"__"$CHR"
            mapfile -t BEDGRAPHS < <(
                find "$CWORLD_SAMPLE_DIR" -maxdepth 1 -type f \
                    -name "$FILENAME_PREFIX--is*--*.insulation.bedGraph" \
                    -newer "$MARKER" -print
            )
            if [ "${{#BEDGRAPHS[@]}}" -ne 1 ]; then
                echo "[ERROR] Expected one new CWorld bedGraph for $SAMPLE $CHR distance $DIST, found ${{#BEDGRAPHS[@]}}" >> "{log}"
                exit 1
            fi
            BEDGRAPH="${{BEDGRAPHS[0]}}"
            printf '%s\t%s\n' "$DIST" "${{BEDGRAPH##*/}}" >> "$MANIFEST_TMP"
            rm -f "$MARKER"
            MARKER=""
        done

        test "$(wc -l < "$MANIFEST_TMP")" -eq 10
        mv "$MANIFEST_TMP" "$MANIFEST"
        """


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
    params:
        outdir=OUTDIR,
        chroms=" ".join(CHROMS),
        distances_csv=CWORLD_IS_DISTANCES_CSV,
        expected_bin_size=CWORLD_IS.source_resolution_bp,
        script_merge=SCRIPT_DIR + "/scripts/convert_cworld_IS_to_meltron_IS.py",
    conda:
        "../env/cworld.yaml"
    shell:
        r"""
        set -euo pipefail

        SAMPLE="{wildcards.sample}"
        OUTDIR="{params.outdir}"
        mkdir -p "$OUTDIR/IS/$SAMPLE" "$OUTDIR/logs/compute_IS" "$OUTDIR/benchmark/compute_IS"

        python "{params.script_merge}" \
            --sample "$SAMPLE" \
            --chroms {params.chroms} \
            --manifests {input.manifests} \
            --outdir "$OUTDIR" \
            --is-distances "{params.distances_csv}" \
            --expected-bin-size "{params.expected_bin_size}" \
            --log "{log}" \
            >> "{log}" 2>&1

        test -s "{output.is_file}"
        """
