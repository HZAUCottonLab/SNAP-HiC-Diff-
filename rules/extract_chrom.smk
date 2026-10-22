# Extract per-chromosome intra-chromosomal interactions from mcool as .cool
# One job per sample — shell block for-loops over all chromosomes internally.
# Pipeline entry point (no upstream dependencies).
#
# Outputs:
#   - balanced/{sample}/{sample}__{chrom}_iced.cool
#   - raw/{sample}/{sample}__{chrom}.raw.cool

checkpoint extract_and_balance_sample:
    input:
        mcool = lambda w: get_sample_mcool(w.sample),
        manifest_check = "data/manifest_mcool.tsv",
    output:
        manifest = OUTDIR + "/logs/extract_balance/{sample}.done",
        cool_files = expand(
            OUTDIR + "/balanced/{{sample}}/{{sample}}__{chrom}_iced.cool",
            chrom=CHROMS
        ),
        raw_files = expand(
            OUTDIR + "/raw/{{sample}}/{{sample}}__{chrom}.raw.cool",
            chrom=CHROMS
        ),
    log:
        OUTDIR + "/logs/extract_balance/{sample}.log"
    benchmark:
        "{}/benchmark/extract_balance/{{sample}}.tsv".format(OUTDIR)
    threads:
        config["hictk"]["threads"]
    params:
        resolution = RESOLUTION,
        tmpdir = TMP,
        outdir = OUTDIR,
        chroms = " ".join(CHROMS)
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        SAMPLE="{wildcards.sample}"
        INPUT_PATH="{input.mcool}"
        RESOLUTION="{params.resolution}"

        # .cool files are single-resolution: use direct path
        # .mcool files need /resolutions/$RESOLUTION suffix
        if [[ "$INPUT_PATH" == *.cool ]]; then
            URI="$INPUT_PATH"
        else
            URI="$INPUT_PATH::/resolutions/$RESOLUTION"
        fi

        mkdir -p "{params.outdir}/balanced/${{SAMPLE}}" \
                 "{params.outdir}/raw/${{SAMPLE}}" \
                 "{params.outdir}/logs/extract_balance" \
                 "{params.tmpdir}"

        echo "[INFO] Input: $INPUT_PATH"
        echo "[INFO] URI: $URI"
        echo "[INFO] Sample: $SAMPLE"

        for CHR in {params.chroms}; do
            echo
            echo "[INFO] Processing chromosome: $CHR"

            BASE="${{SAMPLE}}__${{CHR}}"
            RAW="{params.tmpdir}/${{BASE}}.raw.cool"
            OUT="{params.outdir}/balanced/${{SAMPLE}}/${{BASE}}_iced.cool"
            CHROMSIZES=$(mktemp "{params.tmpdir}/${{BASE}}.XXXXXX.chrom.sizes")

            hictk dump "$URI" \
                --table chroms \
                | awk -v chr="$CHR" '$1 == chr' \
                > "$CHROMSIZES"

            if [ ! -s "$CHROMSIZES" ]; then
                echo "[ERROR] Chromosome $CHR not found in $URI"
                exit 1
            fi

            mkdir -p "{params.outdir}/raw/${{SAMPLE}}"

            hictk dump "$URI" \
                --range "$CHR" \
                --range2 "$CHR" \
                --join \
                --sorted \
            | hictk load - "$RAW" \
                --chrom-sizes "$CHROMSIZES" \
                --bin-size "$RESOLUTION" \
                --format bg2 \
                --assume-sorted \
                --threads {threads} \
                --tmpdir "{params.tmpdir}" \
                --force

            cp "$RAW" "{params.outdir}/raw/${{SAMPLE}}/${{BASE}}.raw.cool"

            cp "$RAW" "$OUT"

            hictk balance ice "$OUT" \
                --tmpdir "{params.tmpdir}" \
                --threads {threads} \
                --name weight \
                --force \
                --no-create-weight-link

            hictk validate "$OUT"

            rm -f "$RAW" "$CHROMSIZES"

            echo "[INFO] Finished chromosome: $CHR"
        done

        echo
        echo "[INFO] Finished sample: $SAMPLE"

        touch "{output.manifest}"
        """
