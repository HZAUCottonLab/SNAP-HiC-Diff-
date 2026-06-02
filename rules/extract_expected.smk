# Calculate distance-stratified expected values from balanced per-chromosome .cool
# Step 1: cooltools expected-cis -> tsv (per sample, loops over all chroms)
# Step 2: tsv -> npz (per sample, calls extract_expected_batch.py)
#
# Inputs: results/balanced/{sample}/{sample}__{chrom}_iced.cool
# Outputs: results/expected/{sample}/{sample}__{chrom}_cooltools.tsv, _expected.npz

# ---------------------------------------------------------------------------
# Step 1: cooltools expected-cis -> tsv (parent: one job per sample, loops chroms)
# Depends on extract_and_balance_sample checkpoint (waits for all chrom cool files)
# ---------------------------------------------------------------------------

rule extract_expected_tsv:
    input:
        lambda wildcards: [
            "{}/logs/extract_balance/{}.done".format(OUTDIR, wildcards.sample)
        ]
    output:
        manifest = "{}/logs/expected/{{sample}}_cooltools.done".format(OUTDIR)
    log:
        "{}/logs/expected/{{sample}}_cooltools.log".format(OUTDIR)
    benchmark:
        "{}/benchmark/expected/{{sample}}_cooltools.tsv".format(OUTDIR)
    threads:
        config["hictk"]["threads"]
    params:
        nproc = config["hictk"]["threads"],
        chroms = " ".join(CHROMS),
        outdir = OUTDIR
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        SAMPLE="{wildcards.sample}"

        mkdir -p "{params.outdir}/expected/${{SAMPLE}}" "{params.outdir}/logs/expected"

        echo "[INFO] Sample: $SAMPLE"
        echo "[INFO] Chromosomes: {params.chroms}"

        for CHR in {params.chroms}; do
            echo
            echo "[INFO] Processing chromosome: $CHR"

            BASE="${{SAMPLE}}__${{CHR}}"
            COOL="{params.outdir}/balanced/${{SAMPLE}}/${{BASE}}_iced.cool"
            TSV="{params.outdir}/expected/${{SAMPLE}}/${{BASE}}_cooltools.tsv"

            if [ ! -f "$COOL" ]; then
                echo "[ERROR] Input not found: $COOL"
                exit 1
            fi

            cooltools expected-cis "$COOL" \
                --clr-weight-name weight \
                --ignore-diags 2 \
                --nproc {params.nproc} \
                -o "$TSV" \
            2>&1 | tee -a "{log}"

            echo "[INFO] Finished chromosome: $CHR"
        done

        echo
        echo "[INFO] Finished sample: $SAMPLE"

        touch "{output.manifest}"
        """


# ---------------------------------------------------------------------------
# Step 2: tsv -> npz (parent: one job per sample, loops chroms via Python script)
# ---------------------------------------------------------------------------

rule extract_expected_npz:
    input:
        lambda wildcards: [
            "{}/logs/expected/{}_cooltools.done".format(OUTDIR, wildcards.sample)
        ],
        script = SCRIPT_DIR + "/scripts/extract_expected_batch.py"
    output:
        manifest = "{}/logs/expected/{{sample}}_to_npz.done".format(OUTDIR)
    log:
        "{}/logs/expected/{{sample}}_to_npz.log".format(OUTDIR)
    benchmark:
        "{}/benchmark/expected/{{sample}}_to_npz.tsv".format(OUTDIR)
    threads:
        config["parallel"]["threads_per_job"]
    params:
        max_dist = 10000,
        chroms = CHROMS,
        outdir = OUTDIR
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        SAMPLE="{wildcards.sample}"

        mkdir -p "{params.outdir}/expected/${{SAMPLE}}" "{params.outdir}/logs/expected"

        INPUTS=""
        OUTPUTS=""
        for chrom in {params.chroms}; do
            BASE="${{SAMPLE}}__${{chrom}}"
            INPUTS="$INPUTS {params.outdir}/expected/${{SAMPLE}}/${{BASE}}_cooltools.tsv"
            OUTPUTS="$OUTPUTS {params.outdir}/expected/${{SAMPLE}}/${{BASE}}_expected.npz"
        done

        python {input.script} \
            --chroms {params.chroms} \
            --outdir {params.outdir} \
            --sample "$SAMPLE" \
            --max-dist {params.max_dist} \
            --inputs $INPUTS \
            --outputs $OUTPUTS \
            --log "{log}" \
            > "{log}" 2>&1

        touch "{output.manifest}"
        """
