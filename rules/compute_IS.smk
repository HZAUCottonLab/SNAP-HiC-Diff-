# Compute insulation scores with cworld-dekker and convert them to Meltron IS.
# One checkpoint job per sample loops over every chromosome and configured IS distance.

checkpoint compute_IS:
    input:
        lambda w: checkpoints.extract_and_balance_sample.get(
            sample=w.sample
        ).output.manifest
    output:
        touch = OUTDIR + "/logs/compute_IS/{sample}.done",
        is_file = OUTDIR + "/IS/{sample}/{sample}_IS.tsv"
    log:
        "{}/logs/compute_IS/{{sample}}_compute_IS.log".format(OUTDIR)
    benchmark:
        "{}/benchmark/compute_IS/{{sample}}_compute_IS.tsv".format(OUTDIR)
    threads:
        config["compute_is"]["threads"]
    params:
        chroms = " ".join(CHROMS),
        is_dists_str = " ".join(str(d) for d in config["meltron"]["is_distances"]),
        is_dists_csv = ",".join(str(d) for d in config["meltron"]["is_distances"]),
        outdir = OUTDIR,
        script_convert = SCRIPT_DIR + "/scripts/convert_cool_to_cworld.py",
        script_merge = SCRIPT_DIR + "/scripts/convert_cworld_IS_to_meltron_IS.py",
        cworld_dir = SCRIPT_DIR + "/env/cworld-dekker",
    conda:
        "../env/cworld.yaml"
    shell:
        r"""
        set -euo pipefail

        SAMPLE="{wildcards.sample}"
        OUTDIR="{params.outdir}"

        mkdir -p "${{OUTDIR}}/cworld_insulation/${{SAMPLE}}" \
                 "${{OUTDIR}}/IS/${{SAMPLE}}" \
                 "${{OUTDIR}}/logs/compute_IS"

        for CHR in {params.chroms}; do
            COOL="${{OUTDIR}}/balanced/${{SAMPLE}}/${{SAMPLE}}__${{CHR}}_iced.cool"
            MATRIX="${{OUTDIR}}/cworld_insulation/${{SAMPLE}}/${{SAMPLE}}__${{CHR}}.matrix.gz"
            PREFIX="${{OUTDIR}}/cworld_insulation/${{SAMPLE}}/${{SAMPLE}}__${{CHR}}"

            python "{params.script_convert}" --cool "$COOL" --chrom "$CHR" --out "$MATRIX" --log "{log}"

            for DIST in {params.is_dists_str}; do
                PERL5LIB="{params.cworld_dir}/lib" perl -I "{params.cworld_dir}/lib" \
                    "{params.cworld_dir}/scripts/perl/matrix2insulation.pl" \
                    --is "$DIST" --ids "$((DIST / 2))" -i "$MATRIX" -o "$PREFIX" >> "{log}" 2>&1
                compgen -G "${{PREFIX}}--is${{DIST}}--*.insulation.bedGraph" > /dev/null || {{
                    echo "[ERROR] Missing CWorld output for ${{SAMPLE}} ${{CHR}} distance ${{DIST}}" >> "{log}"
                    exit 1
                }}
            done
        done

        python "{params.script_merge}" \
            --sample "$SAMPLE" \
            --chroms "{params.chroms}" \
            --outdir "{params.outdir}" \
            --is-distances "{params.is_dists_csv}" \
            --log "{log}" \
            >> "{log}" 2>&1

        test -s "{output.is_file}"
        touch "{output.touch}"
        """
