#!/usr/bin/env bash
# Build and verify the Conda prefixes required by SNAP-HiC-Diff.
# Create one prefix at a time with Mamba: a monolithic Snakemake environment
# bootstrap can exceed the login node's memory limit on linux-aarch64.
set -euo pipefail

ACTION="${1:-}"
ROOT="${SNAP_HIC_DIFF_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
CONDA_PREFIX_DIR="${CONDA_PREFIX_DIR:-$ROOT/.snakemake/conda}"
READY_MARKER="$CONDA_PREFIX_DIR/offline.ready.tsv"
MAMBA_BIN="${MAMBA_BIN:-mamba}"
SPECS=(env/SNAP-HiC-Diff.yaml env/cworld.yaml env/meltron.yaml)

die() {
    printf '[ERROR] %s\n' "$*" >&2
    printf 'Run: bash scripts/conda_envs.sh bootstrap\n' >&2
    exit 1
}

spec_sha256() {
    sha256sum "$ROOT/$1" | awk '{print $1}'
}

environment_hash() {
    local spec="$1"

    {
        printf '%s' "$(realpath "$CONDA_PREFIX_DIR")"
        cat "$ROOT/$spec"
    } | md5sum | awk '{print $1}'
}

desired_prefix() {
    printf '%s/%s_\n' "$CONDA_PREFIX_DIR" "$(environment_hash "$1")"
}

write_snakemake_cache_metadata() {
    local spec="$1"
    local hash cache_yaml setup_marker

    hash="$(environment_hash "$spec")"
    cache_yaml="$CONDA_PREFIX_DIR/${hash}_.yaml"
    setup_marker="$CONDA_PREFIX_DIR/${hash}_.env_setup_done"

    if [[ ! -f "$cache_yaml" ]] || ! cmp -s "$ROOT/$spec" "$cache_yaml"; then
        cp "$ROOT/$spec" "$cache_yaml"
    fi
    : > "$setup_marker"
}

find_completed_prefix() {
    local spec="$1"
    local hash candidate prefix

    hash="$(environment_hash "$spec")"
    for candidate in "${hash:0:8}" "$hash" "${hash}_"; do
        prefix="$CONDA_PREFIX_DIR/$candidate"
        if [[ -s "$prefix/conda-meta/history" ]]; then
            printf '%s\n' "$prefix"
            return 0
        fi
    done

    return 1
}

validate_python_package() {
    local python="$1"
    local package="$2"

    "$python" -c 'from importlib.metadata import version; import sys; print(version(sys.argv[1]))' "$package"
}

validate_prefix() {
    local spec="$1"
    local prefix="$2"

    case "$spec" in
        env/SNAP-HiC-Diff.yaml)
            local snap_prefix="$prefix"
            [[ -x "$snap_prefix/bin/python" ]] || die "Missing Python in $snap_prefix"
            [[ -x "$snap_prefix/bin/hictk" ]] || die "Missing hictk CLI in $snap_prefix"
            validate_python_package "$snap_prefix/bin/python" cooler
            validate_python_package "$snap_prefix/bin/python" cooltools
            validate_python_package "$snap_prefix/bin/python" numpy
            validate_python_package "$snap_prefix/bin/python" pandas
            "$snap_prefix/bin/hictk" --version
            ;;
        env/cworld.yaml)
            local cworld_prefix="$prefix"
            [[ -x "$cworld_prefix/bin/Rscript" ]] || die "Missing Rscript in $cworld_prefix"
            "$cworld_prefix/bin/Rscript" --version
            [[ -x "$cworld_prefix/bin/python" ]] || die "Missing Python in $cworld_prefix"
            [[ -x "$cworld_prefix/bin/perl" ]] || die "Missing Perl in $cworld_prefix"
            [[ -x "$cworld_prefix/bin/bedtools" ]] || die "Missing bedtools CLI in $cworld_prefix"
            validate_python_package "$cworld_prefix/bin/python" cooler
            validate_python_package "$cworld_prefix/bin/python" numpy
            validate_python_package "$cworld_prefix/bin/python" pandas
            "$cworld_prefix/bin/perl" -I "$ROOT/env/cworld-dekker/lib" -Mcworld::dekker -e 'print $cworld::dekker::VERSION, "\n"'
            "$cworld_prefix/bin/bedtools" --version
            ;;
        env/meltron.yaml)
            local meltron_prefix="$prefix"
            [[ -x "$meltron_prefix/bin/Rscript" ]] || die "Missing Rscript in $meltron_prefix"
            "$meltron_prefix/bin/Rscript" --version
            [[ -f "$meltron_prefix/lib/R/library/argparse/DESCRIPTION" ]] || die "Missing R package argparse in $meltron_prefix"
            [[ -f "$meltron_prefix/lib/R/library/data.table/DESCRIPTION" ]] || die "Missing R package data.table in $meltron_prefix"
            [[ -f "$meltron_prefix/lib/R/library/DescTools/DESCRIPTION" ]] || die "Missing R package DescTools in $meltron_prefix"
            ;;
        *)
            die "Unknown environment specification: $spec"
            ;;
    esac
}

write_ready_marker() {
    local marker_tmp="$READY_MARKER.tmp.$$"
    local spec prefix

    : > "$marker_tmp"
    for spec in "${SPECS[@]}"; do
        prefix="$(find_completed_prefix "$spec")" || die "No completed Conda prefix for $spec"
        validate_prefix "$spec" "$prefix"
        printf '%s\t%s\t%s\n' "$spec" "$(spec_sha256 "$spec")" "$prefix" >> "$marker_tmp"
    done
    mv "$marker_tmp" "$READY_MARKER"
}

bootstrap() {
    mkdir -p "$CONDA_PREFIX_DIR"
    rm -f "$READY_MARKER"

    if ! command -v "$MAMBA_BIN" >/dev/null 2>&1; then
        die "Install Mamba or set MAMBA_BIN to a usable executable."
    fi

    local spec prefix backup
    for spec in "${SPECS[@]}"; do
        prefix="$(desired_prefix "$spec")"
        if [[ -s "$prefix/conda-meta/history" ]]; then
            printf '[INFO] Reusing completed Conda prefix for %s: %s\n' "$spec" "$prefix"
        else
            if [[ -e "$prefix" ]]; then
                backup="${prefix}.incomplete.$$"
                printf '[INFO] Preserving incomplete prefix as %s\n' "$backup" >&2
                mv "$prefix" "$backup"
            fi
            printf '[INFO] Creating %s with Mamba: %s\n' "$spec" "$prefix" >&2
            "$MAMBA_BIN" env create --yes --prefix "$prefix" --file "$ROOT/$spec"
        fi
        write_snakemake_cache_metadata "$spec"
    done

    write_ready_marker
    printf '[INFO] Wrote %s\n' "$READY_MARKER"
}

verify() {
    local spec hash prefix extra actual_hash matches

    [[ -f "$READY_MARKER" ]] || die "Missing ready marker: $READY_MARKER"
    while IFS=$'\t' read -r spec hash prefix extra; do
        case "$spec" in
            env/SNAP-HiC-Diff.yaml|env/cworld.yaml|env/meltron.yaml) ;;
            *) die "Unexpected ready-marker specification: $spec" ;;
        esac
    done < "$READY_MARKER"

    for spec in "${SPECS[@]}"; do
        matches="$(awk -F '\t' -v expected="$spec" '$1 == expected {count += 1; line = $0} END {if (count == 1) print line; else exit 1}' "$READY_MARKER")" || die "Ready marker must contain exactly one record for $spec"
        [[ -n "$matches" ]] || die "Ready marker lacks $spec"
        IFS=$'\t' read -r marker_spec hash prefix extra <<< "$matches"
        [[ "$marker_spec" == "$spec" ]] || die "Ready marker specification mismatch: $marker_spec"
        [[ -n "$hash" && -n "$prefix" && -z "${extra:-}" ]] || die "Malformed ready marker: $READY_MARKER"
        actual_hash="$(spec_sha256 "$spec")"
        [[ "$hash" == "$actual_hash" ]] || die "Environment YAML changed: $spec"
        [[ -s "$prefix/conda-meta/history" ]] || die "Incomplete Conda prefix: $prefix"
        validate_prefix "$spec" "$prefix"
        printf '[INFO] Verified %s -> %s\n' "$spec" "$prefix"
    done
}

case "$ACTION" in
    bootstrap) bootstrap ;;
    verify) verify ;;
    *)
        printf 'Usage: %s {bootstrap|verify}\n' "$0" >&2
        exit 2
        ;;
esac
