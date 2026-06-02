#!/bin/bash
set -euo pipefail

ROOT_DIR="/share/home/xhhuang24/project/PanTE/background_pool"
OUTDIR="$ROOT_DIR/results"
WC_DIR="$OUTDIR/window_compare"
BG_OUT="$WC_DIR/HiC_window_background.tsv"

TMP_DIR="$ROOT_DIR/tmp/backfill_ncc_$$"
mkdir -p "$TMP_DIR"

echo "[INFO] Start NCC backfill and background rebuild"
echo "[INFO] window_compare dir: $WC_DIR"

shopt -s nullglob
for comp_dir in "$WC_DIR"/*; do
  [[ -d "$comp_dir" ]] || continue

  comp_name="$(basename "$comp_dir")"
  merged_file="$comp_dir/${comp_name}_window_pair_metrics.tsv"
  gauss_file="$comp_dir/${comp_name}_gaussian_metrics.tsv"

  [[ -f "$merged_file" ]] || continue
  [[ -f "$gauss_file" ]] || continue

  echo "[INFO] Processing $comp_name"

  out_file="$TMP_DIR/${comp_name}_window_pair_metrics.tsv"

  awk -F'\t' -v OFS='\t' '
    NR==FNR {
      if (FNR==1) {
        for (i=1; i<=NF; i++) {
          if ($i=="chrom") c=i
          else if ($i=="window_id") w=i
          else if ($i=="window_start") s=i
          else if ($i=="window_end") e=i
          else if ($i=="NCC") n=i
        }
        next
      }
      key = $c FS $w FS $s FS $e
      ncc[key] = (n>0 ? $n : "")
      next
    }
    FNR==1 {
      has_ncc = 0
      for (i=1; i<=NF; i++) {
        if ($i=="NCC") {
          has_ncc = 1
          ncc_col = i
        }
      }
      if (!has_ncc) {
        print $0, "NCC"
      } else {
        print $0
      }
      for (i=1; i<=NF; i++) {
        if ($i=="chrom") mc=i
        else if ($i=="window_id") mw=i
        else if ($i=="window_start") ms=i
        else if ($i=="window_end") me=i
      }
      next
    }
    {
      key = $mc FS $mw FS $ms FS $me
      val = ((key in ncc) ? ncc[key] : "")
      if (!has_ncc) {
        print $0, val
      } else {
        $ncc_col = val
        print $0
      }
    }
  ' "$gauss_file" "$merged_file" > "$out_file"

  cp "$out_file" "$merged_file"
done

# rebuild background file from updated per-comparison merged files
{
  first=1
  for f in "$WC_DIR"/*/*_window_pair_metrics.tsv; do
    [[ -f "$f" ]] || continue
    if [[ $first -eq 1 ]]; then
      awk -F'\t' 'NR==1{print;next}1' "$f"
      first=0
    else
      awk -F'\t' 'NR==1{next}1' "$f"
    fi
  done
} > "$BG_OUT"

echo "[INFO] Done. Output: $BG_OUT"
rm -rf "$TMP_DIR"
