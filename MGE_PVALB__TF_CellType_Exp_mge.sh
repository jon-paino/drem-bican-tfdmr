#!/usr/bin/env bash
set -euo pipefail

# --- Inputs ---
SUBSET=~/epigenome/bican/protein_coding/chip-atlas_top3_200bp.txt
DMRS=~/epigenome/bican/protein_coding/MGE_P_all_dms2_hypo_merged.bed # chr start end id (id is col 4)

# --- Work dir ---
WORKDIR=tf_cell_acc_dmrs
mkdir -p "$WORKDIR"

DMRS_SORTED="$WORKDIR/dmrs.sorted.bed"
bedtools sort -i "$DMRS" > "$DMRS_SORTED"

# Pre-existing files assumed available:
#   $WORKDIR/acc_list.txt
#   $WORKDIR/ACC_to_tf_cell.tsv
ACC_LIST="$WORKDIR/acc_list.txt"
ACC_META="$WORKDIR/ACC_to_tf_cell.tsv"

echo "[0/2] Found $(wc -l < "$ACC_LIST") peak datasets in subset" >&2

# 1) Build ACC → DMR overlaps
echo "[1/2] Intersecting peaks with DMRs..." >&2
> "$WORKDIR/acc_dmr_pairs.tsv"

while IFS= read -r bedpath; do
  [[ -z "${bedpath:-}" || "${bedpath:0:1}" == "#" ]] && continue
  [[ -f "$bedpath" ]] || { echo "WARNING: missing peak bed: $bedpath" >&2; continue; }

  acc=$(basename "$bedpath" .bed.gz)  # e.g. DRX013181 or SRX2194270
  echo "    [+] $acc" >&2

  # Peaks are BED; DMRs are BED4 with ID in column 4.
  # Output: acc <tab> dmr_id
  zcat -f "$bedpath" \
    | bedtools intersect -a stdin -b "$DMRS_SORTED" -wa -wb \
    | awk -v ACC="$acc" 'BEGIN{OFS="\t"} {print ACC, $NF}' \
    >> "$WORKDIR/acc_dmr_pairs.tsv"

done < "$SUBSET"

sort -u -k1,1 -k2,2 "$WORKDIR/acc_dmr_pairs.tsv" -o "$WORKDIR/acc_dmr_pairs.tsv"
echo "    -> $(wc -l < "$WORKDIR/acc_dmr_pairs.tsv") unique ACC-DMR pairs" >&2

# 2) Using pre-existing metadata
echo "[2/2] Using pre-existing metadata at $ACC_META" >&2
echo "    -> $(wc -l < "$ACC_META") ACCs mapped to TF+Cell" >&2

# 3) Join ACC→DMR with ACC→(TF,Cell) to write TF|Cell|ACC → DMR
echo "[3/3] Joining and writing outputs..." >&2

OUT1="$WORKDIR/TFcellACC_to_dmrs_mge.tsv"
OUT3="$WORKDIR/TF__dmr_3col_mge.tsv"
: > "$OUT1"
: > "$OUT3"

join -t $'\t' -1 1 -2 1 \
  <(sort -k1,1 "$WORKDIR/acc_dmr_pairs.tsv") \
  <(sort -k1,1 "$ACC_META") \
| awk 'BEGIN{OFS="\t"; IGNORECASE=1}
       function is_control(x){x=tolower(x); return (x ~ /(^|[ _-])(input|igg|control)([ _-]|$)/)}
       {
         acc=$1; dmr=$2; tf=$3; cell=$4;
         if (tf=="" || is_control(tf)) next
         key=tf "|" cell "|" acc
         print key, dmr, 1 >> "'"$OUT1"'"
         print tf, dmr, 1 >> "'"$OUT3"'"
       }'

echo "Done."
echo "Outputs:"
echo "  $OUT1   (TF|Cell|ACC -> DMR_ID, DREM-ready)"
echo "  $OUT3   (TF -> DMR_ID, 3-col)"
