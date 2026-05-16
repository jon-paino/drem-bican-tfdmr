#!/usr/bin/env bash
set -euo pipefail

# --- Inputs ---
SUBSET=~/epigenome/bican/protein_coding/chip-atlas_top3_200bp.txt
DMRS=~/epigenome/bican/protein_coding/OPDC_hypo_dmrs.200-250bp.id.bed   # chr start end id (id is col 4)
EXPLIST=/u/project/ernst/ernst/IGVF/experimentList.tab.gz

# --- Work dir ---
WORKDIR=tf_cell_acc_dmrs
mkdir -p "$WORKDIR"

DMRS_SORTED="$WORKDIR/dmrs.sorted.bed"
bedtools sort -i "$DMRS" > "$DMRS_SORTED"

# 0) Extract accessions in subset (DRX/SRX/ERX etc.) for metadata filtering
ACC_LIST="$WORKDIR/acc_list.txt"
awk '
  NF==0 {next}
  $0 ~ /^#/ {next}
  {
    gsub(/^[ \t]+|[ \t]+$/, "", $0)
    n=split($0, a, "/")
    f=a[n]
    sub(/\.bed\.gz$/, "", f)
    print f
  }' "$SUBSET" | sort -u > "$ACC_LIST"

echo "[0/3] Found $(wc -l < "$ACC_LIST") peak datasets in subset" >&2

# 1) Build ACC → DMR overlaps
echo "[1/3] Intersecting peaks with DMRs..." >&2
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

# 2) Parse metadata → ACC_to_tf_cell.tsv (ACC, TF, CellType) restricted to ACC_LIST
echo "[2/3] Parsing experiment metadata..." >&2
ACC_META="$WORKDIR/ACC_to_tf_cell.tsv"

zcat "$EXPLIST" \
| awk -v keepfile="$ACC_LIST" '
  BEGIN{
    FS="\t"; OFS="\t"; IGNORECASE=1
    while ((getline line < keepfile) > 0) {
      gsub(/^[ \t]+|[ \t]+$/, "", line)
      if (line!="") keep[line]=1
    }
    close(keepfile)
  }
  function clean(x){ gsub(/^[ \t]+|[ \t]+$/,"",x); return x }
  function is_control(x){
    x=tolower(x)
    return (x ~ /(^|[ _-])(input|igg|control)([ _-]|$)/)
  }
  {
    acc=$1
    if (!(acc in keep)) next

    tf=clean($4)
    tissue=clean($5)
    cell6=clean($6)
    tail10=(NF>=10?$10:"")

    # If TF missing/looks like control, try to recover from tail metadata
    if (tf=="" || is_control(tf)) {
      n=split(tail10, kvs, /[ \t]*;[ \t]*|[ \t]+/)
      for (i=1; i<=n; i++) if (kvs[i] ~ /^[^=\t]+=[^=\t]*$/) {
        split(kvs[i], kv, "="); key=tolower(clean(kv[1])); val=clean(kv[2])
        if (key=="chip antibody" || key=="chip_antibody" || key=="antibody") tf=val
      }
    }

    # CellType priority: tail10 cell line/tissue -> source_name -> col6 -> col5
    cell=""
    if (tail10!=""){
      n=split(tail10, kvs2, /[ \t]*;[ \t]*|[ \t]+/)
      for (i=1; i<=n; i++) if (kvs2[i] ~ /^[^=\t]+=[^=\t]*$/) {
        split(kvs2[i], kv, "="); key=tolower(clean(kv[1])); val=clean(kv[2])
        if (cell=="" && key=="cell line/tissue") cell=val
        if (cell=="" && key=="source_name") cell=val
      }
    }
    if (cell=="") cell=cell6
    if (cell=="") cell=tissue
    if (cell=="") cell="NA"

    # Drop obvious non-TF / control targets
    if (tf=="" || is_control(tf)) next
    if (tolower(tf) ~ /^histone$/) next

    print acc, tf, cell
  }' \
| sort -u > "$ACC_META"

echo "    -> $(wc -l < "$ACC_META") ACCs mapped to TF+Cell" >&2

# 3) Join ACC→DMR with ACC→(TF,Cell) to write TF|Cell|ACC → DMR
echo "[3/3] Joining and writing outputs..." >&2

OUT1="$WORKDIR/TFcellACC_to_dmrs.tsv"
OUT3="$WORKDIR/TF__dmr_3col.tsv"
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
