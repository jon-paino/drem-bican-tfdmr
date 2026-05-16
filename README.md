# DREM BICAN — final DMR pipeline scripts

Code accompanying the capstone report *Extending DREM to Methylation
Dynamics: A TF–DMR Prior Framework for Brain Cell Lineages* (J. Paino).
This repository contains **scripts only**; the BICAN-derived input and
output data are not redistributed here.

Scripts used to construct the data and outputs for the BICAN DMR-substrate
DREM work (TF–DMR priors and per-lineage DMR β-matrices). Copied here from
their original locations under `epigenome/`; filenames are lineage-prefixed
(`<LINEAGE>__<original_name>`) so the set is collision-free in one flat
directory (two `TF_CellType_Exp.sh` and a generic
`combine_filtered_methylation.sh` would otherwise clobber).

## Pipeline order

1. **β-matrix per lineage:** `*__apm_*_DMR.sh`  →  `*__combine_filtered_methylation_*.sh`
2. **TF–DMR prior:** `regenerate_tf_dmr.sh` (driven by `run_all_lineages.sh`)
3. **BrainSpan TF filter:** `brainspan_filter_tf_prior.py` (mean RPKM > 1.0)
4. (DREM itself is run separately from the Java tool on the β-matrix +
   the prior; not a script in this set.)

## TF–DMR prior construction (canonical / final, all five lineages)

| Script | Origin | Purpose |
|---|---|---|
| `regenerate_tf_dmr.sh` | `IGVF/tf_dmr_alltf_new_v2/` | Final prior builder. For one lineage (`mge\|bach2\|epha4\|opc\|ecge`): sorts the merged lineage DMR BED, intersects ChIP-Atlas top-1000 peak BEDs with it, dedups (accession, DMR) pairs, joins to the curated (TF, cell) label table, and emits the annotated `TF\|cell\|accession` prior and a collapsed TF-only prior. |
| `run_all_lineages.sh` | `IGVF/tf_dmr_alltf_new_v2/` | Driver that calls `regenerate_tf_dmr.sh` sequentially for all five lineages. Expects `regenerate_tf_dmr.sh` adjacent (it is, here). |

## BrainSpan expression filter

| Script | Origin | Purpose |
|---|---|---|
| `brainspan_filter_tf_prior.py` | `~/Downloads/` | Restricts a TF–gene/TF–DMR prior to TFs with BrainSpan bulk-RNA mean RPKM above a threshold (paper uses `--mean-thresh 1.0`). Extracts the TF as the substring before the first `\|`, so it works on annotated `TF\|cell\|accession` priors. Writes the filtered prior plus the kept / below-threshold / not-in-BrainSpan diagnostic lists. Cluster inputs (per its docstring): `…/brainspan/expression_matrix.csv`, `…/brainspan/rows_metadata.csv`. |

The OPC–ODC outputs of this step (mean > 1.0 run) are bundled as provenance
for the §3.5 numbers in `../brainspan_filter_outputs_OPC_ODC_tf_filters_1/`:
`TFs_meanRPKM_gt_1.0.txt` (58 kept), `tfs_missing_from_brainspan.txt`
(3 absent: Epitope, GFP, KMT2D), `brainspan_gene_rpkm_summary.tsv`
(per-gene mean/max RPKM).

## Per-lineage DMR β-matrix construction

| Script | Origin | Purpose |
|---|---|---|
| `<L>__apm_<L>_DMR.sh` | `protein_coding/<lineage>/` | Per developmental-stage cluster: extracts CpG-context calls from the allc methylation file, `bedtools map`s them over the sorted lineage DMR BED, and computes the coverage-weighted methylation β per DMR for that cluster. One per lineage: MGE_PVALB, MSN (BACH2, BACH2_2Te, EPHA4), eCGE_LAMP, glial OPC_ODC. |
| `<L>__combine_filtered_methylation_*.sh` | `protein_coding/<lineage>/` | Pastes the per-cluster β columns side-by-side in developmental-stage order into the final DMR × cluster β-matrix for that lineage. |

## Earlier per-lineage prior-construction scripts (kept for provenance)

| Script | Origin | Purpose |
|---|---|---|
| `<L>__TF_CellType_Exp*.sh` | `protein_coding/<lineage>/` | Original per-lineage TF–DMR prior construction (intersect ChIP-Atlas peaks with that lineage's DMR BED, join to TF/cell metadata). Superseded by the consolidated `regenerate_tf_dmr.sh`; included so the per-lineage history is complete, not required to reproduce the final outputs. |

