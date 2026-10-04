# Nakdong River analysis reproducibility candidate

Status: INTERNAL CANDIDATE, NOT A COMPLETE OR APPROVED REVIEWER RELEASE.

This package contains the numerical definitions and fixed input tables needed for the audited core analyses. It reproduces the current analysis branch; it does not resolve which sampling-date linkage is authoritative. No original source or submitted document is modified.

## Running
Python dependencies are pinned in environment/requirements.txt. R package versions from the audited run are recorded in environment/audited_R_sessionInfo.txt; DESeq2 used R 4.4.3 and Bioconductor 3.20.
Run from any directory:
    python /path/to/package/run_all.py --rscript /path/to/Rscript

The default is core analyses. Add --full-fig5 for the complete GCM nuisance refits and g-computation bootstrap intervals. The process records status.json and per-stage logs. Outputs go to run/ by default; NAKDONG_OUTPUT may specify another directory. No network, upload, or source-data modification is performed.

## Modules and measurements
- 01: reconstruct the 18-class MAG-origin candidate table from RGI, fit adjusted log2(count+1) associations, verify grouped folds, and reconstruct GCM sign-flip tests.
- 02: repeated-site/round PERMANOVA and the submitted DESeq2 filtering/model definitions.
- 03: 1,064 adjusted functional models, MAG size/candidate association, and explicitly labelled multilabel sensitivity.
- 04: independently implemented DML and fitted contrasts, including a labelled coherent-toxin diagnostic.
- 05: the descriptive site-held-out predictive mirror. This is not a formal model-X FDR procedure.
- 06: site-mean spatial tests, sequential PERMANOVA, network summaries and global Procrustes.
- 07: retained source functions for complete GCM and historical g-computation refits; source code is supplied, not just stored scores.

The RGI outcome comprises broad candidates: 573,931 Loose and 312 Strict rows. The 18 class totals and ARG_total are two aggregations of the same MAG-origin records, not independent read-mapped sample abundance measurements. Thirty-two retained samples lack a representative-origin RGI profile and were assigned zero in the original 107-row branch. This is not evidence of biological absence.

## Boundaries that must be resolved before final delivery
1. Sampling-date/season/high-condition input authority is not yet finalized.
2. Some submitted panels differ from the corresponding saved analysis definition/output.
3. The original assembly-process labels used a taxonomy-derived synthetic tree; this package does not validate those labels or silently replace their analysis.
4. The historical g-computation binary-indicator contrast can create combinations inconsistent with MC_ppb. Historical reproduction and the labelled diagnostic must not be described as a corrected final model.
5. Comm_PC1–10 are supplied as fixed derived inputs; their exact generation source remains unresolved. Upstream FASTQ processing and database reconstruction are not performed here.
6. MaAsLin3, beta-diversity partitioning, all supplementary panels, and production figure rendering are not yet included in this core candidate. It is not a complete manuscript reproduction.
7. Final input authority, figure/output reconciliation and remaining source coverage must be checked before releasing this package.

Expected files are comparison targets from the preceding numerical audit, not outputs read in place of model fitting. Source lineage and hashes are recorded in documentation/source_manifest.json. No credentials, account tokens, unrelated projects or raw institutional documents are included.
