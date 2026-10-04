# Nakdong River: scoped analysis code and fixed inputs

This repository provides the numerical code, fixed analysis tables, comparison
outputs and software records contained in the Water Research revision code
package. It contains the audited core snapshot and the approved Table S5
presentation exporter. It does not reconstruct the full manuscript from raw
sequence reads.

## Snapshot status

The **core remains an analysis candidate**: final sampling-date/season input
authority and some figure/output alignment remain under author review. The
original status and all scientific boundaries are preserved in
[`core/README.md`](core/README.md). Uploading this snapshot does not resolve those
items or designate a different input branch as the final analysis.

The Table S5 exporter uses the approved original-model OLS fields. It transforms
stored coefficients and confidence limits into the table's ratio columns; it
does not fit another model or recalculate confidence intervals, FDR or power.

## Contents

| Location | What it provides |
| --- | --- |
| `core/code/` | Seven numbered numerical modules and two Figure 5 support modules |
| `core/data/` | Fixed annotated/derived inputs and retained comparison tables |
| `core/expected/` | Twelve output comparison targets from the preceding audit |
| `core/environment/` | Pinned Python dependencies and audited R software records |
| `Table_S5/` | Approved table inputs, portable exporter, expected values and original MDE helper |
| `validation/prior_audit/` | Existing numerical validation records, distinguished from today's packaging checks |
| `provenance.json` | Original package hashes, included source files and explicitly excluded office-formatting tools |
| `SHA256SUMS.json` | Exact file checksums for this snapshot |

No development archive, unrelated project, raw FASTQ set, portal credential or
manuscript document is included. The RGI table is approximately 96 MB uncompressed.

## Check the files and software without fitting models

Use Python 3.12 and the package versions in `core/environment/requirements.txt`.
The audited R environment used R 4.4.3, Bioconductor 3.20 and DESeq2 1.46.0;
additional version details are in `core/environment/audited_R_sessionInfo.txt`.
No dependencies are installed automatically.

```bash
python tools/check_package.py
python tools/check_environment.py --rscript /path/to/Rscript
```

`check_package.py` verifies hashes and compiles Python source in memory. It does
not execute analysis modules. `check_environment.py` checks installed versions and
required R packages without installing or updating anything.

For a separate Python environment, install the pinned dependencies explicitly:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -r core/environment/requirements.txt
```

The R installation commands in `core/environment/setup_R.R` are documentation;
use a compatible R/Bioconductor installation and check the recorded versions.

## Recreate the approved Table S5 presentation

This command requires only the Python standard library and refuses to overwrite
an existing output:

```bash
python Table_S5/export_table_s5.py --output results_local/Table_S5.csv
```

It verifies all 180 numerical cells against the supplied approved table and
writes a separate validation JSON. The MDE helper is supplied for source
inspection; this presentation export does not recompute MDE values.

## Run the preserved core snapshot

Run only when you intend to execute the supplied historical analysis branch:

```bash
python core/run_all.py --rscript /path/to/Rscript
```

Outputs are written beneath `core/run/`. To choose a separate destination:

```bash
NAKDONG_OUTPUT="$PWD/results_local/core" python core/run_all.py --rscript /path/to/Rscript
```

Add `--full-fig5` for the retained GCM nuisance refits and historical g-computation
bootstrap analysis. Add `--resume` only with the same package and output directory.
The runner writes `status.json`, logs and `qc/package_validation.json`; inspect
these files because a command finishing is not itself evidence that every stage
or comparison passed.

The core includes clearly labelled diagnostic calculations. Those diagnostics
are retained for transparency and do not replace the specified primary models.
The historical binary-indicator g-computation and its separate coherent-toxin
diagnostic must not be described as the same model.

## What the existing validation establishes

The preceding audit recorded seven completed modules, twelve passing output
comparisons and unchanged package inputs. Full Figure 5 refit checks recorded no
numerical mismatches for 506 GCM rows and 16 g-computation rows. These are prior
execution records, not a claim that all models were rerun while preparing this
repository. Current packaging checks verify source integrity, required files,
software availability and the deterministic Table S5 presentation export.

The package starts from supplied community coordinates and annotated tables.
Upstream read processing, construction of those coordinates, MaAsLin3, all
supplementary analyses and production figure rendering are outside its runnable
coverage. See the original core README for the exact remaining boundaries.
