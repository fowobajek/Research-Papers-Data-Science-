# Educational attainment and age at first birth in Nigeria (SAP3)

**Reproducibility code for:** *Educational attainment and the timing of first childbirth across Nigerian birth cohorts* (submitted to *Demographic Research*).

This folder contains the R analysis code for a repeated cross-sectional study of women in the Nigeria Demographic and Health Surveys (DHS) of **1990, 2003, 2008, 2013, 2018 and 2023–24** (the last survey is labelled **2024** in the R code). The primary comparison is the 1950–1959 versus 1990–1999 birth cohorts. The analysis reports a descriptive, model-based decomposition of the difference in restricted mean age at first birth through age 30, using observed educational attainment and an exploratory joint schooling/early-union specification.

## Repository contents

| File | Purpose |
| --- | --- |
| `R/00_config.R` | Original shared configuration, source-file detection, DHS harmonisation, weights, strata and bootstrap helpers. |
| `R/01_build_sap3_only.R` | SAP3-only preprocessing extracted from the original `01_build_data.R` (other SAPs excluded). |
| `R/04_sap3_schooling_mediation.R` | Original SAP3 mediator and survival models, standardisation, bootstrap, CSV exports and figures. |
| `R/99_session_info.R` | Save R version, package versions and session information after an analysis run. |
| `.gitignore` | Prevent accidental commits of microdata, identified files and local generated outputs. |

## Data access

**Individual-level DHS data are not included.** Obtain authorised access from the [DHS Program](https://dhsprogram.com/data/). This analysis expects DHS **Individual Recode (IR)** SPSS `.sav` files from each round and the **1990 Wealth Index (WI)** file when required by the harmonisation code. The extraction paths are discovered recursively under `DHS_DIR`; `R/00_config.R` expects the survey year in the directory hierarchy. Do not upload or redistribute licensed DHS microdata.

The analysis inputs include standard DHS variables for survey design (`v005`, `v021`, `v022`), birth timing (`v010`, `v011`, `v201`, `v211`), education (`v106`, `v133`), first union (`v511`), residence (`v025`), religion (`v130`), and geographic region (`v024`/`sstate`). See the scripts for the complete variable list.

## Software and packages

Use an R installation supporting `rstpm2` and a fork-capable operating system (typically **Linux or macOS** for the original parallel bootstrap implementation). Required R packages: `haven`, `dplyr`, `tidyr`, `purrr`, `readr`, `ggplot2`, `survey`, `data.table`, `broom`, `rstpm2`, `survival`, and `nnet`. The original code uses `parallel::mclapply()`; this approach does **not** provide multicore parallelism on Windows. This public code has not been rerun against licensed microdata in the publishing environment.

## Running the analysis

Place the authorised DHS files in a local directory that is **not tracked in Git**. From a terminal in this folder, set the project directory and the relevant paths. Example (Linux/macOS):

```bash
export SAP_ROOT="$(pwd)"
export DHS_DIR="/absolute/path/to/authorised-DHS-SAV-files"
export SAP_OUT="$(pwd)/output"
export BOOT_B_SAP3=1000
export CORES=3
export OPENBLAS_NUM_THREADS=1

Rscript R/01_build_sap3_only.R
Rscript R/04_sap3_schooling_mediation.R
Rscript R/99_session_info.R
```

**Key outputs:** `output/data/sap3_women.rds` (local individual-level derived data; do not publish), `output/sap3/sap3_table1_cumulative_first_birth.csv`, `sap3_table2_mediation_results.csv`, `sap3_table3_mediation_by_zone.csv`, `sap3_table4_sensitivity.csv`, the coefficient tables and figures, and `output/reproducibility/R_session_info.txt`.

The code intentionally requires at least 200 bootstrap replicates and defaults to 1,000. The original analytical code and its reported outputs were supplied by the study authors. The publication of this code does not imply independent reproduction of its estimates on restricted-access DHS data.

## Interpretation and important caveats

* **Educational attainment is recorded at the survey interview**, not necessarily before the respondent's first birth; the schooling-related component is a statistical decomposition, **not an identified causal effect of schooling**.
* Historical birth cohort and survey round overlap is limited. The analytical adjustment for survey round does not eliminate all cohort-period confounding.
* The analysis is limited to models, eligible observations, covariates, sample selection and censoring defined in the scripts. It does **not** reconstruct a longitudinal schooling trajectory.
* The 2023–24 survey is labelled `year = 2024` in the scripts.
* The bootstrap uses primary-sampling-unit resampling within strata. Running the original `mclapply()` implementation on Windows requires adaptation; no modified results are claimed here.

## Citation

Cite the associated manuscript using the final bibliographic details when published. Until then, cite the descriptive study title above and link to this code folder.

## Contact

Research code maintained by the manuscript author. For access to survey microdata, apply directly to the DHS Program; this repository does not grant data access.
