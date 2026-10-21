# ED triage prediction modelling

My R code contributions to an emergency department triage study: preparing analysis
inputs, selecting predictors with LASSO, fitting CART, random forest, support vector
machine and XGBoost models, and summarising predictive performance.

This repository contains the code for these tasks, rather than the full study
workflow or all analyses and models described in the wider study report.

## Analysis scripts

The numeric prefixes indicate the analysis stages. The four `05_train_*` scripts
are alternative modelling branches after variable selection, not consecutive steps.
The LASSO export in stage 06 also depends on stage 04, rather than on stage 05.

| Script | Purpose |
| --- | --- |
| `01_import_sas_data.R` | Read the SAS cohort and formats, recode variables and save the analysis dataset. |
| `02_create_los48h_partition.R` | Add an outcome-stratified training/validation/test partition for the 48-hour length-of-stay outcome, using probabilities of 60%/30%/10%. |
| `03_clean_predictors.R` | Check predictor levels across partitions, correct sparse levels and remove selected predictors. |
| `04_select_lasso_predictors.R` | Run cross-validated LASSO and save predictor sets at `lambda.min` and `lambda.1se`. |
| `05_train_cart.R` | Tune and fit classification trees, then export AUC summaries. |
| `05_train_random_forest.R` | Random forest tuning and fitting code, including separate historical tuning and refitting sections. |
| `05_train_svm.R` | Fit radial-kernel support vector machines and export AUC summaries. |
| `05_train_xgboost.R` | Tune and fit gradient-boosted trees and export AUC summaries. |
| `06_export_lasso_results.R` | Export selected predictors, coefficients and LASSO AUC summaries to an Excel workbook. |
| `07_generate_auc_tables.Rnw` | Format exported model AUC summaries as LaTeX tables using knitr. |

## Plotting utilities and examples

- `utils/calibration_roc_functions.R`: calibration and ROC plotting functions.
- `examples/plot_calibration_roc.R`: exploratory examples using those functions;
  these retain dependencies on objects from the original analysis session.

## Analysis scope

The modelling scripts address eight outcomes: mortality within 7 and 28 days,
in-hospital mortality, hospital admission, ED re-presentation within 7 and 28 days,
ICU admission and the 48-hour length-of-stay outcome. They use the existing
`part_m1`, `part_m2` and outcome-specific `part_m3_*` partition columns. Only the
additional `part_m3_los48h` partition is created by the scripts in this repository.

## Data and execution

Clinical data, fitted models and generated results are not included. The scripts
retain their original dataset object names, saved-file names and output directory
conventions so existing analysis files remain compatible.

Before running an analysis script, set the working directory to the repository
root. Start each script from that directory again, since several scripts change
the working directory internally. Source the plotting functions before using the
exploratory plotting examples.

The original analysis expects:

- A SAS cohort named `cohort_29nov2022.sas7bdat` and `formats.sas7bcat`.
- Prepared datasets under `RData/Analysis Datasets/`.
- LASSO results under `RData/Trained Models/Variable Selection/`.
- Model outputs under `RData/Trained Models/CART/`, `RF/`, `SVM/` and `XGB/`.

These are historical research scripts and require the original data and some
manual setup. In particular, the SAS import derives its input location from the
original `0.2 Code` directory layout; the random forest refitting section uses a
`Y:/` path; and the AUC report reads tables from a `Z:/` path. The import also
expects a `cohort_29nov2022` object when saving, although its read step assigns
`cohort_eligible`. The plotting example contains session-dependent objects and
unfinished expressions. Review these sections before execution. The AUC report
retains its original file-selection assumptions and is not a general collector
for every model output.

## Requirements

R with `pacman` installed. Individual scripts declare their packages through
`pacman::p_load()`. Principal dependencies include `haven`, `glmnet`, `pROC`,
`rpart`, `ranger`, `caret`, `xgboost`, `mlr`, `parallelMap`, `doParallel` and `xlsx`.
The SVM uses `doMC`, and some tuning scripts use multicore processing; these
require a compatible operating system. Excel export uses the Java-based `xlsx`
package. The `.Rnw` report also requires `knitr`, `xtable`, `kableExtra` and LaTeX.

## Previous filenames

| Previous name | Current path |
| --- | --- |
| `readSASData.R` | `01_import_sas_data.R` |
| `part_m3_LOS48h.R` | `02_create_los48h_partition.R` |
| `variableCheck.R` | `03_clean_predictors.R` |
| `variableSelection.R` | `04_select_lasso_predictors.R` |
| `fitCART.R` | `05_train_cart.R` |
| `fitRF.R` | `05_train_random_forest.R` |
| `fitSVM.R` | `05_train_svm.R` |
| `fitXGBOOST.R` | `05_train_xgboost.R` |
| `saveVariables.R` | `06_export_lasso_results.R` |
| `aucTables.Rnw` | `07_generate_auc_tables.Rnw` |
| `plotFunctions.R` | `utils/calibration_roc_functions.R` |
| `plotExamples.R` | `examples/plot_calibration_roc.R` |

## License

MIT. See `LICENSE`.
