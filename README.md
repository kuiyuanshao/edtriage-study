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

## Requirements

R with `pacman` installed. Individual scripts declare their packages through
`pacman::p_load()`. Principal dependencies include `haven`, `glmnet`, `pROC`,
`rpart`, `ranger`, `caret`, `xgboost`, `mlr`, `parallelMap`, `doParallel` and `xlsx`.
The SVM uses `doMC`, and some tuning scripts use multicore processing; these
require a compatible operating system. Excel export uses the Java-based `xlsx`
package. The `.Rnw` report also requires `knitr`, `xtable`, `kableExtra` and LaTeX.

## License

MIT. See `LICENSE`.
