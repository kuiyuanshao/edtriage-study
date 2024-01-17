### Classification and Regression Trees ###
pacman::p_load(xgboost, pROC, mlr, parallelMap, Matrix)


tuning <- function(data, outcome_name){
  parallelStartMulticore(10, level = "mlr.tuneParams")
  set.seed(100)
  dat.task <- makeClassifTask(id = "xgboost", data = data, target = outcome_name)
  # Define the model
  resamp <- makeResampleDesc("CV", iters = 5)
  # Create the learner
  lrn <- makeLearner("classif.xgboost", predict.type = "prob")
  # Create the grid params
  control.grid <- makeTuneControlGrid() 
  ps <- makeParamSet(
    makeDiscreteParam("max_depth", values = c(3, 6, 7, 10)),
    makeDiscreteParam("min_child_weight", values = c(1, 5, 9)),
    makeDiscreteParam("subsample", values = c(0.5, 0.7, 1)),
    makeDiscreteParam("colsample_bytree", values = c(0.5, 0.7, 1)),
    makeDiscreteParam("eta", values = c(0.1, 0.3, 0.7, 1))
  )
  
  tuned <- tuneParams(lrn, task = dat.task, 
                      resampling = resamp, 
                      control = control.grid, 
                      par.set = ps, 
                      measures = auc)
  parallelStop()
  return (tuned)
}


fitxgb <- function(data, subset_name, outcome_name){
  load(paste0(subset_name, "_", outcome_name, "_cvglmnet_variables.RData"))
  insidexgb <- function(data, subset_name, outcome_name, variables, type = 1){
    data_train <- data[data[[subset_name]] == 1, ]
    formula <- as.formula(paste0(outcome_name, " ~ ", 
                                 paste(variables, 
                                       collapse = ' + ')))
    data_train2 <- data_train[, c(outcome_name, variables)]
    dat <- sparse.model.matrix(formula, data = data_train2)
    dat <- dat[, -1]
    label <- as.numeric(data_train2[[outcome_name]]) - 1
    
    data_train3 <- createDummyFeatures(data_train2, target = outcome_name)
    params <- tuning(data_train3, outcome_name)
    set.seed(469)
    xgbparams <- list(max_depth = params[3]$x$max_depth, 
                      min_child_weight = params[3]$x$min_child_weight,
                      subsample = params[3]$x$subsample, 
                      colsample_bytree = params[3]$x$colsample_bytree,
                      eta = params[3]$x$eta)
    
    bst <- xgb.cv(data = dat, label = label, nrounds = 100, nfold = 5, 
                  params = xgbparams, early_stopping_rounds = 3,
                  metrics = "auc", objective = "binary:logistic")

    mod.xgb <- xgboost(data = dat, label = label,
                       nrounds = bst$best_iteration, 
                       params = xgbparams, objective = "binary:logistic")
    
    data_valid <- data[data[[subset_name]] == 2, ]
    data_test <- data[data[[subset_name]] == 3, ]
    
    X_valid <- sparse.model.matrix(formula,
                                   data = data_valid)
    X_valid <- X_valid[, -1]
    
    X_test <- sparse.model.matrix(formula,
                                  data = data_test)
    X_test <- X_test[, -1]
    
    pred_train_full <- predict(mod.xgb, newdata = dat)
    pred_valid_full <- predict(mod.xgb, newdata = X_valid)
    pred_test_full <- predict(mod.xgb, newdata = X_test)
    aucs <- c(auc(data_train[[outcome_name]], pred_train_full), 
              auc(data_valid[[outcome_name]], pred_valid_full), 
              auc(data_test[[outcome_name]], pred_test_full))
    
    if (type == 1){
      save(mod.xgb,
           pred_valid_full, pred_train_full, pred_test_full, 
           aucs, params, 
           file = paste0(current_dir, "/RData/Trained Models/XGB/", 
                         subset_name, "_", 
                         outcome_name, "_maxVariables",
                         "_xgb", ".RData"))
    }else{
      save(mod.xgb,
           pred_valid_full, pred_train_full, pred_test_full, 
           aucs, params, 
           file = paste0(current_dir, "/RData/Trained Models/XGB/", 
                         subset_name, "_", 
                         outcome_name, "_1seVariables",
                         "_xgb", ".RData"))
    }
  }
  insidexgb(data, subset_name, outcome_name, min_variables, type = 1)
  insidexgb(data, subset_name, outcome_name, se1_variables, type = 2)
}

#load the data 
current_dir <- getwd()
load(paste0(current_dir, "/RData/Analysis Datasets/cohort_29nov2022_corrected.RData"))
setwd(file.path(current_dir, "/RData/Trained Models/Variable Selection"))
dir.create(file.path(current_dir, "/RData/Trained Models/XGB"), 
           showWarnings = FALSE)

#CART For part_m1
fitxgb(cohort_29nov2022_corrected, "part_m1", "outcome_d7")
fitxgb(cohort_29nov2022_corrected, "part_m1", "outcome_d28")
fitxgb(cohort_29nov2022_corrected, "part_m1", "outcome_ICU")
fitxgb(cohort_29nov2022_corrected, "part_m1", "outcome_hd")
fitxgb(cohort_29nov2022_corrected, "part_m1", "outcome_ha")
fitxgb(cohort_29nov2022_corrected, "part_m1", "outcome_tr7")
fitxgb(cohort_29nov2022_corrected, "part_m1", "outcome_tr28")
fitxgb(cohort_29nov2022_corrected, "part_m1", "outcome_LOS48h")

#CART For part_m2
fitxgb(cohort_29nov2022_corrected, "part_m2", "outcome_d7")
fitxgb(cohort_29nov2022_corrected, "part_m2", "outcome_d28")
fitxgb(cohort_29nov2022_corrected, "part_m2", "outcome_ICU")
fitxgb(cohort_29nov2022_corrected, "part_m2", "outcome_hd")
fitxgb(cohort_29nov2022_corrected, "part_m2", "outcome_ha")
fitxgb(cohort_29nov2022_corrected, "part_m2", "outcome_tr7")
fitxgb(cohort_29nov2022_corrected, "part_m2", "outcome_tr28")
fitxgb(cohort_29nov2022_corrected, "part_m2", "outcome_LOS48h")

#CART For part_m3
fitxgb(cohort_29nov2022_corrected, "part_m3_d7", "outcome_d7")
fitxgb(cohort_29nov2022_corrected, "part_m3_d28", "outcome_d28")
fitxgb(cohort_29nov2022_corrected, "part_m3_icu", "outcome_ICU")
fitxgb(cohort_29nov2022_corrected, "part_m3_hd", "outcome_hd")
fitxgb(cohort_29nov2022_corrected, "part_m3_ha", "outcome_ha")
fitxgb(cohort_29nov2022_corrected, "part_m3_tr7", "outcome_tr7")
fitxgb(cohort_29nov2022_corrected, "part_m3_tr28", "outcome_tr28")
fitxgb(cohort_29nov2022_corrected, "part_m3_los48h", "outcome_LOS48h")

#Save AUC results into xlsx
#Note that CART will drop some variables in training
setwd(file.path(current_dir, "/RData/Trained Models/XGB"))
files <- list.files()
aucMat <- function(files, variable, partition){
  saved_files <- sort(files[grep(paste0("_", variable, "_xgb.RData$"), files)])
  if (partition == 1){
    saved_files <- saved_files[1:(partition*8)]
  }else if (partition == 2){
    saved_files <- saved_files[9:(partition*8)]
  }else{
    saved_files <- saved_files[17:(partition*8)]
  }
  auc_full <- NULL
  for (i in 1:length(saved_files)){
    load(saved_files[i])
    auc_full <- rbind(auc_full, round(aucs, 5))
  }
  auc_full <- as.matrix(auc_full)
  auc_full <- cbind(c("outcome_d28", "outcome_d7",       
                      "outcome_ha", "outcome_hd",          
                      "outcome_ICU", "outcome_LOS48h",       
                      "outcome_tr28", "outcome_tr7"), auc_full)
  colnames(auc_full) <- c("Outcomes", "Training AUC", "Validation AUC", "Testing AUC")
  
  write.table(auc_full, file = paste0(variable, "_xgb_", partition,".txt"),
              row.names=FALSE)
}
###Base on 1se variables with hyper-parameter tuned
aucMat(files, "1seVariables", 1)
aucMat(files, "1seVariables", 2)
aucMat(files, "1seVariables", 3)

###Base on max variables with hyper-paramter tuned
aucMat(files, "maxVariables", 1)
aucMat(files, "maxVariables", 2)
aucMat(files, "maxVariables", 3)



