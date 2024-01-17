### Classification and Regression Trees ###
pacman::p_load(ranger, pROC, 
               tuneRanger, parallelMap, mlr)


tuning <- function(data, outcome_name, p){
  parallelStartMulticore(10, level = "mlr.tuneParams")
  set.seed(100)
  dat.task <- makeClassifTask(id = "ranger", data = data, target = outcome_name)
  # Define the model
  resamp <- makeResampleDesc("CV", iters = 3)
  # Create the learner
  lrn <- makeLearner("classif.ranger", predict.type = "prob")
  # Create the grid params
  control.grid <- makeTuneControlGrid() 
  ps <- makeParamSet(
    makeDiscreteParam("mtry", values = c(as.integer(sqrt(p)), 2 * as.integer(sqrt(p)))),
    makeDiscreteParam("min.node.size", values = c(1, 2, 4)),
    makeDiscreteParam("sample.fraction", values = c(0.8, 0.9, 1)),
    makeDiscreteParam("num.trees", values = 100)
  )
  
  tuned <- tuneParams(lrn, task = dat.task, 
                      resampling = resamp, 
                      control = control.grid, 
                      par.set = ps, 
                      measures = auc)
  parallelStop()
  return (tuned)
}



fitRF <- function(data, subset_name, outcome_name){
  load(paste0(subset_name, "_", outcome_name, "_cvglmnet_variables.RData"))
  insideRF <- function(data, subset_name, outcome_name, variables, type = 1){
    data_train <- data[data[[subset_name]] == 1, ]
    formula <- as.formula(paste0(outcome_name, " ~ ", 
                                 paste(variables, 
                                       collapse = ' + ')))
    data_train2 <- data_train[, c(outcome_name, variables)]
    params <- tuning(data_train2, outcome_name, length(variables))
    save(params, file = paste0(current_dir, "/RData/Trained Models/RF/", 
                               subset_name, "_", 
                               outcome_name, "_", type, ".RData"))
    return (NULL)
    set.seed(469)
    mod.rf <- ranger(formula, data = data_train2, num.trees = 500, probability = TRUE,
                     mtry = params[3]$x$mtry, min.node.size = params[3]$x$min.node.size,
                     sample.fraction = params[3]$x$sample.fraction)
    
    data_valid <- data[data[[subset_name]] == 2, ]
    data_test <- data[data[[subset_name]] == 3, ]
    pred_valid_full <- predict(mod.rf, data = data_valid, type = "response")
    pred_train_full <- predict(mod.rf, data = data_train, type = "response")
    pred_test_full <- predict(mod.rf, data = data_test, type = "response")
    aucs <- c(auc(data_train[[outcome_name]], pred_train_full$predictions[, 2]), 
              auc(data_valid[[outcome_name]], pred_valid_full$predictions[, 2]), 
              auc(data_test[[outcome_name]], pred_test_full$predictions[, 2]))

    if (type == 1){
      save(mod.rf,
           pred_valid_full, pred_train_full, pred_test_full, 
           aucs, params, 
           file = paste0(current_dir, "/RData/Trained Models/RF/", 
                         subset_name, "_", 
                         outcome_name, "_maxVariables",
                         "_randomforest", ".RData"))
    }else{
      save(mod.rf,
           pred_valid_full, pred_train_full, pred_test_full, 
           aucs, params, 
           file = paste0(current_dir, "/RData/Trained Models/RF/", 
                         subset_name, "_", 
                         outcome_name, "_1seVariables",
                         "_randomforest", ".RData"))
    }
  }
  insideRF(data, subset_name, outcome_name, min_variables, type = 1)
  insideRF(data, subset_name, outcome_name, se1_variables, type = 2)
}

#load the data 
current_dir <- getwd()
load(paste0(current_dir, "/RData/Analysis Datasets/cohort_29nov2022_corrected.RData"))
setwd(file.path(current_dir, "/RData/Trained Models/Variable Selection"))
dir.create(file.path(current_dir, "/RData/Trained Models/RF"), 
           showWarnings = FALSE)

#CART For part_m1
fitRF(cohort_29nov2022_corrected, "part_m1", "outcome_d7")
fitRF(cohort_29nov2022_corrected, "part_m1", "outcome_d28")
fitRF(cohort_29nov2022_corrected, "part_m1", "outcome_ICU")
fitRF(cohort_29nov2022_corrected, "part_m1", "outcome_hd")
fitRF(cohort_29nov2022_corrected, "part_m1", "outcome_ha")
fitRF(cohort_29nov2022_corrected, "part_m1", "outcome_tr7")
fitRF(cohort_29nov2022_corrected, "part_m1", "outcome_tr28")
fitRF(cohort_29nov2022_corrected, "part_m1", "outcome_LOS48h")

#CART For part_m2
fitRF(cohort_29nov2022_corrected, "part_m2", "outcome_d7")
fitRF(cohort_29nov2022_corrected, "part_m2", "outcome_d28")
fitRF(cohort_29nov2022_corrected, "part_m2", "outcome_ICU")
fitRF(cohort_29nov2022_corrected, "part_m2", "outcome_hd")
fitRF(cohort_29nov2022_corrected, "part_m2", "outcome_ha")
fitRF(cohort_29nov2022_corrected, "part_m2", "outcome_tr7")
fitRF(cohort_29nov2022_corrected, "part_m2", "outcome_tr28")
fitRF(cohort_29nov2022_corrected, "part_m2", "outcome_LOS48h")

#CART For part_m3
fitRF(cohort_29nov2022_corrected, "part_m3_d7", "outcome_d7")
fitRF(cohort_29nov2022_corrected, "part_m3_d28", "outcome_d28")
fitRF(cohort_29nov2022_corrected, "part_m3_icu", "outcome_ICU")
fitRF(cohort_29nov2022_corrected, "part_m3_hd", "outcome_hd")
fitRF(cohort_29nov2022_corrected, "part_m3_ha", "outcome_ha")
fitRF(cohort_29nov2022_corrected, "part_m3_tr7", "outcome_tr7")
fitRF(cohort_29nov2022_corrected, "part_m3_tr28", "outcome_tr28")
fitRF(cohort_29nov2022_corrected, "part_m3_los48h", "outcome_LOS48h")

#Save AUC results into xlsx




#Note that CART will drop some variables in training
setwd(file.path(current_dir, "/RData/Trained Models/RF"))
files <- list.files()

aucMat <- function(files, variable, partition){
  saved_files <- sort(files[grep(paste0("_", variable, "_randomforest", ".RData$"), files)])
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
  
  write.table(auc_full, file = paste0(variable, "_randomforest_", partition,".txt"),
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
































### Classification and Regression Trees ###
pacman::p_load(ranger, pROC, 
               tuneRanger, parallelMap, mlr)


tuning <- function(data, outcome_name, p){
  parallelStartMulticore(10, level = "mlr.tuneParams")
  set.seed(100)
  dat.task <- makeClassifTask(id = "ranger", data = data, target = outcome_name)
  # Define the model
  resamp <- makeResampleDesc("CV", iters = 3)
  # Create the learner
  lrn <- makeLearner("classif.ranger", predict.type = "prob")
  # Create the grid params
  control.grid <- makeTuneControlGrid() 
  ps <- makeParamSet(
    makeDiscreteParam("mtry", values = c(as.integer(sqrt(p)), 2 * as.integer(sqrt(p)))),
    makeDiscreteParam("min.node.size", values = c(1, 2, 4)),
    makeDiscreteParam("sample.fraction", values = c(0.8, 0.9, 1)),
    makeDiscreteParam("num.trees", values = 100)
  )
  
  tuned <- tuneParams(lrn, task = dat.task, 
                      resampling = resamp, 
                      control = control.grid, 
                      par.set = ps, 
                      measures = auc)
  parallelStop()
  return (tuned)
}


fitRF <- function(data, subset_name, outcome_name){
  load(paste0(subset_name, "_", outcome_name, "_cvglmnet_variables.RData"))
  insideRF <- function(data, subset_name, outcome_name, variables, type = 1){
    data_train <- data[data[[subset_name]] == 1, ]
    formula <- as.formula(paste0(outcome_name, " ~ ", 
                                 paste(variables, 
                                       collapse = ' + ')))
    data_train2 <- data_train[, c(outcome_name, variables)]
    #params <- tuning(data_train2, outcome_name, length(variables))
    if (type == 1){
      load(paste0("Y:/R/0.2 Code/RData/Trained Models/RF/", 
                  subset_name, "_", outcome_name, "_maxVariables_randomforest.RData"))
    }else{
      load(paste0("Y:/R/0.2 Code/RData/Trained Models/RF/", 
                  subset_name, "_", outcome_name, "_1seVariables_randomforest.RData"))
    }
    set.seed(469)
    mod.rf <- ranger(formula, data = data_train2, num.trees = 1000, probability = TRUE,
                     mtry = params[3]$x$mtry, min.node.size = params[3]$x$min.node.size,
                     sample.fraction = params[3]$x$sample.fraction)
    
    data_valid <- data[data[[subset_name]] == 2, ]
    data_test <- data[data[[subset_name]] == 3, ]
    pred_valid_full <- predict(mod.rf, data = data_valid, type = "response")
    pred_train_full <- predict(mod.rf, data = data_train, type = "response")
    pred_test_full <- predict(mod.rf, data = data_test, type = "response")
    aucs <- c(auc(data_train[[outcome_name]], pred_train_full$predictions[, 2]), 
              auc(data_valid[[outcome_name]], pred_valid_full$predictions[, 2]), 
              auc(data_test[[outcome_name]], pred_test_full$predictions[, 2]))
    
    if (type == 1){
      save(mod.rf,
           pred_valid_full, pred_train_full, pred_test_full, 
           aucs, params, 
           file = paste0(current_dir, "/RData/Trained Models/RF/", 
                         subset_name, "_", 
                         outcome_name, "_maxVariables",
                         "_randomforest", ".RData"))
    }else{
      save(mod.rf,
           pred_valid_full, pred_train_full, pred_test_full, 
           aucs, params, 
           file = paste0(current_dir, "/RData/Trained Models/RF/", 
                         subset_name, "_", 
                         outcome_name, "_1seVariables",
                         "_randomforest", ".RData"))
    }
  }
  insideRF(data, subset_name, outcome_name, min_variables, type = 1)
  insideRF(data, subset_name, outcome_name, se1_variables, type = 2)
}

#load the data 
current_dir <- getwd()
load(paste0(current_dir, "/RData/Analysis Datasets/cohort_29nov2022_corrected.RData"))
setwd(file.path(current_dir, "/RData/Trained Models/Variable Selection"))
dir.create(file.path(current_dir, "/RData/Trained Models/RF"), 
           showWarnings = FALSE)

#CART For part_m1
fitRF(cohort_29nov2022_corrected, "part_m1", "outcome_d7")
fitRF(cohort_29nov2022_corrected, "part_m1", "outcome_d28")
fitRF(cohort_29nov2022_corrected, "part_m1", "outcome_ICU")
fitRF(cohort_29nov2022_corrected, "part_m1", "outcome_hd")
fitRF(cohort_29nov2022_corrected, "part_m1", "outcome_ha")
fitRF(cohort_29nov2022_corrected, "part_m1", "outcome_tr7")
fitRF(cohort_29nov2022_corrected, "part_m1", "outcome_tr28")
fitRF(cohort_29nov2022_corrected, "part_m1", "outcome_LOS48h")

#CART For part_m2
fitRF(cohort_29nov2022_corrected, "part_m2", "outcome_d7")
fitRF(cohort_29nov2022_corrected, "part_m2", "outcome_d28")
fitRF(cohort_29nov2022_corrected, "part_m2", "outcome_ICU")
fitRF(cohort_29nov2022_corrected, "part_m2", "outcome_hd")
fitRF(cohort_29nov2022_corrected, "part_m2", "outcome_ha")
fitRF(cohort_29nov2022_corrected, "part_m2", "outcome_tr7")
fitRF(cohort_29nov2022_corrected, "part_m2", "outcome_tr28")
fitRF(cohort_29nov2022_corrected, "part_m2", "outcome_LOS48h")

#CART For part_m3
fitRF(cohort_29nov2022_corrected, "part_m3_d7", "outcome_d7")
fitRF(cohort_29nov2022_corrected, "part_m3_d28", "outcome_d28")
fitRF(cohort_29nov2022_corrected, "part_m3_icu", "outcome_ICU")
fitRF(cohort_29nov2022_corrected, "part_m3_hd", "outcome_hd")
fitRF(cohort_29nov2022_corrected, "part_m3_ha", "outcome_ha")
fitRF(cohort_29nov2022_corrected, "part_m3_tr7", "outcome_tr7")
fitRF(cohort_29nov2022_corrected, "part_m3_tr28", "outcome_tr28")
fitRF(cohort_29nov2022_corrected, "part_m3_los48h", "outcome_LOS48h")

#Save AUC results into xlsx
#Note that CART will drop some variables in training
setwd(file.path(current_dir, "/RData/Trained Models/RF"))
files <- list.files()
aucMat <- function(files, variable, tune, partition){
  saved_files <- sort(files[grep(paste0("_", variable, "_randomforest_", tune, ".RData$"), files)])
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
  colnames(auc_full) <- c("Outcomes", "Training AUC", "Validation AUC")
  
  write.table(auc_full, file = paste0(variable, "_randomforest_", partition, "_", tune,".txt"),
              row.names=FALSE)
}
###Base on 1se variables with hyper-parameter tuned
aucMat(files, "1seVariables", "1", 1)
aucMat(files, "1seVariables", "1", 2)
aucMat(files, "1seVariables", "1", 3)

###Base on max variables with hyper-paramter tuned
aucMat(files, "maxVariables", "1", 1)
aucMat(files, "maxVariables", "1", 2)
aucMat(files, "maxVariables", "1", 3)

