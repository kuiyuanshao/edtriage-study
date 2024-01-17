### Classification and Regression Trees ###
pacman::p_load(rpart, pROC, 
               mlr, parallelMap)


cp.select <- function(big.tree){
  min.x <- which.min(big.tree$cptable[, 4]) #column 4 is xerror
  for(i in 1:nrow(big.tree$cptable)) {
    if(big.tree$cptable[i, 4] < big.tree$cptable[min.x, 4] + big.tree$cptable[min.x, 5]) 
      return(big.tree$cptable[i, 1]) #column 5: xstd, column 1: cp 
  }
}

tuning <- function(data, outcome_name){
  parallelStartMulticore(10, level = "mlr.tuneParams")
  set.seed(100)
  dat.task <- makeClassifTask(id = "tree", data = data, target = outcome_name)
  # Define the model
  resamp <- makeResampleDesc("CV", iters = 10)
  # Create the learner
  lrn <- makeLearner("classif.rpart", predict.type = "prob")
  # Create the grid params
  control.grid <- makeTuneControlGrid() 
  ps <- makeParamSet(
    makeDiscreteParam("cp", values = c(0, 1e-5)),
    makeDiscreteParam("minsplit", values = c(5, 25, 100, 300)),
    makeDiscreteParam("maxdepth", values = c(5, 10, 20, 30))
    )
  
  tuned <- tuneParams(lrn, task = dat.task, 
                      resampling = resamp, 
                      control = control.grid, 
                      par.set = ps, 
                      measures = auc)
  parallelStop()
  return (tuned)
}


fitCART <- function(data, subset_name, outcome_name){
  load(paste0(subset_name, "_", outcome_name, "_cvglmnet_variables.RData"))
  insideCART <- function(data, subset_name, outcome_name, variables, type = 1){
    data_train <- data[data[[subset_name]] == 1, ]
    formula <- as.formula(paste0(outcome_name, " ~ ", 
                                 paste(variables, 
                                       collapse = ' + ')))
    
    data_train2 <- data_train[, c(outcome_name, variables)]
    params <- tuning(data_train2, outcome_name)
    set.seed(469)
    mod.rpart <- rpart(formula, data = data_train,
                       control = rpart.control(minsplit = params[3]$x$minsplit,
                                               maxdepth = params[3]$x$maxdepth,
                                               cp = params[3]$x$cp))
    data_valid <- data[data[[subset_name]] == 2, ]
    data_test <- data[data[[subset_name]] == 3, ]
    pred_valid_full <- predict(mod.rpart, newdata = data_valid, type = "prob")
    pred_train_full <- predict(mod.rpart, newdata = data_train, type = "prob")
    pred_test_full <- predict(mod.rpart, newdata = data_test, type = "prob")
    aucs_full <- c(auc(data_train[[outcome_name]], pred_train_full[, 2]), 
                   auc(data_valid[[outcome_name]], pred_valid_full[, 2]),
                   auc(data_test[[outcome_name]], pred_test_full[, 2]))
    
    if (type == 1){
      save(mod.rpart, 
           pred_valid_full, pred_train_full,
           aucs_full, params, 
           file = paste0(current_dir, "/RData/Trained Models/CART/", 
                         subset_name, "_", 
                         outcome_name, "_maxVariables",
                         "_rpart", ".RData"))
    }else{
      save(mod.rpart,
           pred_valid_full, pred_train_full, 
           aucs_full, params, 
           file = paste0(current_dir, "/RData/Trained Models/CART/", 
                         subset_name, "_", 
                         outcome_name, "_1seVariables",
                         "_rpart", ".RData"))
    }
  }
  insideCART(data, subset_name, outcome_name, min_variables, type = 1)
  insideCART(data, subset_name, outcome_name, se1_variables, type = 2)
}

#load the data 
current_dir <- getwd()
load(paste0(current_dir, "/RData/Analysis Datasets/cohort_29nov2022_corrected.RData"))
setwd(file.path(current_dir, "/RData/Trained Models/Variable Selection"))
dir.create(file.path(current_dir, "/RData/Trained Models/CART"), 
           showWarnings = FALSE)

#CART For part_m1
fitCART(cohort_29nov2022_corrected, "part_m1", "outcome_d7")
fitCART(cohort_29nov2022_corrected, "part_m1", "outcome_d28")
fitCART(cohort_29nov2022_corrected, "part_m1", "outcome_ICU")
fitCART(cohort_29nov2022_corrected, "part_m1", "outcome_hd")
fitCART(cohort_29nov2022_corrected, "part_m1", "outcome_ha")
fitCART(cohort_29nov2022_corrected, "part_m1", "outcome_tr7")
fitCART(cohort_29nov2022_corrected, "part_m1", "outcome_tr28")
fitCART(cohort_29nov2022_corrected, "part_m1", "outcome_LOS48h")

#CART For part_m2
fitCART(cohort_29nov2022_corrected, "part_m2", "outcome_d7")
fitCART(cohort_29nov2022_corrected, "part_m2", "outcome_d28")
fitCART(cohort_29nov2022_corrected, "part_m2", "outcome_ICU")
fitCART(cohort_29nov2022_corrected, "part_m2", "outcome_hd")
fitCART(cohort_29nov2022_corrected, "part_m2", "outcome_ha")
fitCART(cohort_29nov2022_corrected, "part_m2", "outcome_tr7")
fitCART(cohort_29nov2022_corrected, "part_m2", "outcome_tr28")
fitCART(cohort_29nov2022_corrected, "part_m2", "outcome_LOS48h")

#CART For part_m3
fitCART(cohort_29nov2022_corrected, "part_m3_d7", "outcome_d7")
fitCART(cohort_29nov2022_corrected, "part_m3_d28", "outcome_d28")
fitCART(cohort_29nov2022_corrected, "part_m3_icu", "outcome_ICU")
fitCART(cohort_29nov2022_corrected, "part_m3_hd", "outcome_hd")
fitCART(cohort_29nov2022_corrected, "part_m3_ha", "outcome_ha")
fitCART(cohort_29nov2022_corrected, "part_m3_tr7", "outcome_tr7")
fitCART(cohort_29nov2022_corrected, "part_m3_tr28", "outcome_tr28")
fitCART(cohort_29nov2022_corrected, "part_m3_los48h", "outcome_LOS48h")

#Save AUC results into xlsx
#Note that CART will drop some variables in training
setwd(file.path(current_dir, "/RData/Trained Models/CART"))
files <- list.files()
aucMat <- function(files, variable, partition){
  saved_files <- sort(files[grep(paste0("_", variable, "_rpart", ".RData$"), files)])
  if (partition == 1){
    saved_files <- saved_files[1:(partition*8)]
  }else if (partition == 2){
    saved_files <- saved_files[9:(partition*8)]
  }else{
    saved_files <- saved_files[17:(partition*8)]
  }
  auc_tune <- NULL
  for (i in 1:length(saved_files)){
    load(saved_files[i])
    auc_tune <- rbind(auc_tune, round(aucs_full, 5))
  }
  auc_tune <- as.matrix(auc_tune)
  auc_tune <- cbind(c("outcome_d28", "outcome_d7",       
                      "outcome_ha", "outcome_hd",          
                      "outcome_ICU", "outcome_LOS48h",       
                      "outcome_tr28", "outcome_tr7"), auc_tune)
  colnames(auc_tune) <- c("Outcomes", "Training AUC", "Validation AUC", "Testing AUC")
  
  write.table(auc_tune, file = paste0(variable, "_rpart_", partition,".txt"),
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



