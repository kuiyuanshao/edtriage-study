### Classification and Regression Trees ###
pacman::p_load(caret, pROC, Matrix, doMC)
registerDoMC(cores = 10)
fitsvm <- function(data, subset_name, outcome_name){
  load(paste0(subset_name, "_", outcome_name, "_cvglmnet_variables.RData"))
  insideSVM <- function(data, subset_name, outcome_name, variables, type = 1){
    data_train <- data[data[[subset_name]] == 1, ]
    formula <- as.formula(paste0("make.names(", outcome_name, ")", " ~ ", 
                                 paste(variables, 
                                       collapse = ' + ')))
    data_train2 <- data_train[, c(outcome_name, variables)]
    
    ctrl <- trainControl(method = "cv",
                         number = 5,
                         summaryFunction = twoClassSummary,
                         classProbs = T, 
                         verboseIter = TRUE,
                         allowParallel = TRUE)
  
    mod.svm <- train(formula, 
                     data = data_train2,
                     method = "svmRadial",
                     tuneLength = 10, 
                     preProc = c("center","scale"),
                     metric = "ROC",
                     trControl = ctrl)
    
    data_valid <- data[data[[subset_name]] == 2, ]
    data_test <- data[data[[subset_name]] == 3, ]
    
    pred_train_full <- predict(mod.svm, newdata = data_train2, type = "prob")
    pred_valid_full <- predict(mod.svm, newdata = data_valid, type = "prob")
    pred_test_full <- predict(mod.svm, newdata = data_test, type = "prob")
    aucs <- c(auc(data_train[[outcome_name]], pred_train_full[, 2]), 
              auc(data_valid[[outcome_name]], pred_valid_full[, 2]), 
              auc(data_test[[outcome_name]], pred_test_full[, 2]))
    
    if (type == 1){
      save(mod.svm,
           pred_valid_full, pred_train_full, pred_test_full, 
           aucs, params, 
           file = paste0(current_dir, "/RData/Trained Models/SVM/", 
                         subset_name, "_", 
                         outcome_name, "_maxVariables",
                         "_SVM", ".RData"))
    }else{
      save(mod.svm,
           pred_valid_full, pred_train_full, pred_test_full, 
           aucs, params, 
           file = paste0(current_dir, "/RData/Trained Models/SVM/", 
                         subset_name, "_", 
                         outcome_name, "_1seVariables",
                         "_SVM", ".RData"))
    }
  }
  insideSVM(data, subset_name, outcome_name, min_variables, type = 1)
  insideSVM(data, subset_name, outcome_name, se1_variables, type = 2)
}

#load the data 
current_dir <- getwd()
load(paste0(current_dir, "/RData/Analysis Datasets/cohort_29nov2022_corrected.RData"))
setwd(file.path(current_dir, "/RData/Trained Models/Variable Selection"))
dir.create(file.path(current_dir, "/RData/Trained Models/SVM"), 
           showWarnings = FALSE)

#CART For part_m1
fitsvm(cohort_29nov2022_corrected, "part_m1", "outcome_d7")
fitsvm(cohort_29nov2022_corrected, "part_m1", "outcome_d28")
fitsvm(cohort_29nov2022_corrected, "part_m1", "outcome_ICU")
fitsvm(cohort_29nov2022_corrected, "part_m1", "outcome_hd")
fitsvm(cohort_29nov2022_corrected, "part_m1", "outcome_ha")
fitsvm(cohort_29nov2022_corrected, "part_m1", "outcome_tr7")
fitsvm(cohort_29nov2022_corrected, "part_m1", "outcome_tr28")
fitsvm(cohort_29nov2022_corrected, "part_m1", "outcome_LOS48h")

#CART For part_m2
fitsvm(cohort_29nov2022_corrected, "part_m2", "outcome_d7")
fitsvm(cohort_29nov2022_corrected, "part_m2", "outcome_d28")
fitsvm(cohort_29nov2022_corrected, "part_m2", "outcome_ICU")
fitsvm(cohort_29nov2022_corrected, "part_m2", "outcome_hd")
fitsvm(cohort_29nov2022_corrected, "part_m2", "outcome_ha")
fitsvm(cohort_29nov2022_corrected, "part_m2", "outcome_tr7")
fitsvm(cohort_29nov2022_corrected, "part_m2", "outcome_tr28")
fitsvm(cohort_29nov2022_corrected, "part_m2", "outcome_LOS48h")

#CART For part_m3
fitsvm(cohort_29nov2022_corrected, "part_m3_d7", "outcome_d7")
fitsvm(cohort_29nov2022_corrected, "part_m3_d28", "outcome_d28")
fitsvm(cohort_29nov2022_corrected, "part_m3_icu", "outcome_ICU")
fitsvm(cohort_29nov2022_corrected, "part_m3_hd", "outcome_hd")
fitsvm(cohort_29nov2022_corrected, "part_m3_ha", "outcome_ha")
fitsvm(cohort_29nov2022_corrected, "part_m3_tr7", "outcome_tr7")
fitsvm(cohort_29nov2022_corrected, "part_m3_tr28", "outcome_tr28")
fitsvm(cohort_29nov2022_corrected, "part_m3_los48h", "outcome_LOS48h")

#Save AUC results into xlsx
#Note that CART will drop some variables in training
setwd(file.path(current_dir, "/RData/Trained Models/SVM"))
files <- list.files()
aucMat <- function(files, variable, partition){
  saved_files <- sort(files[grep(paste0("_", variable, "_SVM.RData$"), files)])
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
  
  write.table(auc_full, file = paste0(variable, "_SVM_", partition,".txt"),
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



