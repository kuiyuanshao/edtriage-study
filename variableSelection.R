pacman::p_load(glmnet, tidyverse, Matrix, hash, Hmisc, MASS, doParallel, pROC, xlsx)

varSel <- function(data, subset_name, outcome_name){
  outcomes <- c("outcome_d7", "outcome_d28", "outcome_ICU",
                "outcome_hd", "outcome_ha", "outcome_tr7",
                "outcome_tr28", "outcome_LOS48h", "EDADU_IP_los", "LOS_total")
  partitions <- c("part_m1", "part_m2", "part_m3_ha", 
                  "part_m3_d7", "part_m3_d28", "part_m3_hd",                     
                  "part_m3_icu", "part_m3_tr7", "part_m3_tr28",
                  "part_m3_los48h")
  data_train <- data[data[[subset_name]] == 1, ]
  outcome_ind <- which(outcomes == outcome_name)
  data_train[, c(partitions, outcomes[-outcome_ind])] <- NULL

  X_train <- sparse.model.matrix(as.formula(paste(outcome_name, "~ .")),
                                 data = data_train)
  X_train <- X_train[, -1]
  
  library(doParallel)
  registerDoParallel(5)
  #Parallel computing cross validation 
  set.seed(1708)
  obj_cvglmnet <- cv.glmnet(x = X_train,
                            y = data_train[[outcome_name]],
                            type.measure = "auc",
                            family = "binomial", 
                            nfolds = 10,
                            alpha = 1, 
                            trace.it = T,
                            parallel = TRUE,
                            maxit = 100000,
                            nlambda = 100)
  setwd(paste0(current_dir, "/RData/Trained Models/Variable Selection"))
  save(obj_cvglmnet, file = paste0(subset_name, "_", 
                                   outcome_name, "_cvglmnet.RData"))
}

matchVars <- function(data, subset_name, outcome_name){
  outcomes <- c("outcome_d7", "outcome_d28", "outcome_ICU",
                "outcome_hd", "outcome_ha", "outcome_tr7",
                "outcome_tr28", "outcome_LOS48h", "EDADU_IP_los", 
                "LOS_total")
  partitions <- c("part_m1", "part_m2", "part_m3_ha", 
                  "part_m3_d7", "part_m3_d28", "part_m3_hd",                     
                  "part_m3_icu", "part_m3_tr7", "part_m3_tr28",
                  "part_m3_los48h")
  data_train <- data[data[[subset_name]] == 1, ]
  outcome_ind <- which(outcomes == outcome_name)
  
  data_valid <- data[data[[subset_name]] == 2, ]
  data_test <- data[data[[subset_name]] == 3, ]
  data_train[, c(partitions, outcomes[-outcome_ind])] <- NULL
  data_valid[, c(partitions, outcomes[-outcome_ind])] <- NULL
  data_test[, c(partitions, outcomes[-outcome_ind])] <- NULL
  
  load(paste0(subset_name, "_", outcome_name, "_cvglmnet.RData"))
  
  optimal_min <- coef(obj_cvglmnet, 
                      s = obj_cvglmnet$lambda.min)
  min_inds <- which(optimal_min != 0)
  min_variables <- row.names(optimal_min)[min_inds]
  min_variables <- min_variables[!(min_variables %in% '(Intercept)')]
  #Based on 1SE
  optimal_1se <- coef(obj_cvglmnet, 
                      s = obj_cvglmnet$lambda.1se)
  se1_inds <- which(optimal_1se != 0)
  se1_variables <- row.names(optimal_1se)[se1_inds]
  se1_variables <- se1_variables[!(se1_variables %in% '(Intercept)')]
  
  
  possibleNames <- hash()
  data_train[[outcome_name]] <- NULL
  for (i in names(data_train)){
    if (length(unique(data_train[[i]])) < 20){
      possibleNames[[i]] <- paste0(i, unique(data_train[[i]]))
    }else{
      possibleNames[[i]] <- i
    }
  }
  inds <- NULL
  for (k in 1:length(possibleNames)){
    if (length(which(values(possibleNames)[[k]] %in% min_variables)) >= 1){
      inds <- c(inds, k)
    }
  }
  
  min_variables <- keys(possibleNames)[inds]
  
  inds <- NULL
  for (k in 1:length(possibleNames)){
    if (length(which(values(possibleNames)[[k]] %in% se1_variables)) >= 1){
      inds <- c(inds, k)
    }
  }
  se1_variables <- keys(possibleNames)[inds]
  
  auc_min <- obj_cvglmnet$cvm[which(obj_cvglmnet$lambda == obj_cvglmnet$lambda.min)]
  auc_1se <- obj_cvglmnet$cvm[which(obj_cvglmnet$lambda == obj_cvglmnet$lambda.1se)]
  X_valid <- sparse.model.matrix(as.formula(paste(outcome_name, "~ .")),
                                 data = data_valid)
  X_valid <- X_valid[, -1]
  pred_auc_min <- predict(obj_cvglmnet, newx = X_valid,
                          type = 'response', s ="lambda.min",
                          gamma = "lambda.min")
  pred_auc_1se <- predict(obj_cvglmnet, newx = X_valid, 
                          type = 'response', s ="lambda.1se",
                          gamma = "lambda.1se")
  
  pred_auc_min <- as.vector(pred_auc_min)
  pred_auc_1se <- as.vector(pred_auc_1se)
  
  auc_df <- data.frame(mod_auc = c(auc_min, auc_1se),
                       est_auc = c(auc(data_valid[[outcome_name]], pred_auc_min), 
                                   auc(data_valid[[outcome_name]], pred_auc_1se)))
  
  save(min_variables, se1_variables, auc_df,
       file = paste0(subset_name, "_", outcome_name, "_cvglmnet_variables.RData"))
}

#load the data 
current_dir <- getwd()
load(paste0(current_dir, "/RData/Analysis Datasets/cohort_29nov2022_corrected.RData"))
dir.create(file.path(current_dir, "/RData/Trained Models/Variable Selection"), 
           showWarnings = FALSE)

#Remove other EWS predictors, only keep ews_g
ews_vars <- names(cohort_29nov2022_corrected)[grep("_ews$", names(cohort_29nov2022_corrected))]
cohort_29nov2022_clean <- cohort_29nov2022_corrected
for (i in ews_vars){
  cohort_29nov2022_clean[[i]] <- NULL
}

#Variable Selection For part_m1
varSel(cohort_29nov2022_clean, "part_m1", "outcome_d7")
varSel(cohort_29nov2022_clean, "part_m1", "outcome_d28")
varSel(cohort_29nov2022_clean, "part_m1", "outcome_ICU")
varSel(cohort_29nov2022_clean, "part_m1", "outcome_hd")
varSel(cohort_29nov2022_clean, "part_m1", "outcome_ha")
varSel(cohort_29nov2022_clean, "part_m1", "outcome_tr7")
varSel(cohort_29nov2022_clean, "part_m1", "outcome_tr28")
varSel(cohort_29nov2022_clean, "part_m1", "outcome_LOS48h")

#Variable Selection For part_m2
varSel(cohort_29nov2022_clean, "part_m2", "outcome_d7")
varSel(cohort_29nov2022_clean, "part_m2", "outcome_d28")
varSel(cohort_29nov2022_clean, "part_m2", "outcome_ICU")
varSel(cohort_29nov2022_clean, "part_m2", "outcome_hd")
varSel(cohort_29nov2022_clean, "part_m2", "outcome_ha")
varSel(cohort_29nov2022_clean, "part_m2", "outcome_tr7")
varSel(cohort_29nov2022_clean, "part_m2", "outcome_tr28")
varSel(cohort_29nov2022_clean, "part_m2", "outcome_LOS48h")

#Variable Selection For part_m3
varSel(cohort_29nov2022_clean, "part_m3_d7", "outcome_d7")
varSel(cohort_29nov2022_clean, "part_m3_d28", "outcome_d28")
varSel(cohort_29nov2022_clean, "part_m3_icu", "outcome_ICU")
varSel(cohort_29nov2022_clean, "part_m3_hd", "outcome_hd")
varSel(cohort_29nov2022_clean, "part_m3_ha", "outcome_ha")
varSel(cohort_29nov2022_clean, "part_m3_tr7", "outcome_tr7")
varSel(cohort_29nov2022_clean, "part_m3_tr28", "outcome_tr28")
varSel(cohort_29nov2022_clean, "part_m3_los48h", "outcome_LOS48h")

# All Done!

#Save the variables selected
#For part_m1
matchVars(cohort_29nov2022_clean, "part_m1", "outcome_d7")
matchVars(cohort_29nov2022_clean, "part_m1", "outcome_d28")
matchVars(cohort_29nov2022_clean, "part_m1", "outcome_ICU")
matchVars(cohort_29nov2022_clean, "part_m1", "outcome_hd")
matchVars(cohort_29nov2022_clean, "part_m1", "outcome_ha")
matchVars(cohort_29nov2022_clean, "part_m1", "outcome_tr7")
matchVars(cohort_29nov2022_clean, "part_m1", "outcome_tr28")
matchVars(cohort_29nov2022_clean, "part_m1", "outcome_LOS48h")

#For part_m2
matchVars(cohort_29nov2022_clean, "part_m2", "outcome_d7")
matchVars(cohort_29nov2022_clean, "part_m2", "outcome_d28")
matchVars(cohort_29nov2022_clean, "part_m2", "outcome_ICU")
matchVars(cohort_29nov2022_clean, "part_m2", "outcome_hd")
matchVars(cohort_29nov2022_clean, "part_m2", "outcome_ha")
matchVars(cohort_29nov2022_clean, "part_m2", "outcome_tr7")
matchVars(cohort_29nov2022_clean, "part_m2", "outcome_tr28")
matchVars(cohort_29nov2022_clean, "part_m2", "outcome_LOS48h")

#For part_m3
matchVars(cohort_29nov2022_clean, "part_m3_d7", "outcome_d7")
matchVars(cohort_29nov2022_clean, "part_m3_d28", "outcome_d28")
matchVars(cohort_29nov2022_clean, "part_m3_icu", "outcome_ICU")
matchVars(cohort_29nov2022_clean, "part_m3_hd", "outcome_hd")
matchVars(cohort_29nov2022_clean, "part_m3_ha", "outcome_ha")
matchVars(cohort_29nov2022_clean, "part_m3_tr7", "outcome_tr7")
matchVars(cohort_29nov2022_clean, "part_m3_tr28", "outcome_tr28")
matchVars(cohort_29nov2022_clean, "part_m3_los48h", "outcome_LOS48h")

