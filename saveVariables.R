pacman::p_load(glmnet, tidyverse, Matrix, hash, Hmisc, MASS, pROC, xlsx)

current_dir <- getwd()
load(paste0(current_dir, "/RData/Analysis Datasets/cohort_29nov2022_corrected.RData"))

#Save to a file
setwd(paste0(current_dir, "/RData/Trained Models/Variable Selection"))

possibleNames <- function(data, names){
  possibleNames <- NULL
  uniqueid <- NULL
  num <- NULL
  for (i in 1:length(names)){
    if (length(unique(data[[names[i]]])) < 20){
      possibleNames <- c(possibleNames, names[i], 
                         paste0(names[i], unique(data[[names[i]]])))
      uniqueid <- c(uniqueid, i, 
                    paste(i, letters[1:length(unique(data[[names[i]]]))]))
    }else{
      possibleNames <- c(possibleNames, names[i])
      uniqueid <- c(uniqueid, i)
      num <- c(num, i)
    }
  }
  df <- cbind(possibleNames, uniqueid)
  return (list(df, num))
}

saveVars <- function(i, wb = NULL){
  files <- list.files()
  glmnet_files <- sort(files[grep("_cvglmnet.RData$", files)])
  variable_files <- sort(files[grep("_variables.RData$", files)])
  
  load(glmnet_files[i])
  load(variable_files[i])
  
  # Get the coefficients related to current outcome variable and partition
  curr_coef_min <- coef(obj_cvglmnet, s = obj_cvglmnet$lambda.min)
  curr_coef_1se <- coef(obj_cvglmnet, s = obj_cvglmnet$lambda.1se)
  
  # Retrieve the indices and variable name with not zero coefficients
  min_inds <- which(curr_coef_min != 0)
  min_variables_levels <- row.names(curr_coef_min)[min_inds]
  
  se1_inds <- which(curr_coef_1se != 0)
  se1_variables_levels <- row.names(curr_coef_1se)[se1_inds]
  
  # Get all the possible combinations for the variables selected by Lasso
  names_min <- possibleNames(cohort_29nov2022_corrected, min_variables)
  num_var_min <- names_min[[2]] # The variable indices that it is a numeric variable
  names_min <- names_min[[1]] # All possible variables with unique ids.
  
  #Similar things, do again for 1 standard error results
  names_1se <- possibleNames(cohort_29nov2022_corrected, se1_variables)
  num_var_1se <- names_1se[[2]]
  names_1se <- names_1se[[1]]
  
  # bind the intercept name to the matrices
  names_min <- rbind(c("(Intercept)", 0), names_min)
  names_1se <- rbind(c("(Intercept)", 0), names_1se)
  
  #Get all the indices that the variable id can be directly converted to a numeric
  #Can be converted to numeric: contains all Variable names, including numeric and categorical variables
  min_id_num <- which(sapply(names_min[, 2], can_convert_to_numeric))[-1]
  se1_id_num <- which(sapply(names_1se[, 2], can_convert_to_numeric))[-1]
  
  # Retrieve the variable indices that are only existing in the non-zero coefficients.
  min_id <- names_min[which(names_min[, 1] %in% min_variables_levels), 2]
  # Some variables will have exactly the same name as the other variable levels
  # E.g. icd10 is variable "icd1" with level "0" in the coefficients, but it can have the same name as variable "icd10"
  # If the unique id can be converted to zero, it means that it is a TRUE variable name (not variable levels)
  min_id_F <- which(sapply(min_id, can_convert_to_numeric))[-1]
  # But it contains some numeric variable names as well, filtering out these numeric variables
  min_id_F <- min_id_F[which(!(names(min_id_F) %in% num_var_min))]
  if (length(min_id_F) > 0){
    min_id <- min_id[-min_id_F]
  }
  
  
  se1_id <- names_1se[which(names_1se[, 1] %in% se1_variables_levels), 2]
  se1_id_F <- which(sapply(se1_id, can_convert_to_numeric))[-1]
  se1_id_F <- se1_id_F[which(!(names(se1_id_F) %in% num_var_1se))]
  if (length(se1_id_F) > 0){
    se1_id <- se1_id[-se1_id_F]
  }
  
  # Only get non-zero coefficients
  coef_min <- curr_coef_min[min_inds, ]
  coef_1se <- curr_coef_1se[se1_inds, ]
  
  # min_id_num would contain all true variable names
  # finding baseline variable names.
  # baseline variable names should not inside all true variable names and should not inside all selected variables by lasso
  baseline_min_ind <- which(names_min[, 1] %in% row.names(curr_coef_min)) #All possible levels selected by lasso, excluding baseline
  baseline_min_id_F <- which(sapply(names_min[baseline_min_ind, 2], can_convert_to_numeric)) #Remove the numeric variables as they are not levels.
  baseline_min_id_F <- baseline_min_id_F[which(!(names(baseline_min_id_F) %in% num_var_min))] #Filter it!
  baseline_min_ind <- baseline_min_ind[-baseline_min_id_F] #All possible levels except baselines
  baseline_min_ind_un <- which(!(names_min[, 2] %in% c(names_min[baseline_min_ind, 2], names_min[min_id_num, 2])))[-1]
  
  df_min <- data.frame(`Selected Variables` = names_min[, 1], 
                       `Variable Levels` = names_min[, 1],
                       `Coefficients` = 0)
  df_min[c(1, min_id_num), 2] <- ""
  df_min[-c(1, min_id_num), 1] <- ""
  #### GET the true baseline indices
  df_min[-c(1, min_id_num), 3] <- 0
  df_min[baseline_min_ind_un, 3] <- "1.0 (ref)"
  df_min[c(1, min_id_num), 3] <- ""
  
  df_min[which(names_min[, 2] %in% min_id), 3] <- round(coef_min[order(names(coef_min))], 5)
  
  baseline_1se_ind <- which(names_1se[, 1] %in% row.names(curr_coef_1se)) #All possible levels selected by lasso, excluding baseline
  baseline_1se_id_F <- which(sapply(names_1se[baseline_1se_ind, 2], can_convert_to_numeric)) #Remove the numeric variables as they are not levels.
  baseline_1se_id_F <- baseline_1se_id_F[which(!(names(baseline_1se_id_F) %in% num_var_1se))] #Filter it!
  baseline_1se_ind <- baseline_1se_ind[-baseline_1se_id_F] #All possible levels except baselines
  baseline_1se_ind_un <- which(!(names_1se[, 2] %in% c(names_1se[baseline_1se_ind, 2], names_1se[se1_id_num, 2])))[-1]
  
  df_1se <- data.frame(`Selected Variables` = names_1se[, 1], 
                       `Variable Levels` = names_1se[, 1],
                       `Coefficients` = 0)
  df_1se[c(1, se1_id_num), 2] <- ""
  df_1se[-c(1, se1_id_num), 1] <- ""
  #### GET the true baseline indices
  df_1se[-c(1, se1_id_num), 3] <- 0
  df_1se[baseline_1se_ind_un, 3] <- "1.0 (ref)"
  df_1se[c(1, se1_id_num), 3] <- ""
  
  df_1se[which(names_1se[, 2] %in% se1_id), 3] <- round(coef_1se[order(names(coef_1se))], 5)
  
  row.names(auc_df) <- c("Maximum AUC", "1 Standard Error AUC")
  colnames(auc_df) <- c("Estimated AUC on Training Data", 
                        "Estimated AUC on Validation Data")
  auc_df <- round(auc_df, 4)
  
  if (i == 1){
    wb <- createWorkbook()
    sh <- createSheet(wb, glmnet_files[i])
    
    rcs <- CellStyle(wb) +
      Font(wb, heightInPoints=15, isBold=TRUE, isItalic=TRUE,
           name="Courier New", color="red") +
      Alignment(h="ALIGN_RIGHT")
    
    bcs <- CellStyle(wb) +
      Font(wb, heightInPoints=15, isBold=TRUE, isItalic=TRUE,
           name="Courier New") +
      Alignment(h="ALIGN_RIGHT")
    
    bbcs <- CellStyle(wb) +
      Font(wb, heightInPoints=12, isBold=TRUE,
           name="Courier New") +
      Alignment(h="ALIGN_RIGHT")
    
    
    addDataFrame(data.frame("No. Variables Selected by Maximum AUC" = length(min_variables),
                            "No. Variables Selected by 1 Standard Error AUC" = length(se1_variables)),
                 sheet = sh, startRow = 1, row.names = FALSE,
                 colnamesStyle = rcs)
    
    addDataFrame(auc_df, sheet = sh, startRow = 5,
                 colnamesStyle = bcs, rownamesStyle = bcs)
    
    addDataFrame(data.frame("Variables Selected by 1 Standard Error AUC" = double()),
                 sheet = sh, startRow = 10, row.names = FALSE,
                 colnamesStyle = rcs)
    
    addDataFrame(df_1se,
                 sheet = sh, startRow = 12, row.names = FALSE,
                 colnamesStyle = bcs)
    
    addDataFrame(data.frame("Variables Selected by Maximum AUC" = double()),
                 sheet = sh, startRow = 12 + nrow(df_1se) + 3, row.names = FALSE,
                 colnamesStyle = rcs)
    
    addDataFrame(df_min,
                 sheet = sh, startRow = 12 + nrow(df_1se) + 5, row.names = FALSE,
                 colnamesStyle = bcs)
    
    rows <- c(12:(12 + nrow(df_1se)), 
              (12 + nrow(df_1se) + 5):(12 + nrow(df_1se) + 5 + nrow(df_min)))
    
    index <- paste(rows, 1, sep = ".")
    
    rows <- getRows(sh, rowIndex = rows)  # get rows
    cells <- getCells(rows, colIndex = 1)  
    lapply(index, function(ii) setCellStyle(cells[[ii]], bbcs))
    
  }else{
    sh <- createSheet(wb, glmnet_files[i])
    
    rcs <- CellStyle(wb) +
      Font(wb, heightInPoints=15, isBold=TRUE, isItalic=TRUE,
           name="Courier New", color="red") +
      Alignment(h="ALIGN_RIGHT")
    
    bcs <- CellStyle(wb) +
      Font(wb, heightInPoints=15, isBold=TRUE, isItalic=TRUE,
           name="Courier New") +
      Alignment(h="ALIGN_RIGHT")
    
    bbcs <- CellStyle(wb) +
      Font(wb, heightInPoints=12, isBold=TRUE,
           name="Courier New") +
      Alignment(h="ALIGN_RIGHT")
    
    addDataFrame(data.frame("No. Variables Selected by Maximum AUC" = length(min_variables),
                            "No. Variables Selected by 1 Standard Error AUC" = length(se1_variables)),
                 sheet = sh, startRow = 1, row.names = FALSE,
                 colnamesStyle = rcs)
    
    addDataFrame(auc_df, sheet = sh, startRow = 5, colnamesStyle = bcs, rownamesStyle = bcs)
    
    addDataFrame(data.frame("Variables Selected by 1 Standard Error AUC" = double()),
                 sheet = sh, startRow = 10, row.names = FALSE,
                 colnamesStyle = rcs)
    
    addDataFrame(df_1se,
                 sheet = sh, startRow = 12, row.names = FALSE, colnamesStyle = bcs)
    
    addDataFrame(data.frame("Variables Selected by Maximum AUC" = double()),
                 sheet = sh, startRow = 12 + nrow(df_1se) + 3, row.names = FALSE,
                 colnamesStyle = rcs)
    
    addDataFrame(df_min,
                 sheet = sh, startRow = 12 + nrow(df_1se) + 5, row.names = FALSE,
                 colnamesStyle = bcs)
    
    rows <- c(12:(12 + nrow(df_1se)), 
              (12 + nrow(df_1se) + 5):(12 + nrow(df_1se) + 5 + nrow(df_min)))
    
    index <- paste(rows, 1, sep = ".")
    
    rows <- getRows(sh, rowIndex = rows)  # get rows
    cells <- getCells(rows, colIndex = 1)  
    lapply(index, function(ii) setCellStyle(cells[[ii]], bbcs))
  }
  return (wb)
  
}

can_convert_to_numeric <- function(x) {
  all(grepl('^(?=.)([+-]?([0-9]*)(\\.([0-9]+))?)$', x, perl = TRUE))  
}

files <- list.files()
glmnet_files <- sort(files[grep("_cvglmnet.RData$", files)])
for (i in 1:length(glmnet_files)){
  # Load current files
  if (i == 1){
    wb <- saveVars(i)
  }else{
    wb <- saveVars(i, wb)
  }
} 

sh <- createSheet(wb, sheetName = "AUC")

bcs <- CellStyle(wb) +
  Font(wb, heightInPoints=15, isBold=TRUE, isItalic=TRUE,
       name="Courier New") +
  Alignment(h="ALIGN_RIGHT")

variable_files <- sort(files[grep("_cvglmnet_variables.RData$", files)])
auc.df <- NULL
for (i in 1:length(variable_files)){
  if (i %% 8 == 1) {
    auc.df <- rbind(auc.df, c(1, 1))
    colnames(auc.df) <- c("mod_auc", "est_auc")
  }
  load(variable_files[i])
  auc.df <- rbind(auc.df, round(auc_df, 5))
}
auc.df <- as.matrix(auc.df)
auc.df <- cbind(rep(c(0, rep(1, 16)), 3),
                rep(c(1, "outcome_d28", 1, "outcome_d7", 1,       
                      "outcome_ha", 1, "outcome_hd", 1,          
                      "outcome_ICU", 1, "outcome_LOS48h", 1,       
                      "outcome_tr28", 1, "outcome_tr7", 1), 3), auc.df)
colnames(auc.df) <- c("Partitions", "Outcomes", "Training AUC", "Validation AUC")

auc.df[which(auc.df[, 1] == 0), 1] <- c("part_m1", "part_m2", "part_m3")
auc.df[which(auc.df == 1)] <- ""

addDataFrame(auc.df, sheet = sh, startRow = 1, colnamesStyle = bcs, rownamesStyle = bcs)



sheets <- getSheets(wb)

for (i in 1:length(sheets)){
  autoSizeColumn(sheets[[i]], colIndex = 1:3)
}


saveWorkbook(wb, file = "variableSelectionLasso.xlsx")


# Merge AUC to One Table
# Continuous los48h outcome
# Final model
# visualization




