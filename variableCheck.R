pacman::p_load(glmnet, tidyverse, Matrix, hash, Hmisc, MASS)

#load the data 
current_dir <- getwd()
load(paste0(current_dir, "/RData/Analysis Datasets/cohort_29nov2022.RData"))

#rm(cohort_29nov2022, cohort_29nov2022_cp)

varCheck <- function(data, subset_name){
  outcomes <- c("outcome_d7", "outcome_d28", "outcome_ICU",
                "outcome_hd", "outcome_ha", "outcome_tr7",
                "outcome_tr28", "outcome_LOS48h")
  partitions <- c("part_m1", "part_m2", "part_m3_ha", 
                  "part_m3_d7", "part_m3_d28", "part_m3_hd",                     
                  "part_m3_icu", "part_m3_tr7", "part_m3_tr28",
                  "part_m3_los48h")
  impvar <- Cs(cr_a, Potassium_blood_a, Sodium_blood_a, 
               Haemoglobin_a, White_cell_count_a, Platelet_Count_a, 
               Neutrophil_a, Lymphocytes_a, Monocytes_a, Eosinophils_a, 
               RBC_a, Albumin_a, Basophils_a, ALT_a, GGT_a, HCT_a, RDW_a, 
               CRP_a, MCV_a, eGFR_blood_a, BILIRUBIN_a, alp_a, Glucose_a,
               PROTEIN_a, Immature_Granulocytes_a, Phosphate_blood_a, 
               Calcium_a, Globulin_a, HbA1c_a, Ferritin_a, TSH_a, cholesterol_a,
               Chloride_blood_a, INR_a, urea_blood_a, TnI_a, aptt_a, 
               Fibrinogen_a, PR_a, Anion_gap_a, Magnesium_blood_a, vitB12_folate_a, 
               ph_blood_a, Bicarbonate_a, base_excess_a, PO2_a, lactate_a, 
               triglyceride_a, Lipase_a, Urine_Culture_a, CK_a, O2sat_a, 
               AST_a, Blood_Culture_a)
  data_train <- data[data[[subset_name]] == 1, ]
  
  data_valid <- data[data[[subset_name]] == 2, ]
  data_test <- data[data[[subset_name]] == 3, ]
  data_train[, c(partitions, outcomes)] <- NULL
  data_valid[, c(partitions, outcomes)] <- NULL
  data_test[, c(partitions, outcomes)] <- NULL
  checkFun <- function(data_train, data_valid, data_test){
    probNames <- list(NULL, NULL)
    for (i in names(data_train)){
      if (length(unique(data_train[[i]])) < 10){
        
        if (!setequal(sort(unique(data_valid[[i]])), sort(unique(data_train[[i]]))) ||
            !setequal(sort(unique(data_valid[[i]])), sort(unique(data_test[[i]]))) ||
            !setequal(sort(unique(data_train[[i]])), sort(unique(data_test[[i]])))){
          probNames[[1]] <- c(probNames[[1]], i)
        }
      }
      if (length(unique(data_train[[i]])) == 1){
        probNames[[2]] <- c(probNames[[2]], i)
      }
    }
    return (probNames)
  }
  
  probNames <- checkFun(data_train, data_valid, data_test)
  realprob <- NULL
  for (i in probNames[[1]]){
    if (!all(c(unique(data_valid[[i]]), 
               unique(data_test[[i]])) %in% unique(data_train[[i]]))){
      realprob <- c(realprob, i)
    }
  }
  labtest <- realprob[grep("_a$", realprob)]
  labtest <- labtest[which(labtest %in% impvar)]
  
  realprob <- c(realprob[-grep("_a$", realprob)], labtest)
  return (realprob)
}



# Running all partitions to find variables with problems

realprob_list <- NULL
partitions <- c("part_m1", "part_m2", "part_m3_ha", 
                "part_m3_d7", "part_m3_d28", "part_m3_hd",                     
                "part_m3_icu", "part_m3_tr7", "part_m3_tr28", 
                "part_m3_los48h")
for (i in partitions){
  realprob_list <- c(realprob_list, varCheck(cohort_29nov2022_cp, i))
}
realprob_list <- unique(realprob_list)
realprob_list 


# Check what's going on

for (i in realprob_list){
  print(paste0("#### ", i))
  print(table(cohort_29nov2022_cp[[i]]))
}


# Correct these variables
## Many of these are too few to be significantly affect the variable selection
### Temporary set to NA for proper implementation in the variable selection.
cohort_29nov2022_cp$temperature_ews[cohort_29nov2022_cp$temperature_ews == 3] <- 99
cohort_29nov2022_cp$HR_ews[cohort_29nov2022_cp$HR_ews == 4] <- 99
cohort_29nov2022_cp$respr_ews[cohort_29nov2022_cp$respr_ews == 6] <- 99
cohort_29nov2022_cp$PR_a[cohort_29nov2022_cp$PR_a == "L"] <- "M"
cohort_29nov2022_cp$Haemoglobin_a[cohort_29nov2022_cp$Haemoglobin_a == "A"] <- "M"
cohort_29nov2022_cp$lactate_a[cohort_29nov2022_cp$lactate_a == "L"] <- "M"

# Drop the corresponding factor levels
cohort_29nov2022_cp$temperature_ews <- droplevels(cohort_29nov2022_cp$temperature_ews)
cohort_29nov2022_cp$HR_ews <- droplevels(cohort_29nov2022_cp$HR_ews)
cohort_29nov2022_cp$respr_ews <- droplevels(cohort_29nov2022_cp$respr_ews)
cohort_29nov2022_cp$PR_a <- droplevels(cohort_29nov2022_cp$PR_a)
cohort_29nov2022_cp$Haemoglobin_a <- droplevels(cohort_29nov2022_cp$Haemoglobin_a)
cohort_29nov2022_cp$lactate_a <- droplevels(cohort_29nov2022_cp$lactate_a)

# One more check
for (i in realprob_list){
  print(paste0("#### ", i))
  print(table(cohort_29nov2022_cp[[i]]))
}

# Drop more useless variables and Save to a new RData file
cohort_29nov2022_corrected <- cohort_29nov2022_cp

#si_766877008 just has one 1, all the rest are 0s
cohort_29nov2022_corrected$si_766877008 <- NULL
cohort_29nov2022_corrected$si_71186008 <- NULL

# Remove some labtest predictors
impvar <- Cs(cr_a, Potassium_blood_a, Sodium_blood_a, 
             Haemoglobin_a, White_cell_count_a, Platelet_Count_a, 
             Neutrophil_a, Lymphocytes_a, Monocytes_a, Eosinophils_a, 
             RBC_a, Albumin_a, Basophils_a, ALT_a, GGT_a, HCT_a, RDW_a, 
             CRP_a, MCV_a, eGFR_blood_a, BILIRUBIN_a, alp_a, Glucose_a,
             PROTEIN_a, Immature_Granulocytes_a, Phosphate_blood_a, 
             Calcium_a, Globulin_a, HbA1c_a, Ferritin_a, TSH_a, cholesterol_a,
             Chloride_blood_a, INR_a, urea_blood_a, TnI_a, aptt_a, 
             Fibrinogen_a, PR_a, Anion_gap_a, Magnesium_blood_a, vitB12_folate_a, 
             ph_blood_a, Bicarbonate_a, base_excess_a, PO2_a, lactate_a, 
             triglyceride_a, Lipase_a, Urine_Culture_a, CK_a, O2sat_a, 
             AST_a, Blood_Culture_a)

labtest_vars <- names(cohort_29nov2022_corrected)[grep("_a$", names(cohort_29nov2022_corrected))]
vars_notinc <- labtest_vars[which(!(labtest_vars  %in% impvar))]

for (i in vars_notinc){
  cohort_29nov2022_corrected[[i]] <- NULL
}

save(cohort_29nov2022_corrected, file = paste0(current_dir, "/RData/Analysis Datasets/cohort_29nov2022_corrected.RData"))
