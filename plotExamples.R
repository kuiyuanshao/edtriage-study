pacman::p_load(caret, tidyverse, rpart, pROC, 
               haven, Hmisc, MASS, predtools, 
               randomForest)
cols <- c("age_g", "gender", "ethnic_gp", "BMI_g", "resident_flag", "smoking_status", "NZdep2018_q",
          "dim_death7_key", "dim_arrival_transport_mode_key", "ed_atte_cate", "attendance_arrived_y",
          "seasons", "TimeofDay2", "ED_number_g", "dim_ED_alcohol_assessment_key", "injury_component_indicator",
          "triage_priority")

#load the data 
current_dir <- getwd()
load(paste0(current_dir, "/RData/Analysis Datasets/cohort_29nov2022.RData"))


##### Create a new category for the variables contains too many NAs


rm(cohort_29nov2022, cohort_29nov2022_cp)

#part_m3_d7
m3_d7_train <- cohort_29nov2022_cp %>% 
  filter(part_m3_d7 == 1) %>% 
  mutate(outcome_d7 = as.factor(outcome_d7))
m3_d7_valid <- cohort_29nov2022_cp %>% 
  filter(part_m3_d7 == 2) %>%
  mutate(outcome_d7 = as.factor(outcome_d7))
m3_d7_test <- cohort_29nov2022_cp %>% 
  filter(part_m3_d7 == 3) %>%
  mutate(outcome_d7 = as.factor(outcome_d7))


##### Calibration plot and ROC curve for outcome_d7
temp_fm_d7 <- reformulate(cols[-c(8, 14)], "outcome_d7")
temp_d7 <- glm(temp_fm_tr7, family = binomial(), data = m3_d7_train, na.action = "na.fail")

m3_d7_valid$attendance_arrived_y[which(m3_d7_valid$attendance_arrived_y %in% c("2020", "2021"))] <- 2019
  
d7_pred <- predict(temp_d7, newdata = m3_d7_valid, type = "response")
m3_d7_valid$d7_pred <- d7_pred
m3_d7_valid$outcome_d7 <- as.numeric(m3_d7_valid) - 1

d7_cali_data <- calibration_plot(data = as.data.frame(m3_d7_valid), obs = "outcome_d7", nTiles = 10, 
                                  pred = "d7_pred", y_lim = c(0, 1), x_lim = c(0, 1))
d7_cali_data <- d7_cali_data$calibration_plot$data
calibration.plot(d7_cali_data, lab = "Death in 7d", title = "d7", 
                 xlim = c(0, 0.1), ylim = c(-0.04, 0.1), 
                 min = -0.035, interval = 0.001)
roc.curve(m3_d7_valid$outcome_d7, d7_pred, varname = "outcome_d7", multi = F)


##### Multiple ROC curves on one plot
rf_d7 <- randomForest(temp_fm_d7, data = m3_d7_train)
rf_d7_pred <- predict(rf_d7, newdata = m3_d7_valid, type = "prob")
predictors <- list(d7_pred, rf_d7_pred[, 2])
roc <- roc.curve(m3_d7_valid$outcome_d7, predictors, 
          varname = "outcome_d7", multi = T, modelnames = c("logistic regression", "random forest"))
roc


m1_train <- cohort_29nov2022_cp_withoutNA %>% 
  filter(part_m1 == 1) 
m1_valid <- cohort_29nov2022_cp_withoutNA %>% 
  filter(part_m1 == 2)
m1_test <- cohort_29nov2022_cp_withoutNA %>% 
  filter(part_m1 == 3)

unique(m1_train$cr_clearance_a)
unique(m1_valid$cr_clearance_a)











