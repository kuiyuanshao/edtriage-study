pacman::p_load(glmnet, tidyverse, Matrix, hash, Hmisc, MASS)

#load the data 
current_dir <- getwd()
load(paste0(current_dir, "/RData/Analysis Datasets/cohort_29nov2022.RData"))

set.seed(12345)
outcome_LOS48h <- cohort_29nov2022_cp$outcome_LOS48h
LOS48h_T <- which(outcome_LOS48h == 1)
LOS48h_F <- which(outcome_LOS48h == 0)

ind_T <- sample(1:3, size = length(LOS48h_T), replace = T,
                prob = c(0.6, 0.3, 0.1))
ind_F <- sample(1:3, size = length(LOS48h_F), replace = T,
                prob = c(0.6, 0.3, 0.1))

train <- c(LOS48h_T[ind_T == 1], LOS48h_F[ind_F == 1])
valid <- c(LOS48h_T[ind_T == 2], LOS48h_F[ind_F == 2])
test <- c(LOS48h_T[ind_T == 3], LOS48h_F[ind_F == 3])

cohort_29nov2022_cp$part_m3_los48h <- 0
cohort_29nov2022_cp$part_m3_los48h[train] <- 1
cohort_29nov2022_cp$part_m3_los48h[valid] <- 2
cohort_29nov2022_cp$part_m3_los48h[test] <- 3

save(cohort_29nov2022, cohort_29nov2022_cp, 
     file = paste0(current_dir, "/RData/Analysis Datasets/cohort_29nov2022.RData"))