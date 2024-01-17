## load libraries
## install.packages("pacman")
pacman::p_load(haven, dplyr, Hmisc, MASS, ggplot2, predtools)

### get current file directory (for R studio only)
datname <- "cohort_29nov2022"

current_dir <- getwd()
path_sasData <- gsub("0.2 Code", "0.1 Database/2. Analysis Dataset", current_dir)
setwd(gsub("0.2 Code", "0.1 Database/2. Analysis Dataset", current_dir))

#For cohort_eligible
cohort_eligible <- read_sas(paste0(datname, ".sas7bdat"), "formats.sas7bcat")
## head(cohort_29nov2022)

## Create a copy
ce_cp <- cohort_eligible

## Transform necessary categorical variables to factors
ce_cp$age_g <- cut(ce_cp$age_g, breaks = c(18, 45, 65, 85, 110), 
                             right = F, labels = c("18-<45", "45-<65", "65-<85", ">=85"))
ce_cp$age_g <- as.character(ce_cp$age_g)

## Match the formats information
do_fmt <- function(x, fmt) {
  lbl <- if (!missing(fmt))
    unlist(unname(fmt)) else attr(x, 'labels')
  
  if (sum(as.vector(attr(x, "labels")) %in% unique(x)) != 0)
    tryCatch(names(lbl[match(unlist(x), lbl)]),
                error = function(e) {
                message(sprintf('formatting failed for %s', attr(x, 'label')),
                        domain = NA)
                  x
            }) else x
}
ce_cp[] <- lapply(ce_cp, do_fmt)
remove_attr <- function(x){
  attributes(x) <- NULL
  x
}
ce_cp[] <- lapply(ce_cp, remove_attr)

#Setting "" to NA
ce_cp[ce_cp == ""] <- NA
NAnames <- names(which(colSums(is.na(ce_cp)) > 0))
for (i in 1:length(NAnames)){
  ce_cp[[NAnames[i]]] <- as.character(ce_cp[[NAnames[i]]])
}
ce_cp[is.na(ce_cp)] <- "NAN"

#Specify the baselines for some variables
cols <- c("age_g", "gender", "ethnic_gp", "BMI_g", "resident_flag", "smoking_status", "NZdep2018_q",
               "dim_death7_key", "dim_arrival_transport_mode_key", "ed_atte_cate", "attendance_arrived_y",
               "seasons", "TimeofDay2", "ED_number_g", "dim_ED_alcohol_assessment_key", "injury_component_indicator",
               "triage_priority")
defaultlevels <- c("18-<45", "Female", "European and Other", "18.5-<25.0", "Y", "Never Smoked", "5-6", "NSH", "Not Specified/Other", 
                   "Not specified", 2016, "Summer", "0-<8", 1, "Not Known", "N", 3)
for (i in 1:length(cols)){
  if (cols[i] %in% names(ce_cp)){
    ce_cp[[cols[i]]] <- relevel(as.factor(ce_cp[[cols[i]]]), ref = defaultlevels[i])
  }
}


cohort_29nov2022 <- cohort_eligible

ce_cp <- subset(as.data.frame(ce_cp), 
                select = -c(patient_reference_no,
                            attendance_encounter_id,
                            attendance_arrived_date,
                            attendance_arrived_y, 
                            attendance_arrived_datepart))
ce_cp$SBP_ews[ce_cp$SBP_ews == 4] <- 0
for (i in (names(ce_cp))){
  if (length(unique(ce_cp[[i]])) < 20){
    ce_cp[[i]] <- as.factor(ce_cp[[i]])
  }
}

cohort_29nov2022_cp <- ce_cp

#Filtering some covariates that shouldn't be included into the Model
setwd(paste0(current_dir, "/RData/Analysis Datasets"))
save(cohort_29nov2022, cohort_29nov2022_cp, file = paste0(datname, ".RData"))





#Changing variables to factors.
covarname <- Cs(si_39579001, si_128053003, si_125605004, si_128069005, si_21522001,
                si_125643001,	si_128477000,	si_25702006, si_3006004, si_52684005,
                si_195967001,	si_82313006, si_712893003, si_279039007, si_283682007,	
                si_271807003,	si_267036007,	si_125666000,	si_13645005, si_62914000,
                si_410429000,	si_cardiacoth, si_128045006, si_40733004, si_262525000,	
                si_29857009, si_309585006, si_88111009, si_49727002, si_766877008,	
                si_3415004, si_248062006,	si_65759007, si_80394007, si_2919008,
                si_108367008,	si_404640003,	si_185389009,	si_18949003, si_162356005,	
                si_371708003, si_81723002, si_249366005, si_371704001, si_371405004,	
                si_125593007,	si_161898004,	si_386661006, si_78164000, si_125670008,	
                si_74474003, si_300479008, si_8765009, si_34093004, si_66857006,
                si_82271004, si_25064002, si_444673007, si_302866003, si_417746004,	
                si_173300003,	si_249307003,	si_212962007,	si_90460009, si_81680005,
                si_55680006, si_12063002, si_289530006, si_80313002, si_180300007,	
                si_385486001,	si_173300003,	si_413307004,	si_262519004,	si_91175000,
                si_127278005,	si_32937002, si_162397003, si_262521009, si_13791008,
                si_249230006,	si_27355003, si_252041008, si_56018004, si_48422000,
                si_315642008,	si_157265008, si_68566005, si_267064002, si_127279002,	
                si_417746004,	si_125600009,	si_34095006, si_10601006, si_102556003,	
                si_406547006,	si_81448000, si_65966004, si_263225007,	si_182888003,	
                si_297217002,	si_278528006,	si_309557009,	si_102570003,	si_77880009,
                si_95673003, si_280816001, si_shorecare, si_38341003, si_34436003,
                si_8510008, si_49436004, si_34801009, si_48694002, si_289195008,	
                si_82991003, si_233604007, si_45007003, si_397706001,	si_23924001,
                si_370380004,	si_90708001, si_95677002, si_22325002, si_48867003,
                si_301120008,	si_418272005,	si_409668002,	si_r40,
                #ED history
                ADU_direct_h, WTBS_period_g_h, SSED_compliant_h, SSED_eligible_h,			
                TTDIP_compliance_h, TTDEM_compliance_h, discharge_from_ED_h,		
                inappropriate_ED_spaces_h, self_discharge_h, LWBS_h,					
                ed_location_path_h_c, ed_specialty_path_h_c, triage_priority_h,
                multiple_EM_clinician_seen_bys_h, six_hr_clock_stop_location_h,		
                ED_planned_return_h, EM_presentation_h, obs_admit_h, ADU_stay_h,						
                bed_booking_compliance_h, return_within_48_hours_h, ward_admission_h,					
                dim_discharge_reason_key_h,
                #Inpatient history
                ED_first_stay_flag_i, inpatient_num_g_i, ip_admission_type_i, 
                event_stays_ward_path_i_c, into_theatre_flag_i, NPF_exception_flag_i, CCL_i,
                #Medchart
                antibiotic_flag_m, opioid_flag_m, infusion_flag_m, diluent_flag_m, medication_count_m,	
                c01_b, c02_b,	c03_b, c04_b,	c05_b, c06_b,	c07_b, c08_b, c09_b, c10_b,				
                c11_b, c12_b,	c13_b, c14_b,	c15_b, c16_b, c17_b, c18_b,	c19_b, c20_b,				
                c21_b, c22_b, c23_b,
                #Testsafe
                Med_Anticoagulant_t, Med_Diuretics_t, Med_Glycaemic_t, Med_Statin_t,				
                Med_BronchodilatorCOPD_t, Med_Opioid_t, Med_Abx_t, Med_Antidepressant_t,		
                Med_Antipsychotic_t, Med_Parkinson_t, Med_Osteoporosis_t, Med_BBlocker_t,
                #ICU
                icu_flag_h,
                #eVitals,
                ews_g, EWS_alert, DBP_ews, sepsis_risk_score_high, pain_at_rest_result,
                sepsis_risk_category, temperature_ews, SBP_ews, HR_ews, respr_ews, O2sat_ews,
                PV_bleeding_result, O2flow_ews, CNS_ews, urine_ews, SBP_ews_contributor, 			
                HR_ews_contributor, respr_ews_contributor, O2sat_ews_contributor,			
                O2flow_ews_contributor, CNS_ews_contributor, urine_ews_contributor, 			
                temperature_ews_contributor, DBP_ews_contributor, pain_at_movement_result,
                dim_eVitals_regime_mod_reason_ke)

