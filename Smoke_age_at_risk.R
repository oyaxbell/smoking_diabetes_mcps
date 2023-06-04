# Smoking and cause-specific mortality in individuals with diabetes in Mexico: an analysis of the Mexico City Prospective Study 
# Data Analysis: Daniel Ramírez García
# Latest version of Analysis June 2023
# For any question regarding analysis contact Omar Yaxmehen Bello-Chavolla at oyaxbell@yahoo.com.mx

#### Libraries ####

library(readr);library(tidyverse);library(lubridate);library(survival)
library(Epi);library(flextable);library(na.tools);library(forestplot)
library(ggpubr);library(ggplotify);library(ggalluvial)

#### Data####
#setwd("~/Proyectos/Data/MCPS_2023")
setwd("~/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/MCPS Projects")
mcps_baseline <- read_csv("Data/2021-004 MCPS BASELINE.csv")
mcps_mort <- read_csv("Data/2022-012 MCPS MORTALITY.csv")
mcps_resampling <- read_csv("Data/2022-012 MCPS RESURVEY.csv")

#Baseline
mcps <- full_join(mcps_baseline, mcps_mort, by = "PATID")

#Resurvey
resurv <- full_join(mcps_baseline, mcps_mort, by = "PATID") %>%
  inner_join(mcps_resampling, by = "PATID")

#### Custom function####
icdselect <- function(..., y, group = FALSE) {
  codigos <- c(...)
  new_var <- 0
  if (group == TRUE) {
    for(i in seq_along(codigos)) {
      new_var[str_detect(y, codigos[i])] <- 1
    }
    new_var
    new_var <- na.replace(new_var, 0)
    return(new_var)
  }
  for(i in seq_along(codigos)) {
    new_var[str_detect(y, codigos[i])] <- i
  }
  new_var
  na.replace(new_var, 0)
}

#### Filters####
#Individuals with diabetes
diab <- filter(mcps, (BASE_DIABETES == 1 |
                        DRUG_D1 == 1 |
                        DRUG_D2 == 1 |
                        DRUG_D3 == 1 |
                        DRUG_D4 == 1 |
                        BASE_HBA1C >= 6.5))

#Exclusion of individuals with comorbidities
diab <- filter(diab, (BASE_EMPHYSEMA == 0 &
                        BASE_HEARTATTACK == 0 &
                        BASE_ANGINA == 0 &
                        BASE_STROKE == 0 &
                        BASE_CKD == 0 &
                        BASE_CIRR == 0 &
                        BASE_LUNGCANCER == 0 &
                        BASE_OTHCANCER == 0 &
                        BASE_PROSTATECANCER == 0 &
                        BASE_CERVCANCER == 0 &
                        BASE_BREASTCANCER == 0 &
                        BASE_STOMCANCER == 0 &
                        BASE_ORALCANCER == 0 &
                        BASE_PAD == 0))

##Excludes those with "U" status and those with accidental or non-defined death.
diab <- diab %>% 
  filter(STATUS == "A" | STATUS == "D") %>%
  mutate(STATUS = ifelse(STATUS == "D", 1, 0)) %>%  #Changes character values for numeric ones.
  filter(D054 != 1) %>% filter(D055 != 1)

#Resurvey filter
diab_resurv <- filter(resurv, ((BASE_DIABETES == 1 |
                                  DRUG_D1 == 1 |
                                  DRUG_D2 == 1 |
                                  DRUG_D3 == 1 |
                                  DRUG_D4 == 1 |
                                  BASE_HBA1C >= 6.5) & #¿Quitar aquellos con diabetes no diagnosticada?
                                 (BASE_EMPHYSEMA == 0 &
                                    BASE_HEARTATTACK == 0 &
                                    BASE_ANGINA == 0 &
                                    BASE_STROKE == 0 &
                                    BASE_CKD == 0 &
                                    BASE_CIRR == 0 &
                                    BASE_LUNGCANCER == 0 &
                                    BASE_OTHCANCER == 0 &
                                    BASE_PROSTATECANCER == 0 &
                                    BASE_CERVCANCER == 0 &
                                    BASE_BREASTCANCER == 0 &
                                    BASE_STOMCANCER == 0 &
                                    BASE_ORALCANCER == 0 &
                                    BASE_PAD == 0)))

#### Variables####
#BMI
diab <- mutate(diab, IMC = (WEIGHT/(HEIGHT^2)*10000))

#Age
diab$edad <- cut(diab$AGE, c(35, 45, 55, 65, 75, 80), right = FALSE)

#Smoking status
diab$cig <- NULL
diab$cig[diab$EVER_SMOK == 0] <- 1
diab$cig[(diab$EVER_SMOK == 1) & (diab$CURR_SMOK == 0)] <- 2
diab$cig[(diab$EVER_SMOK == 1) & (diab$CURR_SMOK == 1)] <- 3

#Smoking status and no. of cigarettes
diab$current <- NULL
diab$current[diab$EVER_SMOK == 0] <- 1
diab$current[diab$CURR_SMOK == 0] <- 2
diab$current[(diab$CURR_SMOK == 1) & (diab$CIGS_PER_DAY >= 1) & (diab$CIGS_PER_DAY < 5)] <- 3
diab$current[(diab$CURR_SMOK == 1) & (diab$CIGS_PER_DAY >= 5) & (diab$CIGS_PER_DAY < 10)] <- 4
diab$current[(diab$CURR_SMOK == 1) & (diab$CIGS_PER_DAY >= 10)] <- 5

diab$edad <- NULL
diab$edad[(diab$AGE >= 35) & (diab$AGE <45)] <- 1
diab$edad[(diab$AGE >= 45) & (diab$AGE <55)] <- 2
diab$edad[(diab$AGE >= 55) & (diab$AGE <65)] <- 3
diab$edad[(diab$AGE >= 65) & (diab$AGE <75)] <- 4
diab$edad[(diab$AGE >= 75)] <- 5

diab$edad2 <- NULL
diab$edad2[(diab$AGE >= 35) & (diab$AGE <55)] <- 1
diab$edad2[(diab$AGE >= 55) & (diab$AGE <65)] <- 2
diab$edad2[(diab$AGE >= 65)] <- 3

#Daily smoking
diab$daily <- NULL
diab$daily[diab$EVER_SMOK == 0] <- 1
diab$daily[(diab$FREQ_SMOK == 1) & (diab$CIGS_PER_DAY >= 1) & (diab$CIGS_PER_DAY < 5)] <- 2
diab$daily[(diab$FREQ_SMOK == 1) & (diab$CIGS_PER_DAY >= 5) & (diab$CIGS_PER_DAY < 10)] <- 3
diab$daily[(diab$FREQ_SMOK == 1) & (diab$CIGS_PER_DAY >= 10)] <- 4

diab$daily1[diab$EVER_SMOK == 0] <- 1
diab$daily1[(diab$FREQ_SMOK == 1)] <- 2

#Non-daily smoking
diab$non_daily[diab$EVER_SMOK == 0] <- 1
diab$non_daily[(diab$FREQ_SMOK == 2 | diab$FREQ_SMOK == 3 | diab$FREQ_SMOK == 4) & (diab$CIGS_PER_DAY >= 1) & (diab$CIGS_PER_DAY < 5)] <- 2
diab$non_daily[(diab$FREQ_SMOK == 2 | diab$FREQ_SMOK == 3 | diab$FREQ_SMOK == 4) & (diab$CIGS_PER_DAY >= 5) & (diab$CIGS_PER_DAY < 10)] <- 3
diab$non_daily[(diab$FREQ_SMOK == 2 | diab$FREQ_SMOK == 3 | diab$FREQ_SMOK == 4) & (diab$CIGS_PER_DAY >= 10)] <- 4

diab$non_daily1[diab$EVER_SMOK == 0] <- 1
diab$non_daily1[(diab$FREQ_SMOK == 2 | diab$FREQ_SMOK == 3 | diab$FREQ_SMOK == 4)] <- 2

#BMI by categories              
diab$imc_1[diab$IMC < 25] <- 1
diab$imc_1[diab$IMC >= 25 & diab$IMC <30] <- 2
diab$imc_1[diab$IMC >= 30 & diab$IMC <35] <- 3
diab$imc_1[diab$IMC >= 35] <- 4

diab$presion <- NULL
diab$presion[(diab$SBP1 <140) & (diab$DBP1 < 90)] <- 1
diab$presion[(diab$SBP1 >= 140) & (diab$DBP1 >= 90)] <- 2

#Glycemic control
diab$control <- NA
diab$control[diab$BASE_HBA1C < 7] <- 1
diab$control[(diab$BASE_HBA1C >= 7) & (diab$BASE_HBA1C <10)] <- 2
diab$control[diab$BASE_HBA1C >= 10] <- 3

#Decade of diagnosis 
diab$DECADE_DX<-NULL
diab$DECADE_DX[diab$BASE_DIABETES_DX==1]<-1950
diab$DECADE_DX[diab$BASE_DIABETES_DX==2]<-1960
diab$DECADE_DX[diab$BASE_DIABETES_DX==3]<-1970
diab$DECADE_DX[diab$BASE_DIABETES_DX==4]<-1980
diab$DECADE_DX[diab$BASE_DIABETES_DX==5]<-1990
diab$DECADE_DX[diab$BASE_DIABETES_DX==6]<-2000

#Time since diagnosis
diab$TM_EVOLUCION <- diab$YEAR_RECRUITED - diab$DECADE_DX
diab$TM_EVOLUCION <- ifelse(is.na(diab$TM_EVOLUCION), 0, diab$TM_EVOLUCION) #Converts to 0 those with undiagnosed diabetes

#Approximate age of diagnosis
diab$EDAD_DX_DM <- (diab$AGE - diab$TM_EVOLUCION) + 5
diab$EDAD_DX_DM <- ifelse(diab$EDAD_DX_DM < 0, 0, diab$EDAD_DX_DM) #Converts to 0 those with negative values on time since dignosis
summary(diab$TM_EVOLUCION)
hist(diab$EDAD_DX_DM)

#Time since diagnosis (categorical)
diab$tiempo_evol[diab$TM_EVOLUCION < 10] <- 1
diab$tiempo_evol[diab$TM_EVOLUCION >= 10] <- 2

diab$tiempo_evol[diab$TM_EVOLUCION < 5] <- 1
diab$tiempo_evol[diab$TM_EVOLUCION >= 5 & diab$TM_EVOLUCION <10] <- 2
diab$tiempo_evol[diab$TM_EVOLUCION >= 10 & diab$TM_EVOLUCION <15] <- 3
diab$tiempo_evol[diab$TM_EVOLUCION >= 15] <- 4

#Undiagnosed diabetes
diab$nodiag <- NULL
diab$nodiag[(diab$BASE_DIABETES == 1)] <- 1
diab$nodiag[(diab$BASE_HBA1C >= 6.5) & (diab$BASE_DIABETES == 0) & (diab$DRUG_D1 == 0) & (diab$DRUG_D2 == 0) & (diab$DRUG_D3 == 0) & (diab$DRUG_D4 == 0)] <- 0

#Age started smoking
diab$edad_inicio <- NULL
diab$edad_inicio[(diab$EVER_SMOK == 0)] <- 1
diab$edad_inicio[diab$AGE_STARTED_SMOKING >= 1 & diab$AGE_STARTED_SMOKING < 15] <- 2
diab$edad_inicio[diab$AGE_STARTED_SMOKING >= 15 & diab$AGE_STARTED_SMOKING < 24] <- 3
diab$edad_inicio[diab$AGE_STARTED_SMOKING >= 24] <- 4

#No. of cigarettes in former smokers
diab$former_cigs <- NULL
diab$former_cigs[diab$EVER_SMOK == 0] <- 1
diab$former_cigs[(diab$CURR_SMOK == 0) & (diab$QUIT_NO_CIGS >= 1) & (diab$QUIT_NO_CIGS < 5)] <- 2
diab$former_cigs[(diab$CURR_SMOK == 0) & (diab$QUIT_NO_CIGS >= 5) & (diab$QUIT_NO_CIGS < 10)] <- 3
diab$former_cigs[(diab$CURR_SMOK == 0) & (diab$QUIT_NO_CIGS >= 10)] <- 4

#Exposure time to tobacco
diab <- mutate(diab, time_exposed = AGE - AGE_STARTED_SMOKING)
diab$time_exposed <- ifelse(diab$CURR_SMOK == 0 & diab$LAST_QUIT_SMOK == 1, diab$time_exposed - 1, diab$time_exposed)
diab$time_exposed <- ifelse(diab$CURR_SMOK == 0 & diab$LAST_QUIT_SMOK == 2, diab$time_exposed - 2, diab$time_exposed)
diab$time_exposed <- ifelse(diab$CURR_SMOK == 0 & diab$LAST_QUIT_SMOK == 3, diab$time_exposed - 5, diab$time_exposed)
diab$time_exposed[diab$EVER_SMOK == 0] <- 0
diab$time_exposed <- ifelse(diab$time_exposed < 0, diab$time_exposed == 0, diab$time_exposed)

#Variables used in resurvey
diab_resurv <- mutate(diab_resurv, IMC = (WEIGHT/(HEIGHT^2)*10000))

diab_resurv$current <- NULL
diab_resurv$current[diab_resurv$EVER_SMOK == 0] <- 1
diab_resurv$current[diab_resurv$CURR_SMOK == 0] <- 2
diab_resurv$current[(diab_resurv$CURR_SMOK == 1) & (diab_resurv$CIGS_PER_DAY >= 1) & (diab_resurv$CIGS_PER_DAY < 5)] <- 3
diab_resurv$current[(diab_resurv$CURR_SMOK == 1) & (diab_resurv$CIGS_PER_DAY >= 5) & (diab_resurv$CIGS_PER_DAY < 10)] <- 4
diab_resurv$current[(diab_resurv$CURR_SMOK == 1) & (diab_resurv$CIGS_PER_DAY >= 10)] <- 5

diab_resurv$cig <- NULL
diab_resurv$cig[diab_resurv$EVER_SMOK == 0] <- 1
diab_resurv$cig[(diab_resurv$EVER_SMOK == 1) & (diab_resurv$CURR_SMOK == 0)] <- 2
diab_resurv$cig[(diab_resurv$EVER_SMOK == 1) & (diab_resurv$CURR_SMOK == 1)] <- 3

diab_resurv$cig_resurv <- NULL
diab_resurv$cig_resurv[diab_resurv$R_EVERSMOK == 0] <- 1
diab_resurv$cig_resurv[diab_resurv$R_EVERSMOK == 1 & diab_resurv$R_CURRSMOK == 0] <- 2
diab_resurv$cig_resurv[diab_resurv$R_EVERSMOK == 1 & diab_resurv$R_CURRSMOK == 1] <- 3

diab_resurv$nodiag <- NULL
diab_resurv$nodiag[(diab_resurv$BASE_DIABETES == 1)] <- 1
diab_resurv$nodiag[(diab_resurv$BASE_HBA1C >= 6.5) & (diab_resurv$BASE_DIABETES == 0) & (diab_resurv$DRUG_D1 == 0) & (diab_resurv$DRUG_D2 == 0) & (diab_resurv$DRUG_D3 == 0) & (diab_resurv$DRUG_D4 == 0)] <- 0

#Individuals with complete covariate data at baseline
diab <- diab %>% filter(if_all(c(AGE, MALE, cig, current, BASE_HBA1C, EDU_LEVEL, COYOACAN, SBP1, IMC, TM_EVOLUCION, ALCGP), ~ !is.na(.x)))

#Individuals with complete covariate data in resurvey
diab_resurv <- diab_resurv %>% filter(if_all(c(AGE, MALE, cig, current, BASE_HBA1C, EDU_LEVEL, COYOACAN, SBP1, IMC, ALCGP), ~ !is.na(.x)))

#Variables used in Cox model for age as time-scale
diab$a0 <- diab$AGE
diab$a1 <- diab$AGE + diab$PERSON_YEARS

#Datasets by sex
men <- diab %>% filter(MALE == 1)
women <- diab %>% filter(MALE == 0)

#### Descriptive analysis####
mean <- list(
  mean = ~ round(mean(.x, na.rm = TRUE),2),
  sd = ~ round(sd(.x, na.rm = TRUE),2))

median <- list(
  median = ~ round(median(.x, na.rm = TRUE),2),
  q1 = ~ round(quantile(.x, probs = c(.25), na.rm = TRUE),2),
  q3 = ~ round(quantile(.x, probs = c(.75), na.rm = TRUE),2)
)

#Men
num_men <- men %>% count(cig) %>% 
  mutate(prct = round(n / sum(n)*100, 2))

fu_men <- median(men$PERSON_YEARS)

var_continuas_mean_m <- men %>% 
  group_by(cig) %>% 
  summarize(across(c(AGE, WEIGHT, HEIGHT, IMC, SBP1, DBP1), mean))

var_continuas_median_m <- men %>% 
  group_by(cig) %>% 
  summarize(across(c(PERSON_YEARS, AGE_STARTED_SMOKING, CIGS_PER_DAY, BASE_HBA1C, QUIT_NO_CIGS), median)) %>% 
  data.frame()

coyoacan_m <- men %>% 
  group_by(cig, COYOACAN) %>% 
  summarize(n = n()) %>% 
  mutate(freq = round((n / sum(n)*100),2)) %>% 
  pivot_wider(names_from = COYOACAN, values_from = c(n, freq), names_prefix = "COYOACAN") %>% 
  relocate(freq_COYOACAN0, .after = n_COYOACAN0)

educacion_m <- men %>% 
  group_by(cig, EDU_LEVEL) %>% 
  summarize(n = n()) %>% 
  mutate(freq = round((n / sum(n)*100),2)) %>% 
  pivot_wider(names_from = EDU_LEVEL, values_from = c(n, freq), names_prefix = "EDU_LEVEL") %>% 
  relocate(freq_EDU_LEVEL1, .after = n_EDU_LEVEL1) %>% 
  relocate(freq_EDU_LEVEL2, .after = n_EDU_LEVEL2) %>% 
  relocate(freq_EDU_LEVEL3, .after = n_EDU_LEVEL3) %>% 
  relocate(freq_EDU_LEVEL4, .after = n_EDU_LEVEL4)

hba1c_m <- men %>% 
  group_by(cig, control) %>% 
  summarize(n = n()) %>% 
  mutate(freq = round((n / sum(n)*100),2)) %>%
  pivot_wider(names_from = control, values_from = c(n, freq), names_prefix = "control") %>% 
  relocate(freq_control1, .after = n_control1) %>% 
  relocate(freq_control2, .after = n_control2) %>% 
  relocate(freq_control3, .after = n_control3)

biguanides_m <- men %>% 
  group_by(cig, DRUG_D1) %>% 
  summarize(n = n()) %>% 
  mutate(freq = round((n / sum(n)*100),2)) %>% 
  pivot_wider(names_from = DRUG_D1, values_from = c(n, freq), names_prefix = "Biguanide") %>% 
  relocate(freq_Biguanide0, .after = n_Biguanide0)

sulf_m <- men %>% 
  group_by(cig, DRUG_D2) %>% 
  summarize(n = n()) %>% 
  mutate(freq = round((n / sum(n)*100),2)) %>% 
  pivot_wider(names_from = DRUG_D2, values_from = c(n, freq), names_prefix = "Sulf") %>% 
  relocate(freq_Sulf0, .after = n_Sulf0)

insulin_m <- men %>% 
  group_by(cig, DRUG_D3) %>% 
  summarize(n = n()) %>% 
  mutate(freq = round((n / sum(n)*100),2)) %>% 
  pivot_wider(names_from = DRUG_D3, values_from = c(n, freq), names_prefix = "Insulin") %>% 
  relocate(freq_Insulin0, .after = n_Insulin0)

#Women
num_women <- women %>% count(cig) %>% 
  mutate(prct = round(n / sum(n)*100, 2))

fu_women <- median(women$PERSON_YEARS)

var_continuas_mean_w <- women %>% 
  group_by(cig) %>% 
  summarize(across(c(AGE, WEIGHT, HEIGHT, IMC, SBP1, DBP1), mean))

var_continuas_median_w <- women %>% 
  group_by(cig) %>% 
  summarize(across(c(PERSON_YEARS, AGE_STARTED_SMOKING, CIGS_PER_DAY, BASE_HBA1C, QUIT_NO_CIGS), median)) %>% 
  data.frame()

coyoacan_w <- women %>% 
  group_by(cig, COYOACAN) %>% 
  summarize(n = n()) %>% 
  mutate(freq = round((n / sum(n)*100),2)) %>% 
  pivot_wider(names_from = COYOACAN, values_from = c(n, freq), names_prefix = "COYOACAN") %>% 
  relocate(freq_COYOACAN0, .after = n_COYOACAN0)

educacion_w <- women %>% 
  group_by(cig, EDU_LEVEL) %>% 
  summarize(n = n()) %>% 
  mutate(freq = round((n / sum(n)*100),2)) %>% 
  pivot_wider(names_from = EDU_LEVEL, values_from = c(n, freq), names_prefix = "EDU_LEVEL") %>% 
  relocate(freq_EDU_LEVEL1, .after = n_EDU_LEVEL1) %>% 
  relocate(freq_EDU_LEVEL2, .after = n_EDU_LEVEL2) %>% 
  relocate(freq_EDU_LEVEL3, .after = n_EDU_LEVEL3) %>% 
  relocate(freq_EDU_LEVEL4, .after = n_EDU_LEVEL4)

hba1c_w <- women %>% 
  group_by(cig, control) %>% 
  summarize(n = n()) %>% 
  mutate(freq = round((n / sum(n)*100),2)) %>%
  pivot_wider(names_from = control, values_from = c(n, freq), names_prefix = "control") %>% 
  relocate(freq_control1, .after = n_control1) %>% 
  relocate(freq_control2, .after = n_control2) %>% 
  relocate(freq_control3, .after = n_control3)

biguanides_w <- women %>% 
  group_by(cig, DRUG_D1) %>% 
  summarize(n = n()) %>% 
  mutate(freq = round((n / sum(n)*100),2)) %>% 
  pivot_wider(names_from = DRUG_D1, values_from = c(n, freq), names_prefix = "Biguanide") %>% 
  relocate(freq_Biguanide0, .after = n_Biguanide0)

sulf_w <- women %>% 
  group_by(cig, DRUG_D2) %>% 
  summarize(n = n()) %>% 
  mutate(freq = round((n / sum(n)*100),2)) %>% 
  pivot_wider(names_from = DRUG_D2, values_from = c(n, freq), names_prefix = "Sulf") %>% 
  relocate(freq_Sulf0, .after = n_Sulf0)

insulin_w <- women %>% 
  group_by(cig, DRUG_D3) %>% 
  summarize(n = n()) %>% 
  mutate(freq = round((n / sum(n)*100),2)) %>% 
  pivot_wider(names_from = DRUG_D3, values_from = c(n, freq), names_prefix = "Insulin") %>% 
  relocate(freq_Insulin0, .after = n_Insulin0)

#Descriptive table
table <- data.frame("Characteristic" = c("No. of participants", "Mean age", "Coyoacán (%)", "Univerity (%)", "High school (%)", 
                                         "Elementary (%)", "Other (%)", "Age started smoking", "No. of cigarettes/d", "Median HbA1C (%)",
                                         "<7 HbA1c (%)", "7-10 HbA1C (%)", ">10 HbA1c (%)", "Weight (kg)", "Height (cm)","BMI (kg/m2)", 
                                         "Systolic Pressure (mmHg)", "Diastolic Pressure (mmHg)", "Biguanides (%)", "Sulfonylureas (%)", "Insulin (%)"),
                    "W_Never smokers" = c(paste0(num_women[1, 2], " (", num_women[1, 3], "%)"), 
                                          paste0(var_continuas_mean_w[1,2], " (±", var_continuas_mean_w[1,3], ")"), 
                                          paste0(coyoacan_w[1,4], " (", coyoacan_w[1,5], "%)"),
                                          paste0(educacion_w[1, 2], " (", educacion_w[1, 3], "%)"),
                                          paste0(educacion_w[1, 4], " (", educacion_w[1, 5], "%)"),
                                          paste0(educacion_w[1, 6], " (", educacion_w[1, 7], "%)"),
                                          paste0(educacion_w[1, 8], " (", educacion_w[1, 9], "%)"),
                                          paste0(var_continuas_median_w[1,5], " (", var_continuas_median_w[1,6], ", ", var_continuas_median_w[1,7], ")"),
                                          paste0(var_continuas_median_w[1,8], " (", var_continuas_median_w[1,9], ", ", var_continuas_median_w[1,10], ")"),
                                          paste0(var_continuas_median_w[1,11], " (", var_continuas_median_w[1,12], ", ", var_continuas_median_w[1,13], ")"),
                                          paste0(hba1c_w[1,2], " (", hba1c_w[1,3], "%)"),
                                          paste0(hba1c_w[1,4], " (", hba1c_w[1,5], "%)"),
                                          paste0(hba1c_w[1,6], " (", hba1c_w[1,7], "%)"),
                                          paste0(var_continuas_mean_w[1,4], " (±", var_continuas_mean_w[1,5], ")"),
                                          paste0(var_continuas_mean_w[1,6], " (±", var_continuas_mean_w[1,7], ")"),
                                          paste0(var_continuas_mean_w[1,8], " (±", var_continuas_mean_w[1,9], ")"),
                                          paste0(var_continuas_mean_w[1,10], " (±", var_continuas_mean_w[1,11], ")"),
                                          paste0(var_continuas_mean_w[1,12], " (±", var_continuas_mean_w[1,13], ")"),
                                          paste0(biguanides_w[1,4], " (", biguanides_w[1,5], "%)"),
                                          paste0(sulf_w[1,4], " (", sulf_w[1,5], "%)"), 
                                          paste0(insulin_w[1,4], " (", insulin_w[1,5], "%)")),
                    "W_Former smokers" = c(paste0(num_women[2, 2], " (", num_women[2, 3], "%)"), 
                                           paste0(var_continuas_mean_w[2,2], " (±", var_continuas_mean_w[2,3], ")"), 
                                           paste0(coyoacan_w[2,4], " (", coyoacan_w[2,5], "%)"),
                                           paste0(educacion_w[2, 2], " (", educacion_w[2, 3], "%)"),
                                           paste0(educacion_w[2, 4], " (", educacion_w[2, 5], "%)"),
                                           paste0(educacion_w[2, 6], " (", educacion_w[2, 7], "%)"),
                                           paste0(educacion_w[2, 8], " (", educacion_w[2, 9], "%)"),
                                           paste0(var_continuas_median_w[2,5], " (", var_continuas_median_w[2,6], ", ", var_continuas_median_w[2,7], ")"),
                                           paste0(var_continuas_median_w[2,8], " (", var_continuas_median_w[2,9], ", ", var_continuas_median_w[2,10], ")"),
                                           paste0(var_continuas_median_w[2,11], " (", var_continuas_median_w[2,12], ", ", var_continuas_median_w[2,13], ")"),
                                           paste0(hba1c_w[2,2], " (", hba1c_w[2,3], "%)"),
                                           paste0(hba1c_w[2,4], " (", hba1c_w[2,5], "%)"),
                                           paste0(hba1c_w[2,6], " (", hba1c_w[2,7], "%)"),
                                           paste0(var_continuas_mean_w[2,4], " (±", var_continuas_mean_w[2,5], ")"),
                                           paste0(var_continuas_mean_w[2,6], " (±", var_continuas_mean_w[2,7], ")"),
                                           paste0(var_continuas_mean_w[2,8], " (±", var_continuas_mean_w[2,9], ")"),
                                           paste0(var_continuas_mean_w[2,10], " (±", var_continuas_mean_w[2,11], ")"),
                                           paste0(var_continuas_mean_w[2,12], " (±", var_continuas_mean_w[2,13], ")"),
                                           paste0(biguanides_w[2,4], " (", biguanides_w[2,5], "%)"),
                                           paste0(sulf_w[2,4], " (", sulf_w[2,5], "%)"), 
                                           paste0(insulin_w[2,4], " (", insulin_w[2,5], "%)")),
                    "W_Current smokers" = c(paste0(num_women[3, 2], " (", num_women[3, 3], "%)"), 
                                            paste0(var_continuas_mean_w[3,2], " (±", var_continuas_mean_w[3,3], ")"), 
                                            paste0(coyoacan_w[3,4], " (", coyoacan_w[3,5], "%)"),
                                            paste0(educacion_w[3, 2], " (", educacion_w[3, 3], "%)"),
                                            paste0(educacion_w[3, 4], " (", educacion_w[3, 5], "%)"),
                                            paste0(educacion_w[3, 6], " (", educacion_w[3, 7], "%)"),
                                            paste0(educacion_w[3, 8], " (", educacion_w[3, 9], "%)"),
                                            paste0(var_continuas_median_w[3,5], " (", var_continuas_median_w[3,6], ", ", var_continuas_median_w[3,7], ")"),
                                            paste0(var_continuas_median_w[3,8], " (", var_continuas_median_w[3,9], ", ", var_continuas_median_w[3,10], ")"),
                                            paste0(var_continuas_median_w[3,11], " (", var_continuas_median_w[3,12], ", ", var_continuas_median_w[3,13], ")"),
                                            paste0(hba1c_w[3,2], " (", hba1c_w[3,3], "%)"),
                                            paste0(hba1c_w[3,4], " (", hba1c_w[3,5], "%)"),
                                            paste0(hba1c_w[3,6], " (", hba1c_w[3,7], "%)"),
                                            paste0(var_continuas_mean_w[3,4], " (±", var_continuas_mean_w[3,5], ")"),
                                            paste0(var_continuas_mean_w[3,6], " (±", var_continuas_mean_w[3,7], ")"),
                                            paste0(var_continuas_mean_w[3,8], " (±", var_continuas_mean_w[3,9], ")"),
                                            paste0(var_continuas_mean_w[3,10], " (±", var_continuas_mean_w[3,11], ")"),
                                            paste0(var_continuas_mean_w[3,12], " (±", var_continuas_mean_w[3,13], ")"),
                                            paste0(biguanides_w[3,4], " (", biguanides_w[3,5], "%)"),
                                            paste0(sulf_w[3,4], " (", sulf_w[3,5], "%)"), 
                                            paste0(insulin_w[3,4], " (", insulin_w[3,5], "%)")),
                    "M_Never smokers" = c(paste0(num_men[1, 2], " (", num_men[1, 3], "%)"), 
                                          paste0(var_continuas_mean_m[1,2], " (±", var_continuas_mean_m[1,3], ")"), 
                                          paste0(coyoacan_m[1,4], " (", coyoacan_m[1,5], "%)"),
                                          paste0(educacion_m[1, 2], " (", educacion_m[1, 3], "%)"),
                                          paste0(educacion_m[1, 4], " (", educacion_m[1, 5], "%)"),
                                          paste0(educacion_m[1, 6], " (", educacion_m[1, 7], "%)"),
                                          paste0(educacion_m[1, 8], " (", educacion_m[1, 9], "%)"),
                                          paste0(var_continuas_median_m[1,5], " (", var_continuas_median_m[1,6], ", ", var_continuas_median_m[1,7], ")"),
                                          paste0(var_continuas_median_m[1,8], " (", var_continuas_median_m[1,9], ", ", var_continuas_median_m[1,10], ")"),
                                          paste0(var_continuas_median_m[1,11], " (", var_continuas_median_m[1,12], ", ", var_continuas_median_m[1,13], ")"),
                                          paste0(hba1c_m[1,2], " (", hba1c_m[1,3], "%)"),
                                          paste0(hba1c_m[1,4], " (", hba1c_m[1,5], "%)"),
                                          paste0(hba1c_m[1,6], " (", hba1c_m[1,7], "%)"),
                                          paste0(var_continuas_mean_m[1,4], " (±", var_continuas_mean_m[1,5], ")"),
                                          paste0(var_continuas_mean_m[1,6], " (±", var_continuas_mean_m[1,7], ")"),
                                          paste0(var_continuas_mean_m[1,8], " (±", var_continuas_mean_m[1,9], ")"),
                                          paste0(var_continuas_mean_m[1,10], " (±", var_continuas_mean_m[1,11], ")"),
                                          paste0(var_continuas_mean_m[1,12], " (±", var_continuas_mean_m[1,13], ")"),
                                          paste0(biguanides_m[1,4], " (", biguanides_m[1,5], "%)"),
                                          paste0(sulf_m[1,4], " (", sulf_m[1,5], "%)"), 
                                          paste0(insulin_m[1,4], " (", insulin_m[1,5], "%)")),
                    "M_Former smokers" = c(paste0(num_men[2, 2], " (", num_men[2, 3], "%)"), 
                                           paste0(var_continuas_mean_m[2,2], " (±", var_continuas_mean_m[2,3], ")"), 
                                           paste0(coyoacan_m[2,4], " (", coyoacan_m[2,5], "%)"),
                                           paste0(educacion_m[2, 2], " (", educacion_m[2, 3], "%)"),
                                           paste0(educacion_m[2, 4], " (", educacion_m[2, 5], "%)"),
                                           paste0(educacion_m[2, 6], " (", educacion_m[2, 7], "%)"),
                                           paste0(educacion_m[2, 8], " (", educacion_m[2, 9], "%)"),
                                           paste0(var_continuas_median_m[2,5], " (", var_continuas_median_m[2,6], ", ", var_continuas_median_m[2,7], ")"),
                                           paste0(var_continuas_median_m[2,8], " (", var_continuas_median_m[2,9], ", ", var_continuas_median_m[2,10], ")"),
                                           paste0(var_continuas_median_m[2,11], " (", var_continuas_median_m[2,12], ", ", var_continuas_median_m[2,13], ")"),
                                           paste0(hba1c_m[2,2], " (", hba1c_m[2,3], "%)"),
                                           paste0(hba1c_m[2,4], " (", hba1c_m[2,5], "%)"),
                                           paste0(hba1c_m[2,6], " (", hba1c_m[2,7], "%)"),
                                           paste0(var_continuas_mean_m[2,4], " (±", var_continuas_mean_m[2,5], ")"),
                                           paste0(var_continuas_mean_m[2,6], " (±", var_continuas_mean_m[2,7], ")"),
                                           paste0(var_continuas_mean_m[2,8], " (±", var_continuas_mean_m[2,9], ")"),
                                           paste0(var_continuas_mean_m[2,10], " (±", var_continuas_mean_m[2,11], ")"),
                                           paste0(var_continuas_mean_m[2,12], " (±", var_continuas_mean_m[2,13], ")"),
                                           paste0(biguanides_m[2,4], " (", biguanides_m[2,5], "%)"),
                                           paste0(sulf_m[2,4], " (", sulf_m[2,5], "%)"), 
                                           paste0(insulin_m[2,4], " (", insulin_m[2,5], "%)")),
                    "M_Current smokers" = c(paste0(num_men[3, 2], " (", num_men[3, 3], "%)"), 
                                            paste0(var_continuas_mean_m[3,2], " (±", var_continuas_mean_m[3,3], ")"), 
                                            paste0(coyoacan_m[3,4], " (", coyoacan_m[3,5], "%)"),
                                            paste0(educacion_m[3, 2], " (", educacion_m[3, 3], "%)"),
                                            paste0(educacion_m[3, 4], " (", educacion_m[3, 5], "%)"),
                                            paste0(educacion_m[3, 6], " (", educacion_m[3, 7], "%)"),
                                            paste0(educacion_m[3, 8], " (", educacion_m[3, 9], "%)"),
                                            paste0(var_continuas_median_m[3,5], " (", var_continuas_median_m[3,6], ", ", var_continuas_median_m[3,7], ")"),
                                            paste0(var_continuas_median_m[3,8], " (", var_continuas_median_m[3,9], ", ", var_continuas_median_m[3,10], ")"),
                                            paste0(var_continuas_median_m[3,11], " (", var_continuas_median_m[3,12], ", ", var_continuas_median_m[3,13], ")"),
                                            paste0(hba1c_m[3,2], " (", hba1c_m[3,3], "%)"),
                                            paste0(hba1c_m[3,4], " (", hba1c_m[3,5], "%)"),
                                            paste0(hba1c_m[3,6], " (", hba1c_m[3,7], "%)"),
                                            paste0(var_continuas_mean_m[3,4], " (±", var_continuas_mean_m[3,5], ")"),
                                            paste0(var_continuas_mean_m[3,6], " (±", var_continuas_mean_m[3,7], ")"),
                                            paste0(var_continuas_mean_m[3,8], " (±", var_continuas_mean_m[3,9], ")"),
                                            paste0(var_continuas_mean_m[3,10], " (±", var_continuas_mean_m[3,11], ")"),
                                            paste0(var_continuas_mean_m[3,12], " (±", var_continuas_mean_m[3,13], ")"),
                                            paste0(biguanides_m[3,4], " (", biguanides_m[3,5], "%)"),
                                            paste0(sulf_m[3,4], " (", sulf_m[3,5], "%)"), 
                                            paste0(insulin_m[3,4], " (", insulin_m[3,5], "%)")))

tab<-`names<-`(table,c("Characteristics","Never smokers","Former smokers","Current smokers","Never smokers ","Former smokers ","Current smokers "))
tab <- align(flextable(tab, cwidth = 1.5), align = "center", part = "body") 
tab <- add_header_row(tab, values = c(NA, "Women", "Men"), colwidths = c(2, 3, 2))
save_as_docx(tab, path = "~/Proyectos/Smoking and diabetes/Figures strata/table1.docx")

#### Lexis expansion for all-cause mortality####
diab_lexis <- diab %>% 
  mutate(date_recruited = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))

diab_lexis <- Lexis(entry = list("period" = date_recruited,
                                   "age" = AGE,
                                   "time_at_entry" = 0),
                      exit = list("period" = date_recruited + PERSON_YEARS),
                      exit.status = STATUS,
                      data = diab_lexis)

diab_lexis <- splitLexis(diab_lexis, "age", breaks = seq(0, 85, 5))
diab_lexis$age_band <- timeBand(diab_lexis, "age", type = "factor")

diab_lexis <- diab_lexis %>% 
  rename(id = lex.id,
         fu_period = lex.dur,
         entry_status = lex.Cst,
         exit_status = lex.Xst) %>% 
  mutate(time_at_exit = time_at_entry + fu_period,
         age_at_exit = age + fu_period) %>% 
  relocate(time_at_exit, .after = time_at_entry) %>% 
  relocate(age_band, .after = period) %>% 
  relocate(age_at_exit, .after = age)

diab_lexis$premature <- ifelse(diab_lexis$exit_status == 1, 1, 0)
diab_lexis$premature[diab_lexis$premature == 1 & diab_lexis$age_at_exit >= 75] <- 0

#Number of premature deaths 
diab_lexis %>% count(premature)

#Number of premature deaths per sex
diab_lexis %>% count(MALE, premature)

#### Selection of ICD-10 codes for cause-specific mortality analyses####
#Cardiovascular mortality
diab$cardiovascular <- icdselect("I21", "I22", "I23", "I24", "I25", "I6", "I26", "I80", "I82", y = diab$ICD10_UNDERLYING, group = TRUE)

#Myocardial infarction
diab$infarto <- icdselect("I20", "I21", "I22", "I23", "I24", "I25", y = diab$ICD10_UNDERLYING, group = TRUE)

#Stroke
diab$stroke <- icdselect("I6", y = diab$ICD10_UNDERLYING, group = TRUE)

#Diabetes
diab$diabetes <- icdselect("E10", "E11", "E12", "E13", "E14", y = diab$ICD10_UNDERLYING, group = TRUE)

diab$comp_agudas <- icdselect("E110", "E111", y = diab$ICD10_UNDERLYING, group = TRUE)

diab$comp_cronicas <- icdselect("E112", "E113", "E114", "E115", "E116", "E117", "E118", y = diab$ICD10_UNDERLYING, group = TRUE)

#Lung cancer
diab$ca_pulmon <- icdselect("C34", y = diab$ICD10_UNDERLYING)

#Non-lung cancer
diab$ca_otros <- icdselect("C0", "C1", "C2", "C30", "C31", "C32", "C33", "C37", "C38", "C39", "C4", "C5", "C6", "C7", "C8", "C9", y = diab$ICD10_UNDERLYING, group = TRUE)

#COPD
diab$epoc <- icdselect("J40", "J41", "J42", "J43", "J44", y = diab$ICD10_UNDERLYING, group = TRUE)

#Kidney mortality
diab$renal <- icdselect("E112", "N18", "N19", y = diab$ICD10_UNDERLYING, group = TRUE)

#### Lexis expansion for cause-specific mortality analyses####
#Cardiovascular mortality
cv_df <- diab %>% 
  mutate(date_recruited = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))

cv_df <- Lexis(entry = list("period" = date_recruited,
                                 "age" = AGE,
                                 "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS),
                    exit.status = cardiovascular,
                    data = cv_df)

cv_df <- splitLexis(cv_df, "age", breaks = seq(0, 85, 5))
cv_df$age_band <- timeBand(cv_df, "age", type = "factor")

cv_df <- cv_df %>% 
  rename(id = lex.id,
         fu_period = lex.dur,
         entry_status = lex.Cst,
         exit_status = lex.Xst) %>% 
  mutate(time_at_exit = time_at_entry + fu_period,
         age_at_exit = age + fu_period) %>% 
  relocate(time_at_exit, .after = time_at_entry) %>% 
  relocate(age_band, .after = period) %>% 
  relocate(age_at_exit, .after = age)

cv_df$premature <- ifelse(cv_df$exit_status == 1, 1, 0)
cv_df$premature[cv_df$premature == 1 & cv_df$age_at_exit >= 75] <- 0

#Myocardial infarction
mi_df <- diab %>% 
  mutate(date_recruited = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))

mi_df <- Lexis(entry = list("period" = date_recruited,
                            "age" = AGE,
                            "time_at_entry" = 0),
               exit = list("period" = date_recruited + PERSON_YEARS),
               exit.status = infarto,
               data = mi_df)

mi_df <- splitLexis(mi_df, "age", breaks = seq(0, 85, 5))
mi_df$age_band <- timeBand(mi_df, "age", type = "factor")

mi_df <- mi_df %>% 
  rename(id = lex.id,
         fu_period = lex.dur,
         entry_status = lex.Cst,
         exit_status = lex.Xst) %>% 
  mutate(time_at_exit = time_at_entry + fu_period,
         age_at_exit = age + fu_period) %>% 
  relocate(time_at_exit, .after = time_at_entry) %>% 
  relocate(age_band, .after = period) %>% 
  relocate(age_at_exit, .after = age)

mi_df$premature <- ifelse(mi_df$exit_status == 1, 1, 0)
mi_df$premature[mi_df$premature == 1 & mi_df$age_at_exit >= 75] <- 0

#Stroke
stroke_df <- diab %>% 
  mutate(date_recruited = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))

stroke_df <- Lexis(entry = list("period" = date_recruited,
                            "age" = AGE,
                            "time_at_entry" = 0),
               exit = list("period" = date_recruited + PERSON_YEARS),
               exit.status = stroke,
               data = stroke_df)

stroke_df <- splitLexis(stroke_df, "age", breaks = seq(0, 85, 5))
stroke_df$age_band <- timeBand(stroke_df, "age", type = "factor")

stroke_df <- stroke_df %>% 
  rename(id = lex.id,
         fu_period = lex.dur,
         entry_status = lex.Cst,
         exit_status = lex.Xst) %>% 
  mutate(time_at_exit = time_at_entry + fu_period,
         age_at_exit = age + fu_period) %>% 
  relocate(time_at_exit, .after = time_at_entry) %>% 
  relocate(age_band, .after = period) %>% 
  relocate(age_at_exit, .after = age)

stroke_df$premature <- ifelse(stroke_df$exit_status == 1, 1, 0)
stroke_df$premature[stroke_df$premature == 1 & stroke_df$age_at_exit >= 75] <- 0

#Diabetes
diabetes_df <- diab %>% 
  mutate(date_recruited = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))

diabetes_df <- Lexis(entry = list("period" = date_recruited,
                            "age" = AGE,
                            "time_at_entry" = 0),
               exit = list("period" = date_recruited + PERSON_YEARS),
               exit.status = diabetes,
               data = diabetes_df)

diabetes_df <- splitLexis(diabetes_df, "age", breaks = seq(0, 85, 5))
diabetes_df$age_band <- timeBand(diabetes_df, "age", type = "factor")

diabetes_df <- diabetes_df %>% 
  rename(id = lex.id,
         fu_period = lex.dur,
         entry_status = lex.Cst,
         exit_status = lex.Xst) %>% 
  mutate(time_at_exit = time_at_entry + fu_period,
         age_at_exit = age + fu_period) %>% 
  relocate(time_at_exit, .after = time_at_entry) %>% 
  relocate(age_band, .after = period) %>% 
  relocate(age_at_exit, .after = age)

diabetes_df$premature <- ifelse(diabetes_df$exit_status == 1, 1, 0)
diabetes_df$premature[diabetes_df$premature == 1 & diabetes_df$age_at_exit >= 75] <- 0

#Acute diabetes
ac_diabetes_df <- diab %>% 
  mutate(date_recruited = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))

ac_diabetes_df <- Lexis(entry = list("period" = date_recruited,
                            "age" = AGE,
                            "time_at_entry" = 0),
               exit = list("period" = date_recruited + PERSON_YEARS),
               exit.status = comp_agudas,
               data = ac_diabetes_df)

ac_diabetes_df <- splitLexis(ac_diabetes_df, "age", breaks = seq(0, 85, 5))
ac_diabetes_df$age_band <- timeBand(ac_diabetes_df, "age", type = "factor")

ac_diabetes_df <- ac_diabetes_df %>% 
  rename(id = lex.id,
         fu_period = lex.dur,
         entry_status = lex.Cst,
         exit_status = lex.Xst) %>% 
  mutate(time_at_exit = time_at_entry + fu_period,
         age_at_exit = age + fu_period) %>% 
  relocate(time_at_exit, .after = time_at_entry) %>% 
  relocate(age_band, .after = period) %>% 
  relocate(age_at_exit, .after = age)

ac_diabetes_df$premature <- ifelse(ac_diabetes_df$exit_status == 1, 1, 0)
ac_diabetes_df$premature[ac_diabetes_df$premature == 1 & ac_diabetes_df$age_at_exit >= 75] <- 0

#Chronic diabetes
chr_diabetes_df <- diab %>% 
  mutate(date_recruited = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))

chr_diabetes_df <- Lexis(entry = list("period" = date_recruited,
                            "age" = AGE,
                            "time_at_entry" = 0),
               exit = list("period" = date_recruited + PERSON_YEARS),
               exit.status = comp_cronicas,
               data = chr_diabetes_df)

chr_diabetes_df <- splitLexis(chr_diabetes_df, "age", breaks = seq(0, 85, 5))
chr_diabetes_df$age_band <- timeBand(chr_diabetes_df, "age", type = "factor")

chr_diabetes_df <- chr_diabetes_df %>% 
  rename(id = lex.id,
         fu_period = lex.dur,
         entry_status = lex.Cst,
         exit_status = lex.Xst) %>% 
  mutate(time_at_exit = time_at_entry + fu_period,
         age_at_exit = age + fu_period) %>% 
  relocate(time_at_exit, .after = time_at_entry) %>% 
  relocate(age_band, .after = period) %>% 
  relocate(age_at_exit, .after = age)

chr_diabetes_df$premature <- ifelse(chr_diabetes_df$exit_status == 1, 1, 0)
chr_diabetes_df$premature[chr_diabetes_df$premature == 1 & chr_diabetes_df$age_at_exit >= 75] <- 0

#Lung cancer
lung_df <- diab %>% 
  mutate(date_recruited = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))

lung_df <- Lexis(entry = list("period" = date_recruited,
                            "age" = AGE,
                            "time_at_entry" = 0),
               exit = list("period" = date_recruited + PERSON_YEARS),
               exit.status = ca_pulmon,
               data = lung_df)

lung_df <- splitLexis(lung_df, "age", breaks = seq(0, 85, 5))
lung_df$age_band <- timeBand(lung_df, "age", type = "factor")

lung_df <- lung_df %>% 
  rename(id = lex.id,
         fu_period = lex.dur,
         entry_status = lex.Cst,
         exit_status = lex.Xst) %>% 
  mutate(time_at_exit = time_at_entry + fu_period,
         age_at_exit = age + fu_period) %>% 
  relocate(time_at_exit, .after = time_at_entry) %>% 
  relocate(age_band, .after = period) %>% 
  relocate(age_at_exit, .after = age)

lung_df$premature <- ifelse(lung_df$exit_status == 1, 1, 0)
lung_df$premature[lung_df$premature == 1 & lung_df$age_at_exit >= 75] <- 0

#Non-lung cancer
non_lung_df <- diab %>% 
  mutate(date_recruited = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))

non_lung_df <- Lexis(entry = list("period" = date_recruited,
                            "age" = AGE,
                            "time_at_entry" = 0),
               exit = list("period" = date_recruited + PERSON_YEARS),
               exit.status = ca_otros,
               data = non_lung_df)

non_lung_df <- splitLexis(non_lung_df, "age", breaks = seq(0, 85, 5))
non_lung_df$age_band <- timeBand(non_lung_df, "age", type = "factor")

non_lung_df <- non_lung_df %>% 
  rename(id = lex.id,
         fu_period = lex.dur,
         entry_status = lex.Cst,
         exit_status = lex.Xst) %>% 
  mutate(time_at_exit = time_at_entry + fu_period,
         age_at_exit = age + fu_period) %>% 
  relocate(time_at_exit, .after = time_at_entry) %>% 
  relocate(age_band, .after = period) %>% 
  relocate(age_at_exit, .after = age)

non_lung_df$premature <- ifelse(non_lung_df$exit_status == 1, 1, 0)
non_lung_df$premature[non_lung_df$premature == 1 & non_lung_df$age_at_exit >= 75] <- 0

#COPD
copd_df <- diab %>% 
  mutate(date_recruited = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))

copd_df <- Lexis(entry = list("period" = date_recruited,
                            "age" = AGE,
                            "time_at_entry" = 0),
               exit = list("period" = date_recruited + PERSON_YEARS),
               exit.status = epoc,
               data = copd_df)

copd_df <- splitLexis(copd_df, "age", breaks = seq(0, 85, 5))
copd_df$age_band <- timeBand(copd_df, "age", type = "factor")

copd_df <- copd_df %>% 
  rename(id = lex.id,
         fu_period = lex.dur,
         entry_status = lex.Cst,
         exit_status = lex.Xst) %>% 
  mutate(time_at_exit = time_at_entry + fu_period,
         age_at_exit = age + fu_period) %>% 
  relocate(time_at_exit, .after = time_at_entry) %>% 
  relocate(age_band, .after = period) %>% 
  relocate(age_at_exit, .after = age)

copd_df$premature <- ifelse(copd_df$exit_status == 1, 1, 0)
copd_df$premature[copd_df$premature == 1 & copd_df$age_at_exit >= 75] <- 0

#### All-cause mortality analyses####
#Full model stratified by sex
model1 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(cig) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = diab_lexis)
summary(model1)

test1 <- cox.zph(model1)
plot(test1[1], resid = FALSE)
abline(0,0, lty = 3)

model2 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(current) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = diab_lexis)
summary(model2)

test2 <- cox.zph(model2)
plot(test2[1], resid = FALSE)
abline(0,0, lty = 3)

#Model in men
men <- diab_lexis %>% filter(MALE == 1)

model3 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band) + factor(current) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = men)
summary(model3)

model3_1 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band) + factor(cig) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = men)
summary(model3_1)

#Model in women
women <- diab_lexis %>% filter(MALE == 0)

model4 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band) + factor(current) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = women)
summary(model4)

model4_1 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band) + factor(cig) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = women)
summary(model4_1)

#### All-cause mortality in current smokers####
active <- diab_lexis %>% filter(EVER_SMOK == 0 | CURR_SMOK == 1)

#Number of cigarettes/d
cox5 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(current) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = active) 
summary(cox5)

#Age started smoking
cox6 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(edad_inicio) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = active) 
summary(cox6)

#### All-cause mortality in former smokers####
#Former smokers and never smokers filter for Cox analysis
former <- diab_lexis %>% filter(EVER_SMOK == 0 | CURR_SMOK == 0)
former <- former %>% filter(!is.na(QUIT_NO_CIGS) | EVER_SMOK == 0) #Excludes those with NA in QUIT_NO_CIGS.

#Number of cigarettes/d
cox7 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(former_cigs) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = former) 
summary(cox7)

#Age started smoking
cox8 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(edad_inicio) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = former) 
summary(cox8)

#### All-cause mortality, smoking status and diabetes####
#Glycemic control
model9 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(cig)*factor(control) + factor(EDU_LEVEL) + factor(COYOACAN) + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = diab_lexis)
summary(model9)

#Glycemic control with HbA1C
model10 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(cig)*scale(BASE_HBA1C, scale = FALSE) + factor(EDU_LEVEL) + factor(COYOACAN) + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = diab_lexis)
summary(model10)

#Glycemic control in men
model9_1 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band) + factor(cig)*factor(control) + factor(EDU_LEVEL) + factor(COYOACAN) + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = men)
summary(model9_1)

#Glycemic control in women
model9_2 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band) + factor(cig)*factor(control) + factor(EDU_LEVEL) + factor(COYOACAN) + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = women)
summary(model9_2)

#Diabetes duration
model11 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(cig)*scale(TM_EVOLUCION, scale = FALSE) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + SBP1 + factor(ALCGP), data = diab_lexis)
summary(model11)

model12 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(cig)*factor(tiempo_evol) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + SBP1 + factor(ALCGP), data = diab_lexis)
summary(model12)

#Diabetes diagnosis
diag <- diab_lexis %>% filter(nodiag == 1)
no_diag <- diab_lexis %>% filter(nodiag == 0)

model14 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(current) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = diag)
summary(model14) #Nada con cig

model15 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(current) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = no_diag)
summary(model15)

#### Cause-specific mortality analyses####
#Cardiovascular mortality
cox_cardiovascular <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(current) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = cv_df)
summary(cox_cardiovascular)

#Myocardial infarction
cox_infarto <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(current) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = mi_df)
summary(cox_infarto)

#Stroke
cox_stroke <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(current) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = stroke_df)
summary(cox_stroke)

#Diabetes
cox_diabetes <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(current) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = diabetes_df)
summary(cox_diabetes)

#Acute diabetes
cox_diabetes_agudas <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(current) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = ac_diabetes_df)
summary(cox_diabetes_agudas)

#Chronic diabetes
cox_diabetes_cronicas <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(current) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = chr_diabetes_df)
summary(cox_diabetes_cronicas)

#Lung cancer
cox_ca_pulmon <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(current) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = lung_df)
summary(cox_ca_pulmon)

#Non-lung cancer
cox_ca_otros <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(current) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = non_lung_df)
summary(cox_ca_otros)

#COPD
cox_epoc <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(current) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = copd_df)
summary(cox_epoc)

#### Sensitivity Analysis####
#Mortality in daily smokers
model_s1 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(daily1) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = diab_lexis)
summary(model_s1)

#Mortality in daily smokers stratified by number of cigarettes
model_s2 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(daily) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = diab_lexis)
summary(model_s2)

#Mortality in non-daily smokers
model_s3 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(non_daily1) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = diab_lexis)
summary(model_s3)

#Mortality in non-daily smokers stratified by number of cigarettes
model_s4 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(non_daily) + factor(EDU_LEVEL) + factor(COYOACAN) + BASE_HBA1C + IMC + TM_EVOLUCION + SBP1 + factor(ALCGP), data = diab_lexis)
summary(model_s4)

#### All-cause mortality vs. healthy individuals####
#Filter excluding individuals with comorbidities
mcps1 <- mcps %>% filter(BASE_EMPHYSEMA == 0 & 
                           BASE_HEARTATTACK == 0 &
                           BASE_ANGINA == 0 & 
                           BASE_STROKE == 0 &
                           BASE_CKD == 0 & 
                           BASE_CIRR == 0 &
                           BASE_LUNGCANCER == 0 &
                           BASE_OTHCANCER == 0 &
                           BASE_PROSTATECANCER == 0 &
                           BASE_CERVCANCER == 0 &
                           BASE_BREASTCANCER == 0 &
                           BASE_STOMCANCER == 0 &
                           BASE_ORALCANCER == 0 &
                           BASE_PAD == 0)

#Excludes those with "U" status and those with accidental or non-defined death 
mcps1 <- mcps1 %>% 
  filter(STATUS == "A" | STATUS == "D") %>%
  mutate(STATUS = ifelse(STATUS == "D", 1, 0)) %>%  #Changes character values for numeric ones.
  filter(D054 != 1) %>% filter(D055 != 1)

##Variables
#BMI
mcps1 <- mutate(mcps1, IMC = (WEIGHT/(HEIGHT^2)*10000))

#Diabetes status
mcps1$diab[mcps1$BASE_DIABETES == 1 |
             mcps1$DRUG_D1 == 1 |
             mcps1$DRUG_D2 == 1 |
             mcps1$DRUG_D3 == 1 |
             mcps1$DRUG_D4 == 1 |
             mcps1$BASE_HBA1C >= 6.5] <- 1
mcps1$diab[mcps1$BASE_DIABETES == 0 &
             mcps1$DRUG_D1 == 0 &
             mcps1$DRUG_D2 == 0 &
             mcps1$DRUG_D3 == 0 &
             mcps1$DRUG_D4 == 0 &
             mcps1$BASE_HBA1C < 6.5] <- 0

mcps1 <- mcps1 %>% filter(!is.na(diab))

#Glycemic control
mcps1$control <- NA
mcps1$control[mcps1$BASE_HBA1C < 7] <- 1
mcps1$control[(mcps1$BASE_HBA1C >= 7) & (mcps1$BASE_HBA1C <10)] <- 2
mcps1$control[mcps1$BASE_HBA1C >= 10] <- 3

mcps1 %>% count(diab)
mcps1 %>% count(diab, control)
diab %>% count(control)

#Healthy individuals and individuals with diabetes by glycemic control (<7, 7-9, >10)
mcps1$sanos[mcps1$diab == 0] <- 0
mcps1$sanos[mcps1$diab == 1 & mcps1$control == 1] <- 1
mcps1$sanos[mcps1$diab == 1 & mcps1$control == 2] <- 2
mcps1$sanos[mcps1$diab == 1 & mcps1$control == 3] <- 3

#Healthy individuals and individuals with diabetes by glycemic control (<9, >9)
mcps1$sanos2[mcps1$diab == 0] <- 0
mcps1$sanos2[mcps1$diab == 1 & mcps1$BASE_HBA1C < 9] <- 1
mcps1$sanos2[mcps1$diab == 1 & mcps1$BASE_HBA1C >= 9] <- 2

mcps1 %>% count(sanos)
mcps1 %>% count(sanos2)

#Smoking status
mcps1$cig <- NULL
mcps1$cig[mcps1$EVER_SMOK == 0] <- 1
mcps1$cig[(mcps1$EVER_SMOK == 1) & (mcps1$CURR_SMOK == 0)] <- 2
mcps1$cig[(mcps1$EVER_SMOK == 1) & (mcps1$CURR_SMOK == 1)] <- 3

#Smokking status and diabetes status
mcps1$cig2[mcps1$diab == 0 & mcps1$cig == 1] <- 1
mcps1$cig2[mcps1$diab == 1 & mcps1$cig == 1] <- 2
mcps1$cig2[mcps1$diab == 1 & mcps1$cig == 2] <- 3
mcps1$cig2[mcps1$diab == 1 & mcps1$cig == 3] <- 4

#Smoking status and no. of cigarettes
mcps1$current <- NULL
mcps1$current[mcps1$EVER_SMOK == 0] <- 1
mcps1$current[mcps1$CURR_SMOK == 0] <- 2
mcps1$current[(mcps1$CURR_SMOK == 1) & (mcps1$CIGS_PER_DAY >= 1) & (mcps1$CIGS_PER_DAY < 5)] <- 3
mcps1$current[(mcps1$CURR_SMOK == 1) & (mcps1$CIGS_PER_DAY >= 5) & (mcps1$CIGS_PER_DAY < 10)] <- 4
mcps1$current[(mcps1$CURR_SMOK == 1) & (mcps1$CIGS_PER_DAY >= 10)] <- 5

mcps1$current2 <- NULL
mcps1$current2[mcps1$EVER_SMOK == 0 & mcps1$diab == 0] <- 1
mcps1$current2[mcps1$EVER_SMOK == 0 & mcps1$diab == 1] <- 2
mcps1$current2[mcps1$CURR_SMOK == 0 & mcps1$diab == 1] <- 3
mcps1$current2[(mcps1$CURR_SMOK == 1) & (mcps1$CIGS_PER_DAY >= 1) & (mcps1$CIGS_PER_DAY < 5) & (mcps1$diab == 1)] <- 4
mcps1$current2[(mcps1$CURR_SMOK == 1) & (mcps1$CIGS_PER_DAY >= 5) & (mcps1$CIGS_PER_DAY < 10) & (mcps1$diab == 1)] <- 5
mcps1$current2[(mcps1$CURR_SMOK == 1) & (mcps1$CIGS_PER_DAY >= 10) & (mcps1$diab == 1)] <- 6

mcps1 %>% count(SMOKEGP)
mcps1 %>% count(current2)

mcps1$a0 <- mcps1$AGE
mcps1$a1 <- mcps1$AGE + mcps1$PERSON_YEARS

#Smoking status prevalence
prevalence <- mcps1 %>% filter(!is.na(cig)) %>% 
  filter(diab == 0) %>% 
  group_by(MALE, cig) %>% 
  summarize(n = n()) %>% 
  mutate(pct = n / sum(n))

#Lexis expansion
mcps_lexis <- mcps1 %>% 
  mutate(date_recruited = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))

mcps_lexis <- Lexis(entry = list("period" = date_recruited,
                                 "age" = AGE,
                                 "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS),
                    exit.status = STATUS,
                    data = mcps_lexis)

mcps_lexis <- splitLexis(mcps_lexis, "age", breaks = seq(0, 85, 5))
mcps_lexis$age_band <- timeBand(mcps_lexis, "age", type = "factor")

mcps_lexis <- mcps_lexis %>% 
  rename(id = lex.id,
         fu_period = lex.dur,
         entry_status = lex.Cst,
         exit_status = lex.Xst) %>% 
  mutate(time_at_exit = time_at_entry + fu_period,
         age_at_exit = age + fu_period) %>% 
  relocate(time_at_exit, .after = time_at_entry) %>% 
  relocate(age_band, .after = period) %>% 
  relocate(age_at_exit, .after = age)

mcps_lexis$premature <- ifelse(mcps_lexis$exit_status == 1, 1, 0)
mcps_lexis$premature[mcps_lexis$premature == 1 & mcps_lexis$age_at_exit >= 75] <- 0

#Modelo con interacción
model0 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(cig) + factor(sanos) + factor(EDU_LEVEL) + factor(COYOACAN) + IMC + factor(ALCGP), data = mcps_lexis)
summary(model0)

model0_1 <- coxph(Surv(time_at_entry, time_at_exit, premature) ~ strata(age_band, MALE) + factor(cig)*factor(sanos) + factor(EDU_LEVEL) + factor(COYOACAN) + IMC + factor(ALCGP), data = mcps_lexis)
summary(model0_1)

#### Figure 1: Smoking prevalence by sex ####
smok_w <- diab %>% filter(MALE == 0) %>% 
  group_by(edad, cig) %>% 
  summarize(n = n()) %>% 
  mutate(pct = n / sum(n))
smok_w$cig <- factor(smok_w$cig, labels = c("Never smoker", "Former smoker", "Current smoker"))
smok_w$edad <- factor(smok_w$edad, labels = c("35-44", "45-54", "55-64", "65-74", "≥75"))

smok_m <- diab %>% filter(MALE == 1) %>% 
  group_by(edad, cig) %>% 
  summarize(n = n()) %>% 
  mutate(pct = n / sum(n))
smok_m$cig <- factor(smok_m$cig, labels = c("Never smoker", "Former smoker", "Current smoker"))
smok_m$edad <- factor(smok_m$edad, labels = c("35-44", "45-54", "55-64", "65-74", "≥75"))

fig1_1 <- smok_w %>% ggplot(aes(x = edad, y = pct, color = cig, group = cig)) +
  geom_point() +
  geom_line() +
  xlab("Age") +
  ylab("Proportion") +
  ylim(0, 0.8) +
  theme(legend.title=element_blank()) +
  labs(colour = NULL) +
  theme_classic() +
  scale_color_manual(values=c("#772986", "#AE4532", "#3172B1"))

fig1_2 <- smok_m %>% ggplot(aes(x = edad, y = pct, color = cig, group = cig)) +
  geom_point() +
  geom_line() +
  xlab("Age") +
  ylab("Proportion") +
  ylim(0, 0.8) +
  labs(colour = NULL) +
  theme_classic() +
  scale_color_manual(values=c("#772986", "#AE4532", "#3172B1"))

fig1 <- ggarrange(fig1_1, fig1_2, nrow = 1, ncol = 2, labels = "AUTO", common.legend = TRUE, legend = "bottom")

ggsave(fig1, file="Proyectos/Smoking and diabetes/fig1.jpg", bg="transparent",
       width=20, height=10, units=c("cm"), dpi=600, limitsize = FALSE)

#### Figure 2: All-cause mortality ####
m1 <- summary(model1)
m2 <- summary(model2)

df_cig <- data.frame(status = c("Former smoker", "Current smoker"),
                     estimate = c(m1[["conf.int"]][1,1], m1[["conf.int"]][2,1]), 
                     lower = c(m1[["conf.int"]][1,3], m1[["conf.int"]][2,3]), 
                     upper = c(m1[["conf.int"]][1,4], m1[["conf.int"]][2,4]))
df_cig$status <- factor(df_cig$status, levels = c("Former smoker", "Current smoker"))

df_current <- data.frame(status = c("<5 cig/d", "5-9 cig/d", "≥10 cig/d"),
                         estimate = c(m2[["conf.int"]][2,1], m2[["conf.int"]][3,1], m2[["conf.int"]][4,1]), 
                         lower = c(m2[["conf.int"]][2,3], m2[["conf.int"]][3,3], m2[["conf.int"]][4,3]), 
                         upper = c(m2[["conf.int"]][2,4], m2[["conf.int"]][3,4], m2[["conf.int"]][4,4]))
df_current$status <- factor(df_current$status, levels = c("<5 cig/d", "5-9 cig/d", "≥10 cig/d"))

fig2_1 <- ggplot(df_cig, aes(status, estimate, ymin = lower, ymax = upper, color = status)) +
  geom_pointrange() +
  ylim(0.95, 1.25) +
  labs(x = "",
       y = "HR with 95% CI",
       title = "Smoking status and all-cause mortality") +
  geom_hline(yintercept = 1, lty = "dashed") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values=c("#AE4532", "#3172B1"))

fig2_2 <- ggplot(df_current, aes(status, estimate, ymin = lower, ymax = upper, color = status)) +
  geom_pointrange() +
  ylim(0.9, 1.40) +
  labs(x = "",
       y = "HR with 95% CI",
       title = "Current smokers and no. of cig/d") +
  geom_hline(yintercept = 1, lty = "dashed") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values=c("#3172B1", "#3172B1", "#3172B1"))

#All-cause mortality stratified by sex
m3 <- summary(model3) #Men
m3_1 <- summary(model3_1) 
m4 <- summary(model4) #Women
m4_1 <- summary(model4_1) 

df_men_current <- data.frame(status = c("<5 cig/d", "5-9 cig/d", "≥10 cig/d"),
                             estimate = c(m3[["conf.int"]][2,1], m3[["conf.int"]][3,1], m3[["conf.int"]][4,1]), 
                             lower = c(m3[["conf.int"]][2,3], m3[["conf.int"]][3,3], m3[["conf.int"]][4,3]), 
                             upper = c(m3[["conf.int"]][2,4], m3[["conf.int"]][3,4], m3[["conf.int"]][4,4]))
df_men_current$status <- factor(df_men_current$status, levels = c("<5 cig/d", "5-9 cig/d", "≥10 cig/d"))

df_men_cig <- data.frame(status = c("Former smoker", "Current smoker"),
                         estimate = c(m3_1[["conf.int"]][1,1], m3_1[["conf.int"]][2,1]), 
                         lower = c(m3_1[["conf.int"]][1,3], m3_1[["conf.int"]][2,3]), 
                         upper = c(m3_1[["conf.int"]][1,4], m3_1[["conf.int"]][2,4]))
df_men_cig$status <- factor(df_men_cig$status, levels = c("Former smoker", "Current smoker"))


df_women_current <- data.frame(status = c("<5 cig/d", "5-9 cig/d", "≥10 cig/d"),
                               estimate = c(m4[["conf.int"]][2,1], m4[["conf.int"]][3,1], m4[["conf.int"]][4,1]), 
                               lower = c(m4[["conf.int"]][2,3], m4[["conf.int"]][3,3], m4[["conf.int"]][4,3]), 
                               upper = c(m4[["conf.int"]][2,4], m4[["conf.int"]][3,4], m4[["conf.int"]][4,4]))
df_women_current$status <- factor(df_women_current$status, levels = c("<5 cig/d", "5-9 cig/d", "≥10 cig/d"))

df_women_cig <- data.frame(status = c("Former smoker", "Current smoker"),
                           estimate = c(m4_1[["conf.int"]][1,1], m4_1[["conf.int"]][2,1]), 
                           lower = c(m4_1[["conf.int"]][1,3], m4_1[["conf.int"]][2,3]), 
                           upper = c(m4_1[["conf.int"]][1,4], m4_1[["conf.int"]][2,4]))
df_women_cig$status <- factor(df_women_cig$status, levels = c("Former smoker", "Current smoker"))


#Women
fig2_3 <- ggplot(df_women_current, aes(status, estimate, ymin = lower, ymax = upper, color = status)) +
  geom_pointrange() +
  ylim(0.8, 1.45) +
  labs(x = "",
       y = "HR with 95% CI",
       title = "Women and no. of cig/d") +
  geom_hline(yintercept = 1, lty = "dashed") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values=c("#3172B1", "#3172B1", "#3172B1"))

fig2_3_1 <- ggplot(df_women_cig, aes(status, estimate, ymin = lower, ymax = upper, color = status)) +
  geom_pointrange() +
  ylim(0.95, 1.30) +
  labs(x = "",
       y = "HR with 95% CI",
       title = "Women and smoking status") +
  geom_hline(yintercept = 1, lty = "dashed") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values=c("#AE4532", "#3172B1"))

#Men
fig2_4 <- ggplot(df_men_current, aes(status, estimate, ymin = lower, ymax = upper, color = status)) +
  geom_pointrange() +
  ylim(0.85, 1.45) +
  labs(x = "",
       y = "HR with 95% CI",
       title = "Men and no. of cig/d") +
  geom_hline(yintercept = 1, lty = "dashed") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values=c("#3172B1", "#3172B1", "#3172B1"))

fig2_4_1 <- ggplot(df_men_cig, aes(status, estimate, ymin = lower, ymax = upper, color = status)) +
  geom_pointrange() +
  ylim(0.90, 1.30) +
  labs(x = "",
       y = "HR with 95% CI",
       title = "Men and smoking status") +
  geom_hline(yintercept = 1, lty = "dashed") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values=c("#AE4532", "#3172B1"))

#### Figure 2: Current and former smokers characteristics ####
m5 <- summary(cox5)
m6 <- summary(cox6)
m7 <- summary(cox7)
m8 <- summary(cox8)

df_age_curr <- data.frame(status = c("<15 years", "15-24 years", "≥25 years"),
                          estimate = c(m6[["conf.int"]][1,1], m6[["conf.int"]][2,1], m6[["conf.int"]][3,1]), 
                          lower = c(m6[["conf.int"]][1,3], m6[["conf.int"]][2,3], m6[["conf.int"]][3,3]), 
                          upper = c(m6[["conf.int"]][1,4], m6[["conf.int"]][2,4], m6[["conf.int"]][3,4]))
df_age_curr$status <- factor(df_age_curr$status, levels = c("<15 years", "15-24 years", "≥25 years"))

df_age_for <- data.frame(status = c("<15 years", "15-24 years", "≥25 years"),
                         estimate = c(m8[["conf.int"]][1,1], m8[["conf.int"]][2,1], m8[["conf.int"]][3,1]), 
                         lower = c(m8[["conf.int"]][1,3], m8[["conf.int"]][2,3], m8[["conf.int"]][3,3]), 
                         upper = c(m8[["conf.int"]][1,4], m8[["conf.int"]][2,4], m8[["conf.int"]][3,4]))
df_age_for$status <- factor(df_age_for$status, levels = c("<15 years", "15-24 years", "≥25 years"))

df_num_for <- data.frame(status = c("<5 cig/d", "5-9 cig/d", "≥10 cig/d"),
                         estimate = c(m7[["conf.int"]][1,1], m7[["conf.int"]][2,1], m7[["conf.int"]][3,1]), 
                         lower = c(m7[["conf.int"]][1,3], m7[["conf.int"]][2,3], m7[["conf.int"]][3,3]), 
                         upper = c(m7[["conf.int"]][1,4], m7[["conf.int"]][2,4], m7[["conf.int"]][3,4]))
df_num_for$status <- factor(df_num_for$status, levels = c("<5 cig/d", "5-9 cig/d", "≥10 cig/d"))

fig3_1 <- ggplot(df_age_curr, aes(status, estimate, ymin = lower, ymax = upper, color = status)) +
  geom_pointrange() +
  ylim(0.9, 1.45) +
  labs(x = "",
       y = "HR with 95% CI",
       title = "Current smokers and age started smoking") +
  geom_hline(yintercept = 1, lty = "dashed") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#3172B1", "#3172B1", "#3172B1"))

fig3_2 <- ggplot(df_age_for, aes(status, estimate, ymin = lower, ymax = upper, color = status)) +
  geom_pointrange() +
  ylim(0.9, 1.3) +
  labs(x = "",
       y = "HR with 95% CI",
       title = "Former smokers and age started smoking") +
  geom_hline(yintercept = 1, lty = "dashed") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#AE4532", "#AE4532", "#AE4532"))

fig3_3 <- ggplot(df_num_for, aes(status, estimate, ymin = lower, ymax = upper, color = status)) +
  geom_pointrange() +
  ylim(0.9, 1.4) +
  labs(x = "",
       y = "HR with 95% CI",
       title = "Former smokers and number of cig/d") +
  geom_hline(yintercept = 1, lty = "dashed") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#AE4532", "#AE4532", "#AE4532"))

fig2_1a <- ggarrange("",fig2_1,"", nrow = 1, ncol = 3, labels = c("", "A", ""), widths = c(0.25,0.5, 0.25))
fig2b <- ggarrange(fig2_2, fig3_1, fig3_3, fig3_2, nrow = 2, ncol = 2, labels = c("B", "C", "D", "E"))

fig2 <- ggarrange(fig2_1a, fig2b, nrow = 2, ncol = 1, heights = c(0.33, 0.66))

ggsave(fig2, file="Proyectos/Smoking and diabetes/fig2.jpg", bg="transparent",
       width=20, height=28, units=c("cm"), dpi=600, limitsize = FALSE)

fig2_sex <- ggarrange(fig2_3_1, fig2_4_1, fig2_3, fig2_4, nrow = 2, ncol = 2, labels = "AUTO")

ggsave(fig2_sex, file="Proyectos/Smoking and diabetes/fig2_sex.jpg", bg="transparent",
       width=20, height=20, units=c("cm"), dpi=600, limitsize = FALSE)

#### Figure 3: Smoking and glycemic control diabetes ####

m9 <- summary(model9)
m9_1 <- summary(model9_1) #Men
m9_2 <- summary(model9_2) #Women

df_control <- data.frame(status = c("Former", "Current", "Former", "Current"),
                         control = c(2, 2, 3, 3),
                         estimate = c(m9[["conf.int"]][16,1], m9[["conf.int"]][17,1], m9[["conf.int"]][18,1], m9[["conf.int"]][19,1]),
                         lower = c(m9[["conf.int"]][16,3], m9[["conf.int"]][17,3], m9[["conf.int"]][18,3], m9[["conf.int"]][19,3]), 
                         upper = c(m9[["conf.int"]][16,4], m9[["conf.int"]][17,4], m9[["conf.int"]][18,4], m9[["conf.int"]][19,4]))
df_control$status <- factor(df_control$status, levels = c("Former", "Current"))
df_control$control <- factor(df_control$control, labels = c("HbA1C 7-9%", "HbA1C ≥10%"))

df_control_men <- data.frame(status = c("Former", "Current", "Former", "Current"),
                             control = c(2, 2, 3, 3),
                             estimate = c(m9_1[["conf.int"]][16,1], m9_1[["conf.int"]][17,1], m9_1[["conf.int"]][18,1], m9_1[["conf.int"]][19,1]),
                             lower = c(m9_1[["conf.int"]][16,3], m9_1[["conf.int"]][17,3], m9_1[["conf.int"]][18,3], m9_1[["conf.int"]][19,3]), 
                             upper = c(m9_1[["conf.int"]][16,4], m9_1[["conf.int"]][17,4], m9_1[["conf.int"]][18,4], m9_1[["conf.int"]][19,4]))
df_control_men$status <- factor(df_control_men$status, levels = c("Former", "Current"))
df_control_men$control <- factor(df_control_men$control, labels = c("HbA1C 7-9%", "HbA1C ≥10%"))

df_control_women <- data.frame(status = c("Former", "Current", "Former", "Current"),
                               control = c(2, 2, 3, 3),
                               estimate = c(m9_2[["conf.int"]][16,1], m9_2[["conf.int"]][17,1], m9_2[["conf.int"]][18,1], m9_2[["conf.int"]][19,1]),
                               lower = c(m9_2[["conf.int"]][16,3], m9_2[["conf.int"]][17,3], m9_2[["conf.int"]][18,3], m9_2[["conf.int"]][19,3]), 
                               upper = c(m9_2[["conf.int"]][16,4], m9_2[["conf.int"]][17,4], m9_2[["conf.int"]][18,4], m9_2[["conf.int"]][19,4]))
df_control_women$status <- factor(df_control_women$status, levels = c("Former", "Current"))
df_control_women$control <- factor(df_control_women$control, labels = c("HbA1C 7-9%", "HbA1C ≥10%"))

fig4_0 <- ggplot(df_control, aes(status, estimate, ymin = lower, ymax = upper)) +
  geom_pointrange() +
  ylim(0.75, 1.7) +
  labs(x = "",
       y = "HR with 95% CI",
       title = "Smoking and glycemic control") +
  geom_hline(yintercept = 1, lty = "dashed") +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(~ control)

fig4_1 <- ggplot(df_control_men, aes(status, estimate, ymin = lower, ymax = upper)) +
  geom_pointrange() +
  ylim(0.75, 2.1) +
  labs(x = "",
       y = "HR with 95% CI",
       title = "Smoking and glycemic control in men") +
  geom_hline(yintercept = 1, lty = "dashed") +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(~ control)

fig4_2 <- ggplot(df_control_women, aes(status, estimate, ymin = lower, ymax = upper)) +
  geom_pointrange() +
  ylim(0.70, 2.3) +
  labs(x = "",
       y = "HR with 95% CI",
       title = "Smoking and glycemic control in women") +
  geom_hline(yintercept = 1, lty = "dashed") +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(~ control)

#Diabetes diagnosis and smoking status
m14 <- summary(model14) #Diag
m15 <- summary(model15) #No diag

df_diag <- data.frame(status = c("Former smokers", "<5 cig/d", "5-9 cig/d", "≥10 cig/d"),
                      estimate = c(m14[["conf.int"]][1,1], m14[["conf.int"]][2,1], m14[["conf.int"]][3,1], m14[["conf.int"]][4,1]),
                      lower = c(m14[["conf.int"]][1,3], m14[["conf.int"]][2,3], m14[["conf.int"]][3,3], m14[["conf.int"]][4,3]), 
                      upper = c(m14[["conf.int"]][1,4], m14[["conf.int"]][2,4], m14[["conf.int"]][3,4], m14[["conf.int"]][4,4]))
df_diag$status <- factor(df_diag$status, levels = c("Former smokers", "<5 cig/d", "5-9 cig/d", "≥10 cig/d"))

df_no_diag <- data.frame(status = c("Former smokers", "<5 cig/d", "5-9 cig/d", "≥10 cig/d"),
                         estimate = c(m15[["conf.int"]][1,1], m15[["conf.int"]][2,1], m15[["conf.int"]][3,1], m15[["conf.int"]][4,1]),
                         lower = c(m15[["conf.int"]][1,3], m15[["conf.int"]][2,3], m15[["conf.int"]][3,3], m15[["conf.int"]][4,3]), 
                         upper = c(m15[["conf.int"]][1,4], m15[["conf.int"]][2,4], m15[["conf.int"]][3,4], m15[["conf.int"]][4,4]))
df_no_diag$status <- factor(df_no_diag$status, levels = c("Former smokers", "<5 cig/d", "5-9 cig/d", "≥10 cig/d"))

fig4_3 <- ggplot(df_diag, aes(status, estimate, ymin = lower, ymax = upper)) +
  geom_pointrange() +
  ylim(0.85, 1.40) +
  labs(x = "",
       y = "HR with 95% CI",
       title = "Smoking and diagnosed diabetes") +
  geom_hline(yintercept = 1, lty = "dashed") +
  theme_classic() +
  theme(legend.position = "none")

fig4_4 <- ggplot(df_no_diag, aes(status, estimate, ymin = lower, ymax = upper)) +
  geom_pointrange() +
  ylim(0.85, 1.9) +
  labs(x = "",
       y = "HR with 95% CI",
       title = "Smoking and undiagnosed diabetes") +
  geom_hline(yintercept = 1, lty = "dashed") +
  theme_classic() +
  theme(legend.position = "none")

fig3a <- ggarrange("",fig4_0,"", nrow = 1, ncol = 3, labels = c("","A", ""), widths = c(0.25, 0.5, 0.25))
fig3b <- ggarrange(fig4_3, fig4_4, nrow = 1, ncol = 2, labels =c("B", "C"))
fig3 <- ggarrange(fig3a, fig3b, nrow = 2, ncol = 1)

ggsave(fig3, file="Proyectos/Smoking and diabetes/fig3.jpg", bg="transparent",
       width=20, height=20, units=c("cm"), dpi=600, limitsize = FALSE)

#### Figure 4: Smoking and diabetes vs. healthy individuals ####
m0 <- summary(model0_1)

df_health <- data.frame(status = c("Former", "Current", "Former", "Current", "Former", "Current"),
                        sanos = c(1, 1, 2, 2, 3, 3),
                        estimate = c(m0[["conf.int"]][15,1], m0[["conf.int"]][16,1], m0[["conf.int"]][17,1], m0[["conf.int"]][18,1], m0[["conf.int"]][19,1], m0[["conf.int"]][20,1]),
                        lower = c(m0[["conf.int"]][15,3], m0[["conf.int"]][16,3], m0[["conf.int"]][17,3], m0[["conf.int"]][18,3], m0[["conf.int"]][19,3], m0[["conf.int"]][20,3]), 
                        upper = c(m0[["conf.int"]][15,4], m0[["conf.int"]][16,4], m0[["conf.int"]][17,4], m0[["conf.int"]][18,4], m0[["conf.int"]][19,4], m0[["conf.int"]][20,4]))
df_health$status <- factor(df_health$status, levels = c("Former", "Current"))
df_health$sanos <- factor(df_health$sanos, labels = c("HbA1C <7%", "HbA1C 7-9%", "HbA1C ≥10%"))

fig4 <- ggplot(df_health, aes(status, estimate, ymin = lower, ymax = upper)) +
  geom_pointrange() +
  ylim(0.65, 1.40) +
  labs(x = "",
       y = "HR with 95% CI",
       title = "Smoking and diabetes vs. healthy individuals") +
  geom_hline(yintercept = 1, lty = "dashed") +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(~ sanos)

ggsave(fig4, file="Proyectos/Smoking and diabetes/fig4.jpg", bg="transparent",
       width=15, height=15, units=c("cm"), dpi=600, limitsize = FALSE)

#### Figure 5: Cause-specific mortality in current smokers ####
cv <- summary(cox_cardiovascular)
im <- summary(cox_infarto)
str <- summary(cox_stroke)
db <- summary(cox_diabetes)
db_a <- summary(cox_diabetes_agudas)
db_c <- summary(cox_diabetes_cronicas)
cp <- summary(cox_ca_pulmon)
cnp <- summary(cox_ca_otros)
epc <- summary(cox_epoc)

m_especif_curr <- data.frame(mean = c(cv[["conf.int"]][2,1], cv[["conf.int"]][3,1], cv[["conf.int"]][4,1],
                                      im[["conf.int"]][2,1], im[["conf.int"]][3,1], im[["conf.int"]][4,1],
                                      str[["conf.int"]][2,1], str[["conf.int"]][3,1], str[["conf.int"]][4,1],
                                      db[["conf.int"]][2,1], db[["conf.int"]][3,1], db[["conf.int"]][4,1],
                                      db_a[["conf.int"]][2,1], db_a[["conf.int"]][3,1], db_a[["conf.int"]][4,1],
                                      db_c[["conf.int"]][2,1], db_c[["conf.int"]][3,1], db_c[["conf.int"]][4,1],
                                      cp[["conf.int"]][2,1], cp[["conf.int"]][3,1], cp[["conf.int"]][4,1],
                                      cnp[["conf.int"]][2,1], cnp[["conf.int"]][3,1], cnp[["conf.int"]][4,1],
                                      epc[["conf.int"]][2,1], epc[["conf.int"]][3,1], epc[["conf.int"]][4,1]),
                             lower = c(cv[["conf.int"]][2,3], cv[["conf.int"]][3,3], cv[["conf.int"]][4,3],
                                       im[["conf.int"]][2,3], im[["conf.int"]][3,3], im[["conf.int"]][4,3],
                                       str[["conf.int"]][2,3], str[["conf.int"]][3,3], str[["conf.int"]][4,3],
                                       db[["conf.int"]][2,3], db[["conf.int"]][3,3], db[["conf.int"]][4,3],
                                       db_a[["conf.int"]][2,3], db_a[["conf.int"]][3,3], db_a[["conf.int"]][4,3],
                                       db_c[["conf.int"]][2,3], db_c[["conf.int"]][3,3], db_c[["conf.int"]][4,3],
                                       cp[["conf.int"]][2,3], cp[["conf.int"]][3,3], cp[["conf.int"]][4,3],
                                       cnp[["conf.int"]][2,3], cnp[["conf.int"]][3,3], cnp[["conf.int"]][4,3],
                                       epc[["conf.int"]][2,3], epc[["conf.int"]][3,3], epc[["conf.int"]][4,3]),
                             upper = c(cv[["conf.int"]][2,4], cv[["conf.int"]][3,4], cv[["conf.int"]][4,4],
                                       im[["conf.int"]][2,4], im[["conf.int"]][3,4], im[["conf.int"]][4,4],
                                       str[["conf.int"]][2,4], str[["conf.int"]][3,4], str[["conf.int"]][4,4],
                                       db[["conf.int"]][2,4], db[["conf.int"]][3,4], db[["conf.int"]][4,4],
                                       db_a[["conf.int"]][2,4], db_a[["conf.int"]][3,4], db_a[["conf.int"]][4,4],
                                       db_c[["conf.int"]][2,4], db_c[["conf.int"]][3,4], db_c[["conf.int"]][4,4],
                                       cp[["conf.int"]][2,4], cp[["conf.int"]][3,4], cp[["conf.int"]][4,4],
                                       cnp[["conf.int"]][2,4], cnp[["conf.int"]][3,4], cnp[["conf.int"]][4,4],
                                       epc[["conf.int"]][2,4], epc[["conf.int"]][3,4], epc[["conf.int"]][4,4]),
                             cause = c("Cardiovascular", NA, NA, "Myocardial infarction", NA, NA, "Stroke", NA, NA, "Diabetes", NA, NA, "Acute diabetes", NA, NA,"Chronic diabetes", NA, NA,  "Lung cancer", NA, NA,
                                  "Non-lung cancer", NA, NA, "COPD", NA, NA),
                             smoke_status = rep(c("<5 cig/d", "5-9 cig/d", "≥10 cig/d"), 9),
                        est = c(paste0(round(cv[["conf.int"]][2,1],2), " (", round(cv[["conf.int"]][2,3],2), "-",  round(cv[["conf.int"]][2,4],2), ")"),
                                paste0(round(cv[["conf.int"]][3,1],2), " (", round(cv[["conf.int"]][3,3],2), "-",  round(cv[["conf.int"]][3,4],2), ")"),
                                paste0(round(cv[["conf.int"]][4,1],2), " (", round(cv[["conf.int"]][4,3],2), "-",  round(cv[["conf.int"]][4,4],2), ")"),
                                paste0(round(im[["conf.int"]][2,1],2), " (", round(im[["conf.int"]][2,3],2), "-",  round(im[["conf.int"]][2,4],2), ")"),
                                paste0(round(im[["conf.int"]][3,1],2), " (", round(im[["conf.int"]][3,3],2), "-",  round(im[["conf.int"]][3,4],2), ")"),
                                paste0(round(im[["conf.int"]][4,1],2), " (", round(im[["conf.int"]][4,3],2), "-",  round(im[["conf.int"]][4,4],2), ")"),
                                paste0(round(str[["conf.int"]][2,1],2), " (", round(str[["conf.int"]][2,3],2), "-",  round(str[["conf.int"]][2,4],2), ")"),
                                paste0(round(str[["conf.int"]][3,1],2), " (", round(str[["conf.int"]][3,3],2), "-",  round(str[["conf.int"]][3,4],2), ")"),
                                paste0(round(str[["conf.int"]][4,1],2), " (", round(str[["conf.int"]][4,3],2), "-",  round(str[["conf.int"]][4,4],2), ")"),
                                paste0(round(db[["conf.int"]][2,1],2), " (", round(db[["conf.int"]][2,3],2), "-",  round(db[["conf.int"]][2,4],2), ")"),
                                paste0(round(db[["conf.int"]][3,1],2), " (", round(db[["conf.int"]][3,3],2), "-",  round(db[["conf.int"]][3,4],2), ")"),
                                paste0(round(db[["conf.int"]][4,1],2), " (", round(db[["conf.int"]][4,3],2), "-",  round(db[["conf.int"]][4,4],2), ")"),
                                paste0(round(db_a[["conf.int"]][2,1],2), " (", round(db_a[["conf.int"]][2,3],2), "-",  round(db_a[["conf.int"]][2,4],2), ")"),
                                paste0(round(db_a[["conf.int"]][3,1],2), " (", round(db_a[["conf.int"]][3,3],2), "-",  round(db_a[["conf.int"]][3,4],2), ")"),
                                paste0(round(db_a[["conf.int"]][4,1],2), " (", round(db_a[["conf.int"]][4,3],2), "-",  round(db_a[["conf.int"]][4,4],2), ")"),
                                paste0(round(db_c[["conf.int"]][2,1],2), " (", round(db_c[["conf.int"]][2,3],2), "-",  round(db_c[["conf.int"]][2,4],2), ")"),
                                paste0(round(db_c[["conf.int"]][3,1],2), " (", round(db_c[["conf.int"]][3,3],2), "-",  round(db_c[["conf.int"]][3,4],2), ")"),
                                paste0(round(db_c[["conf.int"]][4,1],2), " (", round(db_c[["conf.int"]][4,3],2), "-",  round(db_c[["conf.int"]][4,4],2), ")"),
                                paste0(round(cp[["conf.int"]][2,1],2), " (", round(cp[["conf.int"]][2,3],2), "-",  round(cp[["conf.int"]][2,4],2), ")"),
                                paste0(round(cp[["conf.int"]][3,1],2), " (", round(cp[["conf.int"]][3,3],2), "-",  round(cp[["conf.int"]][3,4],2), ")"),
                                paste0(round(cp[["conf.int"]][4,1],2), " (", round(cp[["conf.int"]][4,3],2), "-",  round(cp[["conf.int"]][4,4],2), ")"),
                                paste0(round(cnp[["conf.int"]][2,1],2), " (", round(cnp[["conf.int"]][2,3],2), "-",  round(cnp[["conf.int"]][2,4],2), ")"),
                                paste0(round(cnp[["conf.int"]][3,1],2), " (", round(cnp[["conf.int"]][3,3],2), "-",  round(cnp[["conf.int"]][3,4],2), ")"),
                                paste0(round(cnp[["conf.int"]][4,1],2), " (", round(cnp[["conf.int"]][4,3],2), "-",  round(cnp[["conf.int"]][4,4],2), ")"),
                                paste0(round(epc[["conf.int"]][2,1],2), " (", round(epc[["conf.int"]][2,3],2), "-",  round(epc[["conf.int"]][2,4],2), ")"),
                                paste0(round(epc[["conf.int"]][3,1],2), " (", round(epc[["conf.int"]][3,3],2), "-",  round(epc[["conf.int"]][3,4],2), ")"),
                                paste0(round(epc[["conf.int"]][4,1],2), " (", round(epc[["conf.int"]][4,3],2), "-",  round(epc[["conf.int"]][4,4],2), ")")))

m1<-m_especif_curr %>% 
  forestplot(labeltext = c(cause, smoke_status, est),
             graph.pos = 3,
             xlog = TRUE) %>%
  fp_add_header(cause = c("Cause of death"),
                smoke_status = c("No. of cigs/d"),
                est = c("HR (95% CI)")) %>%
  fp_add_lines(h_1 = gpar(lty = 1),
               h_2 = gpar(lty = 1),
               h_5 = gpar(lty = 1),
               h_8 = gpar(lty = 1),
               h_11 = gpar(lty = 1),
               h_14 = gpar(lty = 1),
               h_17 = gpar(lty = 1),
               h_20 = gpar(lty = 1),
               h_23 = gpar(lty = 1),
               h_26 = gpar(lty = 1),
               h_29 = gpar(lty = 1))
k1<-as.ggplot(~plot(m1))

m_especif_for <- data.frame(mean = c(cv[["conf.int"]][1,1], 
                                     im[["conf.int"]][1,1], 
                                     str[["conf.int"]][1,1], 
                                     db[["conf.int"]][1,1], 
                                     db_a[["conf.int"]][1,1],
                                     db_c[["conf.int"]][1,1],
                                     cp[["conf.int"]][1,1],
                                     cnp[["conf.int"]][1,1],
                                     epc[["conf.int"]][1,1]),
                            lower = c(cv[["conf.int"]][1,3],
                                      im[["conf.int"]][1,3],
                                      str[["conf.int"]][1,3],
                                      db[["conf.int"]][1,3],
                                      db_a[["conf.int"]][1,3],
                                      db_c[["conf.int"]][1,3],
                                      cp[["conf.int"]][1,3], 
                                      cnp[["conf.int"]][1,3],
                                      epc[["conf.int"]][1,3]),
                            upper = c(cv[["conf.int"]][1,4], 
                                      im[["conf.int"]][1,4], 
                                      str[["conf.int"]][1,4],
                                      db[["conf.int"]][1,4],
                                      db_a[["conf.int"]][1,4],
                                      db_c[["conf.int"]][1,4],
                                      cp[["conf.int"]][1,4],
                                      cnp[["conf.int"]][1,4],
                                      epc[["conf.int"]][1,4]),
                             cause = c("Cardiovascular", "Myocardial infarction", "Stroke", "Diabetes", "Acute diabetes", "Chronic diabetes", "Lung cancer",
                                       "Non-lung cancer", "COPD"),
                             smoke_status = rep(c("Former smoker"), 9),
                             est = c(paste0(round(cv[["conf.int"]][1,1],2), " (", round(cv[["conf.int"]][1,3],2), "-",  round(cv[["conf.int"]][1,4],2), ")"),
                                     paste0(round(im[["conf.int"]][1,1],2), " (", round(im[["conf.int"]][1,3],2), "-",  round(im[["conf.int"]][1,4],2), ")"),
                                     paste0(round(str[["conf.int"]][1,1],2), " (", round(str[["conf.int"]][1,3],2), "-",  round(str[["conf.int"]][1,4],2), ")"),
                                     paste0(round(db[["conf.int"]][1,1],2), " (", round(db[["conf.int"]][1,3],2), "-",  round(db[["conf.int"]][1,4],2), ")"),
                                     paste0(round(db_a[["conf.int"]][1,1],2), " (", round(db_a[["conf.int"]][1,3],2), "-",  round(db_a[["conf.int"]][1,4],2), ")"),
                                     paste0(round(db_c[["conf.int"]][1,1],2), " (", round(db_c[["conf.int"]][1,3],2), "-",  round(db_c[["conf.int"]][1,4],2), ")"),
                                     paste0(round(cp[["conf.int"]][1,1],2), " (", round(cp[["conf.int"]][1,3],2), "-",  round(cp[["conf.int"]][1,4],2), ")"),
                                     paste0(round(cnp[["conf.int"]][1,1],2), " (", round(cnp[["conf.int"]][1,3],2), "-",  round(cnp[["conf.int"]][1,4],2), ")"),
                                     paste0(round(epc[["conf.int"]][1,1],2), " (", round(epc[["conf.int"]][1,3],2), "-",  round(epc[["conf.int"]][1,4],2), ")")))

m2<-m_especif_for %>% 
  forestplot(labeltext = c(cause, smoke_status, est),
             graph.pos = 3,
             xlog = TRUE) %>%
  fp_add_header(cause = c("Cause of death"),
                smoke_status = c("Smoking status"),
                est = c("HR (95% CI)")) %>%
  fp_add_lines(h_1 = gpar(lty = 1),
               h_2 = gpar(lty = 1),
               h_3 = gpar(lty = 1),
               h_4 = gpar(lty = 1),
               h_5 = gpar(lty = 1),
               h_6 = gpar(lty = 1),
               h_7 = gpar(lty = 1),
               h_8 = gpar(lty = 1),
               h_9 = gpar(lty = 1),
               h_10 = gpar(lty = 1),
               h_11 = gpar(lty = 1))
k2<-as.ggplot(~plot(m2))


ggsave(k1, file="Proyectos/Smoking and diabetes/fig5.jpg", bg="transparent",
       width=20, height=20, units=c("cm"), dpi=600, limitsize = FALSE)

ggsave(k2, file="Proyectos/Smoking and diabetes/suppfig5.jpg", bg="transparent",
       width=20, height=15, units=c("cm"), dpi=600, limitsize = FALSE)

#### Supplementary Figure 1: Alluvial plots ####
#Convert cig and cig_resurv to factor for alluvial plots
diab_resurv$cig <- factor(diab_resurv$cig, labels = c("Never smoker", "Former smoker", "Current smoker"))
diab_resurv$cig_resurv <- factor(diab_resurv$cig_resurv, labels = c("Never smoker", "Former smoker", "Current smoker"))

#Change in smoking status from baseline to resurvey
curr_currresurv <- diab_resurv %>%
  filter(!is.na(cig_resurv)) %>% 
  group_by(cig, cig_resurv) %>% 
  summarize(freq = n())

#Alluvial plot for change in smoking status
ggplot(as.data.frame(curr_currresurv), aes(axis1 = cig, axis2 = cig_resurv, y = freq)) +
  geom_alluvium(aes(fill = cig)) +
  scale_x_discrete(limits = c("Baseline", "Resurvey"),expand = c(.15, .15)) +
  geom_stratum(alpha = .5) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_classic() +
  scale_fill_brewer(type = "qual", palette = "Set1", direction = 1)

#Convert no diag and R_HXDIAB to factor for alluvial plots
diab_resurv$R_HXDIAB <- if_else(is.na(diab_resurv$R_HXDIAB), 0, diab_resurv$R_HXDIAB)
diab_resurv$R_HXDIAB <- factor(diab_resurv$R_HXDIAB, labels = c("Undiagnosed", "Diagnosed"))
diab_resurv$nodiag <- factor(diab_resurv$nodiag, labels = c("Undiagnosed", "Diagnosed"))

#Change in diagnosis status from baseline to resurvey
dig_nodig <- diab_resurv %>% 
  group_by(nodiag, R_HXDIAB) %>% 
  filter(!is.na(nodiag)) %>% 
  summarize(freq = n())

supp1_1 <- ggplot(as.data.frame(dig_nodig), aes(axis1 = nodiag, axis2 = R_HXDIAB, y = freq)) +
  geom_alluvium(aes(fill = nodiag)) +
  scale_x_discrete(limits = c("Baseline", "Resurvey"),expand = c(.15, .15)) +
  geom_stratum(alpha = .5) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_classic() +
  scale_fill_brewer(type = "qual", palette = "Set1", direction = 1)

#Change in smoking status in undignosed current smokers
resurv_undiag <- diab_resurv %>% 
  filter(nodiag == "Undiagnosed" & R_HXDIAB == "Diagnosed") %>% 
  group_by(cig, cig_resurv) %>% 
  summarize(n = n()) %>% 
  filter(!is.na(cig_resurv)) %>% 
  filter(cig == "Current smoker" | cig == "Former smoker")

supp1_2_undiag <- ggplot(as.data.frame(resurv_undiag), aes(axis1 = cig, axis2 = cig_resurv, y = n)) +
  geom_alluvium(aes(fill = cig)) +
  scale_x_discrete(limits = c("Baseline", "Resurvey"),expand = c(.15, .15)) +
  geom_stratum(alpha = .5) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_classic() +
  scale_fill_brewer(type = "qual", palette = "Set1", direction = 1) +
  ylab("") +
  theme(legend.title=element_blank())

#Change in smoking status in diagnosed current smokers
resurv_diag <- diab_resurv %>% 
  filter(nodiag == "Diagnosed" & R_HXDIAB == "Diagnosed") %>%
  group_by(cig, cig_resurv) %>% 
  summarize(n = n()) %>% 
  filter(!is.na(cig_resurv)) %>% 
  filter(cig == "Current smoker" | cig == "Former smoker")

supp1_3_diag <- ggplot(as.data.frame(resurv_diag), aes(axis1 = cig, axis2 = cig_resurv, y = n)) +
  geom_alluvium(aes(fill = cig)) +
  scale_x_discrete(limits = c("Baseline", "Resurvey"),expand = c(.15, .15)) +
  geom_stratum(alpha = .5) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_classic() +
  scale_fill_brewer(type = "qual", palette = "Set1", direction = 1) +
  ylab("") +
  theme(legend.title=element_blank())

supp1 <- ggarrange(supp1_2_undiag, supp1_3_diag, nrow = 1, ncol = 2, labels = "AUTO", common.legend = TRUE, legend = "bottom")

ggsave(supp1, file="Proyectos/Smoking and diabetes/Suppl1.jpg", bg="transparent",
       width=22, height=10, units=c("cm"), dpi=600, limitsize = FALSE)

