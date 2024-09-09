library(readxl)
library(dplyr)
library(tidyr)
library(here)




## Vaccine related paramter
get_ve_from_trial <- function(df_efficacy) {
  ## Filter the VE based on overall and serostatus stratified but serotype stratified
  df_no_sero <- df_efficacy %>% filter(ve_against %in% 
                                         c("vcd_overall", "hosp_overall", "vcd_sp", "vcd_sn", "hosp_sp", "hosp_sn"))
  
  ## create two groups based VCD and Hosp & Overall or serostatus stratified
  
  df_no_sero_strat <- df_no_sero %>% 
    mutate( group1 = factor(case_when(
      ve_against %in% c("vcd_overall", "vcd_sp", "vcd_sn") ~ "VCD",
      ve_against %in% c("hosp_overall","hosp_sp", "hosp_sn") ~ "HOSPITALIZED"), 
      levels = c("VCD", "HOSPITALIZED")),
      
      group2 = factor(case_when(
        ve_against %in% c("vcd_overall", "hosp_overall") ~ "Overall",
        ve_against %in% c("vcd_sp","hosp_sp") ~ "Sero+",
        ve_against %in% c("vcd_sn","hosp_sn") ~ "Sero-"),
        levels = c("Overall", "Sero+", "Sero-"))
    )
  
  
  ## Now first we remove the values where VE == 100 and replace with the previous and next values
  ## Follow the same for the C.I values
  ## There is only one such instance
  # the value of ve which will replace ve==100
  replace_ve_hosp_sn <- df_no_sero_strat %>%
    filter(ve_against == "hosp_sn", month == 36) %>%
    pull(ve)
  
  # the value of lower bound
  replace_velow_hosp_sn <- df_no_sero_strat %>%
    filter(ve_against == "hosp_sn", month == 36) %>%
    pull(ve_low)
  
  # the value of upper bound
  replace_veup_hosp_sn <- df_no_sero_strat %>%
    filter(ve_against == "hosp_sn", month == 36) %>%
    pull(ve_up)
  
  
  # Modify the dataframe now replacing the VE==100 with above values
  df_no_sero_strat <- df_no_sero_strat %>%
    mutate(
      ## Replace where VE value is 100 with the previous or next value
      ve = ifelse(
        group1 == "HOSPITALIZED" & group2 %in% c("Sero-") & ve == 100,
        replace_ve_hosp_sn,
        ve
      ),
      
      ve_low = ifelse(
        group1 == "HOSPITALIZED" &
          group2 %in% c("Sero-") & is.na(ve_low) == TRUE,
        replace_velow_hosp_sn,
        ve_low
      ),
      
      ve_up = ifelse(
        group1 == "HOSPITALIZED" &
          group2 %in% c("Sero-") & is.na(ve_up) == TRUE,
        replace_veup_hosp_sn,
        ve_up
      )
    )
  
  
  ## Filter the VE based on serotypes
  
  df_sero <- df_efficacy %>%  filter(grepl("_dv[1234]$", ve_against))
  
  ## Create 3 groups
  df_sero_strat <- df_sero %>%
    mutate(
      group1 = factor(case_when(
        grepl("_dv1", ve_against, ignore.case = TRUE) ~ "DENV1",
        grepl("_dv2", ve_against, ignore.case = TRUE) ~ "DENV2",
        grepl("_dv3", ve_against, ignore.case = TRUE) ~ "DENV3",
        grepl("_dv4", ve_against, ignore.case = TRUE) ~ "DENV4"),
        levels = c("DENV1", "DENV2", "DENV3", "DENV4")),
      
      group2 = factor(case_when(
        grepl("_sp_", ve_against, ignore.case = TRUE) ~ "Sero+",
        grepl("_sn_", ve_against, ignore.case = TRUE) ~ "Sero-"),
        levels = c("Sero+", "Sero-")),
      
      group3 = factor(case_when(
        grepl("vcd_", ve_against, ignore.case = TRUE) ~ "VCD",
        grepl("hosp_", ve_against, ignore.case = TRUE) ~ "Hospitalized"),
        levels = c("VCD", "Hospitalized"))
      
    )
  
  ### Now from first data frame we will pull out the estimated of VE against HOSPITALIZED both for Sero+ and Sero-
  # Sero+
  pull_from_no_sero_ve_hosp_sp <- df_no_sero_strat %>%
    filter(ve_against == "hosp_sp") %>%
    pull(ve)
  
  pull_from_no_sero_velow_hosp_sp <- df_no_sero_strat %>%
    filter(ve_against == "hosp_sp") %>%
    pull(ve_low)
  
  pull_from_no_sero_veup_hosp_sp <- df_no_sero_strat %>%
    filter(ve_against == "hosp_sp") %>%
    pull(ve_up)
  
  # Sero- 
  pull_from_no_sero_ve_hosp_sn <- df_no_sero_strat %>%
    filter(ve_against == "hosp_sn") %>%
    pull(ve)
  
  pull_from_no_sero_velow_hosp_sn <- df_no_sero_strat %>%
    filter(ve_against == "hosp_sn") %>%
    pull(ve_low)
  
  pull_from_no_sero_veup_hosp_sn <- df_no_sero_strat %>%
    filter(ve_against == "hosp_sn") %>%
    pull(ve_up)
  
  ## replace where the ve is 100 with next or previous values
  # first instance
  replace_ve_vcd_sp_dv4 <- df_sero_strat %>%
    filter(ve_against == "vcd_sp_dv4", month == 36) %>%
    pull(ve)
  
  replace_velow_vcd_sp_dv4 <- df_sero_strat %>%
    filter(ve_against == "vcd_sp_dv4", month == 36) %>%
    pull(ve_low)
  
  replace_veup_vcd_sp_dv4 <- df_sero_strat %>%
    filter(ve_against == "vcd_sp_dv4", month == 36) %>%
    pull(ve_up)
  
  # second instance
  replace_ve_vcd_sn_dv2 <- df_sero_strat %>%
    filter(ve_against == "vcd_sn_dv2", month == 24) %>%
    pull(ve)
  
  replace_velow_vcd_sn_dv2 <- df_sero_strat %>%
    filter(ve_against == "vcd_sn_dv2", month == 24) %>%
    pull(ve_low)
  
  replace_veup_vcd_sn_dv2 <- df_sero_strat %>%
    filter(ve_against == "vcd_sn_dv2", month == 24) %>%
    pull(ve_up)
  
  # Modify df2 for hospitalization of DENV3 and DENV4 for sero+ and sero-
  df_sero_strat <- df_sero_strat %>%
    mutate(
      
      ## Ve against hosp for denv3 and denv4 both sero+ sero- are zero
      ve = ifelse(group3 == "Hospitalized" & group1 %in% c("DENV3", "DENV4") & group2 %in% c("Sero+", "Sero-"), 0, ve),
      ve_low = ifelse(group3 == "Hospitalized" & group1 %in% c("DENV3", "DENV4") & group2 %in% c("Sero+", "Sero-"), 0, ve_low),
      ve_up = ifelse(group3 == "Hospitalized" & group1 %in% c("DENV3", "DENV4") & group2 %in% c("Sero+", "Sero-"), 0, ve_up),
      
      ## VE against VCD for DENV3 and DENV4 for Sero- are zero
      ve = ifelse(group3 == "VCD" & group1 %in% c("DENV3", "DENV4") & group2 %in% c("Sero-"), 0, ve),
      ve_low = ifelse(group3 == "VCD" & group1 %in% c("DENV3", "DENV4") & group2 %in% c("Sero-"), 0, ve_low),
      ve_up = ifelse(group3 == "VCD" & group1 %in% c("DENV3", "DENV4") & group2 %in% c("Sero-"), 0, ve_up),
      
      ## Replace where VE value is 100 with the previous or next value
      # first instance of VE == 100
      ve = ifelse( group3 == "VCD" & group1 %in% c("DENV4") & group2 %in% c("Sero+") & ve == 100,
                   replace_ve_vcd_sp_dv4, ve),
      ve_low = ifelse( group3 == "VCD" & group1 %in% c("DENV4") & group2 %in% c("Sero+") & is.na(ve_low) == TRUE, 
                       replace_velow_vcd_sp_dv4, ve_low),
      ve_up = ifelse( group3 == "VCD" & group1 %in% c("DENV4") & group2 %in% c("Sero+") & is.na(ve_up) == TRUE,
                      replace_veup_vcd_sp_dv4, ve_up),
      
      ## second instance of VE == 100    
      ve = ifelse( group3 == "VCD" & group1 %in% c("DENV2") & group2 %in% c("Sero-") & ve == 100, 
                   replace_ve_vcd_sn_dv2, ve),
      ve_low = ifelse( group3 == "VCD" & group1 %in% c("DENV2") & group2 %in% c("Sero-") & is.na(ve_low) == TRUE,
                       replace_velow_vcd_sn_dv2, ve_low),
      ve_up = ifelse( group3 == "VCD" & group1 %in% c("DENV2") & group2 %in% c("Sero-") & is.na(ve_up) == TRUE, 
                      replace_veup_vcd_sn_dv2, ve_up),
      
      # Replace hospitalized VE from the unstratified ve estimate
      # seropositive
      ve = ifelse( group3 == "Hospitalized" & group1 %in% c("DENV1", "DENV2") & group2 %in% c("Sero+") ,
                   pull_from_no_sero_ve_hosp_sp, ve),
      ve_low = ifelse( group3 == "Hospitalized" & group1 %in% c("DENV1","DENV2") & group2 %in% c("Sero+") ,
                       pull_from_no_sero_velow_hosp_sp, ve_low),
      ve_up = ifelse( group3 == "Hospitalized" & group1 %in% c("DENV1","DENV2") & group2 %in% c("Sero+") ,
                      pull_from_no_sero_veup_hosp_sp, ve_up),
      
      # Seronegative
      ve = ifelse( group3 == "Hospitalized" & group1 %in% c("DENV1", "DENV2") & group2 %in% c("Sero-"),
                   pull_from_no_sero_ve_hosp_sn, ve),
      ve_low = ifelse( group3 == "Hospitalized" & group1 %in% c("DENV1","DENV2") & group2 %in% c("Sero-") ,
                       pull_from_no_sero_velow_hosp_sn, ve_low),
      ve_up = ifelse( group3 == "Hospitalized" & group1 %in% c("DENV1","DENV2") & group2 %in% c("Sero-") , 
                      pull_from_no_sero_veup_hosp_sn, ve_up)
    )
  
  return(df_sero_strat)
}





get_ve_recursive <- function(df_efficacy,
                             serotype,
                             sero_status,
                             ve_against,
                             month) {
  
  ve <- df_efficacy$ve[df_efficacy$group1 == serotype & 
                         df_efficacy$group2 == sero_status &
                         df_efficacy$group3 == ve_against &
                         df_efficacy$month == month]
  ve_low <- df_efficacy$ve_low[df_efficacy$group1 == serotype & 
                                 df_efficacy$group2 == sero_status &
                                 df_efficacy$group3 == ve_against &
                                 df_efficacy$month == month]
  ve_up <- df_efficacy$ve_up[df_efficacy$group1 == serotype & 
                               df_efficacy$group2 == sero_status &
                               df_efficacy$group3 == ve_against &
                               df_efficacy$month == month]
  return(list(
    ve = ve,
    ve_low = ve_low,
    ve_up = ve_up
  ))
  
}



get_conditional_ve <- function(df_efficacy, n_vac_stage,
                               n_sero, lower_bound) {
  
  month_list <- unique(df_efficacy$month)
  sero_list <- c("DENV1","DENV2", "DENV3", "DENV4")
  
  effi_inf_sero_n_array <- rep(0, n_vac_stage)
  effi_inf_sero_p_array <- rep(0, n_vac_stage)
  
  symp_low_sero_n <- array(NA, dim=c(n_vac_stage, n_sero))
  symp_up_sero_n <- array(NA, dim=c(n_vac_stage, n_sero))
  effi_symp_sero_n_array <- array(NA, dim=c(n_vac_stage, n_sero))
  
  symp_low_sero_p <- array(NA, dim=c(n_vac_stage, n_sero))
  symp_up_sero_p <- array(NA, dim=c(n_vac_stage, n_sero))
  effi_symp_sero_p_array <- array(NA, dim=c(n_vac_stage, n_sero))
  
  hosp_low_sero_n <- array(NA, dim=c(n_vac_stage, n_sero))
  hosp_up_sero_n <- array(NA, dim=c(n_vac_stage, n_sero))
  effi_hosp_sero_n_array <- array(NA, dim=c(n_vac_stage, n_sero))
  
  hosp_low_sero_p <- array(NA, dim=c(n_vac_stage, n_sero))
  hosp_up_sero_p <- array(NA, dim=c(n_vac_stage, n_sero))
  effi_hosp_sero_p_array <- array(NA, dim=c(n_vac_stage, n_sero))
  
  
  
  for (i in 1:n_vac_stage) {
    for (j in 1:n_sero) {
      
      # symptomatic
      symp_low_sero_n[i,j] <- (1/100)*get_ve_recursive(df_efficacy = df_efficacy,
                                                       serotype = sero_list[j],
                                                       sero_status = "Sero-",
                                                       ve_against = "VCD",
                                                       month = month_list[i])$ve_low
      
      
      
      symp_up_sero_n[i,j] <- (1/100)*get_ve_recursive(df_efficacy = df_efficacy,
                                                      serotype = sero_list[j],
                                                      sero_status = "Sero-",
                                                      ve_against = "VCD",
                                                      month = month_list[i])$ve_up
      
      symp_low_sero_p[i,j] <- (1/100)*get_ve_recursive(df_efficacy = df_efficacy,
                                                       serotype = sero_list[j],
                                                       sero_status = "Sero+",
                                                       ve_against = "VCD",
                                                       month = month_list[i])$ve_low
      
      symp_up_sero_p[i,j] <- (1/100)*get_ve_recursive(df_efficacy = df_efficacy,
                                                      serotype = sero_list[j],
                                                      sero_status = "Sero+",
                                                      ve_against = "VCD",
                                                      month = month_list[i])$ve_up
      # hospitalized
      hosp_low_sero_n[i,j] <- (1/100)*get_ve_recursive(df_efficacy = df_efficacy,
                                                       serotype = sero_list[j],
                                                       sero_status = "Sero-",
                                                       ve_against = "Hospitalized",
                                                       month = month_list[i])$ve_low
      
      hosp_up_sero_n[i,j] <- (1/100)*get_ve_recursive(df_efficacy = df_efficacy,
                                                      serotype = sero_list[j],
                                                      sero_status = "Sero-",
                                                      ve_against = "Hospitalized",
                                                      month = month_list[i])$ve_up
      
      hosp_low_sero_p[i,j] <- (1/100)*get_ve_recursive(df_efficacy = df_efficacy,
                                                       serotype = sero_list[j],
                                                       sero_status = "Sero+",
                                                       ve_against = "Hospitalized",
                                                       month = month_list[i])$ve_low
      
      hosp_up_sero_p[i,j] <- (1/100)*get_ve_recursive(df_efficacy = df_efficacy,
                                                      serotype = sero_list[j],
                                                      sero_status = "Sero+",
                                                      ve_against = "Hospitalized",
                                                      month = month_list[i])$ve_up
      
      
      
      if (j==3|j==4) {
        effi_symp_sero_n_array[i,j] <- 0
      }
      
      else {
        effi_symp_sero_n_array[i,j] <-  runif(n = 1, min = max(lower_bound,symp_low_sero_n[i,j]) , max = symp_up_sero_n[i,j] )
        # 1 - ((1 - runif(n = 1, 
        #                                                min = symp_low_sero_n[i,j] , 
        #                                                max = symp_up_sero_n[i,j] ) )/(1 - effi_inf_sero_n_array[i]))
      }
      
      
      effi_symp_sero_p_array[i,j] <- runif(n = 1, min = max(lower_bound,symp_low_sero_p[i,j]), max = symp_up_sero_p[i,j] )
      
      
      
      if (j==3|j==4) {
        effi_hosp_sero_n_array[i,j] <- 0
      }
      
      else {
        effi_hosp_sero_n_array[i,j] <- max(lower_bound, 1 - ((1 - runif(n = 1, 
                                                       min = max(lower_bound,hosp_low_sero_n[i,j]) , 
                                                       max = hosp_up_sero_n[i,j] ) )/(1 - effi_symp_sero_n_array[i,j])) )
      }
      
      if (j==3|j==4) {
        effi_hosp_sero_p_array[i,j] <- 0
      }
      else {
        effi_hosp_sero_p_array[i,j] <-  max(lower_bound, 1 - ((1 - runif(n = 1, 
                                                       min = max(lower_bound,hosp_low_sero_p[i,j]) , 
                                                       max = hosp_up_sero_p[i,j] ) )/(1 - effi_symp_sero_p_array[i,j])) )
      }
      
      
    }
    
  }
  
  
  return(list(inf_n = effi_inf_sero_n_array, 
              inf_p = effi_inf_sero_p_array,
              symp_n = effi_symp_sero_n_array, 
              symp_p = effi_symp_sero_p_array,
              hos_n = effi_hosp_sero_n_array, 
              hos_p = effi_hosp_sero_p_array))
  
}

data_efficacy = read_xlsx( here("data","qdenga_efficacy.xlsx"), 
                           sheet = "Sheet1" )
sero_strat_ve <- get_ve_from_trial(df_efficacy = data_efficacy )
n_vac_stage <- length(unique(sero_strat_ve$month))
n_sero <- 4
n_run <- 200



effi_symp_sero_n_sample = array(NA, dim = c(n_run, n_vac_stage, n_sero))
effi_symp_sero_p_sample = array(NA, dim = c(n_run, n_vac_stage, n_sero))
effi_hos_sero_n_sample = array(NA, dim = c(n_run, n_vac_stage, n_sero))
effi_hos_sero_p_sample = array(NA, dim = c(n_run, n_vac_stage, n_sero))


lower_bound = -0.5

for (i in 1:n_run) {



v_eff <- get_conditional_ve(df_efficacy = sero_strat_ve, 
                            n_vac_stage = n_vac_stage , 
                            n_sero = n_sero,
                            lower_bound = lower_bound)
effi_inf_sero_n = v_eff$inf_n
effi_inf_sero_p = v_eff$inf_p
effi_symp_sero_n = v_eff$symp_n
effi_symp_sero_p = v_eff$symp_p
effi_hos_sero_n = v_eff$hos_n
effi_hos_sero_p = v_eff$hos_p




effi_symp_sero_n_sample[i,,] = effi_symp_sero_n
effi_symp_sero_p_sample[i,,] = effi_symp_sero_p
effi_hos_sero_n_sample[i,,] = effi_hos_sero_n
effi_hos_sero_p_sample[i,,] = effi_hos_sero_p

}


## plot
boxplot(effi_hos_sero_n_sample[,,1])
