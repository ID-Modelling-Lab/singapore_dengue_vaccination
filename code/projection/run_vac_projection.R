
### call all the libraries
library(odin.dust)
library(readxl)
library(dplyr)
library(tidyr)
library(here)
library(foreach)
library(doParallel)


# will remove everything from workspace
rm(list = ls(all = TRUE)) 
# collects garbage memory
gc(reset = TRUE)        

## cal all the functions from utils
source(here("utils","utils_demography.R"))
source(here("utils","utils_sero_prevalence.R"))

## call all the functions related to demography
begin_year = "2013"
end_year = "2021"

# for sero-prevalence data
sero_year = 2013

n_year <- (as.numeric(end_year) - as.numeric(begin_year)) + 1

## number of extra year for which the model will be simulated after 2021
extra_year <- 70

sim_year <- n_year + extra_year

## create array of time in day (required for solving ode)
t <- seq(1,(sim_year*365), by = 1)


## population in each age group in "begin_year"
pop_each_age_group <- get_pop_each_group(begin_year = begin_year,
                                         end_year = end_year )

## birth rate over time
lambda_hum_tt = get_birth_rate(begin_year = begin_year, 
                               end_year = end_year, 
                               sim_year = sim_year)

## death rate over time
death_hum_tt = get_death_rate(begin_year = begin_year, 
                              end_year = end_year, 
                              sim_year = sim_year)

## sero-prevalence data
get_sero <- get_sero_p(sero_year = sero_year )

## sero-positive fraction i.e at least one exposure
sero_p <- get_sero$se_p
## exactly one exposure
primary <- get_sero$pri


##### Demographic
## Number of age group
n_age <- length(pop_each_age_group) 
## 
age_width <- 1
## Transition rate between age groups
age_rate <- c(rep(1/(365*age_width), (n_age - 1)),0)
## mosquito dearth rate
death_mos = 1/10
## Adult mosquito recruitment rate
lambda_mos = death_mos 



get_ve_point_estimate_from_trial <- function(df_efficacy) {
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












## Number of serotype
n_sero <- 4
## read the excel of efficacy of qdenga data

data_efficacy = read_xlsx( here("data","qdenga_efficacy.xlsx"), 
                           sheet = "Sheet1" )


# data_efficacy <- read_xlsx("~/codes/dengue_main/data/qdenga_efficacy.xlsx",sheet="Sheet1")
##
sero_strat_ve <- get_ve_from_trial(df_efficacy = data_efficacy )
n_vac_stage <- length(unique(sero_strat_ve$month))


vac_stage_rate <- c( rep(1/(12*30),3), 0)




###### Initial conditions

## (Primary Infection)
I0 = matrix(0, nrow = n_age, ncol = n_sero, byrow = TRUE)
## Susceptible secondary infection
frac_serotype <-  c(.352860, .490160, .111505,  .045475) 

## Cross-immune
C0 = 0.0*t(frac_serotype*t(array(rep(primary*pop_each_age_group,n_sero), 
                                 dim = c(n_age,n_sero))))
S0 = 1*t(frac_serotype*t(array(rep(primary*pop_each_age_group,n_sero),
                               dim = c(n_age,n_sero))))
##(Secondary infection)
I_ij0 = array(0,dim = c(n_age, n_sero))

## Recovered from Secondary infection
R0 =  sero_p*pop_each_age_group - rowSums(S0) - rowSums(C0)

### Vaccinated population
S0_v_all = array(0, dim = c(n_vac_stage, n_age))
I_v0 = array(0, dim = c(n_vac_stage, n_age, n_sero))
C_v0 = array(0, dim = c(n_vac_stage, n_age, n_sero))
S_v0 = array(0, dim = c(n_vac_stage, n_age, n_sero))
I_v_ij0 = array(0, dim = c(n_vac_stage, n_age, n_sero))
R_v0 = array(0, dim = c(n_vac_stage, n_age))
## Susceptible population
S0_all =  pop_each_age_group - rowSums(I0 + C0 + S0) - 
  rowSums(I_ij0) - (R0)

factor <- .01
total_human <- sum(pop_each_age_group)
## Susceptible mosquito
S_m0 = factor*total_human
# ## Exposed Mosquito
E_m0 = c(0,0,0,0) 
# ## infected mosquito
I_m0 = I_m0 = c(2e-2,5e-3,1e-2,1e-6)   #c(1e-1,1e-1,1e-3,1e-5) 

### yearly sero function

## function for calculating yearly infection serotype wise (cumulative over day)
## this is for new cases/vaccinated
calc_year_sero <- function(arr) {
  
  ## this is for primary infection, where input has only 3 dimensions: (age, serotype, time)
  if (length(dim(arr)) == 3) {
    arr_sero <- apply(arr, c(2,3), sum) ## summing over age, i.e., 1st dimension
    arr_sero_yearly <- t(apply(arr_sero, 1, calc_yearly))
  }
  return(arr_sero_yearly)
}


## this is for prevalence of cases/vaccinated/population
calc_year_sero_wo_sum <- function(arr){
  
  if (length(dim(arr)) == 2){
    arr_sero <- arr
    arr_sero_yearly <- t(apply(arr_sero, 1, calc_yearly_wo_sum))
  }
  ## this is for primary infection, where input has only 3 dimensions: (age, serotype, time)
  if (length(dim(arr)) == 3){
    arr_sero <- apply(arr, c(2,3), sum) ## summing over age, i.e., 1st dimension
    arr_sero_yearly <- t(apply(arr_sero, 1, calc_yearly_wo_sum))
  }
  return(arr_sero_yearly)
}


## function for calculating yearly infection age wise
calc_year_age <- function(arr){
  
  if (length(dim(arr)) == 2) {
    arr_age <- arr
    arr_age_yearly <- t(apply(arr_age, 1, calc_yearly))
  }
  
  ## where input has only 3 dimensions: (age, serotype, time)
  if (length(dim(arr)) == 3) {
    arr_age <- apply(arr, c(1,3), sum) ## summing over sero, i.e., 2nd dimension
    arr_age_yearly <- t(apply(arr_age, 1, calc_yearly))
  }
  return(arr_age_yearly)
}

calc_year_age_wo_sum <- function(arr) {
  
  if (length(dim(arr)) == 2) {
    arr_age <- arr
    arr_age_yearly <- t(apply(arr_age, 1, calc_yearly_wo_sum))
  }
  
  ## where input has only 3 dimensions: (age, serotype, time)
  if (length(dim(arr)) == 3) {
    arr_age <- apply(arr, c(1,3), sum) ## summing over sero, i.e., 2nd dimension
    arr_age_yearly <- t(apply(arr_age, 1, calc_yearly_wo_sum))
  }
  return(arr_age_yearly)
}


## function for calculating yearly infection
calc_yearly <- function(mat) {
  tapply(mat, (seq_along(mat) - 1) %/%365, sum)
}

calc_yearly_wo_sum <- function(mat) {
  mat[seq(1,length(t),by = 365)]
}


###### Transmission related parameters
## biting rate per mosquito
b = 15
# Transmission rate from infected human to mosquito
beta_h = c(0.20, 0.20, 0.20, 0.20)   
#
beta_m = c(0.166, 0.176, 0.152, 0.144)
# Infectious period primary infection
gamma_1 = 1/4
# Infectious period secondary infection
gamma_2 = 1/4
# infectiousness of asymptomatic relative to symptomatic
# kai = 1
# Cross immunity period
alpha = 1/(1*365)
# rate of hospitalization for primary infection
xi_1 = 0.30/4
# rate of hospitalization for secondary infection
xi_2 = 0.30
# extrinsic incubation period
sigma_m = 1/10



## sampling from posterior
n_sample <- 20
post <- readRDS(here("model_output", "posterior_rep_rate.rds"))

sampled_row <- sample(1:nrow(post), n_sample)

rho_2_sample <- post[sampled_row,]


get_post_sample <- function(k) {
  
  rho_2_a1 = rho_2_sample[k,1]
  rho_2_a2 = rho_2_sample[k,2]
  rho_2_a3 = rho_2_sample[k,3]
  rho_2_a4 = rho_2_sample[k,4]
  rho_2_a5 = rho_2_sample[k,5]
  rho_2_a6 = rho_2_sample[k,6]
  rho_2_a7 = rho_2_sample[k,7]
  rho_2_a8 = rho_2_sample[k,8]
  
  rho_2 <- rep(NA, n_age)
  
  rho_2[1:5] <- rho_2_a1
  rho_2[6:15] <- rho_2_a2
  rho_2[16:25] <- rho_2_a3
  rho_2[26:35] <- rho_2_a4
  rho_2[36:45] <- rho_2_a5
  rho_2[46:55] <- rho_2_a6
  rho_2[56:65] <- rho_2_a7
  rho_2[66:91] <- rho_2_a8
  
  return(list(rho_2 = rho_2))
}



# compile the model to run in parallel
out_vaccine = function(par) {
  
  path_to_model <-  here::here("projection", "model_vaccine_waning.R") 
  dengue_vaccine_model <- odin.dust::odin_dust(path_to_model) 
  
  #### instance of the model
  model <- dengue_vaccine_model$new(par, time = 1, n_particles = 1,
                                    ode_control = list(max_steps = 10000000,
                                                       step_size_min = 1e-14,
                                                       debug_record_step_times = TRUE))
  
  # print(1)
  ## simulate the model
  out <- model$simulate(t)
  
  
  # for cross-check that model works well
  # there should not be any change in infection due to vaccination
  ind_inf <- model$info()$index$inc_inf
  ind_inf_v <- model$info()$index$inc_inf_v
  
  ind_symp <- model$info()$index$inc_symp
  ind_symp_v <- model$info()$index$inc_symp_v
  ind_hosp <- model$info()$index$inc_hosp
  ind_hosp_v <- model$info()$index$inc_hosp_v
  ind_tot_age <- model$info()$index$tot_age
  ind_N_age_v <- model$info()$index$tot_vac_age
  ind_N_age_nv <- model$info()$index$tot_nonvac_age
  ind_new_vac <- model$info()$index$new_vac_age
  
  tot_pop_age <- array(out[ind_tot_age,,], dim = c(n_age, length(t)))                           ## pop in each age-group
  tot_pop_age_yearly <- calc_year_age_wo_sum(tot_pop_age)      ## yearly pop. in each group
  
  ## vaccinated population
  tot_pop_age_v <- array(out[ind_N_age_v,,], dim = c(n_age, length(t))) ## vacc. pop in each age-group
  tot_pop_age_v_yearly <- calc_year_age_wo_sum(tot_pop_age_v)  ## yearly vacc. pop. in each group
  
  ## unvaccinated population
  tot_pop_age_nv <- array(out[ind_N_age_nv,,], dim = c(n_age, length(t))) ## vacc. pop in each age-group
  tot_pop_age_nv_yearly <- calc_year_age_wo_sum(tot_pop_age_nv)  ## yearly vacc. pop. in each group
  
  ## incidence of infection
  inf_age_sero <- array(out[ind_inf,,], dim = c(n_age, n_sero, length(t)))  ## for unvaccinated
  inf_age_yearly <- calc_year_age(inf_age_sero)                             ## make this yearly
  
  inf_age_sero_v <- array(out[ind_inf_v,,], dim = c(n_age, n_sero, length(t))) ## for vaccinated
  inf_age_v_yearly <- calc_year_age(inf_age_sero_v)                            ## yearly
  
  ## incidence of symptomatic cases
  symp_age_sero <- array(out[ind_symp,,], dim = c(n_age, n_sero,length(t)))
  symp_age_yearly <- calc_year_age(symp_age_sero)
  
  symp_age_sero_v <- array(out[ind_symp_v,,], dim = c(n_age, n_sero,length(t)))
  symp_age_v_yearly <- calc_year_age(symp_age_sero_v)
  
  ## incidence of hospitalized cases
  hosp_age <- array(out[ind_hosp,,], dim = c(n_age, length(t)))
  hosp_age_yearly <- calc_year_age(hosp_age)
  
  hosp_age_v <- array(out[ind_hosp_v,,], dim = c(n_age, length(t)))
  hosp_age_v_yearly <- calc_year_age(hosp_age_v)
  
  new_vac_age <- array(out[ind_new_vac,,], dim = c(n_age, length(t)))
  new_vac_age_yearly <- calc_year_age(new_vac_age)
  
  rm(list = c("out"))
  
  return(list(tot_pop_age = tot_pop_age_yearly,
              tot_pop_age_v = tot_pop_age_v_yearly,
              tot_pop_age_nv = tot_pop_age_nv_yearly,
              inf_age = inf_age_yearly,
              inf_age_v = inf_age_v_yearly,
              symp_age = symp_age_yearly,
              symp_age_v = symp_age_v_yearly,
              hosp_age = hosp_age_yearly,
              hosp_age_v = hosp_age_v_yearly,
              new_vac_age = new_vac_age_yearly
  ))
}



### nested looping starts here

# Define the number of cores to use
numcores <- 5 # detectCores() - 1  

cl <- makeCluster(numcores)  # Adjust the number of cores based on your HPC setup
registerDoParallel(cl)

# 
# registerDoParallel(numCores)


## Coverage values
v_coverage = 0.0


## indices for vaccinated age-groups
v_a = matrix( c(1,5,6,10,
                11,15,16,20,
                21,25,26,30,
                31,35,36,40,
                41,45,46,50,
                51,55,56,60,
                61,65,66,70,
                71,75,76,80,
                81,85
), nrow = 17, ncol = 2, byrow = TRUE)

## number vaccine coverage scenarios
n_vac_coverage = length(v_coverage)
## number of vaccinated age-groups
n_vac_age = dim(v_a)[1]
### Space for Outputs

## Create space for output quantities
inf_sample_age = array(NA, dim = c(n_sample,n_age, sim_year))
inf_sample_age_v = array(NA, dim = c(n_sample,n_age, sim_year))

symp_sample_age = array(NA, dim = c(n_sample,n_age, sim_year))
symp_sample_age_v = array(NA, dim = c(n_sample,n_age, sim_year))

hosp_sample_age =  array(NA, dim = c(n_sample, n_age,sim_year))
hosp_sample_age_v =  array(NA, dim = c(n_sample,n_age, sim_year))

tot_pop_sample_age = array(NA, dim = c(n_sample, n_age, sim_year))
tot_pop_sample_age_v = array(NA, dim = c(n_sample, n_age, sim_year))
tot_pop_sample_age_nv = array(NA, dim = c(n_sample,  n_age, sim_year))
new_vac_sample_age =  array(NA, dim = c(n_sample, n_age,sim_year))


### save all the sampled efficacy estimates

effi_symp_sero_n_sample = array(NA, dim = c(n_sample, n_vac_stage, n_sero))
effi_symp_sero_p_sample = array(NA, dim = c(n_sample, n_vac_stage, n_sero))
effi_hos_sero_n_sample = array(NA, dim = c(n_sample, n_vac_stage, n_sero))
effi_hos_sero_p_sample = array(NA, dim = c(n_sample, n_vac_stage, n_sero))

results_list <- vector("list", length(v_coverage))

## Run for different vaccination coverage and for different targeted age-groups
## Looping over coverage (i), vaccinated age groups (j), sample(k)

v_year <- 10
vac_start <- v_year*365  ## year of vaccination start
lower_bound <- -0.25  ## lower bound of negative efficacy 


### Run in parallel



### make vaccinated age group run in parallel
for (i in 1:length(v_coverage)) {
  set.seed(12345)
  vac_coverage <- v_coverage[i]
  # print(i)
  
  out_parallel <- foreach(j =  1:n_vac_age, .combine = rbind, .packages = c("dplyr", "tidyr", "here", "odin.dust")) %dopar% {
    
    library(dplyr)
    library(tidyr)
    library(here)
    library(odin.dust)
    
    
    vac_switch <- rep(0,n_age)
    vac_switch[v_a[j,1]:v_a[j,2]] <- 1
    
    for (k in 1:n_sample) {
      print(k)
      
      if (vac_coverage == 0.0) {
        effi_inf_sero_n = array(0, dim=c(n_vac_stage))
        effi_inf_sero_p = array(0, dim=c(n_vac_stage))
        effi_symp_sero_n = array(0, dim=c(n_vac_stage, n_sero))
        effi_symp_sero_p = array(0, dim=c(n_vac_stage, n_sero))
        effi_hos_sero_n = array(0, dim=c(n_vac_stage, n_sero))
        effi_hos_sero_p = array(0, dim=c(n_vac_stage, n_sero))
        
        
      }
      
      else {
      
      
      v_eff <- get_conditional_ve(df_efficacy = sero_strat_ve, 
                                  n_vac_stage = n_vac_stage , 
                                  n_sero = n_sero, lower_bound = lower_bound
                                  )
      effi_inf_sero_n = v_eff$inf_n
      effi_inf_sero_p = v_eff$inf_p
      effi_symp_sero_n = v_eff$symp_n
      effi_symp_sero_p = v_eff$symp_p
      effi_hos_sero_n = v_eff$hos_n
      effi_hos_sero_p = v_eff$hos_p
      }
      
      
      
      
      effi_symp_sero_n_sample[k,,] = effi_symp_sero_n
      effi_symp_sero_p_sample[k,,] = effi_symp_sero_p
      effi_hos_sero_n_sample[k,,] = effi_hos_sero_n
      effi_hos_sero_p_sample[k,,] = effi_hos_sero_p
      
      ratio_p_s <- 1/2
      
      rho_2 <- get_post_sample(k)$rho_2
      rho_1 <- ratio_p_s*rho_2
      
      par <- list(
        n_age = n_age,
        n_sero = n_sero,
        beta_h = beta_h,
        beta_m = beta_m,
        b = b,
        vac_coverage = vac_coverage,
        vac_switch = vac_switch,
        vac_start = vac_start,
        n_vac_stage = n_vac_stage,
        vac_stage_rate = vac_stage_rate,
        rho_1 = rho_1,
        rho_2 = rho_2,
        sigma_m = sigma_m,
        alpha = alpha,
        gamma_1 = gamma_1,
        gamma_2 = gamma_2,
        xi_1 = xi_1,
        xi_2 = xi_2,
        
        lambda_hum_tt = lambda_hum_tt,
        death_hum_tt = death_hum_tt,
        
        lambda_mos = lambda_mos,
        death_mos = death_mos,
        age_rate = age_rate,
        
        I0 = I0,
        C0 = C0,
        S0 = S0,
        I_ij0 = I_ij0,
        R0 = R0,
        S0_all = S0_all,
        S_m0 = S_m0,
        E_m0 = E_m0,
        I_m0 = I_m0,
        
        S0_v_all = S0_v_all,
        I_v0 = I_v0,
        C_v0 = C_v0,
        S_v0 = S_v0,
        I_v_ij0 =  I_v_ij0,
        R_v0 = R_v0,
        
        effi_inf_sero_n = effi_inf_sero_n,
        effi_inf_sero_p = effi_inf_sero_p,
        effi_symp_sero_n = effi_symp_sero_n,
        effi_symp_sero_p = effi_symp_sero_p,
        effi_hos_sero_n = effi_hos_sero_n,
        effi_hos_sero_p = effi_hos_sero_p
      )
      
      xx <-  out_vaccine(par)
      
      ## accumulate the output as a matrix
      tot_pop_sample_age[k,,] <- xx$tot_pop_age
      tot_pop_sample_age_v[k,,] <- xx$tot_pop_age_v
      tot_pop_sample_age_nv[k,,] <- xx$tot_pop_age_nv
      
      inf_sample_age[k,,] <- xx$inf_age
      inf_sample_age_v[k,,] <- xx$inf_age_v
      
      symp_sample_age[k,,] <- xx$symp_age
      symp_sample_age_v[k,,] <- xx$symp_age_v
      
      hosp_sample_age[k,,] <- xx$hosp_age
      hosp_sample_age_v[k,,] <- xx$hosp_age_v
      
      new_vac_sample_age[k,,] <- xx$new_vac_age
      
      rm(list=c("xx"))
    }
    return(list(tot_pop_sample_age = tot_pop_sample_age,
                tot_pop_sample_age_v = tot_pop_sample_age_v,
                tot_pop_sample_age_nv = tot_pop_sample_age_nv,
                inf_sample_age = inf_sample_age,
                inf_sample_age_v = inf_sample_age_v,
                symp_sample_age = symp_sample_age,
                symp_sample_age_v = symp_sample_age_v,
                hosp_sample_age = hosp_sample_age,
                hosp_sample_age_v = hosp_sample_age_v,
                new_vac_sample_age = new_vac_sample_age,
                eff = list(
                  effi_symp_sero_n = effi_symp_sero_n_sample,
                  effi_symp_sero_p = effi_symp_sero_p_sample,
                  effi_hosp_sero_n = effi_hos_sero_n_sample,
                  effi_hosp_sero_p = effi_hos_sero_p_sample)
    ))
  }
  
  
  results_list[[i]] <- out_parallel
  
  rm(out_parallel)
  gc()
}



## TO add v_a in the output list also start of vaccination year
saveRDS(list(coverage = v_coverage, 
             n_sample = n_sample,
             n_sero = n_sero,
             n_year = sim_year,
             n_vac_age = n_vac_age,
             n_age = n_age,
             v_year = v_year,
             
             output=results_list), 
        file = here::here("model_output", paste0("test_model_", "vcov_", v_coverage,
                                                 "_nsamp_", n_sample,"_lb_", 
                                                 lower_bound,"_ve_scenario1.rds" )))


stopCluster(cl)

closeAllConnections()


