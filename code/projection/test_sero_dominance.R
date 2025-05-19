
### call all the libraries
library(odin.dust)
library(readxl)
library(dplyr)
library(tidyr)
library(here)
library(foreach)
library(doParallel)
library(sn) ## for skew normal distribution
# library(doRNG)


# will remove everything from workspace
rm(list = ls(all = TRUE)) 
# collects garbage memory
gc(reset = TRUE)        

## cal all the functions from utils
source(here("utils","utils_demography.R"))
source(here("utils","utils_sero_prevalence.R"))



# ### Accept command-line arguments
# args <- commandArgs(trailingOnly = TRUE)
# 
# # Default values for the parameters
# baseline <- "yes"  # Default: "yes"
# v_coverage <- 0.8  # Default: 0.8
# 
# # Assign values from command-line arguments, if provided
# if (length(args) >= 1) {
#   baseline <- args[1]  # First argument: baseline
# }
# if (length(args) >= 2) {
#   v_coverage <- as.numeric(args[2])  # Second argument: coverage
# }
# 
# ### Check for valid inputs
# if (!baseline %in% c("yes", "no")) {
#   stop("Invalid value for baseline. Must be 'yes' or 'no'.")
# }
# if (is.na(v_coverage) || v_coverage < 0 || v_coverage > 1) {
#   stop("Invalid value for v_coverage. Must be a number between 0 and 1.")
# }
# 
# cat("Running simulation with parameters:\n")
# cat("  baseline =", baseline, "\n")
# cat("  v_coverage =", v_coverage, "\n")
# 
# # The rest of your script remains the same




## call all the functions related to demography
begin_year = "2013"
end_year = "2021"

# for sero-prevalence data
sero_year = 2013

n_year <- (as.numeric(end_year) - as.numeric(begin_year)) + 1

## number of extra year for which the model will be simulated after 2021
vac_program_year <- 60

sim_year <- n_year + vac_program_year

## create array of time in day (required for solving ode)
t <- seq(1,(sim_year*365), by = 1)


## Coverage values should not be zero, for baseine the efficacy should be zero
# v_coverage = 0.2
# baseline = "no" ## should 


v_year <- n_year + 1
vac_start <- v_year*365  ## year of vaccination start

v_a = matrix( c(7,17,
                18,31,
                32,41,
                42,51, 
                52,61,
                62,71,
                72,81
                
), nrow = 7, ncol = 2, byrow = TRUE)


## number vaccine coverage scenarios
n_vac_coverage = length(v_coverage)
## number of vaccinated age-groups
n_vac_age = dim(v_a)[1]


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



get_raw_ve_from_trial <- function(excel_efficacy) {
  
  df_sero <- excel_efficacy %>%  filter(grepl("_dv[1234]$", ve_against))
  
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
  
  return(df_sero_strat)
  
}




get_ve_recursive <- function(df_efficacy,
                             serotype,
                             sero_status,
                             ve_against
                             # month
) {
  
  ve <- df_efficacy$ve[df_efficacy$group1 == serotype & 
                         df_efficacy$group2 == sero_status &
                         df_efficacy$group3 == ve_against 
  ]
  ve_low <- df_efficacy$ve_low[df_efficacy$group1 == serotype & 
                                 df_efficacy$group2 == sero_status &
                                 df_efficacy$group3 == ve_against 
  ]
  ve_up <- df_efficacy$ve_up[df_efficacy$group1 == serotype & 
                               df_efficacy$group2 == sero_status &
                               df_efficacy$group3 == ve_against 
  ]
  return(list(
    ve = ve,
    ve_low = ve_low,
    ve_up = ve_up
  ))
  
}



get_pt_lb_ub_from_excel_as_array <- function(df_efficacy){
  
  
  sero_list <- c("DENV1", "DENV2", "DENV3", "DENV4")
  
  symp_low_sero_n <- array(NA, dim = c(1, n_sero))
  symp_up_sero_n <- array(NA, dim = c(1, n_sero))
  symp_pt_sero_n <- array(NA, dim = c(1, n_sero))
  
  symp_low_sero_p <- array(NA, dim = c(1, n_sero))
  symp_up_sero_p <- array(NA, dim = c(1, n_sero))
  symp_pt_sero_p <- array(NA, dim = c(1, n_sero))
  
  hosp_low_sero_n <- array(NA, dim = c(1, n_sero))
  hosp_up_sero_n <- array(NA, dim = c(1, n_sero))
  hosp_pt_sero_n <- array(NA, dim = c(1, n_sero))
  
  hosp_low_sero_p <- array(NA, dim = c(1, n_sero))
  hosp_up_sero_p <- array(NA, dim = c(1, n_sero))
  hosp_pt_sero_p <- array(NA, dim = c(1, n_sero))
  
  
  for (i in 1:n_sero) {
    ## symp sero-
    symp_low_sero_n[i] <- (1 / 100) * get_ve_recursive(
      df_efficacy = df_efficacy,
      serotype = sero_list[i],
      sero_status = "Sero-",
      ve_against = "VCD"
    )$ve_low
    
    
    symp_pt_sero_n[i] <- (1 / 100) * get_ve_recursive(
      df_efficacy = df_efficacy,
      serotype = sero_list[i],
      sero_status = "Sero-",
      ve_against = "VCD"
    )$ve
    # print(symp_pt_sero_n)
    symp_up_sero_n[i] <- (1 / 100) * get_ve_recursive(
      df_efficacy = df_efficacy,
      serotype = sero_list[i],
      sero_status = "Sero-",
      ve_against = "VCD"
    )$ve_up
    
    
    # symp sero+
    symp_low_sero_p[i] <- (1 / 100) * get_ve_recursive(
      df_efficacy = df_efficacy,
      serotype = sero_list[i],
      sero_status = "Sero+",
      ve_against = "VCD"
    )$ve_low
    
    symp_pt_sero_p[i] <- (1 / 100) * get_ve_recursive(
      df_efficacy = df_efficacy,
      serotype = sero_list[i],
      sero_status = "Sero+",
      ve_against = "VCD"
    )$ve
    
    symp_up_sero_p[i] <- (1 / 100) * get_ve_recursive(
      df_efficacy = df_efficacy,
      serotype = sero_list[i],
      sero_status = "Sero+",
      ve_against = "VCD"
    )$ve_up
    
    
    # hospitalized sero-
    hosp_low_sero_n[i] <- (1 / 100) * get_ve_recursive(
      df_efficacy = df_efficacy,
      serotype = sero_list[i],
      sero_status = "Sero-",
      ve_against = "Hospitalized"
    )$ve_low
    
    hosp_pt_sero_n[i] <- (1 / 100) * get_ve_recursive(
      df_efficacy = df_efficacy,
      serotype = sero_list[i],
      sero_status = "Sero-",
      ve_against = "Hospitalized"
    )$ve
    
    hosp_up_sero_n[i] <- (1 / 100) * get_ve_recursive(
      df_efficacy = df_efficacy,
      serotype = sero_list[i],
      sero_status = "Sero-",
      ve_against = "Hospitalized"
    )$ve_up
    
    
    # hospitalized seo+
    hosp_low_sero_p[i] <- (1 / 100) * get_ve_recursive(
      df_efficacy = df_efficacy,
      serotype = sero_list[i],
      sero_status = "Sero+",
      ve_against = "Hospitalized"
    )$ve_low
    
    hosp_pt_sero_p[i] <- (1 / 100) * get_ve_recursive(
      df_efficacy = df_efficacy,
      serotype = sero_list[i],
      sero_status = "Sero+",
      ve_against = "Hospitalized"
    )$ve
    
    hosp_up_sero_p[i] <- (1 / 100) * get_ve_recursive(
      df_efficacy = df_efficacy,
      serotype = sero_list[i],
      sero_status = "Sero+",
      ve_against = "Hospitalized"
    )$ve_up
    
  }
  
  
  return(list( symp_low_sero_n = symp_low_sero_n,
               symp_pt_sero_n = symp_pt_sero_n,
               symp_up_sero_n = symp_up_sero_n,
               symp_low_sero_p = symp_low_sero_p,
               symp_pt_sero_p = symp_pt_sero_p,
               symp_up_sero_p = symp_up_sero_p,
               
               hosp_low_sero_n = hosp_low_sero_n,
               hosp_pt_sero_n = hosp_pt_sero_n,
               hosp_up_sero_n = hosp_up_sero_n,
               hosp_low_sero_p = hosp_low_sero_p,
               hosp_pt_sero_p = hosp_pt_sero_p,
               hosp_up_sero_p = hosp_up_sero_p
  ))
  
} 



get_sample_from_ci <- function(symp_low_sero_n,
                               symp_pt_sero_n,
                               symp_up_sero_n,
                               symp_low_sero_p,
                               symp_pt_sero_p,
                               symp_up_sero_p,
                               
                               hosp_low_sero_n,
                               hosp_pt_sero_n,
                               hosp_up_sero_n,
                               hosp_low_sero_p,
                               hosp_pt_sero_p,
                               hosp_up_sero_p) {
  
  # 
  ## store for samples
  samp_effi_symp_n <- array(NA, dim = c(1, n_sero))
  samp_effi_symp_p <- array(NA, dim = c(1, n_sero))
  samp_effi_hosp_n <- array(NA, dim = c(1, n_sero))
  samp_effi_hosp_p <- array(NA, dim = c(1, n_sero))
  
  
  ### get the samples from skew normal distribution
  
  ### for infection (not sero-stratified)
  if (ve_inf_sero_n_pt == 0 & ve_inf_sero_n_low == 0 & ve_inf_sero_n_up == 0 ) {
    samp_effi_inf_n = 0} else {
      samp_effi_inf_n <- sample_from_skew_normal(pt = ve_inf_sero_n_pt, 
                                                 low_ci = ve_inf_sero_n_low,
                                                 up_ci = ve_inf_sero_n_up,
                                                 n_sample = 1)
    }
  
  ## if efficcay against infection then zero (no need to sample)    
  if (ve_inf_sero_p_pt == 0 & ve_inf_sero_p_low == 0 & ve_inf_sero_p_up == 0 ) {
    samp_effi_inf_p = 0} else {
      samp_effi_inf_p <- sample_from_skew_normal(pt = ve_inf_sero_p_pt, 
                                                 low_ci = ve_inf_sero_p_low,
                                                 up_ci = ve_inf_sero_p_up,
                                                 n_sample = 1)
    }
  
  
  ### for VCD and Hospitalized
  for (i in 1:n_sero){
    
    samp_effi_symp_n[i] <- sample_from_skew_normal(pt = symp_pt_sero_n[i], 
                                                   low_ci = symp_low_sero_n[i],
                                                   up_ci = symp_up_sero_n[i],
                                                   n_sample = 1)
    
    
    samp_effi_symp_p[i] <- sample_from_skew_normal(pt = symp_pt_sero_p[i], 
                                                   low_ci = symp_low_sero_p[i],
                                                   up_ci = symp_up_sero_p[i],
                                                   n_sample = 1)
    
    ## Missing estimates are replaced by corresponding VE against VCD
    if (i == 2 | i == 4) {
      
      samp_effi_hosp_n[i] <- samp_effi_symp_n[i]
      
    } else{ samp_effi_hosp_n[i] <- sample_from_skew_normal(pt = hosp_pt_sero_n[i], 
                                                           low_ci = hosp_low_sero_n[i],
                                                           up_ci = hosp_up_sero_n[i],
                                                           n_sample = 1)
    }
    
    
    ## Missing estimates are replaced by corresponding VE against VCD
    if (i == 4) {
      
      samp_effi_hosp_p[i] = samp_effi_symp_p[i]
      
    } else{ samp_effi_hosp_p[i] <- sample_from_skew_normal(pt = hosp_pt_sero_p[i], 
                                                           low_ci = hosp_low_sero_p[i],
                                                           up_ci = hosp_up_sero_p[i],
                                                           n_sample = 1)
    }
    
  }
  
  
  
  return(list(
    inf_n = samp_effi_inf_n,
    inf_p = samp_effi_inf_p,
    symp_n = samp_effi_symp_n,
    symp_p = samp_effi_symp_p,
    hosp_n = samp_effi_hosp_n,
    hosp_p = samp_effi_hosp_p
  ))
  
}




get_conditional_ve <- function(inf_n,
                               inf_p,
                               symp_n,
                               symp_p,
                               hosp_n,
                               hosp_p) {
  
  
  cond_inf_n = inf_n
  cond_inf_p = inf_p
  cond_symp_n = (symp_n - rep(inf_n,n_sero))/(1 - rep(inf_n, n_sero))
  cond_symp_p = (symp_p - rep(inf_p,n_sero))/(1 - rep(inf_p, n_sero))
  cond_hosp_n = (hosp_n - symp_n)/(1 - symp_n)
  cond_hosp_p = (hosp_p - symp_p)/(1 - symp_p)
  
  return(list(inf_n_cond = cond_inf_n,
              inf_p_cond = cond_inf_p,
              symp_n_cond = cond_symp_n,
              symp_p_cond = cond_symp_p,
              hosp_n_cond = cond_hosp_n,
              hosp_p_cond = cond_hosp_p ))
  
}




sample_from_skew_normal <- function(pt, low_ci, up_ci, n_sample) {
  # pt: Point estimate (mean)
  # low_ci: Lower bound of the confidence interval
  # up_ci: Upper bound of the confidence interval
  # n_sample: number of samples to draw
  
  # Confidence level (fixed at 95%)
  confidence_level <- 0.95
  z_value <- qnorm(1 - (1 - confidence_level) / 2)  # Z-score for 95% CI
  
  # Estimate standard deviation from CI
  sigma <- (up_ci - low_ci) / (2 * z_value)
  
  # Calculate asymmetry ratio and skewness parameter alpha
  ratio <- (up_ci - pt) / (pt - low_ci)
  k <- 1  # Scaling factor for alpha (adjust as needed)
  alpha <- k * (ratio - 1)
  
  # Initialize an empty vector for samples
  sample_skewnorm <- numeric(0)
  
  # Draw samples until the required number is reached
  while (length(sample_skewnorm) < n_sample) {
    # Draw a batch of samples
    new_samples <- rsn(n = n_sample - length(sample_skewnorm), xi = pt, omega = sigma, alpha = alpha)
    # Filter samples within the desired range
    valid_samples <- new_samples[new_samples >= low_ci & new_samples <= up_ci]
    # Append valid samples
    sample_skewnorm <- c(sample_skewnorm, valid_samples)
  }
  
  return(sample_skewnorm[1:n_sample])  # Ensure the output has exactly n_sample samples
}






get_final_stage_conditional_ve <- function(symp_low_sero_n,
                                           symp_pt_sero_n,
                                           symp_up_sero_n,
                                           symp_low_sero_p,
                                           symp_pt_sero_p,
                                           symp_up_sero_p,
                                           
                                           hosp_low_sero_n,
                                           hosp_pt_sero_n,
                                           hosp_up_sero_n,
                                           hosp_low_sero_p,
                                           hosp_pt_sero_p,
                                           hosp_up_sero_p,
                                           waning){
  
  
  
  ## sample the baseline VE estimates
  base_sample_ve_ci <- get_sample_from_ci(symp_low_sero_n,
                                          symp_pt_sero_n,
                                          symp_up_sero_n,
                                          symp_low_sero_p,
                                          symp_pt_sero_p,
                                          symp_up_sero_p,
                                          
                                          hosp_low_sero_n,
                                          hosp_pt_sero_n,
                                          hosp_up_sero_n,
                                          hosp_low_sero_p,
                                          hosp_pt_sero_p,
                                          hosp_up_sero_p)
  
  ## single them out
  base_samp_effi_inf_n = base_sample_ve_ci$inf_n
  base_samp_effi_inf_p = base_sample_ve_ci$inf_p
  base_samp_effi_symp_n =  base_sample_ve_ci$symp_n 
  base_samp_effi_symp_p =  base_sample_ve_ci$symp_p
  base_samp_effi_hosp_n =  base_sample_ve_ci$hosp_n 
  base_samp_effi_hosp_p =  base_sample_ve_ci$hosp_p 
  
  
  ## calculate the reduced rate of VE (not conditional)
  
  waning_inf_n_stage = array(NA, dim = c(n_vac_stage,1))
  waning_inf_p_stage = array(NA, dim = c(n_vac_stage,1))
  waning_symp_n_stage = array(NA, dim = c(n_vac_stage, n_sero))
  waning_symp_p_stage = array(NA, dim = c(n_vac_stage, n_sero))
  waning_hosp_n_stage = array(NA, dim = c(n_vac_stage, n_sero))
  waning_hosp_p_stage = array(NA, dim = c(n_vac_stage, n_sero))
  
  
  
  for (i in 1:n_vac_stage){
    
    
    ## for efficacy gainst infection, we dont consider any waning
    if (ve_inf_sero_n_pt == 0 & ve_inf_sero_n_low == 0 & ve_inf_sero_n_up == 0 ) {
      waning_inf_n_stage[i] = base_samp_effi_inf_n } else {
        waning_inf_n_stage[i] = base_samp_effi_inf_n - ( i - 1)*waning 
      }
    
    ## seropositive
    if (ve_inf_sero_p_pt == 0 & ve_inf_sero_p_low == 0 & ve_inf_sero_p_up == 0 ) {
      waning_inf_p_stage[i] = base_samp_effi_inf_p } else {
        waning_inf_p_stage[i] = base_samp_effi_inf_p - (i - 1)*waning 
      }
    
    
    waning_symp_n_stage[i,] = base_samp_effi_symp_n - (i - 1)*waning 
    waning_symp_p_stage[i,] = base_samp_effi_symp_p - (i - 1)*waning 
    waning_hosp_n_stage[i,] = base_samp_effi_hosp_n - (i - 1)*waning 
    waning_hosp_p_stage[i,] = base_samp_effi_hosp_p - (i - 1)*waning 
    
    
  }
  
  
  waning_stage_conditional_ve <- get_conditional_ve(inf_n = waning_inf_n_stage,
                                                    inf_p = waning_inf_p_stage,
                                                    symp_n = waning_symp_n_stage,
                                                    symp_p = waning_symp_p_stage,
                                                    hosp_n = waning_hosp_n_stage,
                                                    hosp_p = waning_hosp_p_stage)
  
  
  
  effi_inf_sero_n = waning_stage_conditional_ve$inf_n 
  effi_inf_sero_p =  waning_stage_conditional_ve$inf_p 
  
  effi_symp_sero_n = waning_stage_conditional_ve$symp_n 
  effi_symp_sero_p = waning_stage_conditional_ve$symp_p 
  
  effi_hosp_sero_n = waning_stage_conditional_ve$hosp_n 
  effi_hosp_sero_p = waning_stage_conditional_ve$hosp_p 
  
  return(list(cnd_inf_n = effi_inf_sero_n,
              cnd_inf_p = effi_inf_sero_p,
              cnd_symp_n = effi_symp_sero_n,
              cnd_symp_p = effi_symp_sero_p,
              cnd_hosp_n = effi_hosp_sero_n,
              cnd_hosp_p = effi_hosp_sero_p
  ))
  
}





### Steps to get the final efficacy parameters

n_year_waning <- 5

## 4.5 year with the reported efficcay then each year wanes at a rate "waning"
vac_stage_rate <- c( 1/(54*30), rep( (1/(12*30)),n_year_waning), 0)    #c( rep(1/(12*30),3), 0)

n_vac_stage <- length(vac_stage_rate)
n_sero <- 4
## length of waning rate should be less by 1 than of n_vac_stage
waning <- 0.20


##################

ve_inf_sero_n_low = 0
ve_inf_sero_n_up = 0
ve_inf_sero_n_pt = 0
ve_inf_sero_p_low = 0
ve_inf_sero_p_up = 0
ve_inf_sero_p_pt = 0


# Read the excel file
excel_efficacy = read_xlsx( here("data","qdenga_efficacy.xlsx"), 
                            sheet = "cumulative_54" )

# grab serotype stratified VE estimates
df_efficacy <- get_raw_ve_from_trial(excel_efficacy = excel_efficacy)

## now grab lower, point, and upper bound 
pt_lb_ub <- get_pt_lb_ub_from_excel_as_array(df_efficacy = df_efficacy)

symp_low_sero_n =  pt_lb_ub$symp_low_sero_n
symp_pt_sero_n =  pt_lb_ub$symp_pt_sero_n
symp_up_sero_n = pt_lb_ub$symp_up_sero_n
symp_low_sero_p = pt_lb_ub$symp_low_sero_p
symp_pt_sero_p = pt_lb_ub$symp_pt_sero_p
symp_up_sero_p = pt_lb_ub$symp_up_sero_p

hosp_low_sero_n = pt_lb_ub$hosp_low_sero_n
hosp_pt_sero_n = pt_lb_ub$hosp_pt_sero_n
hosp_up_sero_n = pt_lb_ub$hosp_up_sero_n
hosp_low_sero_p = pt_lb_ub$hosp_low_sero_p
hosp_pt_sero_p = pt_lb_ub$hosp_pt_sero_p
hosp_up_sero_p = pt_lb_ub$hosp_up_sero_p


## now pass these to get the final conditional VE




################################################################################



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
I_m0 =  c(.00005,.001,.1,1e-6) #c(.0001,.01,.1,1e-6)     #c(2e-2,5e-3,1e-2,1e-6) (previous values)  c(.00005,.001,.1,1e-6) 

### importation
import_switch <- c(1,0,0,0)
n_import <- 0.5      ## serotypewise: (0.5, 0.1, 5, 5)
import_start <- vac_start 
import_end <- import_start + 10*365

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



#############################  CHECK POINT 2 ################################

## sampling from posterior
n_sample <- 100

post <- readRDS(here("model_output", "posterior_rep_rate.rds"))


get_post_sample <- function() {
  
  
  sampled_row <- sample(1:nrow(post),1)
  
  rho_2_sample <- post[sampled_row,]
  
  
  rho_2 <- array(NA, dim=c(n_age))
  
  rho_2[1:5] <- rho_2_sample[1]
  rho_2[6:15] <- rho_2_sample[2]
  rho_2[16:25] <- rho_2_sample[3]
  rho_2[26:35] <- rho_2_sample[4]
  rho_2[36:45] <- rho_2_sample[5]
  rho_2[46:55] <- rho_2_sample[6]
  rho_2[56:65] <- rho_2_sample[7]
  rho_2[66:91] <- rho_2_sample[8]
  
  return(list(rho_2_samp = rho_2))
}



# compile the model to run in parallel
out_vaccine = function(par) {
  
  path_to_model <-  here::here("projection", "test_model_sero_dominance.R") 
  dengue_vaccine_model <- odin.dust::odin_dust(path_to_model) 
  
  #### instance of the model
  model <- dengue_vaccine_model$new(par, time = 1, n_particles = 1,
                                    ode_control = list(max_steps = 10000000,
                                                       step_size_min = 1e-14,
                                                       debug_record_step_times = TRUE))
  
  
  ## simulate the model
  out <- model$simulate(t)
  
  ## indices 
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
  
  tot_pop_age_yearly <- calc_year_age_wo_sum(array(out[ind_tot_age,,], dim = c(n_age, length(t))))     
  
  ## vaccinated population
  tot_pop_age_v_yearly <- calc_year_age_wo_sum(array(out[ind_N_age_v,,], dim = c(n_age, length(t)))) 
  
  ## unvaccinated population
  tot_pop_age_nv_yearly <- calc_year_age_wo_sum(array(out[ind_N_age_nv,,], dim = c(n_age, length(t))))  
  
  ## incidence of infection
  inf_age_yearly <- calc_year_age(array(out[ind_inf,,], dim = c(n_age, n_sero, length(t))))                             
  inf_sero_yearly <- calc_year_sero(array(out[ind_inf,,], dim = c(n_age, n_sero, length(t))))  
  
  inf_age_v_yearly <- calc_year_age(array(out[ind_inf_v,,], dim = c(n_age, n_sero, length(t))))                           
  inf_sero_v_yearly <- calc_year_sero(array(out[ind_inf_v,,], dim = c(n_age, n_sero, length(t))))  
  
  ## incidence of symptomatic cases
  symp_age_yearly <- calc_year_age(array(out[ind_symp,,], dim = c(n_age, n_sero,length(t))))
  symp_age_v_yearly <- calc_year_age( array(out[ind_symp_v,,], dim = c(n_age, n_sero,length(t))))
  
  ## incidence of hospitalized cases
  hosp_age_yearly <- calc_year_age(array(out[ind_hosp,,], dim = c(n_age, length(t))))
  hosp_age_v_yearly <- calc_year_age(array(out[ind_hosp_v,,], dim = c(n_age, length(t))))
  
  ## new vaccinated each tyear
  new_vac_age_yearly <- calc_year_age(array(out[ind_new_vac,,], dim = c(n_age, length(t))))
  
  rm(list = c("out"))
  
  return(list(tot_pop_age = tot_pop_age_yearly,
              tot_pop_age_v = tot_pop_age_v_yearly,
              tot_pop_age_nv = tot_pop_age_nv_yearly,
              inf_age = inf_age_yearly,
              inf_age_v = inf_age_v_yearly,
              inf_sero = inf_sero_yearly,
              inf_sero_v = inf_sero_v_yearly,
              symp_age = symp_age_yearly,
              symp_age_v = symp_age_v_yearly,
              hosp_age = hosp_age_yearly,
              hosp_age_v = hosp_age_v_yearly,
              new_vac_age = new_vac_age_yearly
              
  ))
}





vac_switch <- rep(0,n_age)
vac_switch[v_a[1,1]:v_a[1,2]] <- 1
vac_coverage <- 0.0

ratio_p_s <- 1/2

rho_2 <- get_post_sample()$rho_2_samp
rho_1 <- ratio_p_s*rho_2

v_eff <- get_final_stage_conditional_ve(symp_low_sero_n = symp_low_sero_n,
                                        symp_pt_sero_n = symp_pt_sero_n,
                                        symp_up_sero_n = symp_up_sero_n,
                                        symp_low_sero_p = symp_low_sero_p,
                                        symp_pt_sero_p = symp_pt_sero_p,
                                        symp_up_sero_p = symp_up_sero_p,
                                        
                                        hosp_low_sero_n = hosp_low_sero_n,
                                        hosp_pt_sero_n = hosp_pt_sero_n,
                                        hosp_up_sero_n = hosp_up_sero_n,
                                        hosp_low_sero_p = hosp_low_sero_p,
                                        hosp_pt_sero_p = hosp_pt_sero_p,
                                        hosp_up_sero_p = hosp_up_sero_p,
                                        waning = waning)



if (baseline == "yes") {
  
  effi_inf_sero_n = as.numeric(0*v_eff$cnd_inf_n)
  effi_inf_sero_p = as.numeric(0*v_eff$cnd_inf_p)
  effi_symp_sero_n = 0*v_eff$cnd_symp_n
  effi_symp_sero_p = 0*v_eff$cnd_symp_p
  effi_hos_sero_n = 0*v_eff$cnd_hosp_n
  effi_hos_sero_p = 0*v_eff$cnd_hosp_p
  
} else if (baseline == "no") {
  
  effi_inf_sero_n = as.numeric(v_eff$cnd_inf_n)
  effi_inf_sero_p = as.numeric(v_eff$cnd_inf_p)
  effi_symp_sero_n = v_eff$cnd_symp_n
  effi_symp_sero_p = v_eff$cnd_symp_p
  effi_hos_sero_n = v_eff$cnd_hosp_n
  effi_hos_sero_p = v_eff$cnd_hosp_p
  
}



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
  
  n_import = n_import,
  import_switch = import_switch,
  import_start = import_start,
  import_end = import_end,
  
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
tot_pop_age <- xx$tot_pop_age
tot_pop_age_v <- xx$tot_pop_age_v
tot_pop_age_nv <- xx$tot_pop_age_nv

inf_age <- xx$inf_age
inf_age_v <- xx$inf_age_v

inf_sero = xx$inf_sero 
inf_sero_v = xx$inf_sero_v 

symp_age <- xx$symp_age
symp_age_v <- xx$symp_age_v

hosp_age <- xx$hosp_age
hosp_age_v <- xx$hosp_age_v

new_vac_age <- xx$new_vac_age




## plot cases by serotypes
matplot(t(inf_sero[,]), type="l")


## Force of infection
tot_pop <- colSums(tot_pop_age)
tot_inf <- colSums(inf_age) + colSums(inf_age_v)
foi <- 100*tot_inf/tot_pop


plot(foi, type="l")+
  lines(rep(sum(foi)/length(foi),length(foi)), col="red", type="l")

## Serotype distribution over time
percentage_matrix <- apply(inf_sero, 2, function(x) x / sum(x) * 100)

# Convert to a dataframe in long format for ggplot2
years <- 1:10
serotypes <- paste0("Serotype", 1:4)
df <- as.data.frame(percentage_matrix[,10:19])
colnames(df) <- years
df$Serotype <- serotypes

# Reshape the data to long format
library(tidyr)
df_long <- pivot_longer(df, cols = -Serotype, names_to = "Year", values_to = "Percentage")
df_long$Year <- as.numeric(df_long$Year)  # Ensure 'Year' is numeric for plotting

# Plot using ggplot2
library(ggplot2)
ggplot(df_long, aes(x = Year, y = Percentage, fill = Serotype)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar chart
  scale_fill_brewer(palette = "Set3") +  # Choose a color palette
  labs(title = "Percentage of Cases by Serotype Each Year",
       x = "Year",
       y = "Percentage",
       fill = "Serotype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## aggegated contribution for 10 years

agg_contri_sero <- apply(inf_sero[,10:19], 1, function(x) sum(x))

perc_agg_contri_sero <- 100*agg_contri_sero/sum(agg_contri_sero)
barplot(perc_agg_contri_sero)
