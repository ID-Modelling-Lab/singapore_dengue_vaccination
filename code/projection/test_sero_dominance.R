
### call all the libraries
library(odin.dust)
library(readxl)
library(dplyr)
library(tidyr)
library(here)
library(foreach)
library(doParallel)
# library(doRNG)


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

###############################################################################
## Number of serotype
n_sero <- 4

### Assumption on hosp and vcd for denv3 and denv4

####################      CHECK POINT 1    ###########################


n_vac_stage <- 4 #length(unique(sero_strat_ve$month))


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

#############################  CHECK POINT 2 ################################

## sampling from posterior
n_sample <- 20

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
  inf_sero_yearly <- calc_year_sero(inf_age_sero)  ## serotype wise infection
  
  inf_age_sero_v <- array(out[ind_inf_v,,], dim = c(n_age, n_sero, length(t))) ## for vaccinated
  inf_age_v_yearly <- calc_year_age(inf_age_sero_v)                            ## yearly
  inf_sero_v_yearly <- calc_year_sero(inf_age_sero_v)  ## serotype wise infection
  
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
              inf_sero = inf_sero_yearly,
              inf_sero_v = inf_sero_v_yearly,
              symp_age = symp_age_yearly,
              symp_age_v = symp_age_v_yearly,
              hosp_age = hosp_age_yearly,
              hosp_age_v = hosp_age_v_yearly,
              new_vac_age = new_vac_age_yearly
  ))
}

##################  CHECK POINT 3 ###########################


### nested looping starts here



## Coverage values should not be zero, for baseine the efficacy should be zero
vac_coverage = 0.0

v_year <- 15
vac_start <- v_year*365  ## year of vaccination start

vac_switch <- rep(0,n_age)

ratio_p_s <- 1/2

set.seed(12345)      
rho_2 <- get_post_sample()$rho_2_samp
rho_1 <- ratio_p_s*rho_2


effi_inf_sero_n = array(0, dim = c(n_vac_stage))
effi_inf_sero_p = array(0, dim = c(n_vac_stage))
effi_symp_sero_n = array(0, dim = c(n_vac_stage, n_sero))
effi_symp_sero_p = array(0, dim = c(n_vac_stage, n_sero))
effi_hos_sero_n = array(0, dim = c(n_vac_stage, n_sero))
effi_hos_sero_p = array(0, dim = c(n_vac_stage, n_sero))


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
        I_v_ij0 = I_v_ij0,
        R_v0 = R_v0,
        
        effi_inf_sero_n = effi_inf_sero_n,
        effi_inf_sero_p = effi_inf_sero_p,
        effi_symp_sero_n = effi_symp_sero_n,
        effi_symp_sero_p = effi_symp_sero_p,
        effi_hos_sero_n = effi_hos_sero_n,
        effi_hos_sero_p = effi_hos_sero_p
)
      
xx <-  out_vaccine(par)
      
# 
#       
# inf_sample_sero[k,,] = xx$inf_sero #array(NA, dim = c(n_sample,n_sero, sim_year))
# inf_sample_sero_v[k,,] = xx$inf_sero_v # array(NA, dim = c(n_sample,n_sero, sim_year))
# 
#       



## TO add v_a in the output list also start of vaccination year

####  CHECK POINT:  Check the file names 

# 
# if (run_type == "interval") {
#   
#   saveRDS(list(coverage = v_coverage, 
#                n_sample = n_sample,
#                n_sero = n_sero,
#                n_year = sim_year,
#                n_vac_age = n_vac_age,
#                n_age = n_age,
#                v_year = v_year,
#                v_age_list = v_a,
#                lower_bound = lower_bound,
#                output=results_list), 
#           file = here::here("model_output", paste0("yng_old_run_", run_type, "_vcov_", v_coverage,
#                                                    "_nsamp_", n_sample,"_baseline_", baseline, "_lb_", 
#                                                    lower_bound,"_scen_dv34_zero.rds" )))
#   
# } else {
#   
#   saveRDS(list(coverage = v_coverage, 
#                n_sample = n_sample,
#                n_sero = n_sero,
#                n_year = sim_year,
#                n_vac_age = n_vac_age,
#                n_age = n_age,
#                v_year = v_year,
#                v_age_list = v_a,
#                output=results_list), 
#           file = here::here("model_output", paste0("yng_old_run_", run_type, "_vcov_", v_coverage,"_baseline_", baseline,
#                                                    "_nsamp_", n_sample,"_scen_dv34_zero.rds" )))
# }



