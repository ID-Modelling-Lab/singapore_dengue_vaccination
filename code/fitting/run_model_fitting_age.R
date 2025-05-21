##########  Run this script for performing model fitting with surveillance data


### call all the libraries
library(dust)
library(mcstate)
library(odin.dust)
library(readxl)
library(ggplot2)
library(here)
# library(ggpubr)
library(dplyr)
# library(imputeTS)
library(latex2exp)
library(ggpubr)

### will remove everything from workspace
rm(list=ls(all=TRUE)) 
### collects garbage memory
gc(reset=TRUE)        

## cal all the functions from utils
source(here("utils","utils_demography.R"))
source(here("utils","utils_sero_prevalence.R"))


### Demographic Information
#### Population in each age group
## Specify the beginning and and end year for demographic information from the data (singstat.gov.sg)
begin_year = "2013"
end_year = "2021"

n_year <- (as.numeric(end_year) - as.numeric(begin_year))+1

# number of extra year for which the model will be simulated
extra_year <- 0

sim_year <- n_year + extra_year

## read the reported incidence data (yearly)
data_reported_age_yearly = read_xlsx(here("data","serotype_share_age_cases.xlsx" ), sheet="age_case_2014_20")

# select the incidence per 100k data
df_age_rep_inc <- select(data_reported_age_yearly, year,contains("_inc_") )

reported_cases = df_age_rep_inc %>%
  mutate(
    day = c(365*(1:(length(df_age_rep_inc$year))))
  )

## Number of serotype
n_sero <- 4


#### Crude birth rate over time
## population in each age group in "begin_year"
pop_each_age_group <- get_pop_each_group(begin_year = begin_year,
                                         end_year = end_year )
## Number of age group
n_age <- length(pop_each_age_group) 

## birth rate over time
lambda_hum_tt = get_birth_rate(begin_year = begin_year, 
                               end_year = end_year, 
                               sim_year = sim_year)

## death rate over time
death_hum_tt = get_death_rate(begin_year = begin_year, 
                              end_year = end_year, 
                              sim_year = sim_year)


# for sero-prevalence data
sero_year = 2013
## sero-prevalence data
get_sero <- get_sero_p(sero_year = sero_year )

## sero-positive fraction i.e at least one exposure
sero_p <- get_sero$se_p
## exactly one exposure
primary <- get_sero$pri


##### Demographic

age_width <- 1 ## 
## Transition rate between age groups
age_rate <- c(rep(1/(365*age_width), (n_age - 1)),0)

#### Mosquito
## mosquito dearth rate
death_mos = 1/10

## Adult mosquito recruitment rate
lambda_mos = death_mos 


###### Initial conditions

data_year <- reported_cases$day

## upper and lower bound bound of "years" in the data
l_t <- c(1,data_year[-length(data_year)])   
u_t <- c(data_year)   

## upper and lower bounds of age
a_lo <-   c(1, 6, 16, 26, 36, 46, 56, 66) 
a_up <-  c(5,15,25, 35, 45, 55, 65, 91)
n_data_year <- length(data_year)
n_data_age <- length(a_lo)

I0 = matrix(0, nrow = n_age, ncol = n_sero, byrow = TRUE)

frac_serotype <- c(.352860, .490160, .111505,  .045475)    

C0 = 0.0*t(frac_serotype*t(array(rep(primary*pop_each_age_group,n_sero), dim = c(n_age,n_sero))))

S0 = 1*t(frac_serotype*t(array(rep(primary*pop_each_age_group,n_sero), dim = c(n_age,n_sero))))

## Symptomatic (Secondary infection)
I_ij0 = array(0,dim = c(n_age, n_sero))
## Recovered from Secondary infection
R0 =  sero_p*pop_each_age_group - rowSums(S0) - rowSums(C0) # secondary*pop_each_age_group 

## Susceptible population
S0_all =  pop_each_age_group - rowSums(I0 + C0 + S0) - rowSums(I_ij0) - R0

total_pop_age0 <- pop_each_age_group  


inc_symp0 <- array(0, dim = c(n_data_year, n_age))
agg_inc_symp0 <- array(0,n_age)

factor <- .01
total_human <- sum(pop_each_age_group)
## Susceptible mosquito
S_m0 = factor*total_human
# ## Exposed Mosquito
E_m0 = c(0,0,0,0)
# ## infected mosquito
I_m0 = c(.00005,.001,.1,1e-6)  

## biting rate per mosquito
b = 15
# Transmission rate from infected human to mosquito
beta_h = c(0.20,0.20,0.20,0.20)   

# Transmission rate from mosquitoes infected (calibrated from FOI estimates)
beta_m1 = 0.166 
beta_m2 = 0.176 
beta_m3 = 0.152 
beta_m4 = 0.144 
beta_m = c(beta_m1, beta_m2, beta_m3, beta_m4)

ratio_p_s = 1/2 
# Infectious period primary infection
gamma_1 = 1/4
# Infectious period secondary infection
gamma_2 = 1/4
alpha = 1/(1*365)
sigma_m = 1/10


# compile the model
path_to_model <-  here::here("fitting", "model_fitting_age.R")
dengue_model <- odin.dust::odin_dust(path_to_model)


# ## for age specific incidence data
dengue_incidence_case_compare <- function(state, observed, pars = NULL) {
  
  
  ## incidence per 100k
  ll_1 = dpois(x = floor(observed$`0-4y_inc_p_100k`), lambda = state["inc1",,drop = TRUE],log = TRUE)
  
  ll_2 = dpois(x = floor(observed$`5-14y_inc_p_100k`), lambda = state["inc2",,drop = TRUE],log = TRUE)
  
  ll_3 = dpois(x = floor(observed$`15-24y_inc_p_100k`), lambda = state["inc3",,drop = TRUE],log = TRUE)
  
  ll_4 = dpois(x = floor(observed$`25-34y_inc_p_100k`), lambda = state["inc4",,drop = TRUE],log = TRUE)
  
  ll_5 = dpois(x = floor(observed$`35-44y_inc_p_100k`), lambda = state["inc5",,drop = TRUE],log = TRUE)
  
  ll_6 = dpois(x = floor(observed$`45-54y_inc_p_100k`), lambda = state["inc6",,drop = TRUE],log = TRUE)
  
  ll_7 = dpois(x = floor(observed$`55-64y_inc_p_100k`), lambda = state["inc7",,drop = TRUE],log = TRUE)
  
  ll_8 = dpois(x = floor(observed$`65+y_inc_p_100k`), lambda = state["inc8",,drop = TRUE],log = TRUE)

  
  ll = ll_1 + ll_2 + ll_3 + ll_4 + ll_5 + ll_6 + ll_7 + ll_8
  
  return(ll)
}


dengue_incidence_data <- mcstate::particle_filter_data(data = reported_cases, time = "day",
                                                       initial_time = 1, rate = NULL)

# plot the data # TODO

index <- function(info) {
  idx <- unlist(info$index)
  list(run = idx, state = idx)
}


# now do the filtering
filter <- mcstate::particle_deterministic$new(data = dengue_incidence_data, 
                                              model = dengue_model, 
                                              compare = dengue_incidence_case_compare, 
                                              index=index)

priors=list(
  
  mcstate::pmcmc_parameter("rho_2_a1", 0.1, min = 0, max=1,prior = function(a){
    dunif(a,min = 0, max = 1, log = TRUE)
  }),
  
  mcstate::pmcmc_parameter("rho_2_a2", 0.1, min = 0, max=1,prior = function(a){
    dunif(a,min = 0, max = 1, log=TRUE)
  }),
  
  mcstate::pmcmc_parameter("rho_2_a3", 0.1, min = 0, max=1,prior = function(a){
    dunif(a,min=0, max=1, log=TRUE)
  }),
  
  mcstate::pmcmc_parameter("rho_2_a4", 0.1, min = 0, max=1,prior = function(a){
    dunif(a,min=0, max= 1, log=TRUE)
  }),
  
  mcstate::pmcmc_parameter("rho_2_a5", 0.1, min = 0, max=1,prior = function(a){
    dunif(a,min=0, max= 1, log=TRUE)
  }),
  
  mcstate::pmcmc_parameter("rho_2_a6", 0.1, min = 0, max=1,prior = function(a){
    dunif(a,min=0, max= 1, log=TRUE)
  }),
  
  mcstate::pmcmc_parameter("rho_2_a7", 0.1, min = 0, max=1,prior = function(a){
    dunif(a,min=0, max= 1, log=TRUE)
  }),
  
  mcstate::pmcmc_parameter("rho_2_a8", 0.1, min = 0, max=1,prior = function(a){
    dunif(a,min=0, max= 1, log=TRUE)
  })
  
  
  
)



make_transform <- function(n_age,
                           n_sero, 
                           beta_h, 
                           beta_m1, 
                           beta_m2,
                           beta_m3, 
                           beta_m4, 
                           b, 
                           l_t,
                           u_t, 
                           a_lo, 
                           a_up,
                           sigma_m, 
                           alpha,                
                           gamma_1, 
                           gamma_2,              
                           lambda_hum_tt, 
                           death_hum_tt,
                           n_data_year,
                           n_data_age,
                           ratio_p_s,
                           lambda_mos, 
                           death_mos, 
                           age_rate, 
                           I0,
                           C0, 
                           S0, 
                           I_ij0, 
                           R0,
                           S0_all, 
                           S_m0, 
                           E_m0, 
                           I_m0, 
                           total_pop_age0,
                           inc_symp0, 
                           agg_inc_symp0) {
  list(
    n_age = n_age,
    n_sero = n_sero,
    beta_h = beta_h,
    beta_m1 = beta_m1,
    beta_m2 = beta_m2,
    beta_m3 = beta_m3,
    beta_m4 = beta_m4,
    b = b,
    l_t = l_t,
    u_t = u_t,
    a_lo = a_lo,
    a_up = a_up,
    n_data_year = n_data_year,
    n_data_age = n_data_age,
    sigma_m = sigma_m,              
    alpha = alpha,                
    gamma_1 = gamma_1,              
    gamma_2 = gamma_2,              
    ratio_p_s = ratio_p_s,
    
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
    total_pop_age0 = total_pop_age0,
    inc_symp0 = inc_symp0,
    agg_inc_symp0 = agg_inc_symp0
  )
  
  function(theta) {
    list(
      n_age = n_age,
      n_sero = n_sero,
      beta_h = beta_h,
      beta_m1 = beta_m1,
      beta_m2 = beta_m2,
      beta_m3 = beta_m3,
      beta_m4 = beta_m4,
      b = b,
      l_t = l_t,
      u_t = u_t,
      a_lo = a_lo,
      a_up = a_up,
      n_data_year = n_data_year,
      n_data_age = n_data_age,
      rho_2_a1 = theta[["rho_2_a1"]],
      rho_2_a2 = theta[["rho_2_a2"]],
      rho_2_a3 = theta[["rho_2_a3"]],
      rho_2_a4 = theta[["rho_2_a4"]],
      rho_2_a5 = theta[["rho_2_a5"]],
      rho_2_a6 = theta[["rho_2_a6"]],
      rho_2_a7 = theta[["rho_2_a7"]],
      rho_2_a8 = theta[["rho_2_a8"]],
      
      sigma_m = sigma_m,              
      alpha = alpha,                
      gamma_1 = gamma_1,              
      gamma_2 = gamma_2,              
      ratio_p_s = ratio_p_s,
      
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
      total_pop_age0 = total_pop_age0,
      inc_symp0 = inc_symp0,
      agg_inc_symp0 = agg_inc_symp0
    )
  }
}



transform <- make_transform(
  n_age = n_age,
  n_sero = n_sero,
  beta_h = beta_h,
  beta_m1 = beta_m1,
  beta_m2 = beta_m2,
  beta_m3 = beta_m3,
  beta_m4 = beta_m4,
  b = b,
  l_t = l_t,
  u_t = u_t,
  a_lo = a_lo,
  a_up = a_up,
  n_data_year = n_data_year,
  n_data_age = n_data_age,
  
  
  sigma_m = sigma_m,              
  alpha = alpha,                
  
  gamma_1 = gamma_1,              
  gamma_2 = gamma_2,              
  ratio_p_s = ratio_p_s,
  
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
  total_pop_age0 = total_pop_age0,
  inc_symp0 = inc_symp0,
  agg_inc_symp0 = agg_inc_symp0
) 


n_param = 8
proposal_matrix <- diag(1e-5, n_param)

mcmc_pars <- mcstate::pmcmc_parameters$new(priors,
                                           proposal = proposal_matrix,
                                           transform = transform
)

control <- mcstate::pmcmc_control(n_steps = 1e4,
                                  n_workers = 8,
                                  n_threads_total = 8,
                                  n_burnin = 500,
                                  n_chains = 8,
                                  adaptive_proposal = TRUE,
                                  save_state = TRUE,
                                  save_trajectories = TRUE,
                                  progress = TRUE )

# run pmcmc
mcmc_run <- mcstate::pmcmc(mcmc_pars, filter, control = control)

# Plot MCMC output
mcmc_out <- coda::as.mcmc(cbind(mcmc_run$probabilities, mcmc_run$pars))
summary(mcmc_out)
# plot(mcmc_out)

# Calculate effective sample size
coda::effectiveSize(mcmc_out)
coda::rejectionRate(mcmc_out)

# # Pairwise correlation plot
# par(mfrow = c(1,1))
# pairs(mcmc_run$pars[])



post_par <- mcmc_run$pars


## Tuning PMCMC (run again)

proposal_matrix1 <- cov(mcmc_run$pars)

mcmc_pars1 <- mcstate::pmcmc_parameters$new(priors,
                                           proposal_matrix1,
                                           transform)

control1 <- mcstate::pmcmc_control(n_steps = 1e4,
                                    n_workers = 8,
                                    n_threads_total = 8,
                                    n_burnin = 500,
                                    n_chains = 8,
                                    adaptive_proposal = TRUE,
                                    save_state = TRUE,
                                    save_trajectories = TRUE,
                                    progress = TRUE )


  # run pmcmc
  mcmc_run <- mcstate::pmcmc(mcmc_pars1, filter, control = control1)

  # Plot MCMC output
  mcmc_out1 <- coda::as.mcmc(cbind(mcmc_run$probabilities, mcmc_run$pars))
  summary(mcmc_out1)
  
  # Calculate effective sample size
  coda::effectiveSize(mcmc_out1)
  coda::rejectionRate(mcmc_out1)


  
  
# Uncomment below for saving the data  
    
## 1. save the posterior chain 

 #  post_par <- mcmc_run$pars
 # 
 #  
 #  
 # saveRDS(post_par,
 #        file = here::here("fitting/fitting_output", "posterior_rep_rate.rds"))


  
  
# 2. save the posterior along with all the trajectories (states) by sampling 
  # mcmc_sample <- mcstate::pmcmc_sample(mcmc_run, n_sample = 200)
  # 
  # # # extract the trajectories for prediction
  # state <- mcmc_sample$trajectories$state
  # 
  # 
  # saveRDS(state,
  #         file = here::here("fitting/fitting_output", "state.rds"))

  





