###########################################################################################################
###########################################################################################################
################# THIS IS THE MAIN SCRIPT FOR ALL THE VACCINE IMPACT SCENARIOS  ###########################

###### Input: 1.Vaccine Coverage (should be between 0 to 1 )
######        2.Specification if the run is for baseline or not

######        If baseline, then baseline = "yes" otherwise baseline = "no"
###### For a specific vaccine coverage it will run for different targeted age groups in parallel

###### To choose a dengue serotype dominance scenario: change the position of 1 in "import_switch"
###### accordingly. For example, for DENV-2 dominance scenario, it should be: import_switch = (0,1,0,0)

###### For scenario with pre-screening "prescreening" = 1 otherwise 0
###### The default model run assumes that the efficacy against infection is zero
##### For running efficacy against scenario the efficacy should be given 
#### For example we assume the efficacy against infection is same as Dengavaxia
#### In this case, the parameter related efficacy against infection should be given as below:
# ve_inf_sero_n_low = -0.359
# ve_inf_sero_n_up = 0.388
# ve_inf_sero_n_pt = 0.093
# ve_inf_sero_p_low = 0.352
# ve_inf_sero_p_up = 0.585
# ve_inf_sero_p_pt = 0.481

#### The default run assume no waning of efficacy (i.e., waning = 0.0) but can create one such scenario
#### by giving a rate ow waning through the parameter waning
###########################################################################################################
###########################################################################################################

### call all the libraries required
library(odin.dust)
library(readxl)
library(dplyr)
library(tidyr)
library(here)
library(foreach)
library(doParallel)
library(sn) ## for skew normal distribution

# remove everything from workspace
rm(list = ls(all = TRUE)) 

# collects garbage memory
gc(reset = TRUE)        

## cal all the functions from utils
source(here("utils","utils_demography.R"))
source(here("utils","utils_sero_prevalence.R"))
source(here("utils","utils_vac_efficacy.R"))
source(here("utils","utils_aggregated.R"))


####### (uncomment chunk below if want to run on HPC ot in commandline) #######
###### ### Accept command-line arguments 

args <- commandArgs(trailingOnly = TRUE)

# # Assign values from command-line arguments, if provided
if (length(args) >= 1) {
  baseline <- args[1]  # First argument: baseline
}
if (length(args) >= 2) {
  v_coverage <- as.numeric(args[2])  # Second argument: coverage
}

### Check for valid inputs
if (!baseline %in% c("yes", "no")) {
  stop("Invalid value for baseline. Must be 'yes' or 'no'.")
}
if (is.na(v_coverage) || v_coverage < 0 || v_coverage > 1) {
  stop("Invalid value for v_coverage. Must be a number between 0 and 1.")
}

cat("Running simulation with parameters:\n")
cat("baseline =", baseline, "\n")
cat("v_coverage =", v_coverage, "\n")
################################################################################
################################################################################

# Comment(Uncomment) below chuck if want to run on HPC(Locally) resp.)
################################################################################
## Vaccination scenario: 
# baseline <- "yes"
# 
# ## Vaccination coverage
# v_coverage <- 0.8
################################################################################
## this for gathering demographic data from singstat
begin_year = "2013"
end_year = "2021"

#total number of year we have the surveilance data for
n_year <- (as.numeric(end_year) - as.numeric(begin_year)) + 1

## number of extra year for which the model will be simulated after 2021 i.e.,
## duration of vaccination program
vac_program_year <- 15

# total simulation year
sim_year <- n_year + vac_program_year

## create array of time in day (required for solving ode)
t <- seq(1,(sim_year*365), by = 1)

### introduction of vaccination we introduce in 2022
v_year <- n_year + 1  

vac_start <- v_year*365  ## year of vaccination start

### importation
## serotypewise: (0.5, 0.1, 5, 5) by parameter exploration in the baseline model
################################################################################
# change this according for creating serotype dominance scenarios
import_switch <- c(1,0,0,0)
################################################################################

## importation switching
n_import <- if (import_switch[1] == 1) 0.5 else if 
(import_switch[2] == 1) 0.1 else if 
(import_switch[3] == 1) 5 else if 
(import_switch[4] == 1) 5

## scenario based on serotype dominance
sero_domin <- if (import_switch[1] == 1) "denv1" else if 
(import_switch[2] == 1) "denv2" else if 
(import_switch[3] == 1) "denv3" else if 
(import_switch[4] == 1) "denv4"

# period of importation (assumed during vaccination program)
import_start <- vac_start 
import_end <- import_start + vac_program_year*365

## IF Prescreening = 1 means vaccinating only seropositive 
## Prescreening == 0 means vaccinate all with any serostatus
prescreening = 0

## Sensitivity and specificity of the test
test_sensiti <- 0.90
test_specifi <- 0.90

# targeted age groups for vaccination
v_a = matrix( c(7,17,
                18,31,
                32,41,
                42,51,
                52,61,
                62,71,
                72,81
), nrow = 7, ncol = 2, byrow = TRUE)


### This os for alternative age groups for supplementary
# v_a = matrix( c(7,17,
#                 18,31,
#                 32,61,
#                 31,90,
#                 41,90,
#                 51,90,
#                 61,90
# ), nrow = 7, ncol = 2, byrow = TRUE)
# 
# 
# ## number vaccine coverage scenarios
# n_vac_coverage = length(v_coverage)
# 
# ## number of vaccinated age-groups
# n_vac_age = dim(v_a)[1]

################################################################################
# rate of waning per year from baseline values of VE (default model run it is 0)
waning <- 0.0

## number of year after which it will stop waning
n_year_waning <- 2

## 4.5 year with the reported efficacy then each year wanes at a rate "waning"
vac_stage_rate <- c( 1/(54*30), rep( (1/(12*30)), (n_year_waning - 1)), 0)  

# number of vaccination efficacy stage as vaccine will have different efficacy in different year
n_vac_stage <- length(vac_stage_rate)

# number of serotype
n_sero <- 4
################################################################################
################################################################################

################################################################################
## efficacy estimates for Dengvaxia (Uncomment the non-zero estimates for running 
### the scenario with efficacy of infection)
# ve_inf_sero_n_low = 0 #-0.359
# ve_inf_sero_n_up = 0 #0.388
# ve_inf_sero_n_pt = 0 #0.093
# ve_inf_sero_p_low = 0 #0.352
# ve_inf_sero_p_up = 0 #0.585
# ve_inf_sero_p_pt = 0 #0.481
################################################################################

## sensitivity analysis with lower bound 0 

ve_inf_sero_n_low = 0
ve_inf_sero_n_up = 0.388
ve_inf_sero_n_pt = 0.093
ve_inf_sero_p_low = 0.352
ve_inf_sero_p_up = 0.585
ve_inf_sero_p_pt = 0.481

##### Demographic
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
sero_year = 2013
get_sero <- get_sero_p(sero_year = sero_year )

## sero-positive fraction i.e at least one exposure
sero_p <- get_sero$se_p

## exactly one exposure
primary <- get_sero$pri

## Number of age group
n_age <- length(pop_each_age_group)

# total human population
total_human <- sum(pop_each_age_group)

## 1 year age band
age_width <- 1

## Transition rate between age groups
age_rate <- c(rep(1/(365*age_width), (n_age - 1)),0)


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
## Susceptible population
S0_all =  pop_each_age_group - rowSums(I0 + C0 + S0) - rowSums(I_ij0) - (R0)

### Vaccinated population
S0_v_all = array(0, dim = c(n_vac_stage, n_age))
I_v0 = array(0, dim = c(n_vac_stage, n_age, n_sero))
C_v0 = array(0, dim = c(n_vac_stage, n_age, n_sero))
S_v0 = array(0, dim = c(n_vac_stage, n_age, n_sero))
I_v_ij0 = array(0, dim = c(n_vac_stage, n_age, n_sero))
R_v0 = array(0, dim = c(n_vac_stage, n_age))

## ratio of mosquito human population
factor <- .01
## Susceptible mosquito
S_m0 = factor*total_human
# ## Exposed Mosquito
E_m0 = c(0,0,0,0) 
# ## infected mosquito
I_m0 = c(.00005,.001,.1,1e-6) 


###### Transmission related parameters
## biting rate per mosquito
b = 15
# Transmission rate from infected human to mosquito
beta_h = c(0.20, 0.20, 0.20, 0.20)   
# Transmission rate from infected mosquito to human 
beta_m = c(0.166, 0.176, 0.152, 0.144)
# Infectious period primary infection
gamma_1 = 1/4
# Infectious period secondary infection
gamma_2 = 1/4
# Cross immunity period
alpha = 1/(1*365)

# extrinsic incubation period
sigma_m = 1/10
## mosquito dearth rate
death_mos = 1/10
## Adult mosquito recruitment rate
lambda_mos = death_mos 


##### calculate  age-specific hospitalization data
df_hosp <- read_xlsx( here("data","serotype_share_age_cases.xlsx"), 
                      sheet = "age_hosp_frac_2003_17" )

# Function to calculate average for a given year range
calculate_average_hosp_rate <- function(data, start_year, end_year) {
  # Filter the data based on the year range
  filtered_data <- subset(data, year >= start_year & year <= end_year)
  
  # Calculate column-wise mean (excluding 'year' column)
  avg_values <- colMeans(filtered_data[ , -1], na.rm = TRUE)
  
  return(avg_values)
}

# Compute the average values and covert into fraction
## We use data from 2007 to 2017 as the previous data gives very high 
## hospitalization rate due to criteria getting hospitalization was different 
average_hosp_rate <- (1/100)*calculate_average_hosp_rate(data = df_hosp,
                                                         start_year = 2007, 
                                                         end_year =  2017)

ratio_p_s_hosp <- 1/4 ## assuming primary & secondary have same risk

# rate of hospitalization for secondary infection
xi_2 <- array(NA, dim = c(n_age))

xi_2[1:15] <- average_hosp_rate[1]
xi_2[16:25] <- average_hosp_rate[2]
xi_2[26:35] <- average_hosp_rate[3]
xi_2[36:45] <- average_hosp_rate[4]
xi_2[46:55] <- average_hosp_rate[5]
xi_2[56:65] <- average_hosp_rate[6]
xi_2[66:91] <- average_hosp_rate[7]

# rate of hospitalization for primary infection
xi_1 = ratio_p_s_hosp*xi_2

#############################  CHECK POINT 2 ################################

## sampling from posterior
n_sample <- 100

post <- readRDS(here("fitting/fitting_output", "posterior_rep_rate.rds"))

# post <- readRDS(here("fitting/fitting_output", "posterior_rep_rate_pri_eq_sec.rds"))


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
  
  path_to_model <-  here::here("projection", "model_vaccine.R") 
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
  ## yearly pop. in each group
  tot_pop_age_yearly <- calc_year_age_wo_sum(array(out[ind_tot_age,,], dim = c(n_age, length(t))))     
  
  ## vaccinated population
  tot_pop_age_v_yearly <- calc_year_age_wo_sum(array(out[ind_N_age_v,,], dim = c(n_age, length(t))))  ## yearly vacc. pop. in each group
  
  ## unvaccinated population
  tot_pop_age_nv_yearly <- calc_year_age_wo_sum(array(out[ind_N_age_nv,,], dim = c(n_age, length(t))))  ## yearly vacc. pop. in each group
  
  ## incidence of infection
  inf_age_yearly <- calc_year_age(array(out[ind_inf,,], dim = c(n_age, n_sero, length(t))))                             ## make this yearly
  inf_sero_yearly <- calc_year_sero(array(out[ind_inf,,], dim = c(n_age, n_sero, length(t))))  ## serotype wise infection
  
  inf_age_v_yearly <- calc_year_age(array(out[ind_inf_v,,], dim = c(n_age, n_sero, length(t))))                            ## yearly
  inf_sero_v_yearly <- calc_year_sero(array(out[ind_inf_v,,], dim = c(n_age, n_sero, length(t))))  ## serotype wise infection
  
  ## incidence of symptomatic cases
  symp_age_yearly <- calc_year_age(array(out[ind_symp,,], dim = c(n_age, n_sero,length(t))))
  symp_age_v_yearly <- calc_year_age( array(out[ind_symp_v,,], dim = c(n_age, n_sero,length(t))))
  
  ## incidence of hospitalized cases
  hosp_age_yearly <- calc_year_age(array(out[ind_hosp,,], dim = c(n_age, length(t))))
  hosp_age_v_yearly <- calc_year_age(array(out[ind_hosp_v,,], dim = c(n_age, length(t))))
  
  ## new vaccinated each tyear
  new_vac_age_yearly <- calc_year_age(array(out[ind_new_vac,,], dim = c(n_age, length(t))))

  ind_inf_primary <- model$info()$index$inc_inf_primary
  ind_inf_secondary <- model$info()$index$inc_inf_secondary
  
  ind_inf_v_primary <- model$info()$index$inc_inf_v_primary
  ind_inf_v_secondary <- model$info()$index$inc_inf_v_secondary
  
  ind_symp_primary <- model$info()$index$inc_symp_primary
  ind_symp_secondary <- model$info()$index$inc_symp_secondary
  
  ind_symp_v_primary <- model$info()$index$inc_symp_v_primary
  ind_symp_v_secondary <- model$info()$index$inc_symp_v_secondary
  
  ind_hosp_primary <- model$info()$index$inc_hosp_primary
  ind_hosp_secondary <- model$info()$index$inc_hosp_secondary
  
  ind_hosp_v_primary <- model$info()$index$inc_hosp_v_primary
  ind_hosp_v_secondary <- model$info()$index$inc_hosp_v_secondary
  
  ind_new_vac_primary <- model$info()$index$new_vac_age_primary
  ind_new_vac_secondary <- model$info()$index$new_vac_age_secondary
  
  inf_age_primary_yearly <- calc_year_age(array(out[ind_inf_primary,,], dim = c(n_age, n_sero, length(t))))
  inf_age_secondary_yearly <- calc_year_age(array(out[ind_inf_secondary,,], dim = c(n_age, n_sero, length(t))))
  
  inf_age_v_primary_yearly <- calc_year_age(array(out[ind_inf_v_primary,,], dim = c(n_age, n_sero, length(t))))
  inf_age_v_secondary_yearly <- calc_year_age(array(out[ind_inf_v_secondary,,], dim = c(n_age, n_sero, length(t))))
  
  
  symp_age_primary_yearly <- calc_year_age(array(out[ind_symp_primary,,], dim = c(n_age, n_sero,length(t))))
  symp_age_secondary_yearly <- calc_year_age(array(out[ind_symp_secondary,,], dim = c(n_age, n_sero,length(t))))
  
  
  symp_age_v_primary_yearly <- calc_year_age(array(out[ind_symp_v_primary,,], dim = c(n_age, n_sero,length(t))))
  symp_age_v_secondary_yearly <- calc_year_age(array(out[ind_symp_v_secondary,,], dim = c(n_age, n_sero,length(t))))
  
  hosp_age_primary_yearly <- calc_year_age(array(out[ind_hosp_primary,,], dim = c(n_age, length(t))))
  hosp_age_secondary_yearly <- calc_year_age( array(out[ind_hosp_secondary,,], dim = c(n_age, length(t))))
  
  
  
  hosp_age_v_primary_yearly <- calc_year_age(array(out[ind_hosp_v_primary,,], dim = c(n_age, length(t))))
  hosp_age_v_secondary_yearly <- calc_year_age(array(out[ind_hosp_v_secondary,,], dim = c(n_age, length(t))))
  
  new_vac_age_primary_yearly <- calc_year_age(array(out[ind_new_vac_primary,,], dim = c(n_age, length(t))))
  new_vac_age_secondary_yearly <- calc_year_age( array(out[ind_new_vac_secondary,,], dim = c(n_age, length(t))))
  
  
  
  
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
              new_vac_age = new_vac_age_yearly,
              
              ## primary secondary
              inf_age_primary = inf_age_primary_yearly,
              inf_age_v_primary = inf_age_v_primary_yearly,
              symp_age_primary = symp_age_primary_yearly,
              symp_age_v_primary = symp_age_v_primary_yearly,
              hosp_age_primary = hosp_age_primary_yearly,
              hosp_age_v_primary = hosp_age_v_primary_yearly,
              
              inf_age_secondary = inf_age_secondary_yearly,
              inf_age_v_secondary = inf_age_v_secondary_yearly,
              symp_age_secondary = symp_age_secondary_yearly,
              symp_age_v_secondary = symp_age_v_secondary_yearly,
              hosp_age_secondary = hosp_age_secondary_yearly,
              hosp_age_v_secondary = hosp_age_v_secondary_yearly,
              
              new_vac_age_primary = new_vac_age_primary_yearly,
              new_vac_age_secondary = new_vac_age_secondary_yearly
              # 
              
              
              
  ))
}

##################  CHECK POINT 3 ###########################


# Define the number of cores to use
numcores <- 7 # detectCores() - 1  
cl <- makeCluster(numcores)  # Adjust the number of cores based on your HPC setup
registerDoParallel(cl)


## Create space for output quantities
inf_sample_age = array(NA, dim = c(n_sample,n_age, sim_year))
inf_sample_age_v = array(NA, dim = c(n_sample,n_age, sim_year))
inf_sample_sero = array(NA, dim = c(n_sample,n_sero, sim_year))
inf_sample_sero_v = array(NA, dim = c(n_sample,n_sero, sim_year))
symp_sample_age = array(NA, dim = c(n_sample,n_age, sim_year))
symp_sample_age_v = array(NA, dim = c(n_sample,n_age, sim_year))
hosp_sample_age =  array(NA, dim = c(n_sample, n_age,sim_year))
hosp_sample_age_v =  array(NA, dim = c(n_sample,n_age, sim_year))
tot_pop_sample_age = array(NA, dim = c(n_sample, n_age, sim_year))
tot_pop_sample_age_v = array(NA, dim = c(n_sample, n_age, sim_year))
tot_pop_sample_age_nv = array(NA, dim = c(n_sample,  n_age, sim_year))
new_vac_sample_age =  array(NA, dim = c(n_sample, n_age,sim_year))

## Primary secondary stratifiaction

inf_sample_age_primary = array(NA, dim = c(n_sample,n_age, sim_year))
inf_sample_age_secondary = array(NA, dim = c(n_sample,n_age, sim_year))

inf_sample_age_v_primary = array(NA, dim = c(n_sample,n_age, sim_year))
inf_sample_age_v_secondary = array(NA, dim = c(n_sample,n_age, sim_year))

symp_sample_age_primary = array(NA, dim = c(n_sample,n_age, sim_year))
symp_sample_age_secondary = array(NA, dim = c(n_sample,n_age, sim_year))

symp_sample_age_v_primary = array(NA, dim = c(n_sample,n_age, sim_year))
symp_sample_age_v_secondary = array(NA, dim = c(n_sample,n_age, sim_year))

hosp_sample_age_primary =  array(NA, dim = c(n_sample, n_age,sim_year))
hosp_sample_age_secondary =  array(NA, dim = c(n_sample, n_age,sim_year))

hosp_sample_age_v_primary =  array(NA, dim = c(n_sample,n_age, sim_year))
hosp_sample_age_v_secondary =  array(NA, dim = c(n_sample,n_age, sim_year))

new_vac_sample_age_primary =  array(NA, dim = c(n_sample, n_age,sim_year))
new_vac_sample_age_secondary =  array(NA, dim = c(n_sample, n_age,sim_year))

### allocate space for output 
results_list <- vector("list", length(v_coverage))

## Run for different vaccination coverage and for different targeted age-groups
## Looping over coverage (i), vaccinated age groups (j), sample(k)
### Run in parallel
### make vaccinated age group run in parallel
for (i in 1:length(v_coverage)) {
  
  vac_coverage <- v_coverage[i]
  # print(i)
  
  out_parallel <- foreach(j =  1:n_vac_age, .combine = rbind, .packages = c("dplyr", "tidyr", "here", "odin.dust", "sn")) %dopar% {
    
    set.seed(12345)
    
    library(dplyr)
    library(tidyr)
    library(here)
    library(odin.dust)
    library(sn)
    
    
    vac_switch <- rep(0,n_age)
    vac_switch[v_a[j,1]:v_a[j,2]] <- 1
    
    for (k in 1:n_sample) {
      
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
        
        S0_v_all = S0_v_all,
        I_v0 = I_v0,
        C_v0 = C_v0,
        S_v0 = S_v0,
        I_v_ij0 =  I_v_ij0,
        R_v0 = R_v0,
        
        n_import = n_import,
        import_switch = import_switch,
        import_start = import_start,
        import_end = import_end,
        
        test_sensiti = test_sensiti,
        test_specifi = test_specifi,
        prescreening = prescreening,
   
        
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
      
      inf_sample_sero[k,,] = xx$inf_sero 
      inf_sample_sero_v[k,,] = xx$inf_sero_v
      
      symp_sample_age[k,,] <- xx$symp_age
      symp_sample_age_v[k,,] <- xx$symp_age_v
      
      hosp_sample_age[k,,] <- xx$hosp_age
      hosp_sample_age_v[k,,] <- xx$hosp_age_v
      
      new_vac_sample_age[k,,] <- xx$new_vac_age
      
      
      ## primary secondary stratification
      inf_sample_age_primary[k,,] <- xx$inf_age_primary
      inf_sample_age_secondary[k,,] <- xx$inf_age_secondary
      
      inf_sample_age_v_primary[k,,] <- xx$inf_age_v_primary
      inf_sample_age_v_secondary[k,,] <- xx$inf_age_v_secondary
      
      symp_sample_age_primary[k,,] <- xx$symp_age_primary
      symp_sample_age_secondary[k,,] <- xx$symp_age_secondary
      
      symp_sample_age_v_primary[k,,] <- xx$symp_age_v_primary
      symp_sample_age_v_secondary[k,,] <- xx$symp_age_v_secondary
      
      hosp_sample_age_primary[k,,] <- xx$hosp_age_primary
      hosp_sample_age_secondary[k,,] <- xx$hosp_age_secondary
      
      hosp_sample_age_v_primary[k,,] <- xx$hosp_age_v_primary
      hosp_sample_age_v_secondary[k,,] <- xx$hosp_age_v_secondary
      
      new_vac_sample_age_primary[k,,] <- xx$new_vac_age_primary
      new_vac_sample_age_secondary[k,,] <- xx$new_vac_age_secondary
      
      rm(list = c("xx"))
    }
    return(list(tot_pop_sample_age = tot_pop_sample_age,
                tot_pop_sample_age_v = tot_pop_sample_age_v,
                tot_pop_sample_age_nv = tot_pop_sample_age_nv,
                
                inf_sample_age = inf_sample_age,
                inf_sample_age_v = inf_sample_age_v,
                
                inf_sample_sero = inf_sample_sero,
                inf_sample_sero_v = inf_sample_sero_v,
                
                symp_sample_age = symp_sample_age,
                symp_sample_age_v = symp_sample_age_v,
                
                hosp_sample_age = hosp_sample_age,
                hosp_sample_age_v = hosp_sample_age_v,
                
                new_vac_sample_age = new_vac_sample_age,
                
                
                ## primary secondary stratification
                
                new_vac_sample_age_primary = new_vac_sample_age_primary,
                new_vac_sample_age_secondary = new_vac_sample_age_secondary,
                
                inf_sample_age_primary = inf_sample_age_primary,
                inf_sample_age_secondary = inf_sample_age_secondary,
                inf_sample_age_v_primary = inf_sample_age_v_primary,
                inf_sample_age_v_secondary = inf_sample_age_v_secondary,
                
                symp_sample_age_primary = symp_sample_age_primary,
                symp_sample_age_secondary = symp_sample_age_secondary,
                symp_sample_age_v_primary = symp_sample_age_v_primary,
                symp_sample_age_v_secondary = symp_sample_age_v_secondary,
                
                hosp_sample_age_primary = hosp_sample_age_primary,
                hosp_sample_age_secondary = hosp_sample_age_secondary,
                hosp_sample_age_v_primary = hosp_sample_age_v_primary,
                hosp_sample_age_v_secondary = hosp_sample_age_v_secondary
                
    ))
  }
  
  
  results_list[[i]] <- out_parallel
  
  rm(out_parallel)
  gc()
}



saveRDS(list(coverage = v_coverage, 
             n_sample = n_sample,
             n_sero = n_sero,
             n_year = sim_year,
             n_vac_age = n_vac_age,
             n_age = n_age,
             v_year = v_year,
             v_age_list = v_a,
             output = results_list), 
        file = here::here("model_output", paste0("effi_inf_lower_0_pri_sec_neq_detailed_a7_sero_dom_", sero_domin, "_vcov_", v_coverage,"_baseline_", baseline,
                                                 "_nsamp_", n_sample,"_waning_", waning,"_prescrng_", prescreening, ".rds" )))



stopCluster(cl)

closeAllConnections()


