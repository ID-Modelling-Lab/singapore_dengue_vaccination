
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


## Vaccination scenario: 
baseline <- "yes"  

## Vaccination coverage
v_coverage <- 0.8

## call all the functions related to demography
begin_year = "2013"
end_year = "2021"

n_year <- (as.numeric(end_year) - as.numeric(begin_year)) + 1

## number of extra year for which the model will be simulated after 2021
vac_program_year <- 25

sim_year <- n_year + vac_program_year

## create array of time in day (required for solving ode)
t <- seq(1,(sim_year*365), by = 1)

### introduction of vaccination 
v_year <- n_year + 1  # 15

vac_start <- v_year*365  ## year of vaccination start





v_a = matrix( c(7,17,
                18,31,
                32,41,
                42,51,
                52,61,
                62,71,
                72,81
                
), nrow = 7, ncol = 2, byrow = TRUE)


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


##### Demographic
## Number of age group
n_age <- length(pop_each_age_group)
n_sero = 4
## 
age_width <- 1
## Transition rate between age groups
age_rate <- c(rep(1/(365*age_width), (n_age - 1)),0)
## mosquito dearth rate
death_mos = 1/10
## Adult mosquito recruitment rate
lambda_mos = death_mos 

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

# ### Vaccinated population
# S0_v_all = array(0, dim = c(n_vac_stage, n_age))
# I_v0 = array(0, dim = c(n_vac_stage, n_age, n_sero))
# C_v0 = array(0, dim = c(n_vac_stage, n_age, n_sero))
# S_v0 = array(0, dim = c(n_vac_stage, n_age, n_sero))
# I_v_ij0 = array(0, dim = c(n_vac_stage, n_age, n_sero))
# R_v0 = array(0, dim = c(n_vac_stage, n_age))
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
I_m0 = c(.00005,.001,.1,1e-6) #c(2e-2,5e-3,1e-2,1e-6)   


#####
sero_negative_age <- 100*S0_all/pop_each_age_group

sero_positive_age <- 100*sero_p


## sero_negative
# sero_negative_age_6_16 <- 100*sum(S0_all[7:17])/sum(pop_each_age_group[7:17])

sero_negative_age_6_16 <- 100*sum(S0_all[7:17])/sum(pop_each_age_group)

# sero_negative_age_17_30 <- 100*sum(S0_all[18:31])/sum(pop_each_age_group[18:31])

sero_negative_age_17_30 <- 100*sum(S0_all[18:31])/sum(pop_each_age_group)

# sero_negative_age_61_70 <- 100*sum(S0_all[62:71])/sum(pop_each_age_group[62:71])

sero_negative_age_61_70 <- 100*sum(S0_all[62:71])/sum(pop_each_age_group)

# sero_negative_age_71_80 <- 100*sum(S0_all[72:81])/sum(pop_each_age_group[72:81])

sero_negative_age_71_80 <- 100*sum(S0_all[72:81])/sum(pop_each_age_group)


# sero_positive
# sero_positive_age_6_16 <- 100*(sum(rowSums(S0)[7:17]) + sum(R0[7:17]))/sum(pop_each_age_group[7:17])

sero_positive_age_6_16 <- 100*(sum(rowSums(S0)[7:17]) + sum(R0[7:17]))/sum(pop_each_age_group)


sero_positive_age_17_30 <- 100*(sum(rowSums(S0)[18:31]) + sum(R0[18:31]))/sum(pop_each_age_group[18:31])

sero_positive_age_61_70 <- 100*(sum(rowSums(S0)[62:71]) + sum(R0[62:71]))/sum(pop_each_age_group[62:71])

sero_positive_age_71_80 <- 100*(sum(rowSums(S0)[72:81]) + sum(R0[72:81]))/sum(pop_each_age_group[72:81])


## Primary
# primary_age_6_16 <- 100*(sum(rowSums(S0)[7:17]))/sum(pop_each_age_group[7:17])

primary_age_6_16 <- 100*(sum(rowSums(S0)[7:17]))/sum(pop_each_age_group)


# primary_age_17_30 <- 100*(sum(rowSums(S0)[18:31]))/sum(pop_each_age_group[18:31])

primary_age_17_30 <- 100*(sum(rowSums(S0)[18:31]))/sum(pop_each_age_group)

# primary_age_61_70 <- 100*(sum(rowSums(S0)[62:71]))/sum(pop_each_age_group[62:71])

primary_age_61_70 <- 100*(sum(rowSums(S0)[62:71]))/sum(pop_each_age_group)


# primary_age_71_80 <- 100*(sum(rowSums(S0)[72:81]))/sum(pop_each_age_group[72:81])

primary_age_71_80 <- 100*(sum(rowSums(S0)[72:81]))/sum(pop_each_age_group)


## secondary
secondary_age_6_16 <- 100*(sum(R0[7:17]))/sum(pop_each_age_group[7:17])

secondary_age_17_30 <- 100*(sum(R0[18:31]))/sum(pop_each_age_group[18:31])

secondary_age_61_70 <- 100*(sum(R0[62:71]))/sum(pop_each_age_group[62:71])

secondary_age_71_80 <- 100*(sum(R0[72:81]))/sum(pop_each_age_group[72:81])

## serotype wise

# serotypewise_age_6_16 <- 100*colSums(S0[7:17,])/sum(pop_each_age_group[7:17])
# 
# serotypewise_age_17_30 <- 100*colSums(S0[18:31,])/sum(pop_each_age_group[18:31])
# 
# serotypewise_age_62_71 <- 100*colSums(S0[62:71,])/sum(pop_each_age_group[62:71])
# 
# serotypewise_age_72_81 <- 100*colSums(S0[72:81,])/sum(pop_each_age_group[72:81])

serotypewise_age_6_16 <- 100*colSums(S0[7:17,])/sum(pop_each_age_group)

serotypewise_age_17_30 <- 100*colSums(S0[18:31,])/sum(pop_each_age_group)

serotypewise_age_62_71 <- 100*colSums(S0[62:71,])/sum(pop_each_age_group)

serotypewise_age_72_81 <- 100*colSums(S0[72:81,])/sum(pop_each_age_group)


## reporting arte
n_sample <- 100

post <- readRDS(here("model_output", "posterior_rep_rate.rds"))


get_post_sample <- function() {
  
  set.seed(12345)
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



rho_2 <- get_post_sample()$rho_2_samp


repo_rate_6_16 <- sum(rho_2[7:17])/11

repo_rate_17_30 <- sum(rho_2[18:31])/14

repo_rate_61_70 <- sum(rho_2[62:71])/10

repo_rate_71_80 <- sum(rho_2[72:81])/10
