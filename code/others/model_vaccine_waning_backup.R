# total population
N <- sum(S_all[]) + sum(I[, ]) + sum(C[, ]) + sum(S[, ]) + sum(I_ij[, ])  + sum(R[]) +
  sum(S_v_all[, ]) + sum(I_v[, , ]) + sum(C_v[, , ]) + sum(S_v[, , ]) + sum(I_v_ij[, , ]) + sum(R_v[, ])

# total population: age-stratified
N_age[] <- S_all[i] + sum(I[i, ]) + sum(C[i, ]) + sum(S[i, ]) + sum(I_ij[i, ]) + R[i] +
  sum(S_v_all[, i]) + sum(I_v[, i, ]) + sum(C_v[, i, ]) + sum(S_v[, i, ]) + sum(I_v_ij[, i, ]) + sum(R_v[, i])

# infected human: serotype stratified
infected_hum[] <- sum(I[, i]) + sum(I_ij[, i]) + sum(I_v[, , i]) + sum(I_v_ij[, , i])

# total vaccinated: age-stratified
N_v_age[] <- sum(S_v_all[, i]) + sum(I_v[, i, ])  + sum(C_v[, i, ]) + sum(S_v[, i, ]) + sum(I_v_ij[, i, ]) + sum(R_v[, i]) 
# total non-vaccinated: age-stratified
N_nv_age[] <-  S_all[i] + sum(I[i,])  + sum(C[i,]) + sum(S[i,]) + sum(I_ij[i,]) + R[i]

## eligible for vaccination: age-stratified
elig_vac_age[] <- S_all[i] + sum(C[i,]) + sum(S[i,]) + R[i]


## calculate the number of vaccinated for different compartments and introduce vaccination timing
vac_Sall[1] <- if (as.integer(t) == vac_start) vac_switch[i]*vac_coverage*N_age[i]*S_all[i]/elig_vac_age[i] else
  vac_switch[i]*(lambda_hum*N*N_v_age[i]/N_age[i])*(S_all[i]/elig_vac_age[i])

vac_Sall[2:n_age] <- if (as.integer(t) == vac_start) vac_switch[i]*vac_coverage*N_age[i]*S_all[i]/elig_vac_age[i] else
  vac_switch[i]*(age_rate[i - 1]*(N_age[i - 1]*N_v_age[i] - N_v_age[i - 1]*N_age[i])/N_age[i])*(S_all[i]/elig_vac_age[i])

vac_C[1,1:n_sero] <- if (as.integer(t) == vac_start) vac_switch[i]*vac_coverage*N_age[i]*C[i,j]/elig_vac_age[i] else
  vac_switch[i]*(lambda_hum*N*N_v_age[i]/N_age[i])*(C[i,j]/elig_vac_age[i])

vac_C[2:n_age,1:n_sero] <- if (as.integer(t) == vac_start) vac_switch[i]*vac_coverage*N_age[i]*C[i,j]/elig_vac_age[i] else
  vac_switch[i]*(age_rate[i - 1]*(N_age[i - 1]*N_v_age[i] - N_v_age[i - 1]*N_age[i])/N_age[i])*(C[i,j]/elig_vac_age[i])

vac_S[1,1:n_sero] <- if (as.integer(t) == vac_start) vac_switch[i]*vac_coverage*N_age[i]*S[i,j]/elig_vac_age[i] else
  vac_switch[i]*(lambda_hum*N*N_v_age[i]/N_age[i])*(S[i,j]/elig_vac_age[i])

vac_S[2:n_age,1:n_sero] <- if (as.integer(t) == vac_start) vac_switch[i]*vac_coverage*N_age[i]*S[i,j]/elig_vac_age[i] else
  vac_switch[i]*(age_rate[i - 1]*(N_age[i - 1]*N_v_age[i] - N_v_age[i - 1]*N_age[i])/N_age[i])*(S[i,j]/elig_vac_age[i])

vac_R[1] <-  if (as.integer(t) == vac_start) vac_switch[i]*vac_coverage*N_age[i]*R[i]/elig_vac_age[i] else
  vac_switch[i]*(lambda_hum*N*N_v_age[i]/N_age[i])*(R[i]/elig_vac_age[i])

vac_R[2:n_age] <- if (as.integer(t) == vac_start) vac_switch[i]*vac_coverage*N_age[i]*R[i]/elig_vac_age[i] else
  vac_switch[i]*(age_rate[i - 1]*(N_age[i - 1]*N_v_age[i] - N_v_age[i - 1]*N_age[i])/N_age[i])*(R[i]/elig_vac_age[i])



#### importation
import[1:n_sero] <- if (as.integer(t) >= import_start && as.integer(t) <= import_end) n_import else 0

# force of infection from infected mosquitoes 
foi_mos[] <- b*beta_m[i]*I_m[i]/N 

# force of infection from infected human
foi_hum[] <- b*beta_h[i]*infected_hum[i]/N

# for secondary infection
S_sec[] <- sum(S[i,])


# unvaccinated primary infection
deriv(S_all[1]) <- lambda_hum*N - sum(foi_mos[])*S_all[i] - death_hum*S_all[i] -
  age_rate[i]*S_all[i] - vac_switch[i]*vac_Sall[i]

deriv(S_all[2:n_age]) <- -sum(foi_mos[])*S_all[i] - death_hum*S_all[i] +
  age_rate[i - 1]*S_all[i - 1] - age_rate[i]*S_all[i] - vac_switch[i]*vac_Sall[i]

# symptomatic
deriv(I[1,1:n_sero]) <-  foi_mos[j]*S_all[i] - (gamma_1 + death_hum)*I[i,j] - 
  age_rate[i]*I[i,j]

deriv(I[2:n_age,1:n_sero]) <- foi_mos[j]*S_all[i] - (gamma_1  + death_hum)*I[i,j] + 
  age_rate[i - 1]*I[i - 1,j] - age_rate[i]*I[i,j]


# cross-immune
deriv(C[1,1:n_sero]) <- gamma_1*I[i,j] - (alpha + death_hum)*C[i,j] - 
  age_rate[i]*C[i,j] - vac_switch[i]*vac_C[i,j]

deriv(C[2:n_age,1:n_sero]) <- gamma_1*I[i,j] - (alpha + death_hum)*C[i,j] + 
  age_rate[i - 1]*C[i - 1,j] - age_rate[i]*C[i,j] - vac_switch[i]*vac_C[i,j]

# susceptible to secondary
deriv(S[1,1:n_sero ]) <- alpha*C[i,j] - (sum(foi_mos[]) - foi_mos[j])*S[i,j] - death_hum*S[i,j] - 
  age_rate[i]*S[i,j] - vac_switch[i]*vac_S[i,j]

deriv(S[2:n_age,1:n_sero ]) <- alpha*C[i,j] -  (sum(foi_mos[]) - foi_mos[j])*S[i,j] - death_hum*S[i,j]  + 
  age_rate[i - 1]*S[i - 1,j] - age_rate[i]*S[i,j] - vac_switch[i]*vac_S[i,j]

# secondary infection

deriv(I_ij[1,1:n_sero]) <-  foi_mos[j]*(S_sec[i] - S[i,j]) - (gamma_2 +  death_hum)*I_ij[i,j] - 
  age_rate[i]*I_ij[i,j]

deriv(I_ij[2:n_age,1:n_sero]) <-  foi_mos[j]*(S_sec[i] - S[i,j]) - (gamma_2 + death_hum)*I_ij[i,j] + 
  age_rate[i - 1]*I_ij[i - 1,j] - age_rate[i]*I_ij[i,j]

# recovered
deriv(R[1]) <-  gamma_2*sum(I_ij[i,]) - death_hum*R[i] - 
  age_rate[i]*R[i] - vac_switch[i]*vac_R[i]

deriv(R[2:n_age]) <-  gamma_2*sum(I_ij[i,]) - death_hum*R[i] + 
  age_rate[i - 1]*R[i - 1] - age_rate[i]*R[i] - vac_switch[i]*vac_R[i]

# mosquito population

deriv(S_m) <- lambda_mos*(S_m + sum(E_m[]) + sum(I_m[])) - S_m*sum(foi_hum[]) - death_mos*S_m

deriv(E_m[1:n_sero]) <-  S_m*foi_hum[i] - (sigma_m + death_mos)*E_m[i]

deriv(I_m[1:n_sero]) <- import_switch[i]*import[i] + sigma_m*E_m[i] - death_mos*I_m[i]

# vaccinated population
## i- stages of vaccinated population (different estimates of VE from trial data in diff time points)
## j- age
## k - serotype
## For vaccinated we have to define 4 times. First vac stage and age groups separately; 
## and then remaining stages and age-groups

## for 1st Vac stage and age-group 1
deriv(S_v_all[1,1]) <-  vac_switch[j]*vac_Sall[j] -
  sum(foi_mos[])*(1 - effi_inf_sero_n[i])*S_v_all[i,j] - 
  death_hum*S_v_all[i,j] -
  age_rate[j]*S_v_all[i,j] - 
  vac_stage_rate[i]*S_v_all[i,j]

# 1st Vac stage and remaining age-groups
deriv(S_v_all[1,2:n_age]) <-  vac_switch[j]*vac_Sall[j] - 
  sum(foi_mos[])*(1 - effi_inf_sero_n[i])*S_v_all[i,j] - 
  death_hum*S_v_all[i,j] + 
  age_rate[j - 1]*S_v_all[i,j - 1] - 
  age_rate[j]*S_v_all[i,j] - 
  vac_stage_rate[i]*S_v_all[i,j]

## remaining vac stages & age-group 1
deriv(S_v_all[2:n_vac_stage,1]) <- -sum(foi_mos[])*(1 - effi_inf_sero_n[i])*S_v_all[i,j] - 
  death_hum*S_v_all[i,j] -
  age_rate[j]*S_v_all[i,j] + 
  vac_stage_rate[i - 1]*S_v_all[i - 1,j] -
  vac_stage_rate[i]*S_v_all[i,j]

## remaining vac stages & remaining age-groups
deriv(S_v_all[2:n_vac_stage,2:n_age]) <-  -sum(foi_mos[])*(1 - effi_inf_sero_n[i])*S_v_all[i,j] - 
  death_hum*S_v_all[i,j] + 
  age_rate[j - 1]*S_v_all[i,j - 1] - 
  age_rate[j]*S_v_all[i,j] + 
  vac_stage_rate[i - 1]*S_v_all[i - 1,j] -
  vac_stage_rate[i]*S_v_all[i,j]
## 
deriv(I_v[1,1,1:n_sero]) <- foi_mos[k]*(1 - effi_inf_sero_n[i])*S_v_all[i,j] - 
  (gamma_1 + death_hum)*I_v[i,j,k] - 
  age_rate[j]*I_v[i,j,k] -
  vac_stage_rate[i]*I_v[i,j,k]

deriv(I_v[1,2:n_age,1:n_sero]) <- foi_mos[k]*(1 - effi_inf_sero_n[i])*S_v_all[i,j] - 
  (gamma_1 + death_hum)*I_v[i,j,k] +
  age_rate[j - 1]*I_v[i,j - 1,k] - 
  age_rate[j]*I_v[i,j,k]  - 
  vac_stage_rate[i]*I_v[i,j,k]

deriv(I_v[2:n_vac_stage,1,1:n_sero]) <- foi_mos[k]*(1 - effi_inf_sero_n[i])*S_v_all[i,j] - 
  (gamma_1 + death_hum)*I_v[i,j,k] - 
  age_rate[j]*I_v[i,j,k] + 
  vac_stage_rate[i - 1]*I_v[i - 1,j,k] -
  vac_stage_rate[i]*I_v[i,j,k]

deriv(I_v[2:n_vac_stage,2:n_age,1:n_sero]) <- foi_mos[k]*(1 - effi_inf_sero_n[i])*S_v_all[i,j] - 
  (gamma_1 + death_hum)*I_v[i,j,k] +
  age_rate[j - 1]*I_v[i,j - 1,k] - 
  age_rate[j]*I_v[i,j,k] + 
  vac_stage_rate[i - 1]*I_v[i - 1,j,k] -
  vac_stage_rate[i]*I_v[i,j,k]

deriv(C_v[1,1,1:n_sero]) <- vac_switch[j]*vac_C[j,k] +
  gamma_1*I_v[i,j,k]  - 
  (alpha + death_hum)*C_v[i,j,k] -
  age_rate[j]*C_v[i,j,k] -
  vac_stage_rate[i]*C_v[i,j,k]

deriv(C_v[1,2:n_age,1:n_sero]) <- vac_switch[j]*vac_C[j,k] +
  gamma_1*I_v[i,j,k]  - 
  (alpha + death_hum)*C_v[i,j,k] + 
  age_rate[j - 1]*C_v[i,j - 1,k] - 
  age_rate[j]*C_v[i,j,k]  - 
  vac_stage_rate[i]*C_v[i,j,k]

deriv(C_v[2:n_vac_stage,1,1:n_sero]) <- gamma_1*I_v[i,j,k]  - 
  (alpha + death_hum)*C_v[i,j,k] -
  age_rate[j]*C_v[i,j,k] + 
  vac_stage_rate[i - 1]*C_v[i - 1,j,k] -
  vac_stage_rate[i]*C_v[i,j,k]

deriv(C_v[2:n_vac_stage,2:n_age,1:n_sero]) <- gamma_1*I_v[i,j,k]  - 
  (alpha + death_hum)*C_v[i,j,k] + 
  age_rate[j - 1]*C_v[i,j - 1,k] - 
  age_rate[j]*C_v[i,j,k] + 
  vac_stage_rate[i - 1]*C_v[i - 1,j,k] -
  vac_stage_rate[i]*C_v[i,j,k]

deriv(S_v[1,1,1:n_sero ]) <- vac_switch[j]*vac_S[j,k] +
  alpha*C_v[i,j,k] - 
  (sum(foi_mos[]) - foi_mos[k])*(1 - effi_inf_sero_p[i])*S_v[i,j,k] - 
  death_hum*S_v[i,j,k] -
  age_rate[j]*S_v[i,j,k] -
  vac_stage_rate[i]*S_v[i,j,k]


deriv(S_v[1,2:n_age,1:n_sero ]) <-  vac_switch[j]*vac_S[j,k]  +
  alpha*C_v[i,j,k] -  
  (sum(foi_mos[]) - foi_mos[k])*(1 - effi_inf_sero_p[i])*S_v[i,j,k] - 
  death_hum*S_v[i,j,k] + 
  age_rate[j - 1]*S_v[i,j - 1,k] - 
  age_rate[j]*S_v[i,j,k] - 
  vac_stage_rate[i]*S_v[i,j,k]

deriv(S_v[2:n_vac_stage,1,1:n_sero ]) <- alpha*C_v[i,j,k] - 
  (sum(foi_mos[]) - foi_mos[k])*(1 - effi_inf_sero_p[i])*S_v[i,j,k] - 
  death_hum*S_v[i,j,k] -
  age_rate[j]*S_v[i,j,k] + 
  vac_stage_rate[i - 1]*S_v[i - 1,j,k] -
  vac_stage_rate[i]*S_v[i,j,k]

deriv(S_v[2:n_vac_stage,2:n_age,1:n_sero ]) <- alpha*C_v[i,j,k] -  
  (sum(foi_mos[]) - foi_mos[k])*(1 - effi_inf_sero_p[i])*S_v[i,j,k] - 
  death_hum*S_v[i,j,k] + 
  age_rate[j - 1]*S_v[i,j - 1,k] - 
  age_rate[j]*S_v[i,j,k] + 
  vac_stage_rate[i - 1]*S_v[i - 1,j,k] -
  vac_stage_rate[i]*S_v[i,j,k]

# secondary infection for vaccinated
## for secondary infection (vaccinated)
S_sec_v[,] <- sum(S_v[i,j,])

# S_v[i,j,1] + S_v[i,j,2] + S_v[i,j,3] + S_v[i,j,4]

deriv(I_v_ij[1,1,1:n_sero]) <- foi_mos[k]*(1 - effi_inf_sero_p[i])*(S_sec_v[i,j]- S_v[i,j,k]) -
  (gamma_2 + death_hum)*I_v_ij[i,j,k] - 
  age_rate[j]*I_v_ij[i,j,k] -
  vac_stage_rate[i]*I_v_ij[i,j,k]

deriv(I_v_ij[1,2:n_age,1:n_sero]) <- foi_mos[k]*(1 - effi_inf_sero_p[i])*(S_sec_v[i,j]- S_v[i,j,k]) - 
  (gamma_2 + death_hum)*I_v_ij[i,j,k] +
  age_rate[j  - 1]*I_v_ij[i,j - 1,k] - 
  age_rate[j]*I_v_ij[i,j,k] -
  vac_stage_rate[i]*I_v_ij[i,j,k]

deriv(I_v_ij[2:n_vac_stage,1,1:n_sero]) <- foi_mos[k]*(1 - effi_inf_sero_p[i])*(S_sec_v[i,j]- S_v[i,j,k]) -
  (gamma_2 + death_hum)*I_v_ij[i,j,k] - 
  age_rate[j]*I_v_ij[i,j,k] + 
  vac_stage_rate[i - 1]*I_v_ij[i - 1,j,k] -
  vac_stage_rate[i]*I_v_ij[i,j,k]

deriv(I_v_ij[2:n_vac_stage,2:n_age,1:n_sero]) <- foi_mos[k]*(1 - effi_inf_sero_p[i])*(S_sec_v[i,j] - S_v[i,j,k]) - 
  (gamma_2 + death_hum)*I_v_ij[i,j,k] +
  age_rate[j  - 1]*I_v_ij[i,j - 1,k] - 
  age_rate[j]*I_v_ij[i,j,k] + 
  vac_stage_rate[i - 1]*I_v_ij[i - 1,j,k] -
  vac_stage_rate[i]*I_v_ij[i,j,k]

deriv(R_v[1,1]) <- vac_switch[j]*vac_R[j] +
  gamma_2*sum(I_v_ij[i,j,]) - 
  death_hum*R_v[i,j] -
  age_rate[j]*R_v[i,j] -
  vac_stage_rate[i]*R_v[i,j]

deriv(R_v[1,2:n_age]) <-  vac_switch[j]*vac_R[j] +
  gamma_2*sum(I_v_ij[i,j,]) -
  death_hum*R_v[i,j] +
  age_rate[j - 1]*R_v[i,j - 1] - 
  age_rate[j]*R_v[i,j] -
  vac_stage_rate[i]*R_v[i,j]


deriv(R_v[2:n_vac_stage,1]) <-  gamma_2*sum(I_v_ij[i,j,]) - 
  death_hum*R_v[i,j] -
  age_rate[j]*R_v[i,j] + 
  vac_stage_rate[i - 1]*R_v[i - 1,j] -
  vac_stage_rate[i]*R_v[i,j]

deriv(R_v[2:n_vac_stage,2:n_age]) <- gamma_2*sum(I_v_ij[i,j,]) -
  death_hum*R_v[i,j] +
  age_rate[j - 1]*R_v[i,j - 1] - 
  age_rate[j]*R_v[i,j] + 
  vac_stage_rate[i - 1]*R_v[i - 1,j] -
  vac_stage_rate[i]*R_v[i,j]

###### output of interest

####  UNVACCINATED
### Infection
# primary
inc_inf_pri[,] <- foi_mos[j]*S_all[i]
# secondary
inc_inf_sec[,] <- foi_mos[j]*(S_sec[i] - S[i,j])

#### Vaccinated
### Infection
# extra summing over stages of vaccinated
inf_v_pri_vac_stage[,,] <- (1 - effi_inf_sero_n[i])*foi_mos[k]*S_v_all[i,j]
#
inf_v_sec_vac_stage[,,] <- (1 - effi_inf_sero_p[i])*foi_mos[k]*(S_sec_v[i,j] - S_v[i,j,k])
# primary
inc_inf_v_pri[,] <- sum(inf_v_pri_vac_stage[,i,j])
# secondary
inc_inf_v_sec[,] <- sum(inf_v_sec_vac_stage[,i,j])

## symptomatic
symp_v_pri_vac_stage[,,] <- (1 - effi_symp_sero_n[i,k])*inf_v_pri_vac_stage[i,j,k]
inc_symp_v_pri[,] <- sum(symp_v_pri_vac_stage[,i,j])

symp_v_sec_vac_stage[,,] <- (1 - effi_symp_sero_p[i,k])*inf_v_sec_vac_stage[i,j,k]
inc_symp_v_sec[,] <- sum(symp_v_sec_vac_stage[,i,j])

## hospitalized 
hosp_v_pri_vac_stage[,,] <- (1 - effi_hos_sero_n[i,k])*symp_v_pri_vac_stage[i,j,k]
hosp_v_sec_vac_stage[,,] <- (1 - effi_hos_sero_p[i,k])*symp_v_sec_vac_stage[i,j,k]

inc_hosp_v_pri[,] <- sum(hosp_v_pri_vac_stage[,i,j])
inc_hosp_v_sec[,] <- sum(hosp_v_sec_vac_stage[,i,j])


### population in each age group
output(tot_age[]) <- N_age[i]
### Vaccinated 
output(tot_vac_age[]) <- N_v_age[i]
### unVaccinated 
output(tot_nonvac_age[]) <- N_nv_age[i]

##incidence of infection
output(inc_inf[,]) <- inc_inf_pri[i,j] + inc_inf_sec[i,j]
# # incidence of infection 
output(inc_inf_v[,]) <- inc_inf_v_pri[i,j] + inc_inf_v_sec[i,j]

# output(inc_inf_primary[,]) <- inc_inf_pri[i,j] 
# output(inc_inf_secondary[,]) <- inc_inf_sec[i,j]
# # incidence of infection 
# output(inc_inf_v_primary[,]) <- inc_inf_v_pri[i,j] 
# output(inc_inf_v_secondary[,]) <- inc_inf_v_sec[i,j]

## incidence of symptomatic
output(inc_symp[,]) <- rho_1[i]*inc_inf_pri[i,j] + rho_2[i]*inc_inf_sec[i,j]
# # incidence of symptomatic 
output(inc_symp_v[,]) <- rho_1[i]*inc_symp_v_pri[i,j] + rho_2[i]*inc_symp_v_sec[i,j]



# output(inc_symp_primary[,]) <- rho_1[i]*inc_inf_pri[i,j] 
# output(inc_symp_secondary[,]) <- rho_2[i]*inc_inf_sec[i,j]
# 
# 
# # incidence of symptomatic 
# output(inc_symp_v_primary[,]) <- rho_1[i]*inc_symp_v_pri[i,j] 
# output(inc_symp_v_secondary[,]) <- rho_2[i]*inc_symp_v_sec[i,j]


## incidence of hospitalized
output(inc_hosp[]) <- xi_1*rho_1[i]*sum(inc_inf_pri[i,]) + xi_2*rho_2[i]*sum(inc_inf_sec[i,])
# # incidence of hospitalization
output(inc_hosp_v[]) <- xi_1*rho_1[i]*sum(inc_hosp_v_pri[i,]) + xi_2*rho_2[i]*sum(inc_hosp_v_sec[i,])


# output(inc_hosp_primary[]) <- xi_1*rho_1[i]*sum(inc_inf_pri[i,]) 
# output(inc_hosp_secondary[]) <- xi_2*rho_2[i]*sum(inc_inf_sec[i,])
# # incidence of hospitalization
# output(inc_hosp_v_primary[]) <- xi_1*rho_1[i]*sum(inc_hosp_v_pri[i,]) 
# output(inc_hosp_v_secondary[]) <-  xi_2*rho_2[i]*sum(inc_hosp_v_sec[i,])


##
output(new_vac_age[]) <- vac_switch[i]*(vac_Sall[i] + sum(vac_C[i,]) + sum(vac_S[i,]) + vac_R[i])

# output(new_vac_age_primary[]) <- vac_switch[i]*(vac_Sall[i] ) 
# output(new_vac_age_secondary[]) <- vac_switch[i]*(sum(vac_C[i,]) + sum(vac_S[i,]) + vac_R[i])


## Parameters
# time dependent

lambda_hum <- lambda_hum_tt[as.integer(t)]  
death_hum <-  death_hum_tt[as.integer(t)] 

# define the corresponding vector contains per day value
lambda_hum_tt[] <- user()
death_hum_tt[] <- user()

# define the corresponding dimensions 
dim(lambda_hum_tt) <- user()
dim(death_hum_tt) <- user()


## importation
import_switch[] <- user()
dim(import_switch) <- n_sero
dim(import) <- n_sero
n_import <- user()
import_start <- user()
import_end <- user()


## fixed parameters
# scalar valued
b <- user()
alpha <- user()
gamma_1 <- user()
gamma_2 <- user()
xi_1 <- user()
xi_2 <- user()
# kai <- user()
n_age <- user() 
n_sero <- user() 
sigma_m <- user()
lambda_mos <- user()
death_mos <- user()

# fixed but vector valued
rho_1[] <- user()
rho_2[] <- user()
beta_m[] <- user()
beta_h[] <- user()  
age_rate[] <- user()

#dims
dim(rho_1) <- n_age
dim(rho_2) <- n_age
dim(beta_m) <- n_sero
dim(beta_h) <- n_sero
dim(age_rate) <- n_age

# vaccine efficacy
vac_switch[] <- user()
vac_coverage <- user()
vac_start <- user()
vac_stage_rate[] <- user()
n_vac_stage <- user()

effi_inf_sero_n[] <- user()
effi_inf_sero_p[] <- user()
effi_symp_sero_n[,] <- user() 
effi_symp_sero_p[,] <- user() 
effi_hos_sero_n[,] <- user()
effi_hos_sero_p[,] <- user()

# dims
dim(vac_switch) <- n_age
dim(effi_inf_sero_n) <- c(n_vac_stage)
dim(effi_inf_sero_p) <- c(n_vac_stage)
dim(effi_symp_sero_n) <- c(n_vac_stage, n_sero) 
dim(effi_symp_sero_p) <- c(n_vac_stage, n_sero) 
dim(effi_hos_sero_n) <- c(n_vac_stage, n_sero) 
dim(effi_hos_sero_p) <- c(n_vac_stage, n_sero)  

##initial conditions
initial(S_all[]) <- S0_all[i]
initial(I[,]) <- I0[i,j]
initial(C[,]) <- C0[i,j]
initial(S[,]) <- S0[i,j]
initial(I_ij[,]) <- I_ij0[i,j]
initial(R[]) <- R0[i]
initial(S_m) <- S_m0 
initial(E_m[]) <- E_m0[i]
initial(I_m[]) <- I_m0[i]
initial(S_v_all[,]) <- S0_v_all[i,j]
initial(I_v[,,]) <- I_v0[i,j,k]
initial(C_v[,,]) <- C_v0[i,j,k]
initial(S_v[,,]) <- S_v0[i,j,k]
initial(I_v_ij[,,]) <- I_v_ij0[i,j,k]
initial(R_v[,]) <- R_v0[i,j]

# initial conditions
S0_all[] <- user()
I0[,] <- user()
C0[,] <- user()
S0[,] <- user()
I_ij0[,] <- user()
R0[] <- user()
S_m0 <- user()
E_m0[] <- user()
I_m0[] <- user()
S0_v_all[,] <- user()
I_v0[,,] <- user()
C_v0[,,] <- user()
S_v0[,,] <- user()
I_v_ij0[,,] <- user()
R_v0[,] <- user()

# dims
dim(S0_all) <- n_age
dim(I0) <- c(n_age, n_sero)
dim(C0) <- c(n_age, n_sero)
dim(S0) <- c(n_age, n_sero)
dim(I_ij0) <- c(n_age, n_sero)
dim(R0) <- n_age
dim(E_m0) <- n_sero
dim(I_m0) <- n_sero
dim(S0_v_all) <- c(n_vac_stage, n_age)
dim(I_v0) <- c(n_vac_stage,n_age, n_sero)
dim(C_v0) <- c(n_vac_stage,n_age, n_sero)
dim(S_v0) <- c(n_vac_stage,n_age, n_sero)
dim(I_v_ij0) <- c(n_vac_stage,n_age, n_sero)
dim(R_v0) <- c(n_vac_stage, n_age)

# state variables
dim(S_all) <- n_age
dim(I) <- c(n_age, n_sero)
dim(C) <- c(n_age, n_sero)
dim(S) <- c(n_age, n_sero)
dim(I_ij) <- c(n_age, n_sero)
dim(R) <- n_age
dim(E_m) <- n_sero
dim(I_m) <- n_sero
dim(S_v_all) <-  c(n_vac_stage, n_age)
dim(I_v) <- c(n_vac_stage,n_age, n_sero)
dim(C_v) <- c(n_vac_stage,n_age, n_sero)
dim(S_v) <- c(n_vac_stage,n_age, n_sero)
dim(I_v_ij) <- c(n_vac_stage,n_age, n_sero)
dim(R_v) <-  c(n_vac_stage, n_age)

# extra quantities
## EXTRA FOR OUTPUTS 
dim(inc_inf_pri) <- c(n_age, n_sero)
dim(inc_inf_sec) <- c(n_age, n_sero)
dim(inc_inf_v_pri) <- c(n_age, n_sero)
dim(inc_inf_v_sec) <- c(n_age, n_sero)
dim(inc_symp_v_pri) <- c(n_age, n_sero)
dim(inc_symp_v_sec) <- c(n_age, n_sero)
dim(inc_hosp_v_pri) <- c(n_age, n_sero)
dim(inc_hosp_v_sec) <- c(n_age, n_sero)
# dim(inc_hosp_pri) <- c(n_age, n_sero)
# dim(inc_hosp_sec) <- c(n_age, n_sero)
dim(inf_v_pri_vac_stage) <- c(n_vac_stage,n_age, n_sero)
dim(inf_v_sec_vac_stage) <- c(n_vac_stage,n_age, n_sero)
dim(symp_v_pri_vac_stage) <- c(n_vac_stage,n_age, n_sero)
dim(symp_v_sec_vac_stage) <- c(n_vac_stage,n_age, n_sero)
dim(hosp_v_pri_vac_stage) <- c(n_vac_stage,n_age, n_sero)
dim(hosp_v_sec_vac_stage) <- c(n_vac_stage,n_age, n_sero)

dim(infected_hum) <- n_sero
dim(foi_mos) <- n_sero
dim(foi_hum) <- n_sero
dim(S_sec) <- n_age
dim(S_sec_v) <- c(n_vac_stage,n_age)
dim(N_age) <- n_age
dim(N_v_age) <- n_age
dim(N_nv_age) <- n_age

dim(elig_vac_age) <- n_age
dim(vac_Sall) <- n_age
dim(vac_C) <- c(n_age, n_sero)
dim(vac_S) <- c(n_age, n_sero)
dim(vac_R) <- n_age
dim(vac_stage_rate) <- n_vac_stage

## for output
# dim(inc_inf) <- c(n_age,n_sero)
# dim(inc_inf_v) <- c(n_age, n_sero)
# dim(inc_symp) <- c(n_age,n_sero)
# dim(inc_symp_v) <- c(n_age, n_sero)
# dim(inc_hosp) <- n_age
# dim(inc_hosp_v) <- n_age



# dim(inc_inf_primary) <- c(n_age,n_sero)
# dim(inc_inf_secondary) <- c(n_age,n_sero)
# dim(inc_inf_v_primary) <- c(n_age, n_sero)
# dim(inc_inf_v_secondary) <- c(n_age, n_sero)
# dim(inc_symp_primary) <- c(n_age,n_sero)
# dim(inc_symp_secondary) <- c(n_age,n_sero)
# dim(inc_symp_v_primary) <- c(n_age, n_sero)
# dim(inc_symp_v_secondary) <- c(n_age, n_sero)
# dim(inc_hosp_primary) <- n_age
# dim(inc_hosp_secondary) <- n_age
# dim(inc_hosp_v_primary) <- n_age
# dim(inc_hosp_v_secondary) <- n_age

dim(inc_inf) <- c(n_age,n_sero)
dim(inc_inf_v) <- c(n_age, n_sero)
dim(inc_symp) <- c(n_age,n_sero)
dim(inc_symp_v) <- c(n_age, n_sero)
dim(inc_hosp) <- n_age
dim(inc_hosp_v) <- n_age

dim(tot_age) <- n_age
dim(tot_vac_age) <- n_age
dim(tot_nonvac_age) <- n_age
dim(new_vac_age) <- n_age
# dim(new_vac_age_primary) <- n_age
# dim(new_vac_age_secondary) <- n_age