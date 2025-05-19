
# total population
N <- sum(S_all[]) + sum(I[,]) + sum(A[,]) + sum(C[,]) + sum(S[,]) +
  sum(I_ij[,]) + sum(A_ij[,]) + sum(H[]) + sum(R[]) +
  sum(S_v_all[]) + sum(I_v[,]) + sum(A_v[,]) + sum(C_v[,]) + sum(S_v[,]) +
  sum(I_v_ij[,]) + sum(A_v_ij[,]) + sum(H_v[]) + sum(R_v[])


infected_hum[] <- sum(I[,i]) + kai*sum(A[,i]) + sum(I_ij[,i]) + kai*sum(A_ij[,i]) + 
  sum(I_v[,i]) + kai_v*sum(A_v[,i]) + sum(I_v_ij[,i]) + kai_v*sum(A_v_ij[,i])

## population each age group
N_age[] <- S_all[i] + sum(I[i,]) + sum(A[i,]) + sum(C[i,]) + sum(S[i,]) +
  sum(I_ij[i,]) + sum(A_ij[i,]) + H[i] + R[i] +
  S_v_all[i] + sum(I_v[i,]) + sum(A_v[i,]) + sum(C_v[i,]) + sum(S_v[i,]) +
  sum(I_v_ij[i,]) + sum(A_v_ij[i,]) + H_v[i] + sum(R_v[i])

beta_m[1] <- beta_m1
beta_m[2] <- beta_m2
beta_m[3] <- beta_m3
beta_m[4] <- beta_m4

# force of infection of infected mosquitoes
foi_mos[] <- b*beta_m[i]*I_m[i]/N 

# force of infection of infected human
foi_hum[] <- b*beta_h[i]*infected_hum[i]/N

## for secondary infection
S_sec[] <- sum(S[i,])

## for secondary infection
S_sec_v[] <- sum(S_v[i,])


vac_rate[] <- vac_coverage*N_age[i]/(S_all[i] + sum(C[i,]) + sum(S[i,]) + R[i])

# calculates the vaccinated hospitalization as odin does not allow multiplication inside "sum"!
hos_vac_1[,] <- (1 - effi_hos_sero_n[i,j])*I_v[i,j]             # primary infection
hos_vac_2[,] <- (1 - effi_hos_sero_p[i,j])*I_v_ij[i,j]    # secondary infection

# unvaccinated primary infection
deriv(S_all[1]) <- lambda_hum*N - sum(foi_mos[])*S_all[i] -death_hum[i]*S_all[i] -
                   age_rate[i]*S_all[i] - vac_switch[i]*vac_rate[i]*S_all[i]

deriv(S_all[2:n_age]) <- -sum(foi_mos[])*S_all[i] -death_hum[i]*S_all[i] +
                          age_rate[i-1]*S_all[i-1] -age_rate[i]*S_all[i] -
                          vac_switch[i]*vac_rate[i]*S_all[i]

deriv(I[1,1:n_sero]) <-  rho_1[i]*foi_mos[j]*S_all[i] - (gamma_1+xi_1+death_hum[i])*I[i,j] - 
                        age_rate[i]*I[i,j]

deriv(I[2:n_age,1:n_sero]) <- rho_1[i]*foi_mos[j]*S_all[i] - (gamma_1+xi_1+death_hum[i])*I[i,j] + 
                              age_rate[i-1]*I[i-1,j] - age_rate[i]*I[i,j]

deriv(A[1,1:n_sero]) <-   (1-rho_1[i])*foi_mos[j]*S_all[i] - (gamma_1+death_hum[i])*A[i,j] - 
                         age_rate[i]*A[i,j]

deriv(A[2:n_age,1:n_sero]) <-  (1-rho_1[i])*foi_mos[j]*S_all[i] - (gamma_1+death_hum[i])*A[i,j] + 
                              age_rate[i-1]*A[i-1,j] - age_rate[i]*A[i,j]

deriv(C[1,1:n_sero]) <- gamma_1*I[i,j]+ gamma_1*A[i,j] - (alpha + death_hum[i])*C[i,j] - 
                        age_rate[i]*C[i,j] - vac_switch[i]*vac_rate[i]*C[i,j]

deriv(C[2:n_age,1:n_sero]) <- gamma_1*I[i,j]+ gamma_1*A[i,j] - (alpha + death_hum[i])*C[i,j] + 
                              age_rate[i-1]*C[i-1,j] - age_rate[i]*C[i,j] - 
                              vac_switch[i]*vac_rate[i]*C[i,j]

deriv(S[1,1:n_sero ]) <- alpha*C[i,j] - (sum(foi_mos[]) - foi_mos[j])*S[i,j] - death_hum[i]*S[i,j] - 
                         age_rate[i]*S[i,j] - vac_switch[i]*vac_rate[i]*S[i,j] 

deriv(S[2:n_age,1:n_sero ]) <- alpha*C[i,j] -  (sum(foi_mos[]) - foi_mos[j])*S[i,j] - death_hum[i]*S[i,j]  + 
                               age_rate[i-1]*S[i-1,j] - age_rate[i]*S[i,j] - vac_switch[i]*vac_rate[i]*S[i,j]

# secondary infection

deriv(I_ij[1,1:n_sero]) <-  rho_2[i]*foi_mos[j]*(S_sec[i] - S[i,j]) - (gamma_2+xi_2+death_hum[i])*I_ij[i,j] - 
  age_rate[i]*I_ij[i,j]

deriv(I_ij[2:n_age,1:n_sero]) <-  rho_2[i]*foi_mos[j]*(S_sec[i] - S[i,j]) - (gamma_2+xi_2+death_hum[i])*I_ij[i,j] + 
  age_rate[i-1]*I_ij[i-1,j] - age_rate[i]*I_ij[i,j]


deriv(A_ij[1,1:n_sero])  <-  (1-rho_2[i])*foi_mos[j]*(S_sec[i] - S[i,j]) - (gamma_2 + death_hum[i])*A_ij[i,j] - 
  age_rate[i]*A_ij[i,j] 

deriv(A_ij[2:n_age,1:n_sero])  <-  (1-rho_2[i])*foi_mos[j]*(S_sec[i] - S[i,j]) - (gamma_2 + death_hum[i])*A_ij[i,j] + 
  age_rate[i-1]*A_ij[i-1,j] - age_rate[i]*A_ij[i,j] 


# hospitalized and recovered

deriv(H[1]) <- xi_1*sum(I[i,]) + xi_2*sum(I_ij[i,]) - (gamma_2 + death_hum[i])*H[i] - 
  age_rate[i]*H[i]

deriv(H[2:n_age]) <- xi_1*sum(I[i,]) + xi_2*sum(I_ij[i,])  - (gamma_2 + death_hum[i])*H[i] + 
  age_rate[i-1]*H[i-1] - age_rate[i]*H[i]

deriv(R[1]) <-  gamma_2*sum(I_ij[i,]) + gamma_2*sum(A_ij[i,])+ gamma_2*H[i] - death_hum[i]*R[i] - 
  age_rate[i]*R[i] - vac_switch[i]*vac_rate[i]*R[i]

deriv(R[2:n_age]) <-  gamma_2*sum(I_ij[i,]) + gamma_2*sum(A_ij[i,]) + gamma_2*H[i] - death_hum[i]*R[i] + 
  age_rate[i-1]*R[i-1] - age_rate[i]*R[i] - vac_switch[i]*vac_rate[i]*R[i]

# mosquito population

deriv(S_m) <- lambda_mos*(S_m + sum(E_m[]) + sum(I_m[])) - S_m*sum(foi_hum[]) - death_mos*S_m

# deriv(S_m) <- lambda_mos - S_m*sum(foi_hum[]) - death_mos*S_m

deriv(E_m[1:n_sero]) <-  S_m*foi_hum[i] - (sigma_m + death_mos)*E_m[i]

deriv(I_m[1:n_sero]) <- import[i] + sigma_m*E_m[i] - death_mos*I_m[i]

# vaccinated population
deriv(S_v_all[1]) <-  -sum(foi_mos[])*(1 - effi_inf_sero_n[i])*S_v_all[i] - death_hum[i]*S_v_all[i] +
                       vac_switch[i]*vac_rate[i]*S_all[i] - age_rate[i]*S_v_all[i]

deriv(S_v_all[2:n_age]) <-  -sum(foi_mos[])*(1 - effi_inf_sero_n[i])*S_v_all[i] - death_hum[i]*S_v_all[i] + 
                            vac_switch[i]*vac_rate[i]*S_all[i] + age_rate[i-1]*S_v_all[i-1] - age_rate[i]*S_v_all[i]

deriv(I_v[1,1:n_sero]) <- rho_v_1*(1-effi_symp_sero_n[i,j])*foi_mos[j]*(1-effi_inf_sero_n[i])*S_v_all[i] - 
                         (gamma_v_1+xi_v_1*(1-effi_hos_sero_n[i,j])+death_hum[i])*I_v[i,j] - 
                          age_rate[i]*I_v[i,j]

deriv(I_v[2:n_age,1:n_sero]) <- rho_v_1*(1-effi_symp_sero_n[i,j])*foi_mos[j]*(1-effi_inf_sero_n[i])*S_v_all[i] - 
                                (gamma_v_1+xi_v_1*(1-effi_hos_sero_n[i,j])+death_hum[i])*I_v[i,j] +
                                 age_rate[i-1]*I_v[i-1,j] - age_rate[i]*I_v[i,j]

deriv(A_v[1,1:n_sero]) <- (1-rho_v_1*(1-effi_symp_sero_n[i,j]))*foi_mos[j]*(1-effi_inf_sero_n[i])*S_v_all[i] - 
                          (gamma_v_1+death_hum[i])*A_v[i,j] - 
                           age_rate[i]*A_v[i,j]

deriv(A_v[2:n_age,1:n_sero]) <- (1-rho_v_1*(1-effi_symp_sero_n[i,j]))*foi_mos[j]*(1-effi_inf_sero_n[i])*S_v_all[i] - 
                                (gamma_v_1+death_hum[i])*A_v[i,j] +
                                age_rate[i-1]*A_v[i-1,j] - age_rate[i]*A_v[i,j]

deriv(C_v[1,1:n_sero]) <- gamma_v_1*I_v[i,j]+ gamma_v_1*A_v[i,j] - (alpha_v + death_hum[i])*C_v[i,j]+
                          vac_switch[i]*vac_rate[i]*C[i,j] -age_rate[i]*C_v[i,j]

deriv(C_v[2:n_age,1:n_sero]) <- gamma_v_1*I_v[i,j]+ gamma_v_1*A_v[i,j] - (alpha_v + death_hum[i])*C_v[i,j] + 
                                vac_switch[i]*vac_rate[i]*C[i,j]+ age_rate[i-1]*C_v[i-1,j]- age_rate[i]*C_v[i,j]

deriv(S_v[1,1:n_sero ]) <- alpha_v*C_v[i,j] - (sum(foi_mos[]) - foi_mos[j])*(1-effi_inf_sero_p[i])*S_v[i,j] - death_hum[i]*S_v[i,j] + 
                           vac_switch[i]*vac_rate[i]*S[i,j]   - age_rate[i]*S_v[i,j]

deriv(S_v[2:n_age,1:n_sero ]) <- alpha_v*C_v[i,j] -  (sum(foi_mos[]) - foi_mos[j])*(1-effi_inf_sero_p[i])*S_v[i,j] - death_hum[i]*S_v[i,j] + 
                                 vac_switch[i]*vac_rate[i]*S[i,j]+age_rate[i-1]*S_v[i-1,j] - age_rate[i]*S_v[i,j]

# secondary infection for vaccinated

deriv(I_v_ij[1,1:n_sero]) <- rho_v_2*(1-effi_symp_sero_p[i,j])*foi_mos[j]*(1-effi_inf_sero_p[i])*(S_sec_v[i] - S_v[i,j]) -
  (gamma_v_2+xi_v_2*(1-effi_hos_sero_p[i,j])+death_hum[i])*I_v_ij[i,j] - 
  age_rate[i]*I_v_ij[i,j]

deriv(I_v_ij[2:n_age,1:n_sero]) <- rho_v_2*(1-effi_symp_sero_p[i,j])*foi_mos[j]*(1-effi_inf_sero_p[i])*(S_sec_v[i] - S_v[i,j]) - 
  (gamma_v_2+xi_v_2*(1-effi_hos_sero_p[i,j])+death_hum[i])*I_v_ij[i,j]+
  age_rate[i-1]*I_v_ij[i-1,j] - age_rate[i]*I_v_ij[i,j]

deriv(A_v_ij[1,1:n_sero])  <- (1-rho_v_2*(1-effi_symp_sero_p[i,j]))*foi_mos[j]*(1-effi_inf_sero_p[i])*(S_sec_v[i] - S_v[i,j]) - 
  (gamma_v_2 + death_hum[i])*A_v_ij[i,j] - 
  age_rate[i]*A_v_ij[i,j]

deriv(A_v_ij[2:n_age,1:n_sero]) <- (1-rho_v_2*(1-effi_symp_sero_p[i,j]))*foi_mos[j]*(1-effi_inf_sero_p[i])*(S_sec_v[i] - S_v[i,j]) - 
  (gamma_v_2 + death_hum[i])*A_v_ij[i,j] +
  age_rate[i-1]*A_v_ij[i-1,j] - age_rate[i]*A_v_ij[i,j]

# hospitalized and recovered for vaccinated
deriv(H_v[1]) <- xi_v_1*sum(hos_vac_1[i,]) + xi_v_2*sum(hos_vac_2[i,]) - (gamma_v_2 + death_hum[i])*H_v[i] - 
  age_rate[i]*H_v[i]

deriv(H_v[2:n_age]) <- xi_v_1*sum(hos_vac_1[i,]) + xi_v_2*sum(hos_vac_2[i,]) - (gamma_v_2 + death_hum[i])*H_v[i] + 
  age_rate[i-1]*H_v[i-1] - age_rate[i]*H_v[i]

deriv(R_v[1]) <-  gamma_v_2*sum(I_v_ij[i,])+ gamma_v_2*sum(A_v_ij[i,])+ gamma_v_2*H_v[i] - death_hum[i]*R_v[i] +
  vac_switch[i]*vac_rate[i]*R[i]- age_rate[i]*R_v[i]

deriv(R_v[2:n_age]) <- gamma_v_2*sum(I_v_ij[i,])+ gamma_v_2*sum(A_v_ij[i,])+ gamma_v_2*H_v[i]  - death_hum[i]*R_v[i] +
  vac_switch[i]*vac_rate[i]*R[i]+ age_rate[i-1]*R_v[i-1] - age_rate[i]*R_v[i]

# total population in age groups
deriv(total_pop_age[1]) <- lambda_hum*N - death_hum[i]*total_pop_age[i] - age_rate[i]*total_pop_age[i]

deriv(total_pop_age[2:n_age]) <- age_rate[i-1]*total_pop_age[i-1] - age_rate[i]*total_pop_age[i] - death_hum[i]*total_pop_age[i]


### output of interest

## incidence of primary symptomatic
output(inc_primary_symp[,]) <- rho_1[i]*foi_mos[j]*S_all[i]

## incidence of secondary symptomatic

output(inc_secondary_symp[,]) <- rho_2[i]*foi_mos[j]*(S_sec[i] - S[i,j])

## incidence of total primary infection
output(inc_primary_inf[,]) <- foi_mos[j]*S_all[i]

## incidence of total symptomatic infection
# output(inc_symp) <- sum(age_symp_inf_pri[]) + sum(age_symp_inf_sec[])

## incidence of total secondary infection
output(inc_secondary_inf[,]) <- foi_mos[j]*(S_sec[i] - S[i,j])

### for vaccinated

## incidence of primary symptomatic
output(inc_primary_symp_v[,]) <- rho_v_1*(1 - effi_symp_sero_n[i,j])*foi_mos[j]*(1 - effi_inf_sero_n[i])*S_v_all[i]

## incidence of secondary symptomatic

output(inc_secondary_symp_v[,]) <- rho_v_2*(1-effi_symp_sero_p[i,j])*foi_mos[j]*(1-effi_inf_sero_p[i])*(S_sec_v[i] - S_v[i,j])

## incidence of total primary infection
output(inc_primary_inf_v[,]) <- (1-effi_symp_sero_n[i,j])*foi_mos[j]*(1-effi_inf_sero_n[i])*S_v_all[i]

## incidence of total secondary infection
output(inc_secondary_inf_v[,]) <- (1-effi_symp_sero_p[i,j])*foi_mos[j]*(1-effi_inf_sero_p[i])*(S_sec_v[i] - S_v[i,j])

### incidence of hosp
output(inc_hosp[]) <- xi_1*sum(I[i,]) + xi_2*sum(I_ij[i,])

output(inc_hosp_v[]) <- xi_v_1*sum(hos_vac_1[i,]) + xi_v_2*sum(hos_vac_2[i,])


output(mos_foi[]) <- foi_mos[i]

output(hum_foi[]) <- foi_hum[i]

output(N_total_age[]) <- N_age[i]


lambda_hum <- lambda_hum_tt[as.integer(t)]  

death_hum[] <-  death_hum_tt[as.integer(t),i] 

# import_pri[,] <- import_tt_pri[as.integer(t),i,j]
# 
# import_sec[,] <- import_tt_sec[as.integer(t),i,j]

import[] <- import_tt[as.integer(t),i]


b <- b_tt[as.integer(t)] 
vac_rate <- vac_rate_tt[as.integer(t)] #user()


# parameters
# b <- user()
kai <- user()
kai_v <- user()
# rho_1 <- user()
# rho_2 <- user()

rho_1[] <- user()
rho_2[] <- user()

rho_v_1 <- user()
rho_v_2 <- user()
sigma_m <- user()

# disease progression related
alpha <- user()
alpha_v <- user()
gamma_1 <- user()
gamma_2 <- user()
gamma_v_1 <- user()
gamma_v_2 <- user()
xi_1 <- user()
xi_2 <- user()
xi_v_1 <- user()
xi_v_2 <- user()

# demography related
lambda_mos <- user()
death_mos <- user()
lambda_hum_tt[] <- user()
death_hum_tt[,] <- user()
b_tt[] <- user()
vac_rate_tt[] <- user()
import_tt[,] <- user()

# import_tt_pri[,,] <- user()
# import_tt_sec[,,] <- user()
# inc_symp <- user()
# age_symp_inf_pri[] <- user() 
# age_symp_inf_sec[] <- user()
# age_sero_symp_inf_sec[,] <- user()

# LIST OF USER DEFINED PARAMETERS

# time dependent parameter for using interpolation
tt[] <- user()
vac_switch[] <- user()

# vaccine efficacy
effi_inf_sero_n[] <- user()
effi_inf_sero_p[] <- user()
effi_symp_sero_n[,] <- user() 
effi_symp_sero_p[,] <- user() 
effi_hos_sero_n[,] <- user()
effi_hos_sero_p[,] <- user()

# age-group, serotype and transmission rates and age rate
n_age <- user() 
n_sero <- user()  

beta_m1 <- user()
beta_m2 <- user()
beta_m3 <- user()
beta_m4 <- user()

beta_h[] <- user()  
age_rate[] <- user()

# initial conditions
S0_all[] <- user()
I0[,] <- user()
A0[,] <- user()
C0[,] <- user()
S0[,] <- user()
I_ij0[,] <- user()
A_ij0[,] <- user()
H0[] <- user()
R0[] <- user()
S_m0 <- user()
E_m0[] <- user()
I_m0[] <- user()
S0_v_all[] <- user()
I_v0[,] <- user()
A_v0[,] <- user()
C_v0[,] <- user()
S_v0[,] <- user()
I_v_ij0[,] <- user()
A_v_ij0[,] <- user()
H_v0[] <- user()
R_v0[] <- user()
total_pop_age0[] <- user()

## parse initial conditions
initial(S_all[]) <- S0_all[i]
initial(I[,]) <- I0[i,j]
initial(A[,]) <- A0[i,j]
initial(C[,]) <- C0[i,j]
initial(S[,]) <- S0[i,j]

initial(I_ij[,]) <- I_ij0[i,j]
initial(A_ij[,]) <- A_ij0[i,j]
initial(H[]) <- H0[i]
initial(R[]) <- R0[i]
initial(S_m) <- S_m0 
initial(E_m[]) <- E_m0[i]
initial(I_m[]) <- I_m0[i]
initial(S_v_all[]) <- S0_v_all[i]
initial(I_v[,]) <- I_v0[i,j]
initial(A_v[,]) <- A_v0[i,j]
initial(C_v[,]) <- C_v0[i,j]
initial(S_v[,]) <- S_v0[i,j]
initial(I_v_ij[,]) <- I_v_ij0[i,j]
initial(A_v_ij[,]) <- A_v_ij0[i,j]
initial(H_v[]) <- H_v0[i]
initial(R_v[]) <- R_v0[i]
initial(total_pop_age[]) <- total_pop_age0[i]


# DECLARE ALL THE DIMENSIONS OF VECTORS USED AS USER_DEFINED PARAMETERS

# time dependent parameter for interpolation
dim(tt) <- user()
dim(b_tt) <- user()
dim(vac_rate_tt) <- user()

# initial conditions
dim(S0_all) <- n_age
dim(I0) <- c(n_age, n_sero)
dim(A0) <- c(n_age, n_sero)
dim(C0) <- c(n_age, n_sero)
dim(S0) <- c(n_age, n_sero)
dim(I_ij0) <- c(n_age, n_sero)
dim(A_ij0) <- c(n_age, n_sero)

dim(H0) <- n_age
dim(R0) <- n_age
# dim(S_m0) <- 1
dim(E_m0) <- n_sero
dim(I_m0) <- n_sero
dim(S0_v_all) <- n_age
dim(I_v0) <- c(n_age, n_sero)
dim(A_v0) <- c(n_age, n_sero)
dim(C_v0) <- c(n_age, n_sero)
dim(S_v0) <- c(n_age, n_sero)
dim(I_v_ij0) <- c(n_age, n_sero)
dim(A_v_ij0) <- c(n_age, n_sero)


dim(H_v0) <- n_age
dim(R_v0) <- n_age
dim(total_pop_age0) <- n_age


# state variables
dim(S_all) <- n_age
dim(I) <- c(n_age, n_sero)
dim(A) <- c(n_age, n_sero)
dim(C) <- c(n_age, n_sero)
dim(S) <- c(n_age, n_sero)
dim(I_ij) <- c(n_age, n_sero)
dim(A_ij) <- c(n_age, n_sero)

dim(H) <- n_age
dim(R) <- n_age
# dim(S_m) <- 1
dim(E_m) <- n_sero
dim(I_m) <- n_sero
dim(S_v_all) <- n_age
dim(I_v) <- c(n_age, n_sero)
dim(A_v) <- c(n_age, n_sero)
dim(C_v) <- c(n_age, n_sero)
dim(S_v) <- c(n_age, n_sero)
dim(I_v_ij) <- c(n_age, n_sero)
dim(A_v_ij) <- c(n_age, n_sero)

dim(H_v) <- n_age
dim(R_v) <- n_age
dim(total_pop_age) <- n_age

dim(N_total_age) <- n_age



# parameters
dim(beta_m) <- n_sero
dim(beta_h) <- n_sero
dim(effi_inf_sero_n) <- n_age
dim(effi_inf_sero_p) <- n_age
dim(effi_symp_sero_n) <- c(n_age, n_sero) 
dim(effi_symp_sero_p) <- c(n_age, n_sero) 
dim(effi_hos_sero_n) <- c(n_age, n_sero) 
dim(effi_hos_sero_p) <- c(n_age, n_sero) 
dim(age_rate) <- n_age
dim(death_hum) <- n_age
dim(vac_switch) <- n_age
dim(lambda_hum_tt) <- user()
dim(death_hum_tt) <- c(length(tt), n_age)
dim(import_tt) <- c(length(tt), n_sero)
dim(import) <- n_sero


dim(rho_1) <- n_age
dim(rho_2) <- n_age

# extra 
dim(infected_hum) <- n_sero
dim(foi_mos) <- n_sero
dim(foi_hum) <- n_sero
dim(hos_vac_1) <- c(n_age, n_sero)
dim(hos_vac_2) <- c(n_age, n_sero)


## incidence for unvaccinated
dim(inc_primary_symp) <- c(n_age,n_sero)
dim(inc_secondary_symp) <- c(n_age, n_sero)
dim(inc_primary_inf) <- c(n_age,n_sero)
dim(inc_secondary_inf) <- c(n_age, n_sero)
dim(inc_hosp) <- n_age

## incidence for vaccinated
dim(inc_primary_symp_v) <- c(n_age,n_sero)
dim(inc_secondary_symp_v) <- c(n_age, n_sero)
dim(inc_primary_inf_v) <- c(n_age,n_sero)
dim(inc_secondary_inf_v) <- c(n_age, n_sero)
dim(inc_hosp_v) <- n_age

# 
# dim(age_symp_inf_pri) <- n_age 
# dim(age_symp_inf_sec) <- n_age
# dim(age_sero_symp_inf_sec) <- c(n_age, n_sero)

dim(S_sec) <- n_age
dim(S_sec_v) <- n_age

dim(mos_foi) <- n_sero
dim(hum_foi) <- n_sero
dim(N_age) <- n_age


