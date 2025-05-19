



# total population
N <- sum(S_all[]) + sum(I[,]) + sum(C[,]) + sum(S[,]) + sum(I_ij[,]) + sum(R[])


infected_hum[] <- sum(I[,i]) + sum(I_ij[,i]) 


## population each age group
N_age[] <- S_all[i] + sum(I[i,]) + sum(C[i,]) + sum(S[i,]) + sum(I_ij[i,]) + R[i]


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


## age-specific symptomatic infection
age_symp_inf_pri[] <- rho_1[i]*sum(foi_mos[])*S_all[i] 
age_sero_symp_inf_sec[,] <- foi_mos[j]*(S_sec[i] - S[i,j]) #(sum(foi_mos[]) - foi_mos[j])*S[i,j]  
age_symp_inf_sec[] <-  rho_2[i]*sum(age_sero_symp_inf_sec[i,])


# unvaccinated primary infection
deriv(S_all[1]) <- lambda_hum*N - sum(foi_mos[])*S_all[i] - death_hum*S_all[i] -
  age_rate[i]*S_all[i] 

deriv(S_all[2:n_age]) <- -sum(foi_mos[])*S_all[i] - death_hum*S_all[i] +
  age_rate[i-1]*S_all[i-1] - age_rate[i]*S_all[i] 

deriv(I[1,1:n_sero]) <- foi_mos[j]*S_all[i] - (gamma_1+death_hum)*I[i,j] - 
  age_rate[i]*I[i,j]

deriv(I[2:n_age,1:n_sero]) <- foi_mos[j]*S_all[i] - (gamma_1+death_hum)*I[i,j] + 
  age_rate[i-1]*I[i-1,j] - age_rate[i]*I[i,j]

deriv(C[1,1:n_sero]) <- gamma_1*I[i,j] - (alpha + death_hum)*C[i,j] - 
  age_rate[i]*C[i,j] 

deriv(C[2:n_age,1:n_sero]) <- gamma_1*I[i,j] - (alpha + death_hum)*C[i,j] + 
  age_rate[i-1]*C[i-1,j] - age_rate[i]*C[i,j] 

deriv(S[1,1:n_sero ]) <- alpha*C[i,j] - (sum(foi_mos[]) - foi_mos[j])*S[i,j] - death_hum*S[i,j] - 
  age_rate[i]*S[i,j] 

deriv(S[2:n_age,1:n_sero ]) <- alpha*C[i,j] -  (sum(foi_mos[]) - foi_mos[j])*S[i,j] - death_hum*S[i,j] + 
  age_rate[i-1]*S[i-1,j] - age_rate[i]*S[i,j] 

# secondary infection
deriv(I_ij[1,1:n_sero]) <- foi_mos[j]*(S_sec[i] - S[i,j]) - (gamma_2+death_hum)*I_ij[i,j] - 
  age_rate[i]*I_ij[i,j]

deriv(I_ij[2:n_age,1:n_sero]) <- foi_mos[j]*(S_sec[i] - S[i,j]) - (gamma_2+death_hum)*I_ij[i,j] + 
  age_rate[i-1]*I_ij[i-1,j] - age_rate[i]*I_ij[i,j]


deriv(R[1]) <-  gamma_2*sum(I_ij[i,]) - death_hum*R[i] - 
  age_rate[i]*R[i]

deriv(R[2:n_age]) <-  gamma_2*sum(I_ij[i,])  - death_hum*R[i] + 
  age_rate[i-1]*R[i-1] - age_rate[i]*R[i] 

# mosquito population

deriv(S_m) <- lambda_mos*(S_m + sum(E_m[])+sum(I_m[])) - S_m*sum(foi_hum[]) - death_mos*S_m

deriv(E_m[1:n_sero]) <-  S_m*foi_hum[i] - (sigma_m + death_mos)*E_m[i]

deriv(I_m[1:n_sero]) <- sigma_m*E_m[i] - death_mos*I_m[i]

# total population in age groups
deriv(total_pop_age[1]) <- lambda_hum*N - death_hum*total_pop_age[i] - age_rate[i]*total_pop_age[i]

deriv(total_pop_age[2:n_age]) <- age_rate[i-1]*total_pop_age[i-1] - age_rate[i]*total_pop_age[i] - death_hum*total_pop_age[i]



deriv(inc_symp[1:n_data_year, 1:n_age]) <- if (as.integer(t) >= lo[i] && as.integer(t) <= u[i] ) (age_symp_inf_pri[j] + age_symp_inf_sec[j]) else 0


anu_inc_symp[1:n_data_year, 1:n_age] <- if (as.integer(t) >= lo[i] && as.integer(t) <= u[i] ) inc_symp[i,j] else 0

inc_symp_p_100k[1:n_age] <- (sum(anu_inc_symp[,i]))


## for absolute cases
# output(inc[1:n_data_age]) <- sum(inc_symp_p_100k[a_l[i]:a_u[i]]) 

## for incidence per 100k
output(inc[1:n_data_age]) <- 100000*sum(inc_symp_p_100k[a_l[i]:a_u[i]])/sum(N_age[a_l[i]:a_u[i]])



### this is extra not required for the model fitting with incidence of annual cases per 100,000
## but with cumulative one
deriv(agg_inc_symp[1:n_age]) <- (age_symp_inf_pri[i] + age_symp_inf_sec[i]) 

output(agg_inc[1:n_data_age]) <- 100000*sum(agg_inc_symp[a_l[i]:a_u[i]])/sum(N_age[a_l[i]:a_u[i]])


output(mos_foi[]) <- foi_mos[i]

output(hum_foi[]) <- foi_hum[i]

output(N_total_age[]) <- N_age[i]


## Parameters

## time dependent
lambda_hum <- lambda_hum_tt[as.integer(t)]  
death_hum <-  death_hum_tt[as.integer(t)] 
b <-  user()

lo[] <- l_t[i]
u[] <- u_t[i]
a_u[] <- a_up[i]
a_l[] <- a_lo[i]


# parameters

ratio_p_s <- user()

rho_1[] <- ratio_p_s*rho_2[i]   #user()


rho_2[1:5] <- rho_2_a1
rho_2[6:15] <- rho_2_a2
rho_2[16:25] <- rho_2_a3
rho_2[26:35] <- rho_2_a4
rho_2[36:45] <- rho_2_a5
rho_2[46:55] <- rho_2_a6
rho_2[56:65] <- rho_2_a7
rho_2[66:91] <- rho_2_a8


rho_2_a1 <- user()
rho_2_a2 <- user()
rho_2_a3 <- user()
rho_2_a4 <- user()
rho_2_a5 <- user()
rho_2_a6 <- user()
rho_2_a7 <- user()
rho_2_a8 <- user()


sigma_m <- user()

# disease progression related
alpha <- user()
gamma_1 <- user()
gamma_2 <- user()



# demography related
lambda_mos <- user()
death_mos <- user()


# LIST OF USER DEFINED PARAMETERS
lambda_hum_tt[] <- user()
death_hum_tt[] <- user()
# b_tt[] <- user()
l_t[] <- user()
u_t[] <- user()
a_up[] <- user()
a_lo[] <- user()
n_data_year <- user()
n_data_age <- user()


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
C0[,] <- user()
S0[,] <- user()
I_ij0[,] <- user()
R0[] <- user()
S_m0 <- user()
E_m0[] <- user()
I_m0[] <- user()
total_pop_age0[] <- user()
inc_symp0[,] <- user()
agg_inc_symp0[] <- user()




## parse initial conditions
initial(S_all[]) <- S0_all[i]
initial(I[,]) <- I0[i,j]
initial(C[,]) <- C0[i,j]
initial(S[,]) <- S0[i,j]
initial(I_ij[,]) <- I_ij0[i,j]
initial(R[]) <- R0[i]
initial(S_m) <- S_m0 
initial(E_m[]) <- E_m0[i]
initial(I_m[]) <- I_m0[i]
initial(total_pop_age[]) <- total_pop_age0[i]
initial(inc_symp[,]) <- inc_symp0[i,j]
initial(agg_inc_symp[]) <- agg_inc_symp0[i]


# DECLARE ALL THE DIMENSIONS OF VECTORS USED AS USER_DEFINED PARAMETERS

# initial conditions
dim(S0_all) <- n_age
dim(I0) <- c(n_age, n_sero)
# dim(A0) <- c(n_age, n_sero)
dim(C0) <- c(n_age, n_sero)
dim(S0) <- c(n_age, n_sero)
dim(I_ij0) <- c(n_age, n_sero)
dim(R0) <- n_age
dim(E_m0) <- n_sero
dim(I_m0) <- n_sero
dim(total_pop_age0) <- n_age
dim(inc_symp0) <- c(n_data_year, n_age)
dim(agg_inc_symp0) <- n_age

# state variables
dim(S_all) <- n_age
dim(I) <- c(n_age, n_sero)
dim(C) <- c(n_age, n_sero)
dim(S) <- c(n_age, n_sero)
dim(I_ij) <- c(n_age, n_sero)

# dim(H) <- n_age
dim(R) <- n_age
dim(E_m) <- n_sero
dim(I_m) <- n_sero
dim(total_pop_age) <- n_age
dim(inc_symp) <- c(n_data_year, n_age)
dim(agg_inc_symp) <- n_age
dim(N_total_age) <- n_age



# parameters
dim(beta_m) <- n_sero
dim(beta_h) <- n_sero
dim(age_rate) <- n_age
dim(rho_1) <- n_age
dim(rho_2) <- n_age

dim(lambda_hum_tt) <- user() 
dim(death_hum_tt) <- user() 
dim(lo) <- n_data_year
dim(u) <- n_data_year
dim(l_t) <- user()
dim(u_t) <- user()
dim(a_u) <- n_data_age
dim(a_l) <- n_data_age
dim(a_up) <- user()
dim(a_lo) <- user()

# extra 
dim(infected_hum) <- n_sero
dim(foi_mos) <- n_sero
dim(foi_hum) <- n_sero
dim(age_symp_inf_pri) <- n_age 
dim(age_symp_inf_sec) <- n_age
dim(age_sero_symp_inf_sec) <- c(n_age, n_sero)
dim(S_sec) <- n_age
dim(anu_inc_symp) <- c(n_data_year, n_age)
dim(inc_symp_p_100k) <- n_age
dim(inc) <- n_data_age
dim(agg_inc) <- n_data_age

dim(mos_foi) <- n_sero
dim(hum_foi) <- n_sero
dim(N_age) <- n_age


