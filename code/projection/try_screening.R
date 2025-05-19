library(odin.dust)

simple_vac_sir <- odin.dust::odin_dust({
  
  N <- S + I + R + S_v + I_v + R_v
  
  total_vac <- S_v+I_v+R_v
  
  # vac_s1 <- if (as.integer(t)==vac_start) vac_coverage*N*S/(S+R) else 
  #   lambda_hum*(total_vac)*S/(S+R)
  # vac_r1 <- if (as.integer(t)==vac_start) vac_coverage*N*R/(S+R) else 
  #   lambda_hum*(total_vac)*R/(S+R)
  
  
  vac_s1 <- if (as.integer(t) >= vac_start) (1/(365))*(1-test_specificity)*vac_coverage*S*N/(((1-test_specificity)*S)+(test_sensitivity*R)) else 0
    
  vac_r1 <- if (as.integer(t)>=vac_start) (1/(365))*test_sensitivity*vac_coverage*R*N/((1-test_specificity)*S+test_sensitivity*R) else 0
  
  
  
  vac_s2 <- if (as.integer(t) >= vac_start) (1/(365))*vac_coverage*S*N/(S+R) else 0 
  #((1-test_specificity)*S*N/(((1-test_specificity)*S)+(test_sensitivity*R)) + test_sensitivity*R*N/(((1-test_specificity)*S)+(test_sensitivity*R))) else 0 
  
  vac_r2 <- if (as.integer(t)>=vac_start) (1/(365))*vac_coverage*R*N/(S+R) else 0 #((1-test_specificity)*S*N/(((1-test_specificity)*S)+(test_sensitivity*R)) + test_sensitivity*R*N/(((1-test_specificity)*S)+(test_sensitivity*R))) else 0 
  
    
  
  # vac_s2 <- if (as.integer(t)==vac_start) 1 else 1
  #   
  # vac_r2 <- if (as.integer(t)==vac_start) 1 else 1
  #  
  
  
  ############### Un-vaccinated ##################
  deriv(S) <- lambda_hum*N - beta*S*(I+I_v)/N - death_hum*S - vac_s
  deriv(I) <- beta*S*(I+I_v)/N- gamma*I - death_hum*I
  deriv(R) <- gamma*I - death_hum*R - vac_r
  
  
  ################# Vaccinated ###################
  deriv(S_v) <- vac_s  - (1 - v_inf)*beta*S_v*(I+I_v)/N - death_hum*S_v
  deriv(I_v) <- (1 - v_inf)*beta*S_v*(I+I_v)/N - gamma*I_v - death_hum*I_v
  deriv(R_v) <- vac_r + gamma*I_v - death_hum*R_v 
  
  output(n_vac) <- (S_v+I_v+R_v)
  output(total) <- N
  output(new_vac) <- vac_s + vac_r

    
  vac_s <- if (prescreening == 1 && prescreening_non_sero == 0) vac_s1 else
      if (prescreening == 1 && prescreening_non_sero == 1) vac_s2 else 
        if (prescreening == 0) 100 else 0
  
  vac_r <- if (prescreening == 1 && prescreening_non_sero == 0) vac_r1 else
    if (prescreening == 1 && prescreening_non_sero == 1) vac_r2 else 
      if (prescreening == 0) 100 else 0
  

  
  ## Initial conditions
  initial(S) <- S0
  initial(I) <- I0
  initial(R) <- R0
  initial(S_v) <- S_v0
  initial(I_v) <- I_v0
  initial(R_v) <- R_v0
  
 prescreening <- user()
 prescreening_non_sero <- user()
 test_specificity <- user()
 test_sensitivity <- user()
 
  
  ## parameters & Initial conditions
  beta  <- user()
  gamma <- user()
  lambda_hum <- user()
  death_hum <- user() 
  vac_coverage <- user()
  v_inf <- user()
  vac_start <- user()
  
  I0 <- user()
  R0 <- user()
  S_v0 <- user()
  I_v0 <- user()
  R_v0 <- user()
  S0 <- user()
  
  
})



calc_yearly <- function(mat){
  tapply(mat, (seq_along(mat) - 1) %/%365, sum)
}

calc_yearly_wo_sum <- function(mat){
  # tapply(mat, (seq_along(mat) - 1) %/%365, sum)
  mat[seq(1,length(t),by=365)]
}


#########
prescreening_non_sero = 1
prescreening = 1
test_specificity = 0.9
test_sensitivity = 0.9




final_time = 10*365
t <- seq(1,final_time, by = 1)  


lambda_hum = (1/(83*365))
death_hum = (1/(83*365)) 
beta  <- 0.6
gamma <- 0.1
vac_coverage <- 0.7
v_inf <- 0.5
vac_start <- 100


total_pop <- 5400000
I0 <- 1
R0 <- 3000000

I_v0 <- 0
R_v0 <- 0 #vac_coverage*total_pop/2

S_v0 <- 0 #vac_coverage*total_pop/2

S0 <- total_pop - (I0 + R0+S_v0+R_v0)

## parameters & Initial conditions
pars <- list(vac_coverage = vac_coverage,
             vac_start = vac_start,
             lambda_hum = lambda_hum,
             death_hum = death_hum, 
             beta = beta, 
             v_inf = v_inf,
             gamma = gamma,
             S0 = S0,
             I0 = I0,
             R0 = R0,
             S_v0 = S_v0,
             I_v0 = I_v0,
             R_v0 = R_v0,
             prescreening_non_sero = prescreening_non_sero,
             prescreening = prescreening,
             test_sensitivity = test_sensitivity,
             test_specificity = test_specificity
             
)

model <- simple_vac_sir$new(pars, time = 1, n_particles = 1,
                            ode_control = list(max_steps = 3000000,
                                               step_size_min = 1e-11))

out <- model$simulate(t)




ind_S <- model$info()$index$S
ind_I = model$info()$index$I
ind_R = model$info()$index$R
ind_Sv <- model$info()$index$S_v
ind_Iv = model$info()$index$I_v
ind_Rv <- model$info()$index$R_v
ind_nvac = model$info()$index$n_vac
ind_total <- model$info()$index$total# 
ind_newvac <- model$info()$index$new_vac


plot((out[ind_newvac,,]), type = "l", col="red")

# plot((out[ind_S,,]), type = "l", col="red", ylim=c(0,total_pop))+
#   lines((out[ind_I,,]), type = "l", col="blue")+
#   lines((out[ind_R,,]), type = "l", col="green")
# 
# plot((out[ind_Sv,,]), type = "l", col="red", ylim=c(0, total_pop))+
#   lines((out[ind_Iv,,]), type = "l", col="blue")+
#   lines((out[ind_Rv,,]), type = "l", col="green")
# 
# 
# plot((out[ind_nvac,,]), type = "l")+
#   lines((out[ind_total,,]))
# # 
# # 
# plot(calc_yearly_wo_sum(out[ind_nvac,,]), ylim=c(0.3*total_pop,0.5*total_pop))
# # 
# # plot(calc_yearly_wo_sum(out[ind_total,,]), type="l")
# 
# plot(100*calc_yearly_wo_sum(out[ind_nvac,,])/calc_yearly_wo_sum(out[ind_total,,]), type = "l", ylim=c(0,100))
