

make_transform <- function(n_age, n_sero, tt,b0_tt,
                           beta_h, age_rate, death_hum,
                           I0,A0,C0,S0,I_ij0,A_ij0,H0,R0,S0_all,
                           S_m0, E_m0,I_m0,total_pop_age0, C_age0, S_age0,
                           vac_rate_tt,vac_switch,kai,kai_v,rho_1,rho_2,rho_v_1,rho_v_2,               
                           sigma_m,alpha,alpha_v,gamma_1,gamma_2,gamma_v_1,gamma_v_2,             
                           xi_1, xi_2,xi_v_1,xi_v_2,lambda_hum,lambda_mos,death_mos,            
                           S0_v_all,I_v0,A_v0,C_v0,S_v0,I_v_ij0,A_v_ij0,H_v0,R_v0,
                           effi_inf_sero_n,effi_inf_sero_p,effi_symp_sero_n,effi_symp_sero_p,
                           effi_hos_sero_n,effi_hos_sero_p) {
  list(
    n_age=n_age,
    n_sero=n_sero,
    
    tt=tt,
    beta_h=beta_h,
    b0_tt=b0_tt,
    vac_rate_tt=vac_rate_tt,
    vac_switch=vac_switch,
    
    kai=kai,                     
    kai_v=kai_v,                   
    rho_1=rho_1,               
    rho_2=rho_2,                
    rho_v_1=rho_v_1,              
    rho_v_2=rho_v_2,               
    sigma_m=sigma_m,              
    alpha=alpha,                
    alpha_v= alpha_v,             
    gamma_1 =gamma_1,              
    gamma_2 =gamma_2,              
    gamma_v_1=gamma_v_1,             
    gamma_v_2=gamma_v_2,             
    xi_1 =xi_1,         
    xi_2 =xi_2,              
    xi_v_1 =xi_v_1,            
    xi_v_2 =xi_v_2,      
    
    lambda_hum=lambda_hum,   
    lambda_mos=lambda_mos,          
    death_mos=death_mos,            
    
    age_rate=age_rate,        # the last element of age_rate should be zero
    death_hum=death_hum,      # TODO: we can use the function "demography_fitting" function to estimate the age-specific mortality
    I0 = I0,
    A0 = A0,
    C0 = C0,
    S0 = S0,
    I_ij0 = I_ij0,
    A_ij0 = A_ij0,
    H0 = H0,
    R0 = R0,
    S0_all = S0_all,
    S_m0 = S_m0,
    E_m0 = E_m0,
    I_m0 = I_m0,
    total_pop_age0 =total_pop_age0,
    C_age0=C_age0,
    S_age0=S_age0,
    
    S0_v_all = S0_v_all,
    I_v0 =I_v0,
    A_v0 = A_v0,
    C_v0 = C_v0,
    S_v0 = S_v0,
    I_v_ij0 =  I_v_ij0,
    A_v_ij0 = A_v_ij0,
    H_v0 = H_v0,
    R_v0 = R_v0,
    
    effi_inf_sero_n=effi_inf_sero_n,
    effi_inf_sero_p=effi_inf_sero_p,
    effi_symp_sero_n=effi_symp_sero_n,
    effi_symp_sero_p=effi_symp_sero_p,
    effi_hos_sero_n=effi_hos_sero_n,
    effi_hos_sero_p=effi_hos_sero_p)
  
  function(theta) {
    list(
      n_age=n_age,
      n_sero=n_sero,
      
      tt=tt,
      beta_h=beta_h,
      b0_tt=b0_tt,
      vac_rate_tt=vac_rate_tt,
      vac_switch=vac_switch,
      
      kai=kai,                     
      kai_v=kai_v,                   
      rho_1=rho_1,               
      rho_2=rho_2,                
      rho_v_1=rho_v_1,              
      rho_v_2=rho_v_2,               
      sigma_m=sigma_m,              
      alpha=alpha,                
      alpha_v= alpha_v,             
      gamma_1 =gamma_1,              
      gamma_2 =gamma_2,              
      gamma_v_1=gamma_v_1,             
      gamma_v_2=gamma_v_2,             
      xi_1 =xi_1,         
      xi_2 =xi_2,              
      xi_v_1 =xi_v_1,            
      xi_v_2 =xi_v_2,      
      
      lambda_hum=lambda_hum,   
      lambda_mos=lambda_mos,          
      death_mos=death_mos,            
      
      age_rate=age_rate,        # the last element of age_rate should be zero
      death_hum=death_hum,      # TODO: we can use the function "demography_fitting" function to estimate the age-specific mortality
      I0 = I0,
      A0 = A0,
      C0 = C0,
      S0 = S0,
      I_ij0 = I_ij0,
      A_ij0 = A_ij0,
      H0 = H0,
      R0 = R0,
      S0_all = S0_all,
      S_m0 = S_m0,
      E_m0 = E_m0,
      I_m0 = I_m0,
      total_pop_age0 =total_pop_age0,
      C_age0=C_age0,
      S_age0=S_age0,
      
      S0_v_all = S0_v_all,
      I_v0 =I_v0,
      A_v0 = A_v0,
      C_v0 = C_v0,
      S_v0 = S_v0,
      I_v_ij0 =  I_v_ij0,
      A_v_ij0 = A_v_ij0,
      H_v0 = H_v0,
      R_v0 = R_v0,
      
      effi_inf_sero_n=effi_inf_sero_n,
      effi_inf_sero_p=effi_inf_sero_p,
      effi_symp_sero_n=effi_symp_sero_n,
      effi_symp_sero_p=effi_symp_sero_p,
      effi_hos_sero_n=effi_hos_sero_n,
      effi_hos_sero_p=effi_hos_sero_p,
      beta_m1=theta[["beta_m1"]],
      beta_m2=theta[["beta_m2"]],
      beta_m3=theta[["beta_m3"]],
      beta_m4=theta[["beta_m4"]])
  }
}




make_transform_incidence <- function(n_age, n_sero, tt,b0_tt,
                           beta_h, age_rate, death_hum,
                           I0,A0,C0,S0,I_ij0,A_ij0,H0,R0,S0_all,
                           S_m0, E_m0,I_m0,cum_inc0,
                           vac_rate_tt,vac_switch,kai,kai_v,rho_1,rho_2,rho_v_1,rho_v_2,               
                           sigma_m,alpha,alpha_v,gamma_1,gamma_2,gamma_v_1,gamma_v_2,             
                           xi_1, xi_2,xi_v_1,xi_v_2,lambda_hum,lambda_mos,death_mos,            
                           S0_v_all,I_v0,A_v0,C_v0,S_v0,I_v_ij0,A_v_ij0,H_v0,R_v0,
                           effi_inf_sero_n,effi_inf_sero_p,effi_symp_sero_n,effi_symp_sero_p,
                           effi_hos_sero_n,effi_hos_sero_p) {
  list(
    n_age=n_age,
    n_sero=n_sero,
    
    tt=tt,
    beta_h=beta_h,
    b0_tt=b0_tt,
    vac_rate_tt=vac_rate_tt,
    vac_switch=vac_switch,
    
    kai=kai,                     
    kai_v=kai_v,                   
    rho_1=rho_1,               
    rho_2=rho_2,                
    rho_v_1=rho_v_1,              
    rho_v_2=rho_v_2,               
    sigma_m=sigma_m,              
    alpha=alpha,                
    alpha_v= alpha_v,             
    gamma_1 =gamma_1,              
    gamma_2 =gamma_2,              
    gamma_v_1=gamma_v_1,             
    gamma_v_2=gamma_v_2,             
    xi_1 =xi_1,         
    xi_2 =xi_2,              
    xi_v_1 =xi_v_1,            
    xi_v_2 =xi_v_2,      
    
    lambda_hum=lambda_hum,   
    lambda_mos=lambda_mos,          
    death_mos=death_mos,            
    
    age_rate=age_rate,        # the last element of age_rate should be zero
    death_hum=death_hum,      # TODO: we can use the function "demography_fitting" function to estimate the age-specific mortality
    I0 = I0,
    A0 = A0,
    C0 = C0,
    S0 = S0,
    I_ij0 = I_ij0,
    A_ij0 = A_ij0,
    H0 = H0,
    R0 = R0,
    S0_all = S0_all,
    S_m0 = S_m0,
    E_m0 = E_m0,
    I_m0 = I_m0,
    cum_inc0=cum_inc0,
    
    S0_v_all = S0_v_all,
    I_v0 =I_v0,
    A_v0 = A_v0,
    C_v0 = C_v0,
    S_v0 = S_v0,
    I_v_ij0 =  I_v_ij0,
    A_v_ij0 = A_v_ij0,
    H_v0 = H_v0,
    R_v0 = R_v0,
    
    effi_inf_sero_n=effi_inf_sero_n,
    effi_inf_sero_p=effi_inf_sero_p,
    effi_symp_sero_n=effi_symp_sero_n,
    effi_symp_sero_p=effi_symp_sero_p,
    effi_hos_sero_n=effi_hos_sero_n,
    effi_hos_sero_p=effi_hos_sero_p)
  
  function(theta) {
    list(
      n_age=n_age,
      n_sero=n_sero,
      
      tt=tt,
      beta_h=beta_h,
      b0_tt=b0_tt,
      vac_rate_tt=vac_rate_tt,
      vac_switch=vac_switch,
      
      kai=kai,                     
      kai_v=kai_v,                   
      rho_1=rho_1,               
      rho_2=rho_2,                
      rho_v_1=rho_v_1,              
      rho_v_2=rho_v_2,               
      sigma_m=sigma_m,              
      alpha=alpha,                
      alpha_v= alpha_v,             
      gamma_1 =gamma_1,              
      gamma_2 =gamma_2,              
      gamma_v_1=gamma_v_1,             
      gamma_v_2=gamma_v_2,             
      xi_1 =xi_1,         
      xi_2 =xi_2,              
      xi_v_1 =xi_v_1,            
      xi_v_2 =xi_v_2,      
      
      lambda_hum=lambda_hum,   
      lambda_mos=lambda_mos,          
      death_mos=death_mos,            
      
      age_rate=age_rate,        # the last element of age_rate should be zero
      death_hum=death_hum,      # TODO: we can use the function "demography_fitting" function to estimate the age-specific mortality
      I0 = I0,
      A0 = A0,
      C0 = C0,
      S0 = S0,
      I_ij0 = I_ij0,
      A_ij0 = A_ij0,
      H0 = H0,
      R0 = R0,
      S0_all = S0_all,
      S_m0 = S_m0,
      E_m0 = E_m0,
      I_m0 = I_m0,
     cum_inc0=cum_inc0,
      
      S0_v_all = S0_v_all,
      I_v0 =I_v0,
      A_v0 = A_v0,
      C_v0 = C_v0,
      S_v0 = S_v0,
      I_v_ij0 =  I_v_ij0,
      A_v_ij0 = A_v_ij0,
      H_v0 = H_v0,
      R_v0 = R_v0,
      
      effi_inf_sero_n=effi_inf_sero_n,
      effi_inf_sero_p=effi_inf_sero_p,
      effi_symp_sero_n=effi_symp_sero_n,
      effi_symp_sero_p=effi_symp_sero_p,
      effi_hos_sero_n=effi_hos_sero_n,
      effi_hos_sero_p=effi_hos_sero_p,
      beta_m1=theta[["beta_m1"]],
      beta_m2=theta[["beta_m2"]],
      beta_m3=theta[["beta_m3"]],
      beta_m4=theta[["beta_m4"]])
  }
}

fit_sigmoid <- function(x, y) {
  
  sigmoidal <- function(x, A, B) {
    1/ (1 + A*exp(-B*x))
  }
  
  # Initial parameter guesses
  initial_guess <- c(A = 10, B = 0.3)
  
  # Fit the sigmoidal function to the data
  fit <- nls(y ~ sigmoidal(x, A, B), start = initial_guess)
  
  return(fit)
}







