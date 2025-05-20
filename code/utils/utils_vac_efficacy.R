

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
  
  ## for infection (not sero-stratified)
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
  
  
  
  ### this is only for creating baseline scenario with efficacy against infection  
  
  #   samp_effi_inf_n <- sample_from_skew_normal(pt = ve_inf_sero_n_pt, 
  #                                              low_ci = ve_inf_sero_n_low,
  #                                              up_ci = ve_inf_sero_n_up,
  #                                              n_sample = 1)
  #   
  #   if (ve_inf_sero_n_pt == 0 & ve_inf_sero_n_low == 0 & ve_inf_sero_n_up == 0 ) {
  #     
  #     samp_effi_inf_n = 0 } else { 
  #       
  #       samp_effi_inf_n <- samp_effi_inf_n
  # }
  # 
  #   
  #   samp_effi_inf_p <- sample_from_skew_normal(pt = ve_inf_sero_p_pt, 
  #                                              low_ci = ve_inf_sero_p_low,
  #                                              up_ci = ve_inf_sero_p_up,
  #                                              n_sample = 1)
  #   
  #   if (ve_inf_sero_p_pt == 0 & ve_inf_sero_p_low == 0 & ve_inf_sero_p_up == 0 ) {
  #     
  #     samp_effi_inf_p = 0 } else { 
  #       
  #       samp_effi_inf_p <- samp_effi_inf_p
  #     }
  #   
  #   
  
  
  ### for VCD and Hospitalized
  for (i in 1:n_sero){
    
    samp_effi_symp_n[i] <- sample_from_skew_normal(pt = symp_pt_sero_n[i],
                                                   low_ci = symp_low_sero_n[i],
                                                   up_ci = symp_up_sero_n[i],
                                                   n_sample = 1)
    
    # samp_effi_symp_n[i] <- runif(n=1, min = symp_low_sero_n[i], max =symp_up_sero_n[i] )
    
    # 
    samp_effi_symp_p[i] <- sample_from_skew_normal(pt = symp_pt_sero_p[i],
                                                   low_ci = symp_low_sero_p[i],
                                                   up_ci = symp_up_sero_p[i],
                                                   n_sample = 1)
    
    
    
    # samp_effi_symp_p[i] <- runif(n=1, min = symp_low_sero_p[i], max =symp_up_sero_p[i] )
    
    
    ## Missing estimates are replaced by corresponding VE against VCD
    if (i == 2 | i == 4) {
      
      samp_effi_hosp_n[i] <- samp_effi_symp_n[i]
      
    } else{ samp_effi_hosp_n[i] <- sample_from_skew_normal(pt = hosp_pt_sero_n[i],
                                                           low_ci = hosp_low_sero_n[i],
                                                           up_ci = hosp_up_sero_n[i],
                                                           n_sample = 1)
    }
    
    
    # if (i == 2 | i == 4) {
    #   
    #   samp_effi_hosp_n[i] <- samp_effi_symp_n[i]
    #   
    # } else{ samp_effi_hosp_n[i] <- runif(n=1, min = hosp_low_sero_n[i], max = hosp_up_sero_n[i] )
    # 
    # }
    # 
    
    ## Missing estimates are replaced by corresponding VE against VCD
    if (i == 4) {
      
      samp_effi_hosp_p[i] = samp_effi_symp_p[i]
      
    } else{ samp_effi_hosp_p[i] <- sample_from_skew_normal(pt = hosp_pt_sero_p[i],
                                                           low_ci = hosp_low_sero_p[i],
                                                           up_ci = hosp_up_sero_p[i],
                                                           n_sample = 1)
    }
    
    
    # if (i == 4) {
    #   
    #   samp_effi_hosp_p[i] = samp_effi_symp_p[i]
    #   
    # } else{ samp_effi_hosp_p[i] <- runif(n=1, min = hosp_low_sero_p[i], max = hosp_up_sero_p[i] )
    # 
    # }
    
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
  
  sample_sn <- rsn(n = n_sample , xi = pt, omega = sigma, alpha = alpha)
  
  final_sample <- if (sample_sn > up_ci) up_ci else if (sample_sn < low_ci) low_ci else sample_sn 
  
  
  return(final_sample)
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
  
  
  ## n_year_waing == no. of years until it will stop wane then extra baseline VE
  ## that's why its n_year_waning +1
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

