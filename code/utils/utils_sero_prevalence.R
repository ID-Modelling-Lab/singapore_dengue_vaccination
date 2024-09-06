#### Sero-prevalence data (Reconstruction of sero-prevalence using catalytic model)

get_sero_p <- function(sero_year){
  # Read the data
  data_sero <- read_xlsx( here("data","dengue_sero_prevalence_singapore.xlsx"),
                          sheet = "dengue_sero_p_2013_17" )

  ## select the sero_prevalence data
  # sero_year = 2013
  #### age that we are modelling i.e 0-89y and 90+y
  age = 0:90
  seroprevalence <- get_seroprevalence(data_sero = data_sero,
                                       sero_year = sero_year,
                                       age = age)
  sero_p <- seroprevalence[1,]
  primary <- seroprevalence[3,]

  return(list(se_p = sero_p, pri = primary))
}




#### Sero-prevalence data (Reconstruction of sero-prevalence using catalytic model)

get_seroprevalence <- function(data_sero, sero_year, age){
  
  # extract the data from the year 2013
  data_sero_year = data_sero[data_sero$`year`==sero_year,]
  
  ### include number of negative in the data
  data_sero_year <- data_sero_year %>%
    mutate(
      negative = data_sero_year$total - data_sero_year$positive) 
  
  age_str <- data_sero_year$age
  # Function to extract numeric values from the age ranges
  extract_numeric <- function(age_range) {
    # Split the string based on the "-"
    age_parts <- strsplit(age_range, "-|y")[[1]]
    
    if (length(age_parts) == 2) {
      return(as.numeric(age_parts))
    } else if (length(age_parts) == 1) {
      if (age_parts == "60+") {
        
        ## 2013 survey contains age-groups from 16-71
        ## 2017 survey contains age-groups from 16-74
        
        if (sero_year==2013) {
          return(c(61,71))
        }
        
        if (sero_year==2017) {
          return(c(61,74))
        }
        
      } }
  }
  
  # Apply the function to each element of the array
  numeric_age <- sapply(age_str, extract_numeric)
  
  
  
  ##### Calculate likelihood
  
  
  
  likelihood <- function(theta, data, npar, ngroup, d_npar) {
    
    ## pull sero-neg from the data
    no_seronegative <- data |> 
      
      pull(negative)
    
    ## pull sero-pos from the data
    no_seropositive <- data |> 
      
      pull(positive)
    
    ## for constant foi
    
    if (npar == 1) {
      
      theta_f=rep(theta, d_npar)  
      
      #likelihood for each group
      lik_grp <- array(,dim=ngroup)
      
      ## for first n-1 age-groups
      for (n in 1:(ngroup-1)){
        
        lik_grp[n] = 5*sum(theta_f[1:(n+2)]) + 
          3.5*theta_f[n+3]
      }
      
      ## Last age-group
      lik_grp[ngroup] = 5*sum(theta_f[1:(ngroup+2)]) + 
        6.5*theta_f[ngroup+3]
    }
    
    
    ## for time-varying foi
    if (npar > 1) {
      
      ## aggregated lambdas for each age group
      lik_grp <- array(,dim=ngroup)
      ## for first n-1 age-groups
      for (n in 1:(ngroup-1)){
        lik_grp[n] = 5*sum(theta[1:(n+2)]) + 3.5*theta[n+3]
      }
      ## Last age-group
      lik_grp[ngroup] = 5*sum(theta[1:(ngroup+2)]) + 6.5*theta[ngroup+3]
    }
    
    prop_seronegative <- exp(-4*lik_grp)
    
    prop_seropositive <- 1 - prop_seronegative
    
    ## binomial likelihood
    seronegative_term <- (no_seronegative) * log(prop_seronegative)
    
    seropositive_term <- (no_seropositive) * log(prop_seropositive)
    
    logterms <- seronegative_term + seropositive_term
    
    minuslogl <- -sum(logterms)
    
    return(minuslogl)
  }
  
  
  
  get_fit <- function(data, npar, Hes = F, ...) {
    
    initial_guess <- rep(0.05,npar) # abs(rnorm(npar, 0.02, 0.01))
    
    parscale <- 1e-5 + initial_guess
    
    optimized_values <- optim(
      
      par = initial_guess,
      
      fn = likelihood, 
      
      data = data_sero_year,
      
      npar = npar,
      
      method = "L-BFGS-B",
      
      lower = rep(0.000001, length(initial_guess)),
      
      control = list(trace = T, parscale = parscale, maxit = 1000000),
      
      hessian = Hes,
      
      ...
      
    )
    
    # optimized_values$init <- initial_guess
    
    return(optimized_values$par)
    
  }
  
  
  
  ##### Estimate the lambda (foi)
  
  ### number of lambda (npar) is: upto 60 with 5 year age gap then one extra for the last age group
  npar <-  floor((numeric_age[length(numeric_age)-2]/5)) + 1
  lambda_timev <- get_fit(data_sero_year,
                          npar = npar, ngroup=dim(data_sero_year)[1], d_npar="")
  
  
  ##### Reconstruction of sero_prevalence ($\pi_{a,X}$) age-group for all age-groups 
  ##### Function to get seroprevalence estimates from fit 
  
  get_seroprevalence_reconstruct <- function(theta_f, age) {
    
    ## if age ==0, no. of lambda is 1 
    if (age == 0) {
      n_lambda <- 1
    }
    ## if age >0 and <=60, no. of lambda is ceil(age/5)
    if (age > 0 && age <=60) {
      n_lambda <- ceiling(age/5)
    }
    
    ### if >60, we only have 13 lambdas
    if (age > 60) {
      n_lambda <- length(theta_f)
    }
    # calculate the extra age due to the "mid-age" calculation 
    extra= (age - 5*(n_lambda-1)) + 0.5
    
    foi <- 5*sum(head(theta_f,(n_lambda-1))) + extra*theta_f[n_lambda]
    
    ## sero-negative (0-infection)
    sero_neg <- exp(-4*foi)
    
    ## only 1-infection
    primary <- 4*(1- exp(-foi))*exp(-3*foi)
    
    ## secondary (2-infections)
    secondary <- 1 - sero_neg - primary
    
    ### at least one infection
    seroprevalence<- 1 - exp(-4*foi)
    
    out<- c(
      sero_prevalence=seroprevalence,
      sero_negative=sero_neg,
      primary=primary, 
      secondary=secondary
    )
    
    return(out)
  }
  
  
  ### plot for time varying
  theta_f <- lambda_timev
  
  ## get all the outputs: sero-prev, sero_neg, primary, secondary
  output_tv <- sapply(age, get_seroprevalence_reconstruct, theta_f = theta_f)
  
  return(output_tv)
}

