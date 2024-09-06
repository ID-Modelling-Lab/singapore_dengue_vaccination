### Process population data
get_processed_population_data <- function(pop_data, begin_year, end_year){
  
  pop_data_both= pop_data[1:(which(pop_data$Data_Series=="Singapore Male Residents")-1),]
  
  ### Next since "pop_data_both" contains some extra row like "65 years and over",
  ## we discard these
  ### We will keep the age group from 0-89 (1 year age band) and 90 years and over
  
  rows_to_discard <- c("65 Years & Over", "70 Years & Over", "75 Years & Over",
                       "80 Years & Over", "85 Years & Over")
  pop_data_both <- pop_data_both[!(pop_data_both$Data_Series %in% rows_to_discard), ]
  
  ## save the info, only for the specified time interval
  pop_data_final<- select(pop_data_both,c(`Data_Series`,             
                                          begin_year:as.character(as.numeric(end_year)+1)))
  
  return(pop_data_final)
}





get_pop_each_group <- function(begin_year, end_year) {
  
  ## Total number of population in each age-group (Only Singapore Residents (Citizen+PR))
  pop_data = read_xlsx( here("data","singapore_population_data.xlsx"), 
                        sheet = "pop_age_group" )
  
  ## get the population data for n_year
  pop_data_final <- get_processed_population_data(pop_data = pop_data,
                                                  begin_year = begin_year,
                                                  end_year = end_year)
  ## Initial population in each age group
  pop_each_age_group <- as.numeric(unlist(pop_data_final[-1,
                                                         which(names(pop_data_final) %in% c(begin_year))]))
  
  return(pop_each_age_group)
}


get_birth_rate = function(begin_year, end_year, sim_year){
  
  ## Read birth rate data from sgstat.gov.sg
  birth_rate_data = read_xlsx( here("data","singapore_population_data.xlsx"), 
                               sheet = "birth_rate_per_1k" )
  
  
  birth_rate_projection = read_xlsx( here("data","singapore_population_data.xlsx"), 
                                     sheet = "birth_death_projection_medium" )
  
  ### 3. beginning with sgstat.gov.sg and the then projection from WPP
  
  lambda_hum_proj1 <- birth_rate_projection[,which(names(birth_rate_projection) ==
                                                     "Crude Birth Rate (births per 1,000 population)" )]
  
  lambda_hum1 <- birth_rate_data[(birth_rate_data$Data_Series %in%
                                    "Crude Birth Rate (Per Thousand Residents)"),
                                 which(names(birth_rate_data) %in% begin_year:end_year)]
  
  lambda_hum <- c(as.numeric(unlist(rev(lambda_hum1))), as.numeric(unlist(lambda_hum_proj1))[2:length(as.numeric(unlist(lambda_hum_proj1)))] )
  
  ## make it per day per 1k (time dependent)
  lambda_hum_tt <- (1/(365*1000))*rep(lambda_hum[1:(min(length(lambda_hum),sim_year))], each = 365)
  
  
  if (length(lambda_hum_tt) < length(t)) {
    lambda_hum_tt = c(lambda_hum_tt,
                      rep(lambda_hum_tt[length(lambda_hum_tt)],
                          (length(t) - length(lambda_hum_tt)))
    )
  }
  
  return(lambda_hum_tt)
  
}


get_death_rate = function(begin_year, end_year, sim_year){
  age_mort_data = read_xlsx( here("data","singapore_population_data.xlsx"), 
                             sheet = "age_specific_death_rate_per_1k" )
  
  ## Projection from WPP
  death_rate_projection =  read_xlsx( here("data","singapore_population_data.xlsx"), 
                                      sheet = "birth_death_projection_medium")
  
  ### 3. Beginning with sgstat.gov.sg and the then projection from WPP
  total_death_rate <- age_mort_data[age_mort_data$Data_Series ==
                                      "Total Age Specific Death Rate",
                                    -which(names(age_mort_data) %in% c("Data_Series"))]
  death_hum1 <- select(total_death_rate,  begin_year:end_year)
  death_hum_proj1 <- death_rate_projection[,which(names(death_rate_projection) == 
                                                    "Crude Death Rate (deaths per 1,000 population)" )]
  death_hum <- c(as.numeric(unlist(death_hum1)), as.numeric(unlist(death_hum_proj1))[2:length(as.numeric(unlist(death_hum_proj1)))] )
  ### make it daily
  death_hum_tt <- (1/(365*1000))*rep(death_hum[1:(min(length(death_hum),sim_year))],each = 365)
  
  ######### for missing info of death rate takes the last data point and continue #############
  ## death rate
  if (length(death_hum_tt) < length(t)) {
    death_hum_tt = c(death_hum_tt, 
                     rep(death_hum_tt[length(death_hum_tt)],
                         (length(t) - length(death_hum_tt)))
    )
  }
  return(death_hum_tt)
  
}



