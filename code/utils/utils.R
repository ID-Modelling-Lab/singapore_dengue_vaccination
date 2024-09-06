



#### Process age-specific mortality rate
get_age_specific_mortality <- function(age_mort_data, pop_data_final,begin_year,
                                       end_year){

### extract relevant info
death_rate = age_mort_data[1:(which(age_mort_data$Data_Series==
                                     "Male Age Specific Death Rate")-1),]

### Replace "-" with NA 
death_rate[death_rate=="-"] <- NA

## rows to discard
rows <- c("70 Years & Over", "75 Years & Over", "80 Years & Over", 
          "85 Years & Over")
death_rate_1 <- death_rate[!(death_rate$Data_Series %in% rows), ]

### death rate for specific timeline
death_data_final<- select(death_rate_1,c(`Data_Series`, begin_year:end_year))

## Average death rate per year per 1000 (same for all age group)
total_death_rate <- as.numeric(death_data_final[death_data_final$Data_Series==
                                                  "Total Age Specific Death Rate",
                                                -which(names(death_data_final) %in% c("Data_Series"))])
 
## The data contains average death rate (same for all age groups)
age_death_rate <- death_data_final[-which(death_data_final$Data_Series==
                                            "Total Age Specific Death Rate"),]

## Now we need to assign number to age group with year band 1 year

## 0 Year
under_1_year <- age_death_rate[which(age_death_rate$Data_Series=="Under 1 Year"),
                               which(names(age_death_rate) %in% begin_year:end_year)]

## 1-4 year
one_four_year <-age_death_rate[which(age_death_rate$Data_Series=="1 - 4 Years"),
                               which(names(age_death_rate) %in% begin_year:end_year)]

## 5-89 year
five_89_year <- age_death_rate[which(age_death_rate$Data_Series=="5 - 9 Years"):
                                 which(age_death_rate$Data_Series=="85 - 89 Years") ,
                               which(names(age_death_rate) 
                                     %in% begin_year:end_year)]

## 90+
ninety_plus_year <-age_death_rate[which(age_death_rate$Data_Series=="90 Years & Over"), 
                                  which(names(age_death_rate) %in% begin_year:end_year)]

## rbind with replicate values for each age group within an age band
replicated <- rbind(as.numeric(under_1_year), 
                    matrix(rep(as.numeric(one_four_year),each=4), ncol=n_year), 
                    apply(t(five_89_year), 1, function(x) rep(x,each=5)), 
                    as.numeric(ninety_plus_year) )

### Remove NA, by using interpolation
age_death_rate_for_each_age <- t(apply(replicated, 1, 
                                       function(x) na_interpolation(as.numeric(x),
                                                                    method = "linear")))

##### Function to calculate fraction of pop in each age groups

### Create a matrix for population in each age group over the time
pop_data_final_matrix <- matrix(unlist(pop_data_final[-1,which(names(pop_data_final) 
                                                               %in% begin_year:end_year) ]), 
                                ncol=n_year)

### array containing the pattern of how to do column sum
arr=c(dim(under_1_year)[1], rep(4,dim(one_four_year)[1]), 
      rep(5,dim(five_89_year)[1]), dim(ninety_plus_year)[1])

# Function to perform the summation according to the pattern
sum_by_pattern <- function(column, pattern) {
  result <- numeric(length(pattern))
  start_idx <- 1
  
  for (i in seq_along(pattern)) {
    end_idx <- start_idx + pattern[i] - 1
    result[i] <- sum(column[start_idx:end_idx])
    start_idx <- end_idx + 1
  }
  
  return(result)
}

# Apply the summation pattern to each column of the matrix
## (will create replication of total population)
sum_pop_matrix <- matrix(apply(pop_data_final_matrix, 2, sum_by_pattern, pattern = arr),
                         nrow = length(arr), byrow = FALSE)

# Replicate each column according to the pattern a
agg_pop_matrix <- matrix(unlist(lapply(seq_len(ncol(sum_pop_matrix)), 
                                       function(i) rep(sum_pop_matrix[, i], arr))), 
                         ncol = ncol(sum_pop_matrix))

## fraction of population each age age group
pop_age_frac=pop_data_final_matrix/agg_pop_matrix

## multiply the age-specific death rate from data to get approximation 
# of age-specific mortality (age adjusted)
death_hum=age_death_rate_for_each_age*pop_age_frac

return(death_hum)
}



