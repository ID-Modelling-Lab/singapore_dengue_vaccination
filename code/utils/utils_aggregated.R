### yearly sero function

calc_year_sero <- function(arr) {
  
  ## this is for primary infection, where input has only 3 dimensions: (age, serotype, time)
  if (length(dim(arr)) == 3) {
    arr_sero <- apply(arr, c(2,3), sum) ## summing over age, i.e., 1st dimension
    arr_sero_yearly <- t(apply(arr_sero, 1, calc_yearly))
  }
  return(arr_sero_yearly)
}

## this is for prevalence of cases/vaccinated/population
calc_year_sero_wo_sum <- function(arr){
  
  if (length(dim(arr)) == 2){
    arr_sero <- arr
    arr_sero_yearly <- t(apply(arr_sero, 1, calc_yearly_wo_sum))
  }
  ## this is for primary infection, where input has only 3 dimensions: (age, serotype, time)
  if (length(dim(arr)) == 3){
    arr_sero <- apply(arr, c(2,3), sum) ## summing over age, i.e., 1st dimension
    arr_sero_yearly <- t(apply(arr_sero, 1, calc_yearly_wo_sum))
  }
  return(arr_sero_yearly)
}


## function for calculating yearly infection age wise
calc_year_age <- function(arr){
  
  if (length(dim(arr)) == 2) {
    arr_age <- arr
    arr_age_yearly <- t(apply(arr_age, 1, calc_yearly))
  }
  
  ## where input has only 3 dimensions: (age, serotype, time)
  if (length(dim(arr)) == 3) {
    arr_age <- apply(arr, c(1,3), sum) ## summing over sero, i.e., 2nd dimension
    arr_age_yearly <- t(apply(arr_age, 1, calc_yearly))
  }
  return(arr_age_yearly)
}

calc_year_age_wo_sum <- function(arr) {
  
  if (length(dim(arr)) == 2) {
    arr_age <- arr
    arr_age_yearly <- t(apply(arr_age, 1, calc_yearly_wo_sum))
  }
  
  ## where input has only 3 dimensions: (age, serotype, time)
  if (length(dim(arr)) == 3) {
    arr_age <- apply(arr, c(1,3), sum) ## summing over sero, i.e., 2nd dimension
    arr_age_yearly <- t(apply(arr_age, 1, calc_yearly_wo_sum))
  }
  return(arr_age_yearly)
}


## function for calculating yearly infection
calc_yearly <- function(mat) {
  tapply(mat, (seq_along(mat) - 1) %/%365, sum)
}

calc_yearly_wo_sum <- function(mat) {
  mat[seq(1,length(t),by = 365)]
}
