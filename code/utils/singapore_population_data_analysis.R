

population_age_extraction_from_data=function(pop_data,year, age_bracket){
  
  
pop_data_both_male_female= pop_data[1:(which(pop_data$`Data_Series`=="Singapore Male Residents")-1),]


pop_data_both_male_female_year= pop_data_both_male_female[as.character(year)]

pop_data_both_male_female_year_0_to_79=as.numeric(t(as.data.frame(pop_data_both_male_female_year)))[2:81]

pop_80_and_above= as.numeric(t(as.data.frame(pop_data_both_male_female_year)))[95]


out=c(pop_data_both_male_female_year_0_to_79,pop_80_and_above)

# Initialize an empty vector to store the sums
pop_each_age_group <- numeric()

# Set the number of elements to sum
age_bracket <- age_bracket  #TODO : age bracket can be a vector in the case the where age width is not uniform

# Iterate through the array and sum every chunk of 'chunk_size' elements
for (i in seq(1, age_bracket*floor(length(out)/age_bracket), by = age_bracket)) {
  pop_age_group <- sum(out[i:(i + age_bracket - 1)])
  pop_each_age_group <- c(pop_each_age_group, pop_age_group)
}

extra_age_group<-sum(out[(age_bracket*floor(length(out)/age_bracket) +1) :length(out)]) # remianing age group
pop_each_age_group <- c(pop_each_age_group,extra_age_group)


return(pop_each_age_group)

}

