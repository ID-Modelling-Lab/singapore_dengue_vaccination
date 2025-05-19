library(here)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(reshape2)
library(cowplot)


############## plot cases and hospitalization averted for different dominating serotype scenarios


## aggregate over all the age groups
get_agg_over_age <- function(arr) {
  arr_agg_age <- apply(arr, c(1,2,4), sum)
  return(arr_agg_age)
}

## cumulative number over year
get_cumulative_over_time <- function(arr, t_begin, t_end) {
  arr_cum_time <- apply(arr[, , t_begin:t_end], c(1, 2), sum)
  return(arr_cum_time)
}

## calculate the averted 
get_averted <- function(nv_output, v_output, t_begin, t_end, output_type){
  
  nv_agg_age <- get_agg_over_age(nv_output)
  v_agg_age <-  get_agg_over_age(v_output)
  
  nv_cum_t <- get_cumulative_over_time( nv_agg_age, t_begin = t_begin, t_end = t_end)
  v_cum_t <- get_cumulative_over_time( v_agg_age, t_begin = t_begin, t_end = t_end)
  
  if (output_type == "absolute") {
    
    averted <- (nv_cum_t - v_cum_t)
  }
  
  if (output_type == "proportion") {
    
    averted <- 100*(nv_cum_t - v_cum_t)/nv_cum_t
  }
  
  return(averted)
}


## extract variable from .rds file
extract_variable <- function(rds_file, variable_name) {
  
  result <- readRDS(here::here("model_output",  rds_file) )
  
  out_v <- as.data.frame(result$output, stringsAsFactors = FALSE)
  
  n_vac_age = result$n_vac_age
  n_sample = result$n_sample
  n_age = result$n_age
  n_year = result$n_year
  
  x <- out_v[,colnames(out_v) == variable_name] 
  
  y <- array(NA, dim = c(n_vac_age,n_sample, n_age, n_year))
  
  for (i in 1:n_vac_age) {
    
    y[i,,,] <- x[[paste0("result.", i)]]
  }
  
  return(y)
}

## make single figure
get_plot_averted <- function(df_v20, df_v50, df_v80, y_limits, y_breaks) {
  
  # Add a column for age groups
  df_v20$age <- 1:nrow(df_v20)
  df_v50$age <- 1:nrow(df_v50)
  df_v80$age <- 1:nrow(df_v80)
  
  # Reshape the dataframes to long format
  df_v20_long <- gather(df_v20, key = "Realization", value = "Output", -age)
  df_v50_long <- gather(df_v50, key = "Realization", value = "Output", -age)
  df_v80_long <- gather(df_v80, key = "Realization", value = "Output", -age)
  
  # Add a column to indicate the vaccine coverage for each dataset
  df_v20_long$VaccineCoverage <- "20% Coverage"
  df_v50_long$VaccineCoverage <- "50% Coverage"
  df_v80_long$VaccineCoverage <- "80% Coverage"
  
  # Combine all datasets into one long dataframe
  df_long <- bind_rows(df_v20_long, df_v50_long, df_v80_long)
  
  # Calculate the mean and 95% confidence interval for each age group and vaccine coverage
  summary_df <- df_long %>%
    group_by(age, VaccineCoverage) %>%
    summarise(
      mean_output = mean(Output),
      lower_bound = quantile(Output, 0.025),
      upper_bound = quantile(Output, 0.975)
    ) %>%
    ungroup()
  
  # Define colors for different vaccine coverages
  vaccine_colors <- c("20% Coverage" = "#1b7a53", "50% Coverage" = "#d94e28", "80% Coverage" = "#4a5a92")
  
  # Define the dodge position for side-by-side plotting
  dodge <- position_dodge(width = 0.3)
  
  # Plot using ggplot2
  p <- ggplot(summary_df, aes(x = factor(age), y = mean_output, color = VaccineCoverage)) +
    geom_point(size = 9, position = dodge, shape = 16) +  # Points for the mean with dodge
    geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), 
                  width = 0.0, size = 1.8, position = dodge) +  # Error bars with dodge
    scale_color_manual(values = vaccine_colors, name = "") + 
    scale_x_discrete(
      breaks = c(1, 2, 3, 4, 5, 6, 7),   # Set the positions of the x-axis ticks
      labels = c("6-16", "17-30", "31-40", "41-50", "51-60", "61-70", "71-80")  # Labels for age groups
    ) +
    # labs(x = "Age at vaccination", y = ylab, color = NULL) +  # Labels for x and y axes
    scale_y_continuous(limits = y_limits, breaks = y_breaks) +  # Define y-axis tick points )+
    theme_bw() +  # Clean theme
    theme(
      axis.ticks.length = unit(0.4, "cm"),        # Increase the length of the ticks
      axis.ticks = element_line(size = 1.5),      # Increase the thickness of the ticks
      panel.border = element_rect(color = "black", size = 2),  # Thicker border
      plot.margin = margin(t=10, r=10, b=20, l=20),
      axis.text.x = element_text(size = 35, color="black", angle= 0, margin = margin(t = 10, b = 20)),
      axis.text.y = element_text(size = 35,color="black", margin = margin(r = 10, l = 10)),
      # axis.title = element_text(size = 35, color = "black"),
      axis.title = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = element_text(size = 40),
      legend.key.size = unit(2, "cm"),  # Increase size of legend keys
      legend.key.height = unit(1.5, "cm"),  # Increase height of legend keys
      legend.key.width = unit(2.5, "cm"),  # Increase width of legend keys
      legend.background = element_blank()
    )
  
}


get_arranged_plot_avert_percent <- function(set_of_rds_files, variable_name, y_lab, t_begin, t_end) {
  
  plot_list <- list()

  for (i in seq_along(set_of_rds_files)) {
    
    rds_files <- set_of_rds_files[[i]]
    
    v_0 <- extract_variable(rds_file = rds_files[1] , 
                            variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[1], 
                       variable_name = paste0(variable_name, "_v"))
    
    v_20 <- extract_variable(rds_file = rds_files[2] , 
                             variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[2], 
                       variable_name = paste0(variable_name, "_v"))
    
    v_50 <- extract_variable(rds_file = rds_files[3] , 
                             variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[3], 
                       variable_name = paste0(variable_name, "_v"))
    
    v_80 <- extract_variable(rds_file = rds_files[4] , 
                             variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[4], 
                       variable_name = paste0(variable_name, "_v"))
    
    
    
    averted_inc_v20 <- get_averted(nv_output = v_0,
                                   v_output = v_20,
                                   t_begin = t_begin,
                                   t_end = t_end,
                                   output_type = "proportion")
    
    averted_inc_v50 <- get_averted(nv_output = v_0,
                                   v_output = v_50,
                                   t_begin = t_begin,
                                   t_end = t_end,
                                   output_type = "proportion")
    
    averted_inc_v80 <- get_averted(nv_output = v_0,
                                   v_output = v_80,
                                   t_begin = t_begin,
                                   t_end = t_end,
                                   output_type = "proportion")
    
    
    df_v20 <- as.data.frame(averted_inc_v20)
    df_v50 <- as.data.frame(averted_inc_v50)
    df_v80 <- as.data.frame(averted_inc_v80)
    
    
    
    # Determine y-axis settings for this figure
    if (i == 4) {  # Example: The second figure requires different limits
      y_limits <- custom_y_limits
      y_breaks <- custom_y_breaks
    } else {
      y_limits <- default_y_limits
      y_breaks <- default_y_breaks
    }
    
    averted_plot <- get_plot_averted(df_v20 = df_v20 , df_v50 = df_v50, df_v80 = df_v80, 
                                     y_limits = y_limits, y_breaks = y_breaks )
    
    
    # Remove x-axis ticks and labels for the first two plots
    # if (i %in% c(1, 2)) {
    #   averted_plot <- averted_plot +
    #     theme(axis.text.x = element_blank(),
    #           axis.ticks.x = element_blank())
    # }
    # 
    ## append the plots
    plot_list[[i]] <- averted_plot
    
  }

## remove legend in each of the figure
  plot_list_wo_legend <- list()
  
  for (i in seq_along(plot_list)){
    plot_list_wo_legend[[i]] <- plot_list[[i]] + theme(legend.position = "none")
  }  
  
## extrct the legend  
legend_plot <- get_legend(averted_plot)  

# Create the grid of plots without individual x and y labels
plot_wo_labels <- plot_grid(
  plotlist = plot_list_wo_legend,
  nrow = 2,
  ncol = 2,
  align = "hv",
  axis = "tblr",
  labels = c("(a)", "(b)", "(c)", "(d)"),  # Add subfigure labels
  label_size = 40,                        # Adjust label size
  label_x = 0.03,                         # Adjust horizontal position of labels
  label_y = 1.08                          # Adjust vertical position of labels
)

plot_with_labels <- ggdraw() +
  draw_plot(plot_wo_labels, x = 0, y = 0, width = 1, height = 1) +   
  draw_label(y_lab, x = 0.01, y = 0.5, angle = 90, size = 45, fontface = "bold") +  # Y-axis label
  draw_label("Age at vaccination (year)", x = 0.5, y = 0.025, angle = 0, size = 45, fontface = "bold")           # X-axis label

# make grid of plots  
final_plot_with_legend <- plot_grid(
  legend_plot,
  plot_with_labels,
  nrow = 2, 
  rel_heights = c(0.1,1)  
)

return(final_plot_with_legend)
  
}  



### Inputs

## List all the datasets that we want to plot together
set_of_rds_files <- list(
  c("pri_sec_neq_detailed_a7_sero_dom_denv1_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds", 
    "pri_sec_neq_detailed_a7_sero_dom_denv1_vcov_0.2_baseline_no_nsamp_100_waning_0.rds",
    "pri_sec_neq_detailed_a7_sero_dom_denv1_vcov_0.5_baseline_no_nsamp_100_waning_0.rds", 
    "pri_sec_neq_detailed_a7_sero_dom_denv1_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 1
  
  c("pri_sec_neq_detailed_a7_sero_dom_denv2_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds", 
    "pri_sec_neq_detailed_a7_sero_dom_denv2_vcov_0.2_baseline_no_nsamp_100_waning_0.rds",
    "pri_sec_neq_detailed_a7_sero_dom_denv2_vcov_0.5_baseline_no_nsamp_100_waning_0.rds", 
    "pri_sec_neq_detailed_a7_sero_dom_denv2_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 2
  
  c("pri_sec_neq_detailed_a7_sero_dom_denv3_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds", 
    "pri_sec_neq_detailed_a7_sero_dom_denv3_vcov_0.2_baseline_no_nsamp_100_waning_0.rds",
    "pri_sec_neq_detailed_a7_sero_dom_denv3_vcov_0.5_baseline_no_nsamp_100_waning_0.rds", 
    "pri_sec_neq_detailed_a7_sero_dom_denv3_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 3
  
  c("pri_sec_neq_detailed_a7_sero_dom_denv4_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds", 
    "pri_sec_neq_detailed_a7_sero_dom_denv4_vcov_0.2_baseline_no_nsamp_100_waning_0.rds",
    "pri_sec_neq_detailed_a7_sero_dom_denv4_vcov_0.5_baseline_no_nsamp_100_waning_0.rds", 
    "pri_sec_neq_detailed_a7_sero_dom_denv4_vcov_0.8_baseline_no_nsamp_100_waning_0.rds")    # Set 4
)

## grab the start of vaccination year from one of the dataset
t_begin = readRDS(here::here("model_output",set_of_rds_files[[1]][1]) )$v_year 

# General y-axis settings
default_y_limits <- c(0, 17)
default_y_breaks <- seq(0, 15, by = 3)

# Custom y-axis settings for a specific figure (for negative impact)(change it accordingly)
custom_y_limits <- c(-18, 15)  #default_y_limits #c(-18, 15)
custom_y_breaks <- seq(-18, 15, by = 6)  #default_y_breaks  #seq(-18, 15, by = 6)

## Generate figure now
plot_impact_hosp_dominating_sero <- get_arranged_plot_avert_percent(set_of_rds_files = set_of_rds_files,
                                                       variable_name = "symp_sample_age",
                                                       y_lab = "Reported cases averted (%)", 
                                                       t_begin = t_begin,
                                                       t_end = t_begin + 9)

figure_name <- "pri_sec_neq_percent_symp_averted_10_y_dom_sero_wo_waning.jpg"

## save
ggsave(
  filename = file.path(here("figure"), figure_name),  # Combine folder path and file name
  plot = plot_impact_hosp_dominating_sero,                                            # The plot object to save
  width = 29*2,                                           # Width in inches
  height = 40,                                          # Height in inches
  dpi = 300,
  units = "cm")                                  # DPI (resolution)







###################   Absolute averted per 1000 vaccinated

get_averted_per_vac <- function(nv_output, v_output, total_vac_output,  t_begin, t_end){
  
  nv_agg_age <- get_agg_over_age(nv_output)
  v_agg_age <-  get_agg_over_age(v_output)
  total_vac_agg_age <- get_agg_over_age(total_vac_output)
  
  nv_cum_t <- get_cumulative_over_time( nv_agg_age, t_begin = t_begin, t_end = t_end)
  v_cum_t <- get_cumulative_over_time( v_agg_age, t_begin = t_begin, t_end = t_end)
  total_vac_cum_t <- get_cumulative_over_time(total_vac_agg_age, t_begin = t_begin, t_end = t_end)
  # total_vac_cum_t <- get_cumulative_over_time(total_vac_agg_age, t_begin = t_end, t_end = t_end) ## this is beacuse we did not calculate from the beginning of the vaccination
  
  
  averted_per_vac <- (1000)*(nv_cum_t - v_cum_t)/total_vac_cum_t  #100*(nv_cum_t - v_cum_t)/nv_cum_t  #(1000)*(nv_cum_t - v_cum_t)/total_vac_cum_t
  
  
  
  return(averted_per_vac)
  
}



get_plot_averted_per_vac <- function(df_total, df_seroneg, df_seropos, y_limits, y_breaks) {
  
  
  # Add a column for age groups
  df_total$age <- 1:nrow(df_total)
  df_seroneg$age <- 1:nrow(df_seroneg)
  df_seropos$age <- 1:nrow(df_seropos)
  
  # Reshape the dataframes to long format
  df_total_long <- gather(df_total, key = "Realization", value = "Output", -age)
  df_seroneg_long <- gather(df_seroneg, key = "Realization", value = "Output", -age)
  df_seropos_long <- gather(df_seropos, key = "Realization", value = "Output", -age)
  
  # Add a column to indicate the vaccine coverage for each dataset
  df_total_long$poptype <- "All"
  df_seroneg_long$poptype <- "Seronegative"
  df_seropos_long$poptype <- "Seropositive"
  
  # Combine all datasets into one long dataframe
  df_long <- bind_rows(df_total_long, df_seroneg_long, df_seropos_long)
  
  # Calculate the mean and 95% confidence interval for each age group and vaccine coverage
  summary_df <- df_long %>%
    group_by(age, poptype) %>%
    summarise(
      mean_output = mean(Output),
      lower_bound = quantile(Output, 0.025),
      upper_bound = quantile(Output, 0.975)
    ) %>%
    ungroup()
  
  
  
  # Define colors for different vaccine coverages
  poptype_barcolors <- c("All" = "#8dd3c7", "Seronegative" = "#fdb462", "Seropositive" = "#b3de69")
  
  poptype_errorbar_colors <- c("All" = "#55a39b", "Seronegative" = "#e08d3e", "Seropositive" = "#85a84e")
  
  # Define the dodge position for side-by-side plotting
  dodge <- position_dodge(width = 0.3)
  
  # Plot using ggplot2
  p <- ggplot(summary_df, aes(x = factor(age), y = mean_output, fill = poptype)) +
    # geom_point(size = 9, position = dodge, shape = 16) +  # Points for the mean with dodge
    geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.65) +
    
    geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound, color = poptype), 
                  width = 0.0, size = 1.8, position = position_dodge(width = 0.75)) +
    # geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), 
    #               width = 0.0, size = 1.8, position = dodge) +  # Error bars with dodge
    scale_fill_manual(values = poptype_barcolors, name = "") + # Custom colors for bars
    scale_color_manual(values = poptype_errorbar_colors,  name = "") + # Custom colors for error bars
    # scale_color_manual(values = poptype_colors, name = "") + 
    scale_x_discrete(
      breaks = c(1, 2, 3, 4, 5, 6, 7),   # Set the positions of the x-axis ticks
      labels = c("6-16", "17-30", "31-40", "41-50", "51-60", "61-70", "71-80")  # Labels for age groups
    ) +
    # labs(x = "Age at vaccination", y = ylab, color = NULL) +  # Labels for x and y axes
    scale_y_continuous(limits = y_limits, breaks = y_breaks) +  # Define y-axis tick points )+
    theme_bw() +  # Clean theme
    theme(
      axis.ticks.length = unit(0.4, "cm"),        # Increase the length of the ticks
      axis.ticks = element_line(size = 1.5),      # Increase the thickness of the ticks
      panel.border = element_rect(color = "black", size = 2),  # Thicker border
      plot.margin = margin(t=10, r=10, b=20, l=25),
      axis.text.x = element_text(size = 35, color="black", angle= 0, margin = margin(t = 10, b = 20)),
      axis.text.y = element_text(size = 35,color="black", margin = margin(r = 10, l = 10)),
      # axis.title = element_text(size = 35, color = "black"),
      axis.title = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = element_text(size = 40),
      legend.key.size = unit(2, "cm"),  # Increase size of legend keys
      legend.key.height = unit(1.5, "cm"),  # Increase height of legend keys
      legend.key.width = unit(2.5, "cm"),  # Increase width of legend keys
      legend.background = element_blank()
    )
  
  
  
}


get_arranged_plot_averted_per_vac <- function(set_of_rds_files, variable_name, y_lab, t_begin, t_end) {
  
  plot_list <- list()
  
  for (i in seq_along(set_of_rds_files)) {
    
    rds_files <- set_of_rds_files[[i]]
    
    ## for total population
    nv_total <- extract_variable(rds_file = rds_files[1] , 
                            variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[1], 
                       variable_name = paste0(variable_name, "_v"))
    
    
    v_total <- extract_variable(rds_file = rds_files[2] , 
                                 variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[2], 
                       variable_name = paste0(variable_name, "_v"))
    
    new_vac_total <- extract_variable(rds_file = rds_files[2] , 
                                      variable_name = "new_vac_sample_age")
    
    ## seronegative
    nv_seroneg <- extract_variable(rds_file = rds_files[1] , 
                                 variable_name = paste0(variable_name, "_primary")) + 
      extract_variable(rds_file =  rds_files[1], 
                       variable_name = paste0(variable_name, "_v_primary"))
    
    
    v_seroneg <- extract_variable(rds_file = rds_files[2] , 
                                variable_name = paste0(variable_name, "_primary")) + 
      extract_variable(rds_file =  rds_files[2], 
                       variable_name = paste0(variable_name, "_v_primary"))
    
    ## seropositive
    nv_seropos <- extract_variable(rds_file = rds_files[1] , 
                                   variable_name = paste0(variable_name, "_secondary")) + 
      extract_variable(rds_file =  rds_files[1], 
                       variable_name = paste0(variable_name, "_v_secondary"))
    
    
    v_seropos <- extract_variable(rds_file = rds_files[2] , 
                                  variable_name = paste0(variable_name, "_secondary")) + 
      extract_variable(rds_file =  rds_files[2], 
                       variable_name = paste0(variable_name, "_v_secondary"))
    
    
    averted_per_vac_total <- get_averted_per_vac(nv_output = nv_total,
                                                 v_output = v_total,
                                                 total_vac_output = new_vac_total,
                                                 t_begin = t_begin,
                                                 t_end = t_end)
    
    averted_per_vac_seroneg <- get_averted_per_vac(nv_output = nv_seroneg,
                                                   v_output = v_seroneg,
                                                   total_vac_output = new_vac_total,
                                                   t_begin = t_begin,
                                                   t_end = t_end)
    
    averted_per_vac_seropos <- get_averted_per_vac(nv_output = nv_seropos,
                                                   v_output = v_seropos,
                                                   total_vac_output = new_vac_total,
                                                   t_begin = t_begin,
                                                   t_end = t_end)
    
    df_averted_per_vac_total <- as.data.frame(averted_per_vac_total)
    
    df_averted_per_vac_seroneg <- as.data.frame(averted_per_vac_seroneg)
    
    df_averted_per_vac_seropos <- as.data.frame(averted_per_vac_seropos)
    
    
    
    
    # Determine y-axis settings for this figure
    if (i == 3 ) {  # Example: The second figure requires different limits
      y_limits <- custom_y_limits1
      y_breaks <- custom_y_breaks1
      
    } else if (i == 4) {
      y_limits <- custom_y_limits2
      y_breaks <- custom_y_breaks2
      
    }else {
      y_limits <- default_y_limits
      y_breaks <- default_y_breaks
    }
    
    
    averted_plot_per_vac <- get_plot_averted_per_vac(
      df_total = df_averted_per_vac_total,
      df_seroneg = df_averted_per_vac_seroneg,
      df_seropos = df_averted_per_vac_seropos,
      y_limits = y_limits,
      y_breaks = y_breaks)
   
    
    
    plot_list[[i]] <-  averted_plot_per_vac
    
  }
  
  ## remove legend in each of the figure
  plot_list_wo_legend <- list()
  
  for (i in seq_along(plot_list)) {
    plot_list_wo_legend[[i]] <- plot_list[[i]] + theme(legend.position = "none")
  }  
  
  ## extrct the legend  
  legend_plot <- get_legend(averted_plot_per_vac)  
  
  # Create the grid of plots without individual x and y labels
  plot_wo_labels <- plot_grid(
    plotlist = plot_list_wo_legend,
    nrow = 2,
    ncol = 2,
    align = "hv",
    axis = "tblr",
    labels = c("(a)", "(b)", "(c)", "(d)"),  # Add subfigure labels
    label_size = 40,                        # Adjust label size
    label_x = 0.03,                         # Adjust horizontal position of labels
    label_y = 1.08                          # Adjust vertical position of labels
  )
  
  plot_with_labels <- ggdraw() +
    draw_plot(plot_wo_labels, x = 0, y = 0, width = 1, height = 1) +   
    draw_label(y_lab, x = 0.01, y = 0.55, angle = 90, size = 45, fontface = "bold") +  # Y-axis label
    draw_label("Age at vaccination (year)", x = 0.5, y = 0.025, angle = 0, size = 45, fontface = "bold")           # X-axis label
  
  # make grid of plots  
  final_plot_with_legend <- plot_grid(
    legend_plot,
    plot_with_labels,
    nrow = 2, 
    rel_heights = c(0.1,1)  
  )
  
  return(final_plot_with_legend)
  
}  



### Inputs

## List all the datasets that we want to plot together
set_of_rds_files <- list(
  c("pri_sec_neq_detailed_a7_sero_dom_denv1_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
    "pri_sec_neq_detailed_a7_sero_dom_denv1_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 1

  c("pri_sec_neq_detailed_a7_sero_dom_denv2_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
    "pri_sec_neq_detailed_a7_sero_dom_denv2_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 2

  c("pri_sec_neq_detailed_a7_sero_dom_denv3_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
    "pri_sec_neq_detailed_a7_sero_dom_denv3_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 3

  c("pri_sec_neq_detailed_a7_sero_dom_denv4_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
    "pri_sec_neq_detailed_a7_sero_dom_denv4_vcov_0.8_baseline_no_nsamp_100_waning_0.rds")    # Set 4
)

 
# set_of_rds_files <- list(
#   c("detailed_a7_sero_dom_denv1_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
#     "fake_effi_inf_detailed_a17_sero_dom_denv1_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 1
# 
#   c("detailed_a7_sero_dom_denv2_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
#     "fake_effi_inf_detailed_a17_sero_dom_denv2_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 2
# 
#   c("detailed_a7_sero_dom_denv3_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
#     "fake_effi_inf_detailed_a17_sero_dom_denv3_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 3
# 
#   c("detailed_a7_sero_dom_denv4_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
#     "fake_effi_inf_detailed_a17_sero_dom_denv4_vcov_0.8_baseline_no_nsamp_100_waning_0.rds")    # Set 4
# )

## grab the start of vaccination year from one of the dataset
t_begin = readRDS(here::here("model_output",set_of_rds_files[[1]][1]) )$v_year 

# # General y-axis settings
# default_y_limits <- c(0, 200)
# default_y_breaks <- seq(0, 200, by = 100)
# 
# # Custom y-axis settings for a specific figure (for negative impact)(change it accordingly)
# custom_y_limits1 <- c(0, 700)
# custom_y_breaks1 <- seq(0, 700, by = 100)
# 
# custom_y_limits2 <- c(0, 700)
# custom_y_breaks2 <- seq(0, 700, by = 100)



default_y_limits <- c(0, 20)
default_y_breaks <- seq(0, 20, by = 5)

# Custom y-axis settings for a specific figure (for negative impact)(change it accordingly)
custom_y_limits1 <- c(-48, 18)
custom_y_breaks1 <- seq(-45, 15, by = 15)

custom_y_limits2 <- c(-63, 18)
custom_y_breaks2 <- seq(-60, 15, by = 15)

## Generate figure now
plot_impact_hosp_per_vac_dominating_sero <- get_arranged_plot_averted_per_vac(set_of_rds_files = set_of_rds_files,
                                                                    variable_name = "symp_sample_age",
                                                                    y_lab = "Hospitalization cases averted(%)", 
                                                                    t_begin = t_begin,
                                                                    t_end = t_begin + 9)

figure_name <- "pri_sec_neq_hosp_sero_p_seron_serop_percent.jpg"

## save
ggsave(
  filename = file.path(here("figure"), figure_name),  # Combine folder path and file name
  plot = plot_impact_hosp_per_vac_dominating_sero,                                            # The plot object to save
  width = 29*2,                                           # Width in inches
  height = 40,                                          # Height in inches
  dpi = 300,
  units = "cm")                                  # DPI (resolution)








#################### sero share plot


extract_variable_sero <- function(rds_file, variable_name) {
  
  result <- readRDS(here::here("model_output",  rds_file) )
  
  out_v <- as.data.frame(result$output, stringsAsFactors = FALSE)
  
  
  
  n_sample = result$n_sample
  n_sero = result$n_sero
  n_year = result$n_year
  n_vac_age = result$n_vac_age
  
  x <- out_v[,colnames(out_v) == variable_name] 
  
  y <- array(NA, dim = c(n_vac_age,n_sample,n_sero, n_year))
  
  for (i in 1:n_vac_age) {
    
    y[i,,,] <- x[[paste0("result.", i)]]
  }
  
  return(y)
}


get_aggregated_sero_share_plot <- function(df_inf, t_begin, t_end){
  
  ## select only one scenario and one sample as the infection will remain same
  inf_sero_select <- df_inf[1,1,,t_begin:t_end]
  
  # Calculate the percentage of infection for each serotype
  # total_infection <- rowSums(inf_sero_select)
  percentage_vector <- 100*rowSums(inf_sero_select)/sum(inf_sero_select)
  
  serotypes <- paste("D", 1:4, sep = "-")
  
  df <- data.frame(
    serotype = serotypes,
    percentage = percentage_vector)
  # Plot stacked area plot using ggplot2 with improved aesthetics
  ggplot(df, aes(x = serotype, y = percentage, fill = serotype)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8, alpha = 0.9) +  # Stacked bar plot
    # labs(title = "",
    #      x = "Serotype", y = "Percentage of infection") +
    # Manually defining colors for the four serotypes
    scale_fill_manual(values = c("D-1" = "#8da0cb",   # Tomato
                                 "D-2" = "#e78ac3",   # SteelBlue
                                 "D-3" = "#a6d854",   # LimeGreen
                                 "D-4" = "#ffd92f"),
                      name="") + # Gold
    scale_y_continuous(limits = y_limits, breaks = y_breaks) + 
    theme_bw() + # Classic theme for a clean look
    theme(
      axis.ticks.length = unit(0.4, "cm"),        # Increase the length of the ticks
      axis.ticks = element_line(size = 1.5),      # Increase the thickness of the ticks
      panel.border = element_rect(color = "black", size = 2),  # Thicker border
      plot.margin = margin(t=10, r=10, b=10, l=40),
      axis.text.x = element_text(size = 35, color="black", angle= 0, margin = margin(t = 10, b = 20)),
      axis.text.y = element_text(size = 35,color="black", margin = margin(r = 10, l = 10)),
      # axis.title = element_text(size = 35, color = "black"),
      axis.title = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = element_text(size = 40),
      legend.key.size = unit(2, "cm"),  # Increase size of legend keys
      legend.key.height = unit(1.5, "cm"),  # Increase height of legend keys
      legend.key.width = unit(2.5, "cm"),  # Increase width of legend keys
      legend.background = element_blank()
      
    )
  
}


get_arranged_plot_aggregated_sero_share <- function(set_of_rds_files, y_lab, t_begin, t_end) {
  
  plot_list <- list()
  
  for (i in seq_along(set_of_rds_files)) {
    
    rds_files <- set_of_rds_files[[i]]
    
    inf_sero <- extract_variable_sero(rds_file = rds_files, variable_name = "inf_sample_sero") +
      extract_variable_sero(rds_file = rds_files, variable_name = "inf_sample_sero_v") 
    
    
    aggregated_sero_share_plot <- get_aggregated_sero_share_plot(df_inf = inf_sero, t_begin = t_begin, t_end = t_end )
    
    plot_list[[i]] <-  aggregated_sero_share_plot
    
  }
  
  ## remove legend in each of the figure
  plot_list_wo_legend <- list()
  
  for (i in seq_along(plot_list)){
    plot_list_wo_legend[[i]] <- plot_list[[i]] + theme(legend.position = "none")
  }  
  
  ## extract the legend  
  legend_plot <- get_legend(aggregated_sero_share_plot)  
  
  # Create the grid of plots without individual x and y labels
  plot_wo_labels <- plot_grid(
    plotlist = plot_list_wo_legend,
    nrow = 2,
    ncol = 2,
    align = "hv",
    axis = "tblr",
    labels = c("(a)", "(b)", "(c)", "(d)"),  # Add subfigure labels
    label_size = 40,                        # Adjust label size
    label_x = 0.03,                         # Adjust horizontal position of labels
    label_y = 1.1                          # Adjust vertical position of labels
  )
  
  plot_with_labels <- ggdraw() +
    draw_plot(plot_wo_labels, x = 0, y = 0, width = 1, height = 1) +   
    draw_label(y_lab, x = 0.015, y = 0.55, angle = 90, size = 45, fontface = "bold") +  # Y-axis label
    draw_label("Serotype", x = 0.55, y = 0.025, angle = 0, size = 45, fontface = "bold")           # X-axis label
  
  # make grid of plots  
  final_plot_with_legend <- plot_grid(
    legend_plot,
    plot_with_labels,
    nrow = 2, 
    rel_heights = c(0.12,1)  
  )
  
  return(final_plot_with_legend)
  
}  


### Inputs

## List all the datasets that we want to plot together
set_of_rds_files <- list(
  "pri_sec_neq_detailed_a7_sero_dom_denv1_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds", 
  "pri_sec_neq_detailed_a7_sero_dom_denv2_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds", 
  "pri_sec_neq_detailed_a7_sero_dom_denv3_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds", 
  "pri_sec_neq_detailed_a7_sero_dom_denv4_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds"
  )   


## grab the start of vaccination year from one of the dataset
t_begin = readRDS(here::here("model_output",set_of_rds_files[[1]]) )$v_year 

y_limits <- c(0, 80)
y_breaks <- seq(0, 80, by = 20)
## Generate figure now
plot_aggregated_sero_share <- get_arranged_plot_aggregated_sero_share(set_of_rds_files = set_of_rds_files,
                                                                              y_lab = "Share of infection(%)",
                                                                              t_begin = t_begin,
                                                                              t_end = t_begin + 9)

plot_aggregated_sero_share

figure_name <- "aggregated_sero_share_cases.jpg"




# # ### save
# ggsave(
#   filename = file.path(here("figure"), figure_name),  # Combine folder path and file name
#   plot = plot_aggregated_sero_share,                                            # The plot object to save
#   width = 20*2,                                           # Width in inches
#   height = 40,                                          # Height in inches
#   dpi = 300,
#   units = "cm")                                  # DPI (resolution)






##########   Impact matrix plot for vaccinated in different age groups


get_plot_averted_matrix <- function(df_averted_mean) {
  
  
  age_group_labels <- c("6-16","17-30",
                        "31-40","41-50",
                        "51-60","61-70", "71-80")  
  df <- melt(df_averted_mean)
  
  # Create heatmap using ggplot2
  p <- ggplot(df, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white", size = 0.1) +  # Add white borders for 
    # geom_text(aes(label = ifelse(floor(value) != 0, round(value, 1), "")), color = "black", size = 3) +  # Add numbers for non-zero values
    # # geom_text(aes(label = round(value, 1)), color = "white", size = 3) +
    geom_text(aes(label = floor(value)), color = "black", size = 10)+
    
    # Use a diverging color scale to handle negative and positive values
    scale_fill_gradient2(
      low = "#a6611a",       # Color for negative values
      mid = "#ffffbf",      # Color for zero
      high = "#018571",       # Color for positive values
      midpoint = 0,       # Define the midpoint as 0 for neutral color
      name = "Hospitalization \naverted (%)" # Label for the legend
    ) +
    theme_minimal()+
  
    # labs(x = "Vaccinated age group", y = "Age group", title = "") +
    theme(
      # axis.ticks.length = unit(0.4, "cm"),        # Increase the length of the ticks
      # axis.ticks = element_line(size = 1.5),      # Increase the thickness of the ticks
      # panel.border = element_rect(color = "black", size = 2),  # Thicker border
      plot.margin = margin(t=0, r=0, b=20, l=0),
      axis.text.x = element_text(margin = margin(t= 20, b = 0),angle = 45, color = "black", size = 30),
      axis.text.y = element_text(margin = margin(t = 0),angle = 45, color = "black", size = 30),#element_text(angle = 45, hjust = 1.0, vjust = 1.0, color = "black", size = 35),
      axis.title = element_blank(),
      # plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      panel.grid.major = element_blank(),  # Remove gridlines
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = element_text(size = 25 ),
      legend.title = element_text(size = 20, face = "bold"),
      legend.key.size = unit(1, "cm"),  # Increase size of legend keys
      legend.key.height = unit(1.0, "cm"),  # Increase height of legend keys
      legend.key.width = unit(1.2, "cm"),  # Increase width of legend keys
      legend.background = element_blank()
      # legend.position = "top",  # Position the legend on the right for clarity
      # legend.direction = "horizontal",
      # legend.key.width = unit(1.0, "cm")  # Make the legend key wider for better visibility
    ) + scale_x_continuous(breaks = 1:7, 
                           labels = age_group_labels) +
    scale_y_continuous(breaks = 1:7, 
                       labels = age_group_labels) +
    coord_fixed()  # Ensure square tiles
  
  
  p
}


get_cumulative_over_time_age_wise <- function(arr, t_begin, t_end, age_list) {
  arr_cum_time_age <- apply(arr[, , ,t_begin:t_end], c(1, 2, 3), sum)
  
  arr_cum_age_gr <- array(NA, dim = c(dim(arr)[1], dim(arr)[2], nrow(age_list) ))
  
  for (i in 1:nrow(age_list)) {
    
    arr_cum_age_gr[,,i] <-  apply(arr_cum_time_age[,,age_list[i,1]:age_list[i,2]], c(1, 2), sum)
    
    
  }
  
  return(arr_cum_age_gr)
}



get_averted_age <- function(nv_output, v_output, t_begin, t_end, age_list, output_type){
  
  
  nv_cum_t <- get_cumulative_over_time_age_wise( nv_output, t_begin = t_begin, t_end = t_end, age_list)
  v_cum_t <- get_cumulative_over_time_age_wise( v_output, t_begin = t_begin, t_end = t_end, age_list)
  
  if (output_type == "absolute") {
    
    averted <- (nv_cum_t - v_cum_t)
  }
  
  if (output_type == "proportion") {
    
    averted <- 100*(nv_cum_t - v_cum_t)/nv_cum_t
  }
  return(averted)
  
}


get_age_list <- function(rds_file) {
  
  result <- readRDS(here::here("model_output",  rds_file) )
  
  age_list <- result$v_age_list
  
  return(age_list)
}





get_arranged_plot_impact_matrix_age_gr <- function(set_of_rds_files, variable_name, y_lab, t_begin, t_end, age_list) {
  
  plot_list <- list()
  
  for (i in seq_along(set_of_rds_files)) {
    
    rds_files <- set_of_rds_files[[i]]
    
    
    
    nv_symp <- extract_variable(rds_file = rds_files[1], 
                                variable_name = variable_name) +
      extract_variable(rds_file = rds_files[1],
                       variable_name =  paste0(variable_name, "_v"))

    v_symp <- extract_variable(rds_file = rds_files[2], 
                                variable_name = variable_name) +
      extract_variable(rds_file = rds_files[2],
                       variable_name = paste0(variable_name, "_v"))

    
    ## to get rid of small number which will create large number when we divide to get the proportion
    # nv_symp[nv_symp < 10] <- 0
    # v_symp[v_symp < 10] <- 0
    
    averted_inc <- get_averted_age(nv_output = nv_symp,
                                   v_output = v_symp,
                                   t_begin = t_begin,
                                   t_end = t_end,
                                   age_list = age_list,
                                   output_type = "proportion")
    
    averted_inc_mean <-  apply(averted_inc, c(1, 3), mean)
    averted_inc_mean[is.nan(averted_inc_mean)] <- 0 
    
    averted_matrix_plot <- get_plot_averted_matrix(df_averted_mean = (averted_inc_mean))
    
    plot_list[[i]] <-  averted_matrix_plot
    
  }
  
  # plot_list <- lapply(plot_list, function(plot) {
  #   plot + theme(plot.margin = margin(t = 5, r = 0, b = 5, l = 0))
  # })
  # 
  # 
  ## remove legend in each of the figure
  # plot_list_wo_legend <- list()
  # 
  # for (i in seq_along(plot_list)) {
  #   plot_list_wo_legend[[i]] <- plot_list[[i]] + theme(legend.position = "none")
  # }  
  
  ## extrct the legend  
  # legend_plot <- get_legend(averted_plot_per_vac)  
  
  # Create the grid of plots without individual x and y labels
  plot_wo_labels <- plot_grid(
    plotlist = plot_list,
    nrow = 2,
    ncol = 2,
    align = "hv",
    axis = "tblr",
    labels = c("(a)", "(b)", "(c)", "(d)"),  # Add subfigure labels
    label_size = 40,                        # Adjust label size
    label_x = 0.08,                         # Adjust horizontal position of labels
    label_y = 1.0                          # Adjust vertical position of labels
  )
  
  plot_with_labels <- ggdraw() +
    draw_plot(plot_wo_labels, x = 0, y = 0, width = 1, height = 1) +   
    draw_label(y_lab, x = 0.01, y = 0.55, angle = 90, size = 45, fontface = "bold") +  # Y-axis label
    draw_label("Vaccinated age group (year)", x = 0.5, y = 0.023, angle = 0, size = 45, fontface = "bold")           # X-axis label
  
  # make grid of plots  
  # final_plot_with_legend <- plot_grid(
  #   legend_plot,
  #   plot_with_labels,
  #   nrow = 2, 
  #   rel_heights = c(0.1,1)  
  # )
  
  return(plot_with_labels)
  
}  



### Inputs

## List all the datasets that we want to plot together
set_of_rds_files <- list(
  c("detailed_a7_sero_dom_denv1_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds", 
    "detailed_a7_sero_dom_denv1_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 1
  
  c("detailed_a7_sero_dom_denv2_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds", 
    "detailed_a7_sero_dom_denv2_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 2
  
  c("detailed_a7_sero_dom_denv3_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds", 
    "detailed_a7_sero_dom_denv3_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 3
  
  c("detailed_a7_sero_dom_denv4_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds", 
    "detailed_a7_sero_dom_denv4_vcov_0.8_baseline_no_nsamp_100_waning_0.rds")    # Set 4
)

## grab the start of vaccination year from one of the dataset
t_begin = readRDS(here::here("model_output",set_of_rds_files[[1]][1]) )$v_year 

age_list <- get_age_list(rds_file = set_of_rds_files[[1]][1])

# y_limits <- c(0, 80)
# y_breaks <- seq(0, 80, by = 20)
## Generate figure now
plot_aggregated_impact_matrix_symp <- get_arranged_plot_impact_matrix_age_gr(set_of_rds_files = set_of_rds_files,
                                                                             variable_name = "hosp_sample_age",
                                                                      y_lab = "Age group (year)",
                                                                      t_begin = t_begin,
                                                                      t_end = t_begin + 9, 
                                                                      age_list = age_list)


plot_aggregated_impact_matrix_symp

figure_name <- "aggregated_impact_matrix_hosp.jpg"


# # ### save
ggsave(
  filename = file.path(here("figure"), figure_name),  # Combine folder path and file name
  plot = plot_aggregated_impact_matrix_symp,                                            # The plot object to save
  width = 20*2,                                           # Width in inches
  height = 40,                                          # Height in inches
  dpi = 300,
  units = "cm")                                  # DPI (resolution)

























###################   Absolute averted per 1000 vaccinated with and without efficacy against infection


get_agg_over_age <- function(arr) {
  arr_agg_age <- apply(arr, c(1,2,4), sum)
  return(arr_agg_age)
}

## cumulative number over year
get_cumulative_over_time <- function(arr, t_begin, t_end) {
  arr_cum_time <- apply(arr[, , t_begin:t_end], c(1, 2), sum)
  return(arr_cum_time)
}


get_averted_per_vac <- function(nv_output, v_output, total_vac_output,  t_begin, t_end){
  
  nv_agg_age <- get_agg_over_age(nv_output)
  v_agg_age <-  get_agg_over_age(v_output)
  total_vac_agg_age <- get_agg_over_age(total_vac_output)
  
  nv_cum_t <- get_cumulative_over_time( nv_agg_age, t_begin = t_begin, t_end = t_end)
  v_cum_t <- get_cumulative_over_time( v_agg_age, t_begin = t_begin, t_end = t_end)
  total_vac_cum_t <- get_cumulative_over_time(total_vac_agg_age, t_begin = t_begin, t_end = t_end)
  # total_vac_cum_t <- get_cumulative_over_time(total_vac_agg_age, t_begin = t_end, t_end = t_end) ## this is beacuse we did not calculate from the beginning of the vaccination
  
  
  averted_per_vac <- 100*(nv_cum_t - v_cum_t)/nv_cum_t #1000*(nv_cum_t - v_cum_t)/total_vac_cum_t
  
  
  
  return(averted_per_vac)
  
}





get_plot_averted_per_vac_with_without_effi_inf <- function(df_with_effi_inf, df_without_effi_inf, y_limits, y_breaks) {
  
  
  # Add a column for age groups
  df_with_effi_inf$age <- 1:nrow(df_with_effi_inf)
  df_without_effi_inf$age <- 1:nrow(df_without_effi_inf)
  
  
  # Reshape the dataframes to long format
  df_with_effi_inf_long <- gather(df_with_effi_inf, key = "Realization", value = "Output", -age)
  df_without_effi_inf_long <- gather(df_without_effi_inf, key = "Realization", value = "Output", -age)
  
  
  # Add a column to indicate the vaccine coverage for each dataset
  df_with_effi_inf_long$efficacy_type <- "With efficacy against infection"
  df_without_effi_inf_long$efficacy_type <- "Without efficacy against infection"
  
  # Combine all datasets into one long dataframe
  df_long <- bind_rows(df_with_effi_inf_long, df_without_effi_inf_long)
  
  # Calculate the mean and 95% confidence interval for each age group and vaccine coverage
  summary_df <- df_long %>%
    group_by(age, efficacy_type) %>%
    summarise(
      mean_output = mean(Output),
      lower_bound = quantile(Output, 0.025),
      upper_bound = quantile(Output, 0.975)
    ) %>%
    ungroup()
  
  # # Define colors for different vaccine coverages
  # pype_barcolors <- c("All" = "#8dd3c7", "Seronegative" = "#fdb462", "Seropositive" = "#b3de69")
  # 
  # poptype_errorbar_colors <- c("All" = "#55a39b", "Seronegative" = "#e08d3e", "Seropositive" = "#85a84e")
  # 
  # 
  # Define colors for different vaccine coverages
  efficacy_type_barcolors <- c("With efficacy against infection" = "#fb8072", "Without efficacy against infection" = "#80b1d3")
  efficacy_type_errorbar_colors <- c("With efficacy against infection" = "#cc5045", "Without efficacy against infection" = "#40759a")
  # Define the dodge position for side-by-side plotting
  dodge <- position_dodge(width = 0.5)
  
  # Plot using ggplot2
  p <- ggplot(summary_df, aes(x = factor(age), y = mean_output, fill = efficacy_type)) +
    # geom_point(size = 9, position = dodge, shape = 16) +  # Points for the mean with dodge
    geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.65) +
    
    geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound, color = efficacy_type), 
                  width = 0.0, size = 1.8, position = position_dodge(width = 0.75)) +
    # geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), 
    #               width = 0.0, size = 1.8, position = dodge) +  # Error bars with dodge
    scale_fill_manual(values = efficacy_type_barcolors, name = "") + # Custom colors for bars
    scale_color_manual(values = efficacy_type_errorbar_colors,  name = "") + # Custom colors for error bars
    # scale_color_manual(values = poptype_colors, name = "") + 
    scale_x_discrete(
      breaks = c(1, 2, 3, 4, 5, 6, 7),   # Set the positions of the x-axis ticks
      labels = c("6-16", "17-30", "31-40", "41-50", "51-60", "61-70", "71-80")  # Labels for age groups
    ) +
    # labs(x = "Age at vaccination", y = ylab, color = NULL) +  # Labels for x and y axes
    scale_y_continuous(limits = y_limits, breaks = y_breaks) +  # Define y-axis tick points )+
    theme_bw() +  # Clean theme
    theme(
      axis.ticks.length = unit(0.4, "cm"),        # Increase the length of the ticks
      axis.ticks = element_line(size = 1.5),      # Increase the thickness of the ticks
      panel.border = element_rect(color = "black", size = 2),  # Thicker border
      plot.margin = margin(t=10, r=10, b=20, l=25),
      axis.text.x = element_text(size = 35, color="black", angle= 0, margin = margin(t = 10, b = 20)),
      axis.text.y = element_text(size = 35,color="black", margin = margin(r = 10, l = 10)),
      # axis.title = element_text(size = 35, color = "black"),
      axis.title = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = element_text(size = 40),
      legend.key.size = unit(2, "cm"),  # Increase size of legend keys
      legend.key.height = unit(1.5, "cm"),  # Increase height of legend keys
      legend.key.width = unit(2.5, "cm"),  # Increase width of legend keys
      legend.background = element_blank()
    )
  
  
  
}



get_arranged_plot_averted_per_vac_with_without_effi_inf <- function(set_of_rds_files, variable_name, y_lab, t_begin, t_end) {
  
  plot_list <- list()
  
  for (i in seq_along(set_of_rds_files)) {
    
    rds_files <- set_of_rds_files[[i]]
    
    ## for total population
    nv_total <- extract_variable(rds_file = rds_files[1] , 
                                 variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[1], 
                       variable_name = paste0(variable_name, "_v"))
    
    
    v_without_effi_inf <- extract_variable(rds_file = rds_files[2] , 
                                variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[2], 
                       variable_name = paste0(variable_name, "_v"))
    
    v_with_effi_inf <- extract_variable(rds_file = rds_files[3] , 
                                        variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[3], 
                       variable_name = paste0(variable_name, "_v"))
    
    new_vac_total <- extract_variable(rds_file = rds_files[2] , 
                                      variable_name = "new_vac_sample_age")
    
    
   
    averted_per_vac_with_effi_inf <- get_averted_per_vac(nv_output = nv_total,
                                                 v_output = v_with_effi_inf,
                                                 total_vac_output = new_vac_total,
                                                 t_begin = t_begin,
                                                 t_end = t_end)
    
    averted_per_vac_without_effi_inf <- get_averted_per_vac(nv_output = nv_total,
                                                   v_output = v_without_effi_inf,
                                                   total_vac_output = new_vac_total,
                                                   t_begin = t_begin,
                                                   t_end = t_end)

    averted_per_vac_with_effi_inf <- as.data.frame( averted_per_vac_with_effi_inf)
    
    averted_per_vac_without_effi_inf <- as.data.frame(averted_per_vac_without_effi_inf)
    
    # df_averted_per_vac_seropos <- as.data.frame(averted_per_vac_seropos)
    
    
    
    
    # Determine y-axis settings for this figure
    if (i == 1 ) {  # Example: The second figure requires different limits
      y_limits <- custom_y_limits1
      y_breaks <- custom_y_breaks1
      
    } else if (i == 2) {
      y_limits <- custom_y_limits2
      y_breaks <- custom_y_breaks2
      
    }else if (i == 3) {
      y_limits <- custom_y_limits3
      y_breaks <- custom_y_breaks3
      
    }else if (i == 4) {
      y_limits <- custom_y_limits4
      y_breaks <- custom_y_breaks4
      
    }
    
    
    averted_plot_per_vac_with_without_effi_inf <- get_plot_averted_per_vac_with_without_effi_inf(
      df_with_effi_inf =  averted_per_vac_with_effi_inf,
      df_without_effi_inf = averted_per_vac_without_effi_inf,
      y_limits = y_limits,
      y_breaks = y_breaks)
    
    
    
    plot_list[[i]] <-  averted_plot_per_vac_with_without_effi_inf
    
  }
  
  ## remove legend in each of the figure
  plot_list_wo_legend <- list()
  
  for (i in seq_along(plot_list)) {
    plot_list_wo_legend[[i]] <- plot_list[[i]] + theme(legend.position = "none")
  }  
  
  ## extrct the legend  
  legend_plot <- get_legend(averted_plot_per_vac_with_without_effi_inf)  
  
  # Create the grid of plots without individual x and y labels
  plot_wo_labels <- plot_grid(
    plotlist = plot_list_wo_legend,
    nrow = 2,
    ncol = 2,
    align = "hv",
    axis = "tblr",
    labels = c("(a)", "(b)", "(c)", "(d)"),  # Add subfigure labels
    label_size = 40,                        # Adjust label size
    label_x = 0.03,                         # Adjust horizontal position of labels
    label_y = 1.08                          # Adjust vertical position of labels
  )
  
  plot_with_labels <- ggdraw() +
    draw_plot(plot_wo_labels, x = 0, y = 0, width = 1, height = 1) +   
    draw_label(y_lab, x = 0.01, y = 0.55, angle = 90, size = 45, fontface = "bold") +  # Y-axis label
    draw_label("Age at vaccination (year)", x = 0.5, y = 0.025, angle = 0, size = 45, fontface = "bold")           # X-axis label
  
  # make grid of plots  
  final_plot_with_legend <- plot_grid(
    legend_plot,
    plot_with_labels,
    nrow = 2, 
    rel_heights = c(0.1,1)  
  )
  
  return(final_plot_with_legend)
  
}  



### Inputs

## List all the datasets that we want to plot together
set_of_rds_files <- list(
  c("pri_sec_neq_detailed_a7_sero_dom_denv1_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds", 
    "pri_sec_neq_detailed_a7_sero_dom_denv1_vcov_0.8_baseline_no_nsamp_100_waning_0.rds",
    "effi_inf_detailed_a7_sero_dom_denv1_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 1
  
  c("pri_sec_neq_detailed_a7_sero_dom_denv2_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds", 
    "pri_sec_neq_detailed_a7_sero_dom_denv2_vcov_0.8_baseline_no_nsamp_100_waning_0.rds",
    "effi_inf_detailed_a7_sero_dom_denv2_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 2
  
  c("pri_sec_neq_detailed_a7_sero_dom_denv3_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds", 
    "pri_sec_neq_detailed_a7_sero_dom_denv3_vcov_0.8_baseline_no_nsamp_100_waning_0.rds",
    "effi_inf_detailed_a7_sero_dom_denv3_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 3
  
  c("pri_sec_neq_detailed_a7_sero_dom_denv4_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds", 
    "pri_sec_neq_detailed_a7_sero_dom_denv4_vcov_0.8_baseline_no_nsamp_100_waning_0.rds",
    "effi_inf_detailed_a7_sero_dom_denv4_vcov_0.8_baseline_no_nsamp_100_waning_0.rds")    # Set 4
)

## grab the start of vaccination year from one of the dataset
t_begin = readRDS(here::here("model_output",set_of_rds_files[[1]][1]) )$v_year 

# General y-axis settings

# custom_y_limits1 <- c(-100, 400)
# custom_y_breaks1 <- seq(-100, 400, by = 100)
# 
# custom_y_limits2 <- c(-300, 700)
# custom_y_breaks2 <- seq(-200, 600, by = 200)
# 
# # Custom y-axis settings for a specific figure (for negative impact)(change it accordingly)
# custom_y_limits3 <- c(0, 200)
# custom_y_breaks3 <- seq(0, 200, by = 50)
# 
# custom_y_limits4 <- c(-200, 250)



custom_y_limits1 <- c(-10, 50)
custom_y_breaks1 <- seq(-10, 50, by = 10)

custom_y_limits2 <- c(-40, 80)
custom_y_breaks2 <- seq(-40, 80, by = 20)

# Custom y-axis settings for a specific figure (for negative impact)(change it accordingly)
custom_y_limits3 <- c(0, 20)
custom_y_breaks3 <- seq(0, 20, by = 5)

custom_y_limits4 <- c(-20, 30)
custom_y_breaks4 <- seq(-20, 30, by = 10)

## Generate figure now
plot_fig <- get_arranged_plot_averted_per_vac_with_without_effi_inf(set_of_rds_files = set_of_rds_files,
                                                                              variable_name = "symp_sample_age",
                                                                              y_lab = "Reported cases averted (%)", 
                                                                              t_begin = t_begin,
                                                                              t_end = t_begin + 9)

figure_name <- "prop_symp_averted_per_vac_10_y_dom_sero_wo_waning_effi_inf.jpg"

### save
ggsave(
  filename = file.path(here("figure"), figure_name),  # Combine folder path and file name
  plot = plot_fig,                                            # The plot object to save
  width = 29*2,                                           # Width in inches
  height = 40,                                          # Height in inches
  dpi = 300,
  units = "cm")                                  # DPI (resolution)
















############################## Waning of efficacy



get_averted_waning <- function(nv_output, v_output,  t_begin, t_end){
  
  nv_agg_age <- get_agg_over_age(nv_output)
  v_agg_age <-  get_agg_over_age(v_output)
  # total_vac_agg_age <- get_agg_over_age(total_vac_output)
  
  nv_cum_t <- get_cumulative_over_time( nv_agg_age, t_begin = t_begin, t_end = t_end)
  v_cum_t <- get_cumulative_over_time( v_agg_age, t_begin = t_begin, t_end = t_end)
  # total_vac_cum_t <- get_cumulative_over_time(total_vac_agg_age, t_begin = t_begin, t_end = t_end)
  # total_vac_cum_t <- get_cumulative_over_time(total_vac_agg_age, t_begin = t_end, t_end = t_end) ## this is beacuse we did not calculate from the beginning of the vaccination
  
  
  averted <- 100*(nv_cum_t - v_cum_t)/nv_cum_t  #(1000)*(nv_cum_t - v_cum_t)/total_vac_cum_t
  
  
  
  return(averted)
  
}



get_plot_averted_waning <- function(df_wan0, df_wan1, df_wan2, 
                                    df_wan3, df_wan4, df_wan5,
                                    y_limits, y_breaks) {
  
  
  # Add a column for age groups
  df_wan0$age <- 1:nrow(df_wan0)
  df_wan1$age <- 1:nrow(df_wan1)
  df_wan2$age <- 1:nrow(df_wan2)
  df_wan3$age <- 1:nrow(df_wan3)
  df_wan4$age <- 1:nrow(df_wan4)
  df_wan5$age <- 1:nrow(df_wan5)
  
  
  # Reshape the dataframes to long format
  df_wan0_long <- gather(df_wan0, key = "Realization", value = "Output", -age)
  df_wan1_long <- gather(df_wan1, key = "Realization", value = "Output", -age)
  df_wan2_long <- gather(df_wan2, key = "Realization", value = "Output", -age)
  df_wan3_long <- gather(df_wan3, key = "Realization", value = "Output", -age)
  df_wan4_long <- gather(df_wan4, key = "Realization", value = "Output", -age)
  df_wan5_long <- gather(df_wan5, key = "Realization", value = "Output", -age)
  

  
  # Add a column to indicate the vaccine coverage for each dataset
  df_wan0_long$waning_rate <- "0%"
  df_wan1_long$waning_rate <- "5%"
  df_wan2_long$waning_rate <- "10%"
  df_wan3_long$waning_rate <- "15%"
  df_wan4_long$waning_rate <- "20%"
  df_wan5_long$waning_rate <- "25%"
  # 
  
  # Combine all datasets into one long dataframe
  df_long <- bind_rows(df_wan0_long, df_wan1_long, df_wan2_long,df_wan3_long,df_wan4_long, df_wan5_long)
  # Ensure waning_rate is a factor with the correct order
  df_long$waning_rate <- factor(
    df_long$waning_rate, 
    levels = c("0%", "5%", "10%", "15%", "20%", "25%")
  )
  
  # Calculate the mean and 95% confidence interval for each age group and vaccine coverage
  summary_df <- df_long %>%
    group_by(age, waning_rate) %>%
    summarise(
      mean_output = mean(Output),
      lower_bound = quantile(Output, 0.025),
      upper_bound = quantile(Output, 0.975)
    ) %>%
    ungroup()
  

  # Define colors for different vaccine coverages
  waning_colors <- c("0%" = "#003f5c", "5%" = "#2f4b7c", "10%" = "#665191", "15%" = "#a05195", "20%" = "#d45087", "25%" = "#f95d6a" )
  
  # Define the dodge position for side-by-side plotting
  dodge <- position_dodge(width = 0.9)
  
  # Plot using ggplot2
  p <- ggplot(summary_df, aes(x = factor(age), y = mean_output, color = waning_rate)) +
    geom_point(size = 9, position = dodge, shape = 16) +  # Points for the mean with dodge
    geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), 
                  width = 0.0, size = 1.8, position = dodge) +  # Error bars with dodge
    scale_color_manual(values = waning_colors, name = "Waning \nrate") + 
    scale_x_discrete(
      breaks = c(1, 2, 3, 4, 5, 6, 7),   # Set the positions of the x-axis ticks
      labels = c("6-16", "17-30", "31-40", "41-50", "51-60", "61-70", "71-80")  # Labels for age groups
    ) +
    # labs(x = "Age at vaccination", y = ylab, color = NULL) +  # Labels for x and y axes
    scale_y_continuous(limits = y_limits, breaks = y_breaks) +  # Define y-axis tick points )+
    theme_bw() +  # Clean theme
    theme(
      axis.ticks.length = unit(0.4, "cm"),        # Increase the length of the ticks
      axis.ticks = element_line(size = 1.5),      # Increase the thickness of the ticks
      panel.border = element_rect(color = "black", size = 2),  # Thicker border
      plot.margin = margin(t = 40, r = 10, b = 20, l = 20),
      axis.text.x = element_text(size = 35, color = "black", angle = 0, margin = margin(t = 10, b = 20)),
      axis.text.y = element_text(size = 35, color = "black", margin = margin(r = 10, l = 10)),
      # axis.title = element_text(size = 35, color = "black"),
      axis.title = element_blank(),
      legend.position = "right",
      legend.direction = "vertical",
      legend.text = element_text(size = 30),
      legend.title = element_text(size = 30),  # Adjust legend title size
      legend.key.size = unit(2, "cm"),  # Increase size of legend keys
      legend.key.height = unit(1.5, "cm"),  # Increase height of legend keys
      legend.key.width = unit(2.5, "cm"),  # Increase width of legend keys
      legend.background = element_blank()
    )
  
  
  
  
  
}


get_arranged_plot_averted_waning <- function(set_of_rds_files, variable_name, y_lab, t_begin, t_end) {
  
  plot_list <- list()
  
  for (i in seq_along(set_of_rds_files)) {
    
    rds_files <- set_of_rds_files[[i]]
    
    ## for total population
    nv <- extract_variable(rds_file = rds_files[1] , 
                                 variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[1], 
                       variable_name = paste0(variable_name, "_v"))
    
    
    v_0 <- extract_variable(rds_file = rds_files[2] , 
                                variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[2], 
                       variable_name = paste0(variable_name, "_v"))
    
    v_1 <- extract_variable(rds_file = rds_files[3] , 
                            variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[3], 
                       variable_name = paste0(variable_name, "_v"))
    
    v_2 <- extract_variable(rds_file = rds_files[4] , 
                            variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[4], 
                       variable_name = paste0(variable_name, "_v"))
    
    v_3 <- extract_variable(rds_file = rds_files[5] , 
                            variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[5], 
                       variable_name = paste0(variable_name, "_v"))
    
    v_4 <- extract_variable(rds_file = rds_files[6] , 
                            variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[6], 
                       variable_name = paste0(variable_name, "_v"))
    
    v_5 <- extract_variable(rds_file = rds_files[7] , 
                            variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[7], 
                       variable_name = paste0(variable_name, "_v"))
    
    
    
    averted_wan0 <- get_averted_waning(nv_output = nv,
                                                 v_output = v_0,
                                                 t_begin = t_begin,
                                                 t_end = t_end)
    
    averted_wan1 <- get_averted_waning(nv_output = nv,
                                       v_output = v_1,
                                       t_begin = t_begin,
                                       t_end = t_end)
    
    averted_wan2 <- get_averted_waning(nv_output = nv,
                                       v_output = v_2,
                                       t_begin = t_begin,
                                       t_end = t_end)
    
    averted_wan3 <- get_averted_waning(nv_output = nv,
                                       v_output = v_3,
                                       t_begin = t_begin,
                                       t_end = t_end)
    
    averted_wan4 <- get_averted_waning(nv_output = nv,
                                       v_output = v_4,
                                       t_begin = t_begin,
                                       t_end = t_end)
    
    
    averted_wan5 <- get_averted_waning(nv_output = nv,
                                       v_output = v_5,
                                       t_begin = t_begin,
                                       t_end = t_end)

    # 
    df_averted_wan0 <- as.data.frame(averted_wan0 )
    df_averted_wan1 <- as.data.frame(averted_wan1 )
    df_averted_wan2 <- as.data.frame(averted_wan2 )
    df_averted_wan3 <- as.data.frame(averted_wan3 )
    df_averted_wan4 <- as.data.frame(averted_wan4 )
    df_averted_wan5 <- as.data.frame(averted_wan5 )
    

    
    
    
    # Determine y-axis settings for this figure
    if (i == 3 ) {  # Example: The second figure requires different limits
      y_limits <- custom_y_limits1
      y_breaks <- custom_y_breaks1
      
    } else if (i == 4) {
      y_limits <- custom_y_limits2
      y_breaks <- custom_y_breaks2
      
    }else {
      y_limits <- default_y_limits
      y_breaks <- default_y_breaks
    }
    
    
    averted_plot_waning <- get_plot_averted_waning(
      df_wan0 = df_averted_wan0,
      df_wan1 = df_averted_wan1,
      df_wan2 = df_averted_wan2,
      df_wan3 = df_averted_wan3,
      df_wan4 = df_averted_wan4,
      df_wan5 = df_averted_wan5,
      y_limits = y_limits,
      y_breaks = y_breaks)
    
    
    
    plot_list[[i]] <-   averted_plot_waning
    
  }
  
  ## remove legend in each of the figure
  plot_list_wo_legend <- list()
  
  for (i in seq_along(plot_list)) {
    plot_list_wo_legend[[i]] <- plot_list[[i]] + theme(legend.position = "none")
  }  
  
  ## extrct the legend  
  legend_plot <- get_legend(averted_plot_waning)  
  
  # Create the grid of plots without individual x and y labels
  plot_wo_labels <- plot_grid(
    plotlist = plot_list_wo_legend,
    nrow = 2,
    ncol = 2,
    align = "hv",
    axis = "tblr",
    labels = c("(a)", "(b)", "(c)", "(d)"),  # Add subfigure labels
    label_size = 40,                        # Adjust label size
    label_x = 0.03,                         # Adjust horizontal position of labels
    label_y = 1.01                          # Adjust vertical position of labels
  )
  
  plot_with_labels <- ggdraw() +
    draw_plot(plot_wo_labels, x = 0, y = 0, width = 1, height = 1) +   
    draw_label(y_lab, x = 0.01, y = 0.55, angle = 90, size = 45, fontface = "bold") +  # Y-axis label
    draw_label("Age at vaccination (year)", x = 0.5, y = 0.025, angle = 0, size = 45, fontface = "bold")           # X-axis label
  
  # make grid of plots  
  final_plot_with_legend <- plot_grid(
    
    plot_with_labels,
    legend_plot,
    nrow = 1,
    ncol = 2,
    rel_widths = c(1,.1)
  )
  
  return(final_plot_with_legend)
  
}  





set_of_rds_files <- list(
  c("detailed_a7_sero_dom_denv1_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
    "detailed_a7_sero_dom_denv1_vcov_0.8_baseline_no_nsamp_100_waning_0.rds",
    "detailed_a7_sero_dom_denv1_vcov_0.8_baseline_no_nsamp_100_waning_0.05.rds",
    "detailed_a7_sero_dom_denv1_vcov_0.8_baseline_no_nsamp_100_waning_0.1.rds",
    "detailed_a7_sero_dom_denv1_vcov_0.8_baseline_no_nsamp_100_waning_0.15.rds",
    "detailed_a7_sero_dom_denv1_vcov_0.8_baseline_no_nsamp_100_waning_0.2.rds",
    "detailed_a7_sero_dom_denv1_vcov_0.8_baseline_no_nsamp_100_waning_0.25.rds"),  # Set 1
  
  c("detailed_a7_sero_dom_denv2_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
    "detailed_a7_sero_dom_denv2_vcov_0.8_baseline_no_nsamp_100_waning_0.rds",
    "detailed_a7_sero_dom_denv2_vcov_0.8_baseline_no_nsamp_100_waning_0.05.rds",
    "detailed_a7_sero_dom_denv2_vcov_0.8_baseline_no_nsamp_100_waning_0.1.rds",
    "detailed_a7_sero_dom_denv2_vcov_0.8_baseline_no_nsamp_100_waning_0.15.rds",
    "detailed_a7_sero_dom_denv2_vcov_0.8_baseline_no_nsamp_100_waning_0.2.rds",
    "detailed_a7_sero_dom_denv2_vcov_0.8_baseline_no_nsamp_100_waning_0.25.rds"),
  
  c("detailed_a7_sero_dom_denv3_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
    "detailed_a7_sero_dom_denv3_vcov_0.8_baseline_no_nsamp_100_waning_0.rds",
    "detailed_a7_sero_dom_denv3_vcov_0.8_baseline_no_nsamp_100_waning_0.05.rds",
    "detailed_a7_sero_dom_denv3_vcov_0.8_baseline_no_nsamp_100_waning_0.1.rds",
    "detailed_a7_sero_dom_denv3_vcov_0.8_baseline_no_nsamp_100_waning_0.15.rds",
    "detailed_a7_sero_dom_denv3_vcov_0.8_baseline_no_nsamp_100_waning_0.2.rds",
    "detailed_a7_sero_dom_denv3_vcov_0.8_baseline_no_nsamp_100_waning_0.25.rds"),
  
  c("detailed_a7_sero_dom_denv4_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
    "detailed_a7_sero_dom_denv4_vcov_0.8_baseline_no_nsamp_100_waning_0.rds",
    "detailed_a7_sero_dom_denv4_vcov_0.8_baseline_no_nsamp_100_waning_0.05.rds",
    "detailed_a7_sero_dom_denv4_vcov_0.8_baseline_no_nsamp_100_waning_0.1.rds",
    "detailed_a7_sero_dom_denv4_vcov_0.8_baseline_no_nsamp_100_waning_0.15.rds",
    "detailed_a7_sero_dom_denv4_vcov_0.8_baseline_no_nsamp_100_waning_0.2.rds",
    "detailed_a7_sero_dom_denv4_vcov_0.8_baseline_no_nsamp_100_waning_0.25.rds")
)
  
  
 

## grab the start of vaccination year from one of the dataset
t_begin = readRDS(here::here("model_output",set_of_rds_files[[1]][1]) )$v_year 

# General y-axis settings
default_y_limits <- c(0, 15)
default_y_breaks <- seq(0, 15, by = 5)

# Custom y-axis settings for a specific figure (for negative impact)(change it accordingly)
custom_y_limits1 <- c(-10,10 )
custom_y_breaks1 <- seq(-10, 10, by = 5)

custom_y_limits2 <- c(-30, 15)
custom_y_breaks2 <- seq(-30, 15, by = 10)



# default_y_limits <- c(-5, 10)
# default_y_breaks <- seq(-5, 10, by = 1)
# 
# # Custom y-axis settings for a specific figure (for negative impact)(change it accordingly)
# custom_y_limits1 <- c(-5, 10)
# custom_y_breaks1 <- seq(-5, 10, by = 5)
# 
# custom_y_limits2 <- c(-5, 15)
# custom_y_breaks2 <- seq(-5, 15, by = 5)

## Generate figure now
plot_fig <- get_arranged_plot_averted_waning(set_of_rds_files = set_of_rds_files,
                                                                              variable_name = "symp_sample_age",
                                                                              y_lab = "Reported cases averted(%)", 
                                                                              t_begin = t_begin,
                                                                              t_end = t_begin + 9)

figure_name <- "cases_waning.jpg"

## save
ggsave(
  filename = file.path(here("figure"), figure_name),  # Combine folder path and file name
  plot = plot_fig,                                            # The plot object to save
  width = 29*2,                                           # Width in inches
  height = 40,                                          # Height in inches
  dpi = 300,
  units = "cm")                                  # DPI (resolution)






############################################ Pre screening scenario





get_agg_over_age <- function(arr) {
  arr_agg_age <- apply(arr, c(1,2,4), sum)
  return(arr_agg_age)
}

## cumulative number over year
get_cumulative_over_time <- function(arr, t_begin, t_end) {
  arr_cum_time <- apply(arr[, , t_begin:t_end], c(1, 2), sum)
  return(arr_cum_time)
}


get_averted_screening <- function(nv_output, v_output,  t_begin, t_end){
  
  nv_agg_age <- get_agg_over_age(nv_output)
  v_agg_age <-  get_agg_over_age(v_output)
  
  nv_cum_t <- get_cumulative_over_time( nv_agg_age, t_begin = t_begin, t_end = t_end)
  v_cum_t <- get_cumulative_over_time( v_agg_age, t_begin = t_begin, t_end = t_end)
  
  averted_screening <- 100*(nv_cum_t - v_cum_t)/nv_cum_t 
  
  return( averted_screening)
}





get_plot_averted_screening <- function(df_with_screening, df_without_screening, y_limits, y_breaks) {
  
  
  # Add a column for age groups
  df_with_screening$age <- 1:nrow(df_with_screening)
  df_without_screening$age <- 1:nrow(df_without_screening)
  
  
  # Reshape the dataframes to long format
  df_with_screening_long <- gather(df_with_screening, key = "Realization", value = "Output", -age)
  df_without_screening_long <- gather(df_without_screening, key = "Realization", value = "Output", -age)
  
  
  # Add a column to indicate the vaccine coverage for each dataset
  df_with_screening_long$screening_type <- "With pre-screening"
  df_without_screening_long$screening_type <- "Without pre-screening"
  
  # Combine all datasets into one long dataframe
  df_long <- bind_rows(df_with_screening_long, df_without_screening_long)
  
  # Calculate the mean and 95% confidence interval for each age group and vaccine coverage
  summary_df <- df_long %>%
    group_by(age, screening_type) %>%
    summarise(
      mean_output = mean(Output),
      lower_bound = quantile(Output, 0.025),
      upper_bound = quantile(Output, 0.975)
    ) %>%
    ungroup()
  
  # # Define colors for different vaccine coverages
  # pype_barcolors <- c("All" = "#8dd3c7", "Seronegative" = "#fdb462", "Seropositive" = "#b3de69")
  # 
  # poptype_errorbar_colors <- c("All" = "#55a39b", "Seronegative" = "#e08d3e", "Seropositive" = "#85a84e")
  # 
  # 
  # Define colors for different vaccine coverages
  screening_type_barcolors <- c("With pre-screening" = "#d8b365", "Without pre-screening" = "#5ab4ac")
  screening_type_errorbar_colors <- c("With pre-screening" = "#8c510a", "Without pre-screening" = "#01665e")
  # Define the dodge position for side-by-side plotting
  dodge <- position_dodge(width = 0.5)
  
  # Plot using ggplot2
  p <- ggplot(summary_df, aes(x = factor(age), y = mean_output, fill = screening_type)) +
    # geom_point(size = 9, position = dodge, shape = 16) +  # Points for the mean with dodge
    geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.65) +
    
    geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound, color = screening_type), 
                  width = 0.0, size = 1.8, position = position_dodge(width = 0.75)) +
    # geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), 
    #               width = 0.0, size = 1.8, position = dodge) +  # Error bars with dodge
    scale_fill_manual(values = screening_type_barcolors, name = "") + # Custom colors for bars
    scale_color_manual(values = screening_type_errorbar_colors,  name = "") + # Custom colors for error bars
    # scale_color_manual(values = poptype_colors, name = "") + 
    scale_x_discrete(
      breaks = c(1, 2, 3, 4, 5, 6, 7),   # Set the positions of the x-axis ticks
      labels = c("6-16", "17-30", "31-40", "41-50", "51-60", "61-70", "71-80")  # Labels for age groups
    ) +
    # labs(x = "Age at vaccination", y = ylab, color = NULL) +  # Labels for x and y axes
    scale_y_continuous(limits = y_limits, breaks = y_breaks) +  # Define y-axis tick points )+
    theme_bw() +  # Clean theme
    theme(
      axis.ticks.length = unit(0.4, "cm"),        # Increase the length of the ticks
      axis.ticks = element_line(size = 1.5),      # Increase the thickness of the ticks
      panel.border = element_rect(color = "black", size = 2),  # Thicker border
      plot.margin = margin(t=10, r=10, b=20, l=25),
      axis.text.x = element_text(size = 35, color="black", angle= 0, margin = margin(t = 10, b = 20)),
      axis.text.y = element_text(size = 35,color="black", margin = margin(r = 10, l = 10)),
      # axis.title = element_text(size = 35, color = "black"),
      axis.title = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = element_text(size = 40),
      legend.key.size = unit(2, "cm"),  # Increase size of legend keys
      legend.key.height = unit(1.5, "cm"),  # Increase height of legend keys
      legend.key.width = unit(2.5, "cm"),  # Increase width of legend keys
      legend.background = element_blank()
    )
  
  
  
}



get_arranged_plot_averted_screening <- function(set_of_rds_files, variable_name, y_lab, t_begin, t_end) {
  
  plot_list <- list()
  
  for (i in seq_along(set_of_rds_files)) {
    
    rds_files <- set_of_rds_files[[i]]
    
    ## for total population
    nv_total <- extract_variable(rds_file = rds_files[1] , 
                                 variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[1], 
                       variable_name = paste0(variable_name, "_v"))
    
    
    v_without_screening <- extract_variable(rds_file = rds_files[2] , 
                                           variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[2], 
                       variable_name = paste0(variable_name, "_v"))
    
    v_with_screening <- extract_variable(rds_file = rds_files[3] , 
                                        variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[3], 
                       variable_name = paste0(variable_name, "_v"))
    # 
    # new_vac_total <- extract_variable(rds_file = rds_files[2] , 
    #                                   variable_name = "new_vac_sample_age")
    # 
    
    
    averted_with_screening <- get_averted_screening(nv_output = nv_total,
                                                         v_output = v_with_screening,
                                                         t_begin = t_begin,
                                                         t_end = t_end)
    
    averted_without_screening <- get_averted_screening(nv_output = nv_total,
                                                            v_output = v_without_screening,
                                                            t_begin = t_begin,
                                                            t_end = t_end)
    
    df_averted_with_screening <- as.data.frame(averted_with_screening )
    
    df_averted_without_screening <- as.data.frame(averted_without_screening)
    
    
    # Determine y-axis settings for this figure
    if (i == 1 ) {  # Example: The second figure requires different limits
      y_limits <- custom_y_limits1
      y_breaks <- custom_y_breaks1
      
    } else if (i == 2) {
      y_limits <- custom_y_limits2
      y_breaks <- custom_y_breaks2
      
    }else if (i == 3) {
      y_limits <- custom_y_limits3
      y_breaks <- custom_y_breaks3
      
    }else if (i == 4) {
      y_limits <- custom_y_limits4
      y_breaks <- custom_y_breaks4
      
    }
    
    
    averted_plot_screening <- get_plot_averted_screening(
      df_with_screening =  df_averted_with_screening ,
      df_without_screening = df_averted_without_screening,
      y_limits = y_limits,
      y_breaks = y_breaks)
    
    
    
    plot_list[[i]] <-  averted_plot_screening 
    
  }
  
  ## remove legend in each of the figure
  plot_list_wo_legend <- list()
  
  for (i in seq_along(plot_list)) {
    plot_list_wo_legend[[i]] <- plot_list[[i]] + theme(legend.position = "none")
  }  
  
  ## extrct the legend  
  legend_plot <- get_legend( averted_plot_screening )  
  
  # Create the grid of plots without individual x and y labels
  plot_wo_labels <- plot_grid(
    plotlist = plot_list_wo_legend,
    nrow = 2,
    ncol = 2,
    align = "hv",
    axis = "tblr",
    labels = c("(a)", "(b)", "(c)", "(d)"),  # Add subfigure labels
    label_size = 40,                        # Adjust label size
    label_x = 0.03,                         # Adjust horizontal position of labels
    label_y = 1.08                          # Adjust vertical position of labels
  )
  
  plot_with_labels <- ggdraw() +
    draw_plot(plot_wo_labels, x = 0, y = 0, width = 1, height = 1) +   
    draw_label(y_lab, x = 0.01, y = 0.55, angle = 90, size = 45, fontface = "bold") +  # Y-axis label
    draw_label("Age at vaccination (year)", x = 0.5, y = 0.025, angle = 0, size = 45, fontface = "bold")           # X-axis label
  
  # make grid of plots  
  final_plot_with_legend <- plot_grid(
    legend_plot,
    plot_with_labels,
    nrow = 2, 
    rel_heights = c(0.1,1)  
  )
  
  return(final_plot_with_legend)
  
}  



### Inputs

## List all the datasets that we want to plot together
set_of_rds_files <- list(
  c("pri_sec_neq_detailed_a7_sero_dom_denv1_vcov_0.8_baseline_yes_nsamp_100_waning_0_prescrng_1.rds", 
    "pri_sec_neq_detailed_a7_sero_dom_denv1_vcov_0.8_baseline_no_nsamp_100_waning_0.rds",
    "pri_sec_neq_detailed_a7_sero_dom_denv1_vcov_0.8_baseline_no_nsamp_100_waning_0_prescrng_1.rds"),  # Set 1
  
  c("pri_sec_neq_detailed_a7_sero_dom_denv2_vcov_0.8_baseline_yes_nsamp_100_waning_0_prescrng_1.rds", 
    "pri_sec_neq_detailed_a7_sero_dom_denv2_vcov_0.8_baseline_no_nsamp_100_waning_0.rds",
    "pri_sec_neq_detailed_a7_sero_dom_denv2_vcov_0.8_baseline_no_nsamp_100_waning_0_prescrng_1.rds"),   # Set 2
  
  c("pri_sec_neq_detailed_a7_sero_dom_denv3_vcov_0.8_baseline_yes_nsamp_100_waning_0_prescrng_1.rds", 
    "pri_sec_neq_detailed_a7_sero_dom_denv3_vcov_0.8_baseline_no_nsamp_100_waning_0.rds",
    "pri_sec_neq_detailed_a7_sero_dom_denv3_vcov_0.8_baseline_no_nsamp_100_waning_0_prescrng_1.rds"),  # Set 3
  
  c("pri_sec_neq_detailed_a7_sero_dom_denv4_vcov_0.8_baseline_yes_nsamp_100_waning_0_prescrng_1_sensiv_0.99.rds", 
    "pri_sec_neq_detailed_a7_sero_dom_denv4_vcov_0.8_baseline_no_nsamp_100_waning_0.rds",
    "pri_sec_neq_detailed_a7_sero_dom_denv4_vcov_0.8_baseline_no_nsamp_100_waning_0_prescrng_1_sensiv_0.99.rds") # set 4 
)

## grab the start of vaccination year from one of the dataset
t_begin = readRDS(here::here("model_output",set_of_rds_files[[1]][1]) )$v_year 

# General y-axis settings

# custom_y_limits1 <- c(-100, 400)
# custom_y_breaks1 <- seq(-100, 400, by = 100)
# 
# custom_y_limits2 <- c(-300, 700)
# custom_y_breaks2 <- seq(-200, 600, by = 200)
# 
# # Custom y-axis settings for a specific figure (for negative impact)(change it accordingly)
# custom_y_limits3 <- c(0, 200)
# custom_y_breaks3 <- seq(0, 200, by = 50)
# 
# custom_y_limits4 <- c(-200, 250)



custom_y_limits1 <- c(0, 20)
custom_y_breaks1 <- seq(0, 20, by = 5)

custom_y_limits2 <- c(0, 20)
custom_y_breaks2 <- seq(0, 20, by = 5)

# Custom y-axis settings for a specific figure (for negative impact)(change it accordingly)
custom_y_limits3 <- c(0, 20)
custom_y_breaks3 <- seq(0, 20, by = 5)

custom_y_limits4 <- c(-20, 20)
custom_y_breaks4 <- seq(-20, 20, by = 5)

## Generate figure now
plot_fig <- get_arranged_plot_averted_screening(set_of_rds_files = set_of_rds_files,
                                                                    variable_name = "symp_sample_age",
                                                                    y_lab = "Cases averted (%)", 
                                                                    t_begin = t_begin,
                                                                    t_end = t_begin + 9)

figure_name <- "prop_symp_averted_screening_try.jpg"

### save
ggsave(
  filename = file.path(here("figure"), figure_name),  # Combine folder path and file name
  plot = plot_fig,                                            # The plot object to save
  width = 29*2,                                           # Width in inches
  height = 40,                                          # Height in inches
  dpi = 300,
  units = "cm")                                  # DPI (resolution)







##############  for poster 


get_averted_per_vac <- function(nv_output, v_output, total_vac_output,  t_begin, t_end){
  
  nv_agg_age <- get_agg_over_age(nv_output)
  v_agg_age <-  get_agg_over_age(v_output)
  total_vac_agg_age <- get_agg_over_age(total_vac_output)
  
  nv_cum_t <- get_cumulative_over_time( nv_agg_age, t_begin = t_begin, t_end = t_end)
  v_cum_t <- get_cumulative_over_time( v_agg_age, t_begin = t_begin, t_end = t_end)
  total_vac_cum_t <- get_cumulative_over_time(total_vac_agg_age, t_begin = t_begin, t_end = t_end)
  # total_vac_cum_t <- get_cumulative_over_time(total_vac_agg_age, t_begin = t_end, t_end = t_end) ## this is beacuse we did not calculate from the beginning of the vaccination
  
  
  averted_per_vac <- 100*(nv_cum_t - v_cum_t)/nv_cum_t  #(1000)*(nv_cum_t - v_cum_t)/total_vac_cum_t
  
  
  
  return(averted_per_vac)
  
}





get_plot_averted_per_vac <- function(df_total, df_seroneg, df_seropos, y_limits, y_breaks) {
  
  
  # Add a column for age groups
  df_total$age <- 1:nrow(df_total)
  df_seroneg$age <- 1:nrow(df_seroneg)
  df_seropos$age <- 1:nrow(df_seropos)
  
  # Reshape the dataframes to long format
  df_total_long <- gather(df_total, key = "Realization", value = "Output", -age)
  df_seroneg_long <- gather(df_seroneg, key = "Realization", value = "Output", -age)
  df_seropos_long <- gather(df_seropos, key = "Realization", value = "Output", -age)
  
  # Add a column to indicate the vaccine coverage for each dataset
  df_total_long$poptype <- "All"
  df_seroneg_long$poptype <- "Seronegative"
  df_seropos_long$poptype <- "Seropositive"
  
  # Combine all datasets into one long dataframe
  df_long <- bind_rows(df_total_long, df_seroneg_long, df_seropos_long)
  
  # Calculate the mean and 95% confidence interval for each age group and vaccine coverage
  summary_df <- df_long %>%
    group_by(age, poptype) %>%
    summarise(
      mean_output = mean(Output),
      lower_bound = quantile(Output, 0.025),
      upper_bound = quantile(Output, 0.975)
    ) %>%
    ungroup()
  
  # 7,11,
  # 12,16,
  # 17,21,
  # 22,26,
  # 27,36,
  # 37,46,
  # 47,56,
  # 57,101
  
  # Define colors for different vaccine coverages
  poptype_barcolors <- c("All" = "#377eb8", "Seronegative" = "#4daf4a", "Seropositive" = "#984ea3")
  
  poptype_errorbar_colors <- c("All" = "#000000", "Seronegative" = "#000000", "Seropositive" = "#000000")
  
  # Define the dodge position for side-by-side plotting
  dodge <- position_dodge(width = 0.3)
  
  # Plot using ggplot2
  p <- ggplot(summary_df, aes(x = factor(age), y = mean_output, fill = poptype)) +
    # geom_point(size = 9, position = dodge, shape = 16) +  # Points for the mean with dodge
    geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.65) +
    
    geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound, color = poptype), 
                  width = 0.0, size = 1.8, position = position_dodge(width = 0.75)) +
    # geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), 
    #               width = 0.0, size = 1.8, position = dodge) +  # Error bars with dodge
    scale_fill_manual(values = poptype_barcolors, name = "") + # Custom colors for bars
    scale_color_manual(values = poptype_errorbar_colors,  name = "") + # Custom colors for error bars
    # scale_color_manual(values = poptype_colors, name = "") + 
    scale_x_discrete(
      breaks = c(1, 2, 3, 4, 5, 6, 7),   # Set the positions of the x-axis ticks
      labels = c("6-16", "17-30", "31-40", "41-50", "51-60", "61-70", "71-80")  
    ) +
    # labs(x = "Age at vaccination", y = ylab, color = NULL) +  # Labels for x and y axes
    scale_y_continuous(limits = y_limits, breaks = y_breaks) +  # Define y-axis tick points )+
    theme_bw() +  # Clean theme
    theme(
      axis.ticks.length = unit(0.4, "cm"),        # Increase the length of the ticks
      axis.ticks = element_line(size = 1.5),      # Increase the thickness of the ticks
      panel.border = element_rect(color = "black", size = 2),  # Thicker border
      plot.margin = margin(t=10, r=10, b=20, l=25),
      axis.text.x = element_text(size = 35, color="black", angle= 30, margin = margin(t = 10, b = 20)),
      axis.text.y = element_text(size = 35,color="black", margin = margin(r = 10, l = 10)),
      # axis.title = element_text(size = 35, color = "black"),
      axis.title = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = element_text(size = 40),
      legend.key.size = unit(2, "cm"),  # Increase size of legend keys
      legend.key.height = unit(1.5, "cm"),  # Increase height of legend keys
      legend.key.width = unit(2.5, "cm"),  # Increase width of legend keys
      legend.background = element_blank()
    )
  
  
  
}


get_arranged_plot_averted_per_vac <- function(set_of_rds_files, variable_name, y_lab, t_begin, t_end) {
  
  plot_list <- list()
  
  for (i in seq_along(set_of_rds_files)) {
    
    rds_files <- set_of_rds_files[[i]]
    
    ## for total population
    nv_total <- extract_variable(rds_file = rds_files[1] , 
                                 variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[1], 
                       variable_name = paste0(variable_name, "_v"))
    
    
    v_total <- extract_variable(rds_file = rds_files[2] , 
                                variable_name = variable_name) + 
      extract_variable(rds_file =  rds_files[2], 
                       variable_name = paste0(variable_name, "_v"))
    
    new_vac_total <- extract_variable(rds_file = rds_files[2] , 
                                      variable_name = "new_vac_sample_age")
    
    ## seronegative
    nv_seroneg <- extract_variable(rds_file = rds_files[1] , 
                                   variable_name = paste0(variable_name, "_primary")) + 
      extract_variable(rds_file =  rds_files[1], 
                       variable_name = paste0(variable_name, "_v_primary"))
    
    
    v_seroneg <- extract_variable(rds_file = rds_files[2] , 
                                  variable_name = paste0(variable_name, "_primary")) + 
      extract_variable(rds_file =  rds_files[2], 
                       variable_name = paste0(variable_name, "_v_primary"))
    
    ## seropositive
    nv_seropos <- extract_variable(rds_file = rds_files[1] , 
                                   variable_name = paste0(variable_name, "_secondary")) + 
      extract_variable(rds_file =  rds_files[1], 
                       variable_name = paste0(variable_name, "_v_secondary"))
    
    
    v_seropos <- extract_variable(rds_file = rds_files[2] , 
                                  variable_name = paste0(variable_name, "_secondary")) + 
      extract_variable(rds_file =  rds_files[2], 
                       variable_name = paste0(variable_name, "_v_secondary"))
    
    
    averted_per_vac_total <- get_averted_per_vac(nv_output = nv_total,
                                                 v_output = v_total,
                                                 total_vac_output = new_vac_total,
                                                 t_begin = t_begin,
                                                 t_end = t_end)
    
    averted_per_vac_seroneg <- get_averted_per_vac(nv_output = nv_seroneg,
                                                   v_output = v_seroneg,
                                                   total_vac_output = new_vac_total,
                                                   t_begin = t_begin,
                                                   t_end = t_end)
    
    averted_per_vac_seropos <- get_averted_per_vac(nv_output = nv_seropos,
                                                   v_output = v_seropos,
                                                   total_vac_output = new_vac_total,
                                                   t_begin = t_begin,
                                                   t_end = t_end)
    
    df_averted_per_vac_total <- as.data.frame(averted_per_vac_total)
    
    df_averted_per_vac_seroneg <- as.data.frame(averted_per_vac_seroneg)
    
    df_averted_per_vac_seropos <- as.data.frame(averted_per_vac_seropos)
    
    
    
    
    # Determine y-axis settings for this figure
    if (i == 3 ) {  # Example: The second figure requires different limits
      y_limits <- custom_y_limits3
      y_breaks <- custom_y_breaks3
      
    } else if (i == 4) {
      y_limits <- custom_y_limits4
      y_breaks <- custom_y_breaks4
      
    } else if (i == 1) {
      y_limits <- default_y_limits1
      y_breaks <- default_y_breaks1
    }  else if (i == 2) {
      y_limits <- default_y_limits2
      y_breaks <- default_y_breaks2
    }
    
    
    averted_plot_per_vac <- get_plot_averted_per_vac(
      df_total = df_averted_per_vac_total,
      df_seroneg = df_averted_per_vac_seroneg,
      df_seropos = df_averted_per_vac_seropos,
      y_limits = y_limits,
      y_breaks = y_breaks)
    
    
    
    plot_list[[i]] <-  averted_plot_per_vac
    
  }
  
  ## remove legend in each of the figure
  plot_list_wo_legend <- list()
  
  for (i in seq_along(plot_list)) {
    plot_list_wo_legend[[i]] <- plot_list[[i]] + theme(legend.position = "none")
  }  
  
  ## extrct the legend  
  legend_plot <- get_legend(averted_plot_per_vac)  
  
  # Create the grid of plots without individual x and y labels
  plot_wo_labels <- plot_grid(
    plotlist = plot_list_wo_legend,
    nrow = 2,
    ncol = 2,
    align = "hv",
    axis = "tblr",
    labels = c("(a)", "(b)", "(c)", "(d)"),  # Add subfigure labels
    label_size = 40,                        # Adjust label size
    label_x = 0.03,                         # Adjust horizontal position of labels
    label_y = 1.08                          # Adjust vertical position of labels
  )
  
  plot_with_labels <- ggdraw() +
    draw_plot(plot_wo_labels, x = 0, y = 0, width = 1, height = 1) +   
    draw_label(y_lab, x = 0.01, y = 0.55, angle = 90, size = 45, fontface = "bold") +  # Y-axis label
    draw_label("Age at vaccination (year)", x = 0.5, y = 0.025, angle = 0, size = 45, fontface = "bold")           # X-axis label
  
  # make grid of plots  
  final_plot_with_legend <- plot_grid(
    legend_plot,
    plot_with_labels,
    nrow = 2, 
    rel_heights = c(0.1,1)  
  )
  
  return(final_plot_with_legend)
  
}  



### Inputs

## List all the datasets that we want to plot together
set_of_rds_files <- list(
  c("pri_sec_neq_detailed_a7_sero_dom_denv1_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
    "pri_sec_neq_detailed_a7_sero_dom_denv1_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 1
  
  c("pri_sec_neq_detailed_a7_sero_dom_denv2_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
    "pri_sec_neq_detailed_a7_sero_dom_denv2_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 2
  
  c("pri_sec_neq_detailed_a7_sero_dom_denv3_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
    "pri_sec_neq_detailed_a7_sero_dom_denv3_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 3
  
  c("pri_sec_neq_detailed_a7_sero_dom_denv4_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
    "pri_sec_neq_detailed_a7_sero_dom_denv4_vcov_0.8_baseline_no_nsamp_100_waning_0.rds")    # Set 4
)


# set_of_rds_files <- list(
#   c("detailed_a7_sero_dom_denv1_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
#     "fake_effi_inf_detailed_a17_sero_dom_denv1_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 1
# 
#   c("detailed_a7_sero_dom_denv2_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
#     "fake_effi_inf_detailed_a17_sero_dom_denv2_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 2
# 
#   c("detailed_a7_sero_dom_denv3_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
#     "fake_effi_inf_detailed_a17_sero_dom_denv3_vcov_0.8_baseline_no_nsamp_100_waning_0.rds"),  # Set 3
# 
#   c("detailed_a7_sero_dom_denv4_vcov_0.8_baseline_yes_nsamp_100_waning_0.rds",
#     "fake_effi_inf_detailed_a17_sero_dom_denv4_vcov_0.8_baseline_no_nsamp_100_waning_0.rds")    # Set 4
# )

## grab the start of vaccination year from one of the dataset
t_begin = readRDS(here::here("model_output",set_of_rds_files[[1]][1]) )$v_year 

# # General y-axis settings
# default_y_limits <- c(0, 200)
# default_y_breaks <- seq(0, 200, by = 100)
# 
# # Custom y-axis settings for a specific figure (for negative impact)(change it accordingly)
# custom_y_limits1 <- c(0, 700)
# custom_y_breaks1 <- seq(0, 700, by = 100)
# 
# custom_y_limits2 <- c(0, 700)
# custom_y_breaks2 <- seq(0, 700, by = 100)




# default_y_limits <- c(0, 20)
# default_y_breaks <- seq(0, 20, by = 5)
# 
# # Custom y-axis settings for a specific figure (for negative impact)(change it accordingly)
# custom_y_limits1 <- c(-48, 18)
# custom_y_breaks1 <- seq(-45, 15, by = 15)
# 
# custom_y_limits2 <- c(-63, 18)
# custom_y_breaks2 <- seq(-60, 15, by = 15)



default_y_limits1 <- c(0, 15)
default_y_breaks1 <- seq(0, 15, by = 5)


default_y_limits2 <- c(0, 17)
default_y_breaks2 <- seq(0, 15, by = 5)

# Custom y-axis settings for a specific figure (for negative impact)(change it accordingly)
custom_y_limits3 <- c(-8, 15)
custom_y_breaks3 <- seq(-5, 15, by = 5)

custom_y_limits4 <- c(-63, 15)
custom_y_breaks4 <- seq(-60, 15, by = 15)

## Generate figure now
plot_impact_hosp_per_vac_dominating_sero <- get_arranged_plot_averted_per_vac(set_of_rds_files = set_of_rds_files,
                                                                              variable_name = "symp_sample_age",
                                                                              y_lab = "Reported cases averted(%)", 
                                                                              t_begin = t_begin,
                                                                              t_end = t_begin + 9)

figure_name <- "singa_pri_sec_neq_symp_sero_p_seron_serop_percent.jpg"

## save
ggsave(
  filename = file.path(here("figure"), figure_name),  # Combine folder path and file name
  plot = plot_impact_hosp_per_vac_dominating_sero,                                            # The plot object to save
  width = 29*2,                                           # Width in inches
  height = 40,                                          # Height in inches
  dpi = 300,
  units = "cm")                                  # DPI (resolution)









