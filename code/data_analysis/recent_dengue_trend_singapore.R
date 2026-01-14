library(readxl)
library(ggplot2)
library(here)
library(dplyr)
library(tidyr)
library(cowplot)



get_plot_dengue_cases_trend <- function(){


# read the data
df_case_age_annual = read_xlsx(here("data","serotype_share_age_cases.xlsx"),
                               sheet = "age_case_2014_20")

# select the incidence per 100k data
df_case_age_inc <- select(df_case_age_annual, year, contains("_inc_"))


df_long_case <- df_case_age_inc %>%
  pivot_longer(
    cols = -year,
    names_to = "age_group",
    values_to = "incidence"
  ) %>%
  mutate(
    age_group = gsub("y_inc_p_100k", "", age_group)  # optional cleaning
  ) %>%
  mutate(
    age_group = factor(
      age_group,
      levels = c("0-4", "5-14", "15-24", "25-34",
                 "35-44", "45-54", "55-64", "65+")
    )
  )


# color_cases <- c(
#   "0-4"  = "#f7fbff",
#   "5-14" = "#deebf7",
#   "15-24" = "#c6dbef",
#   "25-34" = "#9ecae1",
#   "35-44" = "#6baed6",
#   "45-54" = "#4292c6",
#   "55-64" = "#2171b5",
#   "65+"   = "#084594"
# )

color_cases <- c(
  "0-4"   = "#c6dbef",
  "5-14"  = "#9ecae1",
  "15-24" = "#6baed6",
  "25-34" = "#4292c6",
  "35-44" = "#2171b5",
  "45-54" = "#08519c",
  "55-64" = "#08306b",
  "65+"   = "#041f4a"
)




p_cases <- ggplot(df_long_case, aes(x = factor(year), y = incidence, fill = age_group)) +
     geom_col(position = position_dodge(width = 0.9), color = "white", alpha = 0.9) +
     scale_fill_manual(values = color_cases) +
     scale_y_continuous(limits = c(0, 800), breaks = seq(0,800, by = 200)) +
    labs(
    x = "",
    y = "",
    fill = "Age (year)",
    title = ""
  ) +
  theme_bw() +  # Clean theme
  theme(
    axis.ticks.length = unit(0.4, "cm"),        # Increase the length of the ticks
    axis.ticks = element_line(size = 1.5),      # Increase the thickness of the ticks
    panel.border = element_rect(color = "black", size = 2),  # Thicker border
    plot.margin = margin(t=10, r=0, b=10, l=25),
    axis.text.x = element_text(size = 35, color="black", angle = 30, margin = margin(t = 10, b = 0)),
    axis.text.y = element_text(size = 35,color="black", margin = margin(r = 10, l = 20)),
    # axis.title.y = element_text(size = 35, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = c(0.45, 0.85),
    legend.direction = "horizontal",
    legend.text = element_text(size = 35),
    legend.title = element_text(size = 35, face = "bold"),
    legend.key.size = unit(1.5, "cm"),   # Increase size of legend keys
    legend.key.height = unit(1, "cm"),  # Increase height of legend keys
    legend.key.width = unit(1, "cm"),  # Increase width of legend keys
    legend.background = element_blank()
  )+
  guides(fill = guide_legend(
    nrow = 2,
    title.position = "top"    #THIS moves the title above legend keys
  ))


################ Now plot the serotyped incidence of cases


# read the data aggregated cases data
df_case_aggregated_annual = read_xlsx(here("data","serotype_share_age_cases.xlsx"),
                               sheet = "yearly_reported_cases")

df_case_aggregated_annual_selected_years <- df_case_aggregated_annual %>%
  filter(year >= 2014 & year <= 2020) %>%
  select(year, cases_p_100k)
  # select(year, absolute_cases)


## read sero share data
df_sero_share = read_xlsx(here("data","serotype_share_age_cases.xlsx"),
                          sheet = "sero_share")


df_sero_share_selected_years <- df_sero_share %>%
  filter(year >= 2014 & year <= 2020)


cases_vec <- df_case_aggregated_annual_selected_years$cases_p_100k   # extract vector of length 7

# cases_vec <- df_case_aggregated_annual_selected_years$absolute_cases

df_sero_inc_cases <- df_sero_share_selected_years %>%
  mutate(
    across(
      D1:D4,
      ~ (1/100)*.x * cases_vec,
      .names = "incidence_sero_{col}"
    )
  ) %>%
  select(year, starts_with("incidence_sero_")) %>%
  pivot_longer(
    cols = starts_with("incidence_sero"),
    names_to = "serotype",
    values_to = "incidence"
  ) %>%
  mutate(
    serotype = factor(serotype,
                      levels = c("incidence_sero_D1",
                                 "incidence_sero_D2",
                                 "incidence_sero_D3",
                                 "incidence_sero_D4"),
                      labels = c("DENV-1", "DENV-2", "DENV-3", "DENV-4"))
  )



color_sero <- c(
  "DENV-1" = "#7570b3",
  "DENV-2" = "#b2182b",
  "DENV-3" = "#66a61e",
  "DENV-4" = "#e6ab02"
)

p_sero <- ggplot(df_sero_inc_cases, aes(x = factor(year), y = incidence, fill = serotype)) +
  geom_col(position = position_dodge(width = 0.9), color = "white", alpha = 0.9) +
  scale_fill_manual(values = color_sero) +
  scale_y_continuous(limits = c(0, 300), breaks = seq(0,300, by = 50)) +
  labs(
    x = "",
    y = "",
    fill = "Serotype",
    title = ""
  ) +
  theme_bw() +  # Clean theme
  theme(
    axis.ticks.length = unit(0.4, "cm"),        # Increase the length of the ticks
    axis.ticks = element_line(size = 1.5),      # Increase the thickness of the ticks
    panel.border = element_rect(color = "black", size = 2),  # Thicker border
    plot.margin = margin(t=10, r=10, b=10, l=10),
    axis.text.x = element_text(size = 35, color="black", angle = 30, margin = margin(t = 10, b = 0)),
    axis.text.y = element_text(size = 35,color="black", margin = margin(r = 10, l = 20)),
    # axis.title.y = element_text(size = 35, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = c(0.35, 0.85),
    legend.direction = "horizontal",
    legend.text = element_text(size = 35),
    legend.title = element_text(size = 35, face = "bold"),
    legend.key.size = unit(1.5, "cm"),  # Increase size of legend keys
    legend.key.height = unit(1, "cm"),  # Increase height of legend keys
    legend.key.width = unit(1, "cm"),  # Increase width of legend keys
    legend.background = element_blank()
  )+
  guides(fill = guide_legend(
    nrow = 2,
    title.position = "top"    #THIS moves the title above legend keys
  ))

figure_list <- list()
figure_list[[1]] <- p_cases
figure_list[[2]] <- p_sero


plot_wo_labels <- plot_grid(
  plotlist = figure_list,
  nrow = 1,
  ncol = 2,
  align = "hv",
  axis = "tblr",
  labels = c("(a)", "(b)"),  # Add subfigure labels
  label_size = 35,                        # Adjust label size
  label_x = 0.02,                         # Adjust horizontal position of labels
  label_y = 1.02                          # Adjust vertical position of labels
)

plot_with_labels <- ggdraw() +
  draw_plot(plot_wo_labels, x = 0, y = 0, width = 1, height = 1) +   
  draw_label("Incidence per 100000 population", x = 0.01, y = 0.5, angle = 90, size = 35, fontface = "bold") +  # Y-axis label
  draw_label("Year", x = 0.52, y = 0.04, angle = 0, size = 35, fontface = "bold")           # X-axis label

# make grid of plots  


return(plot_with_labels)


}


plot_cases_trend <- get_plot_dengue_cases_trend()


figure_name <- "cases_recent_trend_singapore.jpg"

# ## save
ggsave(
  filename = file.path(here("figure"), figure_name),  # Combine folder path and file name
  plot = plot_cases_trend,                                            # The plot object to save
  width = 29*2,                                           # Width in inches
  height = 20,                                          # Height in inches
  dpi = 300,
  units = "cm") 








#### Now plot weekly dengue cases



df_weekly_cases = read.csv(here("data","dengue-cases-climate_with-2023.csv"))


df <- df_weekly_cases %>%
  filter(year >= 2014 & year <= 2023) %>%
  select(year, eweek, dengue_cases) %>%
  mutate(year = factor(year))

# generate 10 colors from light to dark
# my_cols <- colorRampPalette(c("#f7fcf5", "#00441b"))(length(unique(df$year)))
# names(my_cols) <- levels(df$year)   # important: map colors to year labels

my_cols <- c(
  "2014" = "#1f77b4",  # blue
  "2015" = "#ff7f0e",  # orange
  "2016" = "#2ca02c",  # green
  "2017" = "#d62728",  # red
  "2018" = "#9467bd",  # purple
  "2019" = "#8c564b",  # brown
  "2020" = "#e377c2",  # pink
  "2021" = "#7f7f7f",  # gray
  "2022" = "#bcbd22",  # olive
  "2023" = "#17becf"   # cyan
)

p_weekly_cases <- ggplot(df, aes(x = eweek, y = dengue_cases, 
               color = year, group = year)) +
  geom_line(linewidth = 2) +
  
  scale_color_manual(values = my_cols, name = "Year") +
  
  scale_x_continuous(breaks = seq(1, 52, 7)) +
  scale_y_continuous(limits = c(0, 1800), breaks = seq(0,1800,300)) +
  
  labs(
    x = "Epidemiological week",
    y = "Dengue cases"
  ) +
  
  theme_bw() +
  theme(
    axis.ticks.length = unit(0.4, "cm"),
    axis.ticks = element_line(size = 1.5),
    panel.border = element_rect(color = "black", size = 2),
    plot.margin = margin(t=10, r=10, b=10, l=10),
    axis.text.x = element_text(size = 35, color="black"),
    axis.text.y = element_text(size = 35, color="black"),
    axis.title.y = element_text(size = 35, face="bold"),
    axis.title.x = element_text(size = 35, face="bold"),
    
    legend.position = c(.90,0.65),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25, face = "bold"),
    
    legend.key.size = unit(1, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.key.width = unit(1, "cm"),
    
    legend.background = element_blank(),
    
    # ensure legend shows line, not a box
    legend.key = element_rect(fill = NA)
  )

figure_name <- "weekly_cases_singapore.jpg"

# ## save
ggsave(
  filename = file.path(here("figure"), figure_name),  # Combine folder path and file name
  plot = p_weekly_cases,                                            # The plot object to save
  width = 29*1,                                           # Width in inches
  height = 20,                                          # Height in inches
  dpi = 300,
  units = "cm") 


