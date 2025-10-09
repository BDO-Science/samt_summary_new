# Load required libraries
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(janitor)
library(patchwork)

#setting parameters
runs <- c('WR', 'LF')
year <- year(Sys.Date())

#renaming files to include date at the beginning for easier filtering
files <- list.files(path = 'StarsFiles/', pattern = '^Data.*\\.csv$', ignore.case = TRUE)

file_names <- lapply(files, function(file){
  parts <- strsplit(file, '\\.|_')[[1]]
  month <- parts[3]
  day <- parts[5]
  filedate <- mdy(paste0(month,'_',day,'_', year))
  file.rename(paste0('StarsFiles/',file),
              paste0('StarsFiles/',filedate,'_',file))
})

# pull most recent files for both WR and LF and read in data
data <- lapply(runs, function(run){
  file_list <- list.files(path = 'StarsFiles/', pattern = paste0("Data\\.",run,".*\\.csv$"), full.names = TRUE)
  if (length(file_list) > 0) {
    # Get the most recent file based on modification time
    most_recent_file <- max(file_list)
    
    # Read the most recent file into a data frame
    data <- read_csv(most_recent_file) %>% clean_names()
  } else {
    # Handle the case where no files are found
    data <- NULL
    message("No files found matching the pattern.")
  }
  data <- data %>% mutate(run = run)
})

wr_new <- data[[1]]
lf_new <- data[[2]]
#################
#winter-run STARs
#################

###making dataframes more readable for graphing
wr_survival <- wr_new |>
  select(contains("survival"), date) %>%
  pivot_longer(
    cols = starts_with("survival_"),
    names_to = c("location", "stat"),
    names_pattern = "survival_(.*)_(est|lcl_80|ucl_80)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  )%>%
  mutate(metric = 'Survival')

wr_routing <- wr_new |>
  select(contains("routing"), date) %>%
  pivot_longer(
    cols = starts_with("routing_"),
    names_to = c("location", "stat"),
    names_pattern = "routing_probability_(.*)_(est|lcl_80|ucl_80)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  ) %>%
  mutate(metric = 'Routing Probability')

wr_travel_time <- wr_new |>
  select(contains("travel"), date) %>%
  pivot_longer(
    cols = starts_with("median_"),
    names_to = c("location", "stat"),
    names_pattern = "median_travel_time_(.*)_(est|lcl_80|ucl_80)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  ) %>%
  mutate(metric = 'Travel Time')

wr_long <- bind_rows(wr_survival, wr_routing, wr_travel_time)

####wr plots
wr_survival_plot <- filter(wr_long, metric == 'Survival' &
                             location == 'overall') %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = lcl_80, ymax = ucl_80), fill = 'steelblue3', alpha = 0.2) +
  geom_line(aes(y = est), color = 'steelblue3', linewidth = 1)+
  labs(x = 'Date', y = 'Survival Probability', title = 'Overall Survival Probability') +
  ylim(0,1) +
  scale_x_date(limits = c(as.Date('2025-10-01'), as.Date('2026-04-30')), date_breaks = '1 months', date_labels = '%b') +
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 12))
  #geom_ribbon(aes(ymin = lcl_80, ymax = ucl_80), alpha = 0.2)
wr_survival_plot

wr_routing_plot <- filter(wr_long, metric == 'Routing Probability' &
                             location == 'interior_delta') %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = lcl_80, ymax = ucl_80), fill = 'steelblue3', alpha = 0.2) +
  geom_line(aes(y = est), color = 'steelblue3', linewidth = 1)+
  labs(x = 'Date', y = 'Routing Probability', title = 'Interior Delta Routing Probability') +
  ylim(0,1) +
  scale_x_date(limits = c(as.Date('2025-10-01'), as.Date('2026-04-30')), date_breaks = '1 months', date_labels = '%b') +
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 12))
#geom_ribbon(aes(ymin = lcl_80, ymax = ucl_80), alpha = 0.2)
wr_routing_plot

wr_timing_plot <- filter(wr_long, metric == 'Travel Time' &
                            location == 'overall') %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = lcl_80, ymax = ucl_80), fill = 'steelblue3', alpha = 0.2) +
  geom_line(aes(y = est), color = 'steelblue3', linewidth = 1)+
  labs(x = 'Date', y = 'Travel Time (in days)', title = 'Overall Travel Time') +
  scale_x_date(limits = c(as.Date('2025-10-01'), as.Date('2026-04-30')), date_breaks = '1 months', date_labels = '%b') +
  theme_bw() +
  theme(legend.position = 'bottom',
        plot.title = element_text(size = 12))
#geom_ribbon(aes(ymin = lcl_80, ymax = ucl_80), alpha = 0.2)
wr_timing_plot

wr_plot <- wr_survival_plot/wr_routing_plot/wr_timing_plot
wr_plot <- wr_plot + plot_annotation(
  title = 'Winter Run STARs model predictions',
  caption = 'Data were queried from: https://www.cbr.washington.edu/shiny/STARS/'
) & theme(plot.title = element_text(face = 'bold'))
wr_plot

#################
#late-fall STARs
#################
###making dataframes more readable for graphing
lf_survival <- lf_new |>
  select(contains("survival"), date) %>%
  pivot_longer(
    cols = starts_with("survival_"),
    names_to = c("location", "stat"),
    names_pattern = "survival_(.*)_(est|lcl|ucl)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  ) %>%
  mutate(metric = 'Survival')

lf_routing <- lf_new |>
  select(contains("routing"), date) %>%
  pivot_longer(
    cols = starts_with("routing_probability_"),
    names_to = c("location", "stat"),
    names_pattern = "routing_probability_(.*)_(est|lcl|ucl)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  ) %>%
  mutate(metric = 'Routing')

lf_travel <- lf_new |>
  select(contains("median"), date) %>%
  pivot_longer(
    cols = starts_with("median_travel_time_"),
    names_to = c("location", "stat"),
    names_pattern = "median_travel_time_(.*)_(est|lcl|ucl)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  ) %>%
  mutate(metric = 'Travel Time')

lf_long <- bind_rows(lf_survival,lf_routing, lf_travel)

lf_survival_plot <- filter(lf_long, metric == 'Survival' &
                             location == 'all_reaches') %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = 'steelblue3', alpha = 0.2) +
  geom_line(aes(y = est), color = 'steelblue3', linewidth = 1)+
  labs(x = 'Date', y = 'Survival Probability', title = 'Overall Survival Probability') +
  ylim(0,1) +
  scale_x_date(limits = c(as.Date('2025-10-01'), as.Date('2026-04-30')), date_breaks = '1 months', date_labels = '%b') +
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 12))
#geom_ribbon(aes(ymin = lcl_80, ymax = ucl_80), alpha = 0.2)
lf_survival_plot

lf_routing_plot <- filter(lf_long, metric == 'Routing' &
                            location == 'georgiana_slough') %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = 'steelblue3', alpha = 0.2) +
  geom_line(aes(y = est), color = 'steelblue3', linewidth = 1)+
  labs(x = 'Date', y = 'Routing Probability', title = 'Georgiana Slough Routing Probability') +
  ylim(0,1) +
  scale_x_date(limits = c(as.Date('2025-10-01'), as.Date('2026-04-30')), date_breaks = '1 months', date_labels = '%b') +
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 12))
#geom_ribbon(aes(ymin = lcl_80, ymax = ucl_80), alpha = 0.2)
lf_routing_plot

lf_timing_plot <- filter(lf_long, metric == 'Travel Time' &
                           location == 'all_reaches') %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = 'steelblue3', alpha = 0.2) +
  geom_line(aes(y = est), color = 'steelblue3', linewidth = 1)+
  labs(x = 'Date', y = 'Travel Time (in days)', title = 'Overall Travel Time') +
  scale_x_date(limits = c(as.Date('2025-10-01'), as.Date('2026-04-30')), date_breaks = '1 months', date_labels = '%b') +
  theme_bw() +
  theme(legend.position = 'bottom',
        plot.title = element_text(size = 12))
#geom_ribbon(aes(ymin = lcl_80, ymax = ucl_80), alpha = 0.2)
lf_timing_plot

lf_plot <- lf_survival_plot/lf_routing_plot/lf_timing_plot
lf_plot <- lf_plot + plot_annotation(
  title = 'Late-fall Run STARs model predictions',
  caption = 'Data were queried from: https://www.cbr.washington.edu/shiny/STARS/'
) & theme(plot.title = element_text(face = 'bold'))
lf_plot