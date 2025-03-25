library(tidyverse)
library(busdater)

wy <- get_fy(Sys.Date(), opt_fy_start = '07-01')  #pull the water year based on BY designation in LTO docs
jpe <- 98893

#set source folder, destination folder, and read in filenames of newest files
source_folder <- 'SalvageFiles/' #where the current files for the summary script to process live
file_names <- list.files(path = source_folder, 
                         pattern = "\\.csv$", full.names = TRUE) #read csv files in source folder

salmon_raw <- read_csv(max(file_names[grep('salmon', file_names, ignore.case = TRUE)])) %>%
  mutate(Date = case_when(
    grepl('/', SampleDate) ~ mdy(SampleDate),
    grepl('-', SampleDate) ~ ymd(SampleDate),
    TRUE ~ NA_Date_  # Ensures missing values are handled properly
  ))

wr <- salmon_raw %>%
  filter(AdClip == 'N') %>%
  filter(.[[8]] == 'W' & (.[[9]] == 'W'| is.na(.[[9]]))) %>%
  mutate(facility = if_else(FACILITY == 2, 'CVP', 'SWP')) %>%
  select(Date, facility, loss=LOSS) %>%
  mutate(species = 'Winter-run') %>%
  mutate(cumul_loss = cumsum(loss))

cumul_max <- max(wr$cumul_loss)
max_date <- max(wr$Date)
thresholds <- data.frame(threshold = c('100% Threshold', '75% Threshold', '50% Threshold'),
                         value = c(jpe*.005, jpe*.005*.75, jpe*.005*.5))

wr_cumul_graph <- ggplot() +
  geom_line(wr, mapping = aes(x = Date, y = cumul_loss), linewidth = 1, color = 'steelblue3') +
  geom_hline(thresholds, mapping = aes(yintercept = value), linetype = 'dotted', linewidth = 1, color = 'red') +
  geom_vline(mapping = aes(xintercept = Sys.Date()), linetype = 'dashed', linewidth = 1, color = 'darkgrey') +
  geom_text(thresholds, mapping = aes(x = as.Date('2025-01-10'), y = value + 20, label = threshold, fontface = 'bold')) +
  geom_label(mapping = aes(x = max_date + 1, y = cumul_max *1.8, 
                                    label = paste0(cumul_max, ' (', round((cumul_max/(jpe*.005))*100,0), '%)')), 
             fontface = 'bold', size = 4) +
  scale_x_date(limits = c(as.Date('2025-01-01'), as.Date('2025-06-01'))) +
  labs(y = 'Cumulative Loss', title = 'WY2025 Cumulative Loss of Natural Winter-run') +
  theme_bw() +
  theme(plot.margin = margin(0.5, 0.5, 0.2, 0.2, unit = 'cm'),
        axis.title.y = element_text(margin = margin(r = 15), size = 15),
        axis.title.x = element_text(margin = margin(t = 15), size = 15),
        strip.text = element_text(size = 13),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1),  # Adjust angle
        axis.text.y = element_text(size = 13),
        legend.position = 'bottom')
wr_cumul_graph

ggsave(wr_cumul_graph, file = 'outputs/wr_cumul_graph.png', width = 9, height = 6)
