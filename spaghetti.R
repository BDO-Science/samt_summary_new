library(tidyverse)
library(busdater)

wday <- readRDS('CodeFiles/waterDay.rds')
url <- 'https://www.cbr.washington.edu/sacramento/data/php/rpt/juv_loss_detail.php?sc=1&outputFormat=csv&year=all&species=1%3Af&dnaOnly=no&age=no'
#set source folder, destination folder, and read in filenames of newest files
source_folder <- 'SalvageFiles/' #where the current files for the summary script to process live
file_names <- list.files(path = source_folder, 
                         pattern = "\\.csv$", full.names = TRUE) #read csv files in source folder

######Genetic WR spaghetti Graph
#pull in current loss data
salmon_raw <- read_csv(max(file_names[grep('salmon', file_names, ignore.case = TRUE)])) %>%
  mutate(Date = case_when(
    grepl('/', SampleDate) ~ mdy(SampleDate),
    grepl('-', SampleDate) ~ ymd(SampleDate),
    TRUE ~ NA_Date_  # Ensures missing values are handled properly
  )) %>%
  filter(DNA_Run == 'W')


genetic_wr <- bind_rows(read.csv('CodeFiles/CVP_genetics.csv'), read.csv('CodeFiles/SWP_genetics.csv')) %>%
  filter(Genetic_Assignment == 'Winter') %>%
  mutate(count = 1) %>%
  mutate(Date = as.Date(SampleDateTime)) %>%
  mutate(WY = get_fy(Date, opt_fy_start = '10-01')) %>%
  group_by(WY, Date) %>%
  summarize(count = sum(count)) %>%
  mutate(cumul = cumsum(count)) %>%
  ungroup() %>%
  mutate(wday = wday(Date)) %>%
  mutate(shortDate = format(Date, "%b %d")) %>%
  filter(wday >= wday(Sys.Date()))

realtime_wr <- data.frame(Date = seq(as.Date('2024-11-01'), Sys.Date(), 1)) %>%
  left_join(select(salmon_raw, Date, CATCH)) %>%
  replace(is.na(.), 0) %>%
  mutate(cumul = cumsum(CATCH)) %>%
  mutate(wday = wday(Date))

labels_wr <- genetic_wr %>% group_by(WY) %>%
  summarize(wday = max(wday),
            cumul = max(cumul))

graph_wr <- ggplot() +
  geom_line(genetic_wr, mapping = aes(x = wday, y = cumul, group = factor(WY)), 
            color = 'grey', linewidth = 1, position = position_dodge(width = 0.01)) +
  geom_line(realtime_wr, mapping = aes(x = wday, y = cumul), color = 'steelblue3', linewidth = 1) +
  scale_x_continuous(breaks = c(31,93,153,214), labels = c('Nov', 'Jan', 'Mar', 'May')) +
  geom_text(labels_wr, mapping = aes(x = wday+4, y = cumul+4, label = factor(WY))) +
  geom_vline(xintercept = wday(Sys.Date()), linetype = 'dotted', linewidth = 1) +
  labs(x = 'Date', y = 'Cumulative Count in Salvage') +
  theme_bw() +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, unit = 'cm'),
        axis.title.y = element_text(margin=margin(r=15), size = 15),
        axis.title.x = element_text(margin=margin(t=15), size = 15),
        strip.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))
graph_wr

ggsave(graph_wr, file = 'outputs/genetic_wr.png', width = 7, height = 5)


######steelhead spaghetti graph
steelhead_raw <- read_csv(max(file_names[grep('steelhead', file_names, ignore.case = TRUE)])) %>% 
  select(1,4,10) %>%
  pivot_longer(!Date, names_to = 'facility', values_to = 'loss') %>%
  mutate(Date = mdy(Date)) %>%
  mutate(facility = if_else(facility == 'Loss_SWPN', 'SWP', 'CVP'))

realtime_sh <- data.frame(Date = seq(as.Date('2024-11-01'), Sys.Date(), 1)) %>%
  left_join(select(steelhead_raw, Date, loss)) %>%
  replace(is.na(.), 0) %>%
  mutate(cumul = cumsum(loss)) %>%
  mutate(wday = wday(Date))

historic_sh <- read.csv('https://www.cbr.washington.edu/sacramento/data/php/rpt/juv_loss_detail.php?sc=1&outputFormat=csv&year=all&species=2%3Af&dnaOnly=no&age=no') %>%
  mutate(Date = as.Date(ymd_hms(Sample.Time))) %>%
  mutate(WY = get_fy(Date, opt_fy_start = '10-01')) %>%
  filter(WY < 2025 & WY > 2008) %>%
  group_by(WY, Date) %>%
  summarize(loss = sum(Loss)) %>%
  mutate(cumul = cumsum(loss)) %>%
  ungroup() %>%
  mutate(wday = wday(Date)) %>%
  filter(wday >= wday(Sys.Date()))

labels_sh <- historic_sh %>% group_by(WY) %>%
  summarize(wday = max(wday),
            cumul = max(cumul)) %>%
  filter(cumul >= max(realtime_sh$cumul))

graph_sh <- ggplot() +
  geom_line(historic_sh, mapping = aes(x = wday, y = cumul, group = factor(WY)), 
            color = 'grey', linewidth = 1, position = position_dodge(width = 0.01)) +
  geom_line(realtime_sh, mapping = aes(x = wday, y = cumul), color = 'steelblue3', linewidth = 1) +
  scale_x_continuous(breaks = c(31,93,153,214,244,274), labels = c('Nov', 'Jan', 'Mar', 'May','Jun','Jul')) +
  geom_text(labels_sh, mapping = aes(x = wday+4, y = cumul+4, label = factor(WY))) +
  geom_vline(xintercept = wday(Sys.Date()), linetype = 'dotted', linewidth = 1) +
  labs(x = 'Date', y = 'Cumulative Loss') +
  theme_bw() +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, unit = 'cm'),
        axis.title.y = element_text(margin=margin(r=15), size = 15),
        axis.title.x = element_text(margin=margin(t=15), size = 15),
        strip.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))
graph_sh
