library(tidyverse)
library(busdater)

wday <- readRDS('ReferenceCode/waterDay.rds')
url <- 'https://www.cbr.washington.edu/sacramento/data/php/rpt/juv_loss_detail.php?sc=1&outputFormat=csv&year=all&species=1%3Af&dnaOnly=no&age=no'
#set source folder, destination folder, and read in filenames of newest files
source_folder <- 'SalvageFiles/' #where the current files for the summary script to process live
file_names <- list.files(path = source_folder, 
                         pattern = "\\.csv$", full.names = TRUE) #read csv files in source folder

#pull in current loss data
salmon_raw <- read_csv(max(file_names[grep('salmon', file_names, ignore.case = TRUE)])) %>%
  mutate(Date = ymd(SampleDate)) %>%
  filter(DNA_Run == 'W')

data <- read.csv(url) %>% 
  select(1,2,14,15)

genetic <- bind_rows(read.csv('geneticdata/CVP_genetics.csv'), read.csv('geneticdata/SWP_genetics.csv')) %>%
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

realtime <- data.frame(Date = seq(as.Date('2024-11-01'), Sys.Date(), 1)) %>%
  left_join(select(salmon_raw, Date, CATCH)) %>%
  replace(is.na(.), 0) %>%
  mutate(cumul = cumsum(CATCH)) %>%
  mutate(wday = wday(Date))

labels <- genetic %>% group_by(WY) %>%
  summarize(wday = max(wday),
            cumul = max(cumul))

graph1 <- ggplot() +
  geom_line(genetic, mapping = aes(x = wday, y = cumul, group = factor(WY)), 
            color = 'darkgrey', linewidth = 1, linetype = 'dashed') +
  geom_line(realtime, mapping = aes(x = wday, y = cumul), color = 'steelblue3', linewidth = 1) +
  scale_x_continuous(breaks = c(31,93,153,214), labels = c('Nov', 'Jan', 'Mar', 'May')) +
  #geom_text(aes(x = 153, y = 100, label = Sys.Date())) +
  facet_wrap(~WY) +
  geom_vline(xintercept = wday(Sys.Date()), linetype = 'dotted', linewidth = 1) +
  theme_bw()
graph1

graph2 <- ggplot() +
  geom_line(genetic, mapping = aes(x = wday, y = cumul, group = factor(WY)), 
            color = 'grey', linewidth = 1, position = position_dodge(width = 0.01)) +
  geom_line(realtime, mapping = aes(x = wday, y = cumul), color = 'steelblue3', linewidth = 1) +
  scale_x_continuous(breaks = c(31,93,153,214), labels = c('Nov', 'Jan', 'Mar', 'May')) +
  geom_text(labels, mapping = aes(x = wday+4, y = cumul+4, label = factor(WY))) +
  geom_vline(xintercept = wday(Sys.Date()), linetype = 'dotted', linewidth = 1) +
  labs(x = 'Date', y = 'Cumulative Count in Salvage') +
  theme_bw() +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, unit = 'cm'),
        axis.title.y = element_text(margin=margin(r=15), size = 15),
        axis.title.x = element_text(margin=margin(t=15), size = 15),
        strip.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))
graph2

ggsave(graph2, file = 'outputs/genetic.png', width = 7, height = 5)