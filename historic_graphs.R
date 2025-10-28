library(tidyverse)
library(readxl)
library(busdater)
library(janitor)

#url for importing steelhead data from sacpas
url <- 'https://www.cbr.washington.edu/sacramento/data/php/rpt/juv_loss_detail.php?sc=1&outputFormat=csv&year=all&species=2%3Af&dnaOnly=no&age=no'

###############################
#Histograms of historic loss
###############################

CVP <- read_csv('CodeFiles/CVP_genetics.csv')
SWP <- read_csv('CodeFiles/SWP_genetics.csv')
all <- bind_rows(CVP, SWP) %>%
  mutate(day = yday(SampleDateTime)) %>%
  filter(Genetic_Assignment == 'Winter') %>%
  mutate(Date = format(SampleDateTime, '%m/%d'))

today <- yday(Sys.Date())
#####genetics histogram
genetics <- ggplot(all, aes(day, fill = Facility)) +
  geom_histogram(position = 'identity', color = 'darkgrey') +
  facet_wrap(~Facility, ncol = 1) +
  geom_vline(xintercept = yday(Sys.Date()), linewidth = 1, linetype = 'dotted', color = '#999999') +
  labs(y = 'Count', x = 'Date') +
  scale_x_continuous(breaks = c(1,32,61,92,122,153),
                     labels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun'),
                     limits = c(1,160)) +
  annotate("text", x = yday(Sys.Date()) + 6, y = 20, label = "Today", fontface = 'italic') +
  theme_bw() +
  theme(plot.margin = margin(0.5, 0.5, 0.1, 0.2, unit = 'cm'),
        axis.title.y = element_text(margin = margin(r = 15), size = 15),
        axis.title.x = element_text(margin = margin(t = 10), size = 15),
        strip.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),  # Adjust angle
        axis.text.y = element_text(size = 13),
        legend.position = 'none')
genetics

###steelhead loss histogram
sh_import <- read_csv(url) %>% #import data from sacpas
  clean_names()

sh_all <- sh_import %>% filter(adipose_clip %in% c('Unclipped', 'Clipped')) %>%
  mutate(Date = as.Date(ymd_hms(sample_time))) %>%
  mutate(day = yday(Date))

sh_all_graph <- sh_all %>%
  ggplot() +
  geom_histogram(aes(x = day, fill = adipose_clip), color = 'darkgrey') +
  facet_wrap(~facility) +
  scale_x_continuous(breaks = c(1,32,61,92,122,153,183),
                     labels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul'),
                     limits = c(1,200)) +
  geom_vline(xintercept = yday(Sys.Date()), linewidth = 1, linetype = 'dotted', color = '#555555') +
  annotate("text", x = yday(Sys.Date()) + 15, y = 1000, label = "Today", fontface = 'italic') +
  labs(x = 'Date', y = 'Count') +
  theme_bw() +
  theme(plot.margin = margin(0.5, 0.5, 0.1, 0.2, unit = 'cm'),
        axis.title.y = element_text(margin = margin(r = 15), size = 15),
        axis.title.x = element_text(margin = margin(t = 10), size = 15),
        strip.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),  # Adjust angle
        axis.text.y = element_text(size = 13),
        legend.position = 'bottom',
        legend.title = element_blank())
sh_all_graph

###############################
#Spaghetti Graphs
###############################

wDay <- readRDS('CodeFiles/waterDay.rds')
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

realtime_wr <- data.frame(Date = seq(as.Date('2025-10-21'), Sys.Date(), 1)) %>%
  left_join(select(salmon_raw, Date, CATCH)) %>%
  replace(is.na(.), 0) %>%
  mutate(cumul = cumsum(CATCH)) %>%
  mutate(wday = wday(Date))

labels_wr <- genetic_wr %>% group_by(WY) %>%
  summarize(wday = max(wday),
            cumul = max(cumul))

spaghetti_wr <- ggplot() +
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
spaghetti_wr


######steelhead spaghetti graph
steelhead_raw <- read_csv(max(file_names[grep('steelhead', file_names, ignore.case = TRUE)])) %>% 
  select(1,4,10) %>%
  pivot_longer(!Date, names_to = 'facility', values_to = 'loss') %>%
  mutate(Date = mdy(Date)) %>%
  mutate(facility = if_else(facility == 'Loss_SWPN', 'SWP', 'CVP'))

realtime_sh <- data.frame(Date = seq(as.Date('2025-10-21'), Sys.Date(), 1)) %>%
  left_join(select(steelhead_raw, Date, loss)) %>%
  replace(is.na(.), 0) %>%
  mutate(cumul = cumsum(loss)) %>%
  mutate(wday = wDay(Date))

historic_sh <- read_csv(url) %>%
  clean_names() %>%
  mutate(Date = as.Date(ymd_hms(sample_time))) %>%
  mutate(WY = get_fy(Date, opt_fy_start = '10-01')) %>%
  filter(WY < 2025 & WY > 2008) %>%
  group_by(WY, Date) %>%
  summarize(loss = sum(loss)) %>%
  mutate(cumul = cumsum(loss)) %>%
  ungroup() %>%
  mutate(wday = wDay(Date)) %>%
  filter(wday >= wDay(Sys.Date()) & wday < 300)

labels_sh <- historic_sh %>% group_by(WY) %>%
  summarize(wday = max(wday),
            cumul = max(cumul)) %>%
  filter(cumul >= max(realtime_sh$cumul))

spaghetti_sh <- ggplot() +
  geom_line(historic_sh, mapping = aes(x = wday, y = cumul, group = factor(WY)), 
            color = 'grey', linewidth = 1, position = position_dodge(width = 0.01)) +
  geom_line(realtime_sh, mapping = aes(x = wday, y = cumul), color = 'steelblue3', linewidth = 1) +
  scale_x_continuous(breaks = c(31,93,153,214,274), labels = c('Nov', 'Jan', 'Mar', 'May','Jul'),
                     limits = c(0,310)) +
  geom_text(labels_sh, mapping = aes(x = wday+4, y = cumul+4, label = factor(WY))) +
  geom_vline(xintercept = wDay(Sys.Date()), linetype = 'dotted', linewidth = 1) +
  labs(x = 'Date', y = 'Cumulative Loss') +
  theme_bw() +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, unit = 'cm'),
        axis.title.y = element_text(margin=margin(r=15), size = 15),
        axis.title.x = element_text(margin=margin(t=15), size = 15),
        strip.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))
spaghetti_sh
