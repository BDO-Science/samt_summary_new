library(tidyverse)
library(readxl)
CVP <- read_csv('CodeFiles/CVP_genetics.csv')
SWP <- read_csv('CodeFiles/SWP_genetics.csv')
all <- bind_rows(CVP, SWP) %>%
  mutate(day = yday(SampleDateTime)) %>%
  filter(Genetic_Assignment == 'Winter') %>%
  mutate(Date = format(SampleDateTime, '%m/%d'))

today <- yday(Sys.Date())
#####genetics histogram
genetics <- ggplot(all, aes(day, fill = Facility)) +
  geom_histogram(position = 'identity', alpha = 0.5) +
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
ggsave(genetics, file = 'genetics_graph.png', width = 9, height = 7)

###steelhead loss
sh_import <- read.csv('https://www.cbr.washington.edu/sacramento/data/php/rpt/juv_loss_detail.php?sc=1&outputFormat=csv&year=all&species=2%3Aall&dnaOnly=no&age=no')

sh_all <- sh_import %>% filter(Adipose.Clip %in% c('Unclipped', 'Clipped')) %>%
  mutate(Date = as.Date(ymd_hms(Sample.Time))) %>%
  mutate(day = yday(Date))

sh_all_graph <- sh_all %>%
  ggplot() +
  geom_histogram(aes(x = day, fill = Adipose.Clip), color = 'darkgrey') +
  facet_wrap(~Facility) +
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
