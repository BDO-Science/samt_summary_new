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

#genetics special graph
lad <- read_excel('CodeFiles/LADModel.xlsx') %>% 
  right_join(all, by = 'Date') %>%
  filter(Date <= '06/01') %>%
  mutate(realDate = mdy(paste0(Date,'/',year(Sys.Date()))))


ggplot(lad, aes(x = realDate)) +
  geom_ribbon(aes(ymin = `Min FL`, ymax = `Max FL`), linewidth = 1, alpha = 0.5, fill = 'grey', color = 'darkgrey') + 
  geom_point(aes(y = ForkLength), color = 'darkorange') +
  labs(x = 'Date', y = 'Fork Length') +
  theme_bw() +
  theme(plot.margin = margin(0.5, 0.5, 0.1, 0.2, unit = 'cm'),
        axis.title.y = element_text(margin = margin(r = 15), size = 15),
        axis.title.x = element_text(margin = margin(t = 10), size = 15),
        strip.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),  # Adjust angle
        axis.text.y = element_text(size = 13),
        legend.position = 'none')
