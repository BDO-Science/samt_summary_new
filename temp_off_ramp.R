library(tidyverse)
library(CDECRetrieve)

startDate <- paste0(year(Sys.Date()),'-06-01')

mossdale <- cdec_query('MSD',
                       '25',
                       'H',
                       startDate,
                       Sys.Date()) %>%
  select(station = location_id, 
         date = datetime, 
         temp = parameter_value) %>%
  mutate(date = as.Date(date)) %>%
  group_by(station, date) %>%
  summarize(temp = mean(temp, na.rm = TRUE)) %>%
  mutate(temp = (temp -32)/1.8)

prisoners <- cdec_query('PPT',
                        '25',
                        'H',
                        startDate,
                        Sys.Date()) %>%
  select(station = location_id, 
         date = datetime, 
         temp = parameter_value) %>%
  mutate(date = as.Date(date)) %>%
  group_by(station, date) %>%
  summarize(temp = mean(temp, na.rm = TRUE)) %>%
  mutate(temp = (temp -32)/1.8)

all_temps <- bind_rows(mossdale, prisoners) %>%
  mutate(trigger = if_else(temp >= 22.2, 'YES', 'NO')) %>%
  na.omit() %>%
  mutate(station = factor(station, levels = c('MSD', 'PPT'), labels = c('Mossdale', 'Prisoners Point')))

yes_triggers <- all_temps %>%
  filter(trigger == "YES") %>%
  group_by(station) %>%
  summarize(n_yes = n(), .groups = "drop") %>%
  complete(station = unique(all_temps$station), fill = list(n_yes = 0))

temp_offramp <- ggplot(all_temps, mapping = aes(x = date, y = temp)) +
  geom_line(linewidth = 1, alpha = 0.5) +
  geom_point(aes(fill = trigger), size = 4, shape = 21) +
  geom_label(yes_triggers, mapping = aes(x = max(all_temps$date) - 1, 
                                         y = 19, 
                                         label = paste0('n = ',n_yes)),
             size = 5) +
  geom_hline(yintercept = 22.2, linetype = 'dashed', color = 'red', linewidth = 1) +
  labs(x = 'Date', y = 'Water Temperature (C)') +
  ylim(c(18,25)) +
  scale_fill_manual(values = c('#33cc33', 'red')) +
  facet_wrap(~station, ncol = 1) +
  theme_bw() +
  theme(legend.position = 'none',
        plot.margin = margin(0.2, 0.5, 0.2, 0.2, unit = 'cm'),
        axis.title.y = element_text(margin = margin(r = 15), size = 15),
        axis.title.x = element_text(margin = margin(t = 15), size = 15),
        strip.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))
temp_offramp
ggsave(temp_offramp, file = 'outputs/temp_offramp.png', height = 7, width = 9)
