library(tidyverse)
library(CDECRetrieve)

stations <- c('JER', 'HOL', 'BAC', 'BET')

data_list <- lapply(stations, function(station){
  data <- cdec_query(station, '100', 'D', '2025-09-01', Sys.Date())
})

all_data <- bind_rows(data_list) %>%
  mutate(location_id = factor(location_id, levels = c('JER', 'BET', 'HOL', 'BAC'),
         labels = c('Jersey Point', 'Bethel Island', 'Holland Cut', 'Bacon Island')))

thresholds <- data.frame(location_id = unique(all_data$location_id), threshold = c(1800,800,700,1000))

graph <- ggplot() +
  geom_line(all_data, mapping = aes(x = datetime, y = parameter_value), linewidth = 1, color = 'steelblue3') +
  geom_hline(thresholds, mapping = aes(yintercept = threshold), 
             linewidth = 1.5, linetype = 'dashed', color = '#888888') +
  labs(x = 'Date', y = 'Conductivity') +
  facet_wrap(~location_id, scale = 'free_y') +
  theme_bw()
graph
