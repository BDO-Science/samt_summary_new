library(tidyverse)
library(busdater)
library(flextable)
library(officer)
library(zoo)

wy <- get_fy(Sys.Date(), opt_fy_start = '07-01')  #pull the water year based on BY designation in LTO docs
jpe <- 98893

#set source folder, destination folder, and read in filenames of newest files
source_folder <- 'SalvageFiles/' #where the current files for the summary script to process live
file_names <- list.files(path = source_folder, 
                         pattern = "\\.csv$", full.names = TRUE) #read csv files in source folder

#pull in loss data
salmon_raw <- read_csv(max(file_names[grep('salmon', file_names, ignore.case = TRUE)])) %>%
  mutate(Date = case_when(
    grepl('/', SampleDate) ~ mdy(SampleDate),
    grepl('-', SampleDate) ~ ymd(SampleDate),
    TRUE ~ NA_Date_  # Ensures missing values are handled properly
  ))

steelhead_raw <- read_csv(max(file_names[grep('steelhead', file_names, ignore.case = TRUE)]))

############winter-run weekly threshold
#pull in weekly threshold data and manipulate for easier joining with loss table
wr_thresholds <- read_csv('CodeFiles/weeklyThresholds.csv') %>% #pulling in weekly distributed loss thresholds
  mutate(StartDate = dmy(paste0(StartDate,'-',wy))) %>% #converting to date format with current water year
  mutate(EndDate = dmy(paste0(EndDate,'-',wy))) %>% #ditto
  rowwise() %>%
  mutate(Date = list(seq.Date(StartDate, EndDate, by = "day"))) %>%
  unnest(Date) %>%
  select(Date, HistoricPresent) %>%
  mutate(threshold = ((jpe*.005)*.5)*HistoricPresent)

#processing salmon loss table for WR
temp <- salmon_raw %>%
  filter(AdClip == 'N' & (.[[8]] == 'W' | .[[9]] == 'W')) %>%
  select(14,Size_Race = 8,DNA_Race = 9,13) %>%
  mutate(confirmed = case_when(Size_Race == 'W' & DNA_Race != 'W' ~ 'NO',
                               DNA_Race == 'W' ~ 'YES',
                               Size_Race == 'W' & is.na(DNA_Race) ~ 'PARTIAL')) %>%
  filter(confirmed != 'NO') %>%
  group_by(Date, confirmed) %>%
  summarize(loss = sum(LOSS))

wr_weekly <- data.frame(Date = seq(as.Date('2024-12-01'), as.Date('2025-06-30'), 1)) %>%
  left_join(temp, by = 'Date') %>%
  left_join(wr_thresholds, by = 'Date') %>%
  replace(is.na(.), 0) %>%
  mutate(threshold = round(threshold, 2)) %>%
  mutate(sum_7D_loss = rollsum(loss, k = 7, fill = 'left', align = 'right')) %>%
  filter(Date <= Sys.Date() & Date >= Sys.Date()-6) %>%
  mutate(triggered = if_else(sum_7D_loss < threshold, 'No', 'Yes')) %>%
  mutate(loss = if_else(confirmed == 'PARTIAL', paste0(as.character(loss),'*'), as.character(loss)))

wr_table <- wr_weekly %>%
  mutate(Date = format(Date, "%b %d")) %>%
  select(Date, 'Winter-run Daily Loss' = 3, 'Winter-run 7-day rolling sum loss' = 6, 
         'Winter-run Daily Threshold' = 5, 'Winter-run Daily Trigger' = 7)


############steelhead weekly threshold
salvage <- steelhead_raw %>% select(1,3,9) %>%
  pivot_longer(!Date, names_to = 'facility', values_to = 'salvage') %>%
  mutate(Date = mdy(Date)) %>%
  group_by(Date) %>%
  summarize(salvage = sum(salvage, na.rm = TRUE))

SH_weekly <- data.frame(Date = seq(as.Date('2024-12-01'), as.Date('2025-06-30'), 1)) %>%
  left_join(salvage, by = 'Date') %>%
  replace(is.na(.), 0) %>%
  mutate(sum_7D_salvage = rollsum(salvage, k = 7, fill = 'left', align = 'right')) %>%
  filter(Date <= Sys.Date() & Date >= Sys.Date()-6) %>%
  mutate(triggered = if_else(sum_7D_salvage < 120, 'No', 'Yes'))


thick_border <- fp_border(color = "black", width = 2)

weekly_table <- SH_weekly %>%
  mutate(Date = format(Date, "%b %d")) %>%
  select(Date, 'Steelhead Daily Salvage' = 2, 'Steelhead 7-day rolling sum salvage' = 3, 'Steelhead Daily Trigger' = 4) %>%
  left_join(wr_table, by = 'Date') %>%
  flextable() %>%
  vline() %>%
  hline() %>%
  border_outer() %>%
  align(align = 'center', part = 'all') %>%
  width(width = 1, unit = "in") %>%
  font(fontname = 'Segoe UI', part = "body") %>%
  font(fontname = 'Segoe UI Semibold', part = "header") %>%
  fontsize(size = 12, part = "header") %>%
  fontsize(size = 10, part = "body") %>%
  border(j = 4, border.right = thick_border, part = "all") %>%
  set_caption(caption = 'Summary of weekly salvage of steelhead and loss of winter-run to inform weekly distributed loss thresholds. Steelhead thresholds are triggered when the 7-day rolling sum of expanded salvage exceeds 120 fish.  Winter-run thresholds are triggered when 7-day rolling sum of loss exceeds the daily threshold based on historic proportion of winter-run present in the delta.*All or a portion of Winter-run loss has not been genetically confirmed',
              align_with_table = FALSE)
weekly_table

############cumulative loss of steelhead and winter-run
temp_wr <- salmon_raw %>%
  filter(AdClip == 'N') %>%
  filter(.[[8]] == 'W' & (.[[9]] == 'W'| is.na(.[[9]]))) %>%
  mutate(facility = if_else(FACILITY == 2, 'CVP', 'SWP')) %>%
  select(Date, facility, loss=LOSS) %>%
  mutate(species = 'Winter-run')


cumulative_loss <- steelhead_raw %>% select(1,4,10) %>%
  pivot_longer(!Date, names_to = 'facility', values_to = 'loss') %>%
  mutate(Date = mdy(Date)) %>%
  mutate(facility = if_else(facility == 'Loss_SWPN', 'SWP', 'CVP')) %>%
  mutate(species = 'Steelhead') %>%
  bind_rows(temp_wr) %>%
  na.omit() %>%
  group_by(Date, facility, species) %>%
  summarize(loss = sum(loss, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(species) %>%
  mutate(cumul_loss = cumsum(loss))

last_winter_run_loss <- cumulative_loss %>%
  filter(species == "Winter-run") %>%
  arrange(desc(Date)) %>%
  slice(1) %>%
  pull(cumul_loss)

# Create the new row
new_row <- data.frame(
  Date = max(cumulative_loss$Date),
  facility = "CVP",  # Assuming CVP since past Winter-run entries are CVP
  species = "Winter-run",
  loss = 0,  # No new loss, just carrying forward
  cumul_loss = last_winter_run_loss  # Maintain the last cumulative loss
)

# Bind the new row to the dataframe
cumulative_loss <- cumulative_loss %>%
  bind_rows(new_row) %>%
  arrange(Date, species)

lossmax <- cumulative_loss %>% group_by(species) %>% summarize(Date = max(Date), cumul_loss = max(cumul_loss)) %>%
  mutate(threshold = if_else(species == 'Steelhead', (cumul_loss/3000)*100, (cumul_loss/(jpe*.005))*100)) %>%
  mutate(threshold = round(threshold, 1))

thresholds <- data.frame(species = c('Steelhead', 'Steelhead', 'Steelhead', 
                                    'Winter-run', 'Winter-run', 'Winter-run'),
                         threshold = c('perc100', 'perc75', 'perc50',
                                       'perc100', 'perc75', 'perc50'),
                        thresholds = c(3000, 3000*.75, 3000*.5, 
                                      jpe*0.005, jpe*0.005*.75, jpe*0.005*.5)) %>%
  mutate(thresholds = round(thresholds, 0)) %>%
  pivot_wider(names_from = threshold, values_from = thresholds) %>%
  mutate(label = paste0('Thresholds:\n',perc100, ' (100%)\n',
                        perc75, ' (75%)\n',
                        perc50, ' (50%)')) %>%
  left_join(select(lossmax, species,cumul_loss), by = 'species')

#annual_thresholds <- expand.grid(threshold_perc = c(.5, .75, 1), species = c('steelhead', 'winter-run')) %>%
  #mutate(threshold = if_else(species == 'steelhead', 3000 * threshold_perc, jpe*.005*threshold_perc))

loss_graph <- ggplot(cumulative_loss) +
  geom_line(aes(x = Date, y = cumul_loss, color = "Cumulative Loss"), linetype = 'dashed', linewidth = 1) +
  geom_col(aes(x = Date, y = loss, fill = facility), position = 'dodge') +
  geom_text(lossmax, mapping = aes(x = Date + 1, y = cumul_loss*1.05, 
                                    label = paste0(cumul_loss,' (',threshold,'%)')), fontface = 'bold', size = 4) +
  geom_text(data = thresholds, aes(x = min(cumulative_loss$Date)-15, 
                                   y = cumul_loss*.90, 
                                   label = label), 
            hjust = 0, fontface = "italic", size = 4) + 
  facet_wrap(~species, scales = 'free_y') +
  scale_color_manual(name = "", values = 'black') +
  scale_x_date(date_breaks = '6 weeks', date_labels = "%b %d", limits = c(as.Date('2024-12-01'), Sys.Date() + 14)) +
  guides(color = guide_legend(override.aes = list(linetype = 'dashed'))) +
  labs(y = 'Loss', fill = 'Facility') +
  theme_bw() +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, unit = 'cm'),
        axis.title.y = element_text(margin=margin(r=15), size = 15),
        axis.title.x = element_text(margin=margin(t=15), size = 15),
        strip.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = 'bottom')
loss_graph

ggsave(plot = loss_graph, file = 'outputs/loss_summary.png', width = 9, height = 6)
save_as_docx(weekly_table, loss_graph, path = paste0('outputs/salvage-summary_',Sys.Date(),'.docx'))
