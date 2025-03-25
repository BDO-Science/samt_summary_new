library(tidyverse)
library(busdater)
library(flextable)
library(officer)
library(zoo)
# Load necessary library
library(patchwork)

wy <- get_fy(Sys.Date(), opt_fy_start = '07-01')  #pull the water year based on BY designation in LTO docs
jpe <- 98893
jpe_hatch <- 135342
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
  ), DNA_Run = ifelse(
    Date == as.Date("2025-03-18") & 
      format(Time, "%H:%M:%S") == "11:00:00" & 
      LGT == 107, 
    "UW", 
    DNA_Run  # Ensure to use the correct column name
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
  filter(AdClip == 'N' & 
           (.[[8]] == 'W' | DNA_Run == 'W' | DNA_Run == 'UW')) %>%
  select(Date, Size_Run = 8, DNA_Run, LOSS) %>%
  mutate(confirmed = case_when(
    Size_Run == 'W' & DNA_Run != 'W' ~ 'NO',
    DNA_Run == 'W' ~ 'YES',
    Size_Run == 'W' & is.na(DNA_Run) ~ 'PARTIAL'
  )) %>%
  filter(confirmed != 'NO') %>%
  group_by(Date, confirmed) %>%
  summarize(loss = sum(LOSS), .groups = 'drop')  # Adding .groups = 'drop' to avoid warnings



# Count occurrences of "W" in Size_Run and NA in DNA_Run and CWT_Run for the entire dataset
no_wr_salvage <- salmon_raw %>%
  summarize(
    DNA_Run_W = sum(DNA_Run == "W", na.rm = TRUE),  # Count "W" in DNA_Run
    CWT_Run_W = sum(CWT_Run == "W", na.rm = TRUE),   # Count "W" in CWT_Run
    unconfirmed_W = sum(Size_Run == "W" & is.na(DNA_Run) & is.na(CWT_Run) & AdClip == 'N')  # Count "W" in Size_Run where DNA_Run and CWT_Run are NA
  )

# Get the last week of data
last_week_data <- salmon_raw %>%
  filter(Date >= (Sys.Date() - 7))

# Count occurrences of "W" for the last week
last_week_summary <- last_week_data %>%
  summarize(
    DNA_Run_W = sum(DNA_Run == "W", na.rm = TRUE),  # Count "W" in DNA_Run
    CWT_Run_W = sum(CWT_Run == "W", na.rm = TRUE),   # Count "W" in CWT_Run
    unconfirmed_W = sum(Size_Run == "W" & is.na(DNA_Run) & is.na(CWT_Run) & AdClip == 'N')  # Count "W" in Size_Run where DNA_Run and CWT_Run are NA
  )

# Combine results
final_summary <- list(
  overall_summary = no_wr_salvage,
  last_week_summary = last_week_summary
)

print(final_summary)
  

wr_weekly <- data.frame(Date = seq(as.Date('2024-12-01'), as.Date('2025-06-30'), 1)) %>%
  left_join(temp, by = 'Date') %>%
  left_join(wr_thresholds, by = 'Date') %>%
  replace(is.na(.), 0) %>%
  mutate(threshold = round(threshold, 2)) %>%
  mutate(sum_7D_loss = rollsum(loss, k = 7, fill = NA, align = 'right')) %>%
  filter(Date <= Sys.Date() & Date >= Sys.Date() - 6) %>%
  mutate(triggered = if_else(sum_7D_loss < threshold, 'No', 'Yes')) %>%
  mutate(loss = if_else(confirmed == 'PARTIAL', paste0(as.character(loss), '*'), as.character(loss))) %>%
  group_by(Date) %>%
  summarise(
    confirmed = paste(unique(confirmed[confirmed != 0]), collapse = ", "),
    loss = if_else(any(str_detect(loss, "\\*")), 
                   paste0(sum(as.numeric(str_remove(loss, "\\*")), na.rm = TRUE), "*"), 
                   as.character(sum(as.numeric(loss), na.rm = TRUE))),
    HistoricPresent = first(HistoricPresent),
    threshold = first(threshold),
    sum_7D_loss = max(sum_7D_loss),
    triggered = first(triggered)
  ) %>%
  ungroup()

wr_table <- wr_weekly %>%
  mutate(Date = format(Date, "%b %d")) %>%
  select(Date, 'Winter-run Daily Loss' = 3, 'Winter-run 7-day rolling sum loss' = 6, 
         'Winter-run Daily Threshold' = 5, 'Winter-run Daily Trigger' = 7)

#this is for the entire water year up until this point
wr_weekly_WY <- data.frame(Date = seq(as.Date('2024-12-01'), as.Date(Sys.Date()), 1)) %>%
  left_join(temp, by = 'Date') %>%
  left_join(wr_thresholds, by = 'Date') %>%
  replace(is.na(.), 0) %>%
  mutate(threshold = round(threshold, 2)) %>%
  mutate(sum_7D_loss = rollsum(loss, k = 7, fill = NA, align = 'right')) %>%
  mutate(triggered = if_else(sum_7D_loss < threshold, 'No', 'Yes')) %>%
  mutate(loss = if_else(confirmed == 'PARTIAL', paste0(as.character(loss), '*'), as.character(loss))) 

# Extend the threshold column for one week beyond Sys.Date()
threshold_extension <- data.frame(Date = seq(as.Date(Sys.Date()) + 1, as.Date(Sys.Date()) + 14, 1)) %>%
  left_join(wr_thresholds, by = 'Date') %>%
  mutate(threshold = round(threshold, 2))

# Combine the two data frames
wr_weekly_WY <- bind_rows(wr_weekly_WY %>% filter(Date <= Sys.Date()), threshold_extension) %>%
  arrange(Date) %>%
  group_by(Date) %>%
  mutate(threshold = if_else(is.na(threshold), lag(threshold, order_by = Date), threshold)) %>%
  ungroup()

############steelhead weekly threshold
salvage <- steelhead_raw %>% select(1,3,9) %>%
  pivot_longer(!Date, names_to = 'facility', values_to = 'salvage') %>%
  mutate(Date = mdy(Date)) %>%
  group_by(Date) %>%
  summarize(salvage = sum(salvage, na.rm = TRUE))

SH_weekly <- data.frame(Date = seq(as.Date('2024-12-01'), as.Date('2025-06-30'), 1)) %>%
  left_join(salvage, by = 'Date') %>%
  replace(is.na(.), 0) %>%
  mutate(sum_7D_salvage = rollsum(salvage, k = 7, fill = NA, align = 'right')) %>%
  filter(Date <= Sys.Date() & Date >= Sys.Date() - 6) %>%
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

#weekly steelhead thresholds for the entire water year
SH_weekly_WY <- data.frame(Date = seq(as.Date('2024-12-01'), as.Date(Sys.Date()), 1)) %>%
  left_join(salvage, by = 'Date') %>%
  replace(is.na(.), 0) %>%
  mutate(sum_7D_salvage = rollsum(salvage, k = 7, fill = NA, align = 'right')) %>%
  mutate(triggered = if_else(sum_7D_salvage < 120, 'No', 'Yes'))

############cumulative loss of steelhead and winter-run
temp_wr <- salmon_raw %>%
  filter(AdClip == 'N') %>%
  filter(.[[8]] == 'W' & (.[[9]] == 'W'| .[[9]] == 'UW' | is.na(.[[9]]))) %>%
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
  left_join(select(lossmax, species,cumul_loss), by = 'species') %>%
  mutate(max_thresh = if_else(species == 'Steelhead', 120, max(wr_weekly_WY$threshold))) %>%
  rowwise() %>%
  mutate(max_loss = max(cumul_loss, max_thresh))

#annual_thresholds <- expand.grid(threshold_perc = c(.5, .75, 1), species = c('steelhead', 'winter-run')) %>%
#mutate(threshold = if_else(species == 'steelhead', 3000 * threshold_perc, jpe*.005*threshold_perc))

####need to bind salmon and steelhead together to play nicely with figure code as is
# Add a species column to each data frame
SH_weekly_WY <- SH_weekly_WY %>%
  mutate(species = "Steelhead") %>%
  rename(sum_7D = sum_7D_salvage)

wr_weekly_WY <- wr_weekly_WY %>%
  mutate(species = "Winter-run") %>%
  rename(sum_7D = sum_7D_loss)

# Combine both data frames
combined_weekly <- bind_rows(SH_weekly_WY, wr_weekly_WY)

# Create a horizontal line data frame for Steelhead
hline_data <- data.frame(yintercept = 120, species = "Steelhead")

loss_graph <- ggplot(cumulative_loss) +
  geom_line(aes(x = Date, y = cumul_loss, color = "Annual Cumulative Loss"), linetype = 'dashed', linewidth = 1) +
  #geom_col(aes(x = Date, y = loss, fill = facility), position = 'dodge') +
  geom_label(lossmax, mapping = aes(x = Date + 1, y = cumul_loss, 
                                   label = paste0(cumul_loss, ' (', threshold, '%)')), 
            fontface = 'bold', size = 4, nudge_x = -10) +
  geom_text(data = thresholds, aes(x = min(cumulative_loss$Date) - 15, 
                                   y = max_loss * .90, 
                                   label = label), 
           hjust = 0, fontface = "italic", size = 4) + 
  geom_line(data = combined_weekly, aes(x = Date, y = sum_7D, color = "7 d Rolling Sum"), alpha = 0.9, size = 1) +  # Add lines for SH and WR data
  geom_line(data = combined_weekly, aes(x = Date, y = threshold, color = "Weekly Distributed Loss Threshold"), linetype = "dotted", size = 1) + #winter-run threshold
  # Add a condition to only show the line on the Steelhead facet
  geom_hline(data = hline_data, aes(yintercept = yintercept), 
             color = "red", linetype = "dotted", size = 1) +
  geom_col(aes(x = Date, y = loss, fill = facility), position = 'dodge') +
  facet_wrap(~species, scales = 'free_y') +
  scale_color_manual(name = "", values = c('#999999', 'blue', 'red')) +  # Adjust colors as needed
  scale_x_date(date_breaks = '3 weeks', date_labels = "%b %d", limits = c(as.Date('2024-12-01'), Sys.Date() + 14)) +
  guides(color = guide_legend(override.aes = list(linetype = c('solid', 'dashed', 'dotted')))) +
  scale_fill_manual(values = c('#0066cc', 'orange3')) +
  labs(y = 'Loss', fill = 'Facility', x = NULL) +
  theme_bw() +
  theme(plot.margin = margin(0.2, 0.5, 0.2, 0.2, unit = 'cm'),
        axis.title.y = element_text(margin = margin(r = 15), size = 15),
        axis.title.x = element_text(margin = margin(t = 15), size = 15),
        strip.text = element_text(size = 13),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1),  # Adjust angle
        axis.text.y = element_text(size = 13),
        legend.position = 'bottom')

# Display the graph
print(loss_graph)

#####summarizing hatchery winter-run loss
wr_hatch <- salmon_raw %>%
  filter(.[[10]] == 'W') %>%
  mutate(facility = if_else(FACILITY == 2, 'CVP', 'SWP')) %>%
  select(Date, facility, loss=LOSS) %>%
  mutate(species = 'Winter-run') %>%
  mutate(cumul_loss = cumsum(loss))

cumul_max <- max(wr_hatch$cumul_loss)
max_date <- max(wr_hatch$Date)
thresholds <- data.frame(threshold = c('100% Threshold', '75% Threshold', '50% Threshold'),
                         value = c(jpe_hatch*.0012, jpe_hatch*.0012*.75, jpe_hatch*.0012*.5))

wr_hatch_cumul_graph <- ggplot() +
  geom_col(wr_hatch, mapping = aes(x = Date, y = loss, fill = facility), position = 'dodge') +
  geom_line(wr_hatch, mapping = aes(x = Date, y = cumul_loss), linewidth = 1, color = 'blue', linetype = "dashed") +
  geom_hline(thresholds, mapping = aes(yintercept = value), linetype = 'dotted', linewidth = 1, color = 'red') +
  geom_vline(mapping = aes(xintercept = Sys.Date(), linetype = 'Today', color = 'Today'), linewidth = 1) +
  geom_text(thresholds, mapping = aes(x = as.Date('2025-04-05'), y = value + 10, label = threshold, fontface = 'bold')) +
  geom_label(mapping = aes(x = max_date + 1, y = cumul_max *1.1, 
                           label = paste0(cumul_max, ' (', round((cumul_max/(jpe_hatch*.0012))*100,1), '%)')), 
             fontface = 'bold', size = 4) +
  scale_x_date(date_breaks = '3 weeks', date_labels = "%b %d", limits = c(as.Date('2025-03-01'), as.Date('2025-06-01'))) +
  scale_linetype_manual(name = "", values = c("Today" = "dashed")) +
  scale_color_manual(name = "", values = c('Today' = 'darkgrey')) +
  labs(y = 'Loss', x = NULL, subtitle = 'Cumulative Loss of Hatchery Winter-run Chinook Salmon') +
  scale_fill_manual(values = c('#0066cc', 'orange3')) +
  theme_bw() +
  theme(plot.margin = margin(0.5, 0.5, 0.1, 0.2, unit = 'cm'),
        axis.title.y = element_text(margin = margin(r = 15), size = 15),
        axis.title.x = element_text(margin = margin(t = 10), size = 15),
        strip.text = element_text(size = 13),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1),  # Adjust angle
        axis.text.y = element_text(size = 13),
        legend.position = 'bottom',
        legend.title = element_blank())
wr_hatch_cumul_graph

ggsave(plot = loss_graph, file = 'outputs/loss_summary.png', width = 9, height = 6)
ggsave(wr_hatch_cumul_graph, file = 'outputs/hatchery_summary.png', width = 9, height = 6)
save_as_docx(weekly_table, loss_graph, path = paste0('outputs/salvage-summary_',Sys.Date(),'.docx'))

# Combine the two graphs
combined_graph <- loss_graph / wr_hatch_cumul_graph +
  plot_annotation(title = "WY2025 Loss for Steelhead and Winter-run Chinook Salmon") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

# Display the combined graph
print(combined_graph)
