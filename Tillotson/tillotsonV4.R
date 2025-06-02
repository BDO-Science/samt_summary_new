#This script generates quantile predictions of steelhead and winter-run incidental take based
#on quantile regression forest models publised in Tillotson et al. (2022), but
#updated to be trained with 1999-2021 data. This model version also excludes
#the Delta beach seine index since so few steelhead are caught and the variable 
#was of essentially zero importance.

#Load packages
require(randomForest)
require(quantregForest)
require(tidyverse)
require(stringr)
require(lubridate)
require(zoo)
require(caret)
library(CDECRetrieve)
#Source model setup and training script, requires "ITMData.rda" to be in same folder
source("Tillotson/Steelhead_Model_Setup.r")
source("Tillotson/WR_Model_Setup.r")

waterDay <- readRDS('Tillotson/waterDay.rds') #function for parsing water day from Date

#finding file paths to pull in
files <- list.files(path = 'SalvageFiles/')
shfile <- max(grep('steelhead', files, value = TRUE, ignore.case = TRUE))
#tillfile <- max(grep('tillotson', files, value = TRUE, ignore.case = TRUE))
wrfile <- max(grep('salmon', files, value = TRUE, ignore.case = TRUE))

#pull in most recent Geir steelhead file
steelhead <- read.csv(paste0('SalvageFiles/',shfile)) %>%
  select(1,4,10) %>% #select interested columns
  mutate(Date = as.Date(Date, "%m/%d/%Y")) %>% #change Date to date format
  gather(key = 'Facility', value = 'Loss', 2:3) %>% #more data cleaning
  group_by(Date) %>%
  summarize(Loss = sum(Loss, na.rm = TRUE)) %>% #sum Loss for both facilities
  mutate(wyWeek = week(as.Date(waterDay(Date), origin = as.Date("2024-01-01")))) #calculate week of water year

#pull in most recent Geir winter-run file
wr <- read.csv(paste0('SalvageFiles/', wrfile)) %>%
  mutate(Date = case_when(
    grepl('/', SampleDate) ~ mdy(SampleDate),
    grepl('-', SampleDate) ~ ymd(SampleDate),
    TRUE ~ NA_Date_  # Ensures missing values are handled properly
  )) %>%
  filter(AdClip == 'N' & Size_Run == 'W') %>%
  mutate(wyWeek = week(as.Date(waterDay(Date), origin = as.Date("2024-01-01"))))

#########inputs for Tillotson model
wy_week = week(as.Date(waterDay(Sys.Date()), origin = as.Date("2024-01-01"))) #parse current water week

#exports
cvp <- cdec_query('TRP', '70', 'D', '2025-01-01', Sys.Date()) %>% #combined CVP and SWP exports
  filter(!is.na(parameter_value)) %>% select(Date = 3, TRP = 5)
exports <- cdec_query('HRO', '70', 'D', '2025-01-01', Sys.Date()) %>% 
  filter(!is.na(parameter_value)) %>% select(Date = 3, HRO = 5) %>%
  left_join(cvp, by = 'Date') %>% 
  mutate(Flow = HRO + TRP) %>%
  mutate(Date = as.Date(Date)) %>%
  mutate(wyWeek = week(as.Date((waterDay(Date)-1), origin = as.Date("2024-01-01")))) %>%
  filter(Date == max(Date)) %>%
  pull(Flow)

#freeport and vernalis
sac <- cdec_query('FPT', '20', 'D', "2025-01-01", Sys.Date()) %>% #Sacramento flows at Freeport
  filter(!is.na(parameter_value)) %>% select(Date = 3, FPT = 5) %>%
  mutate(Date = as.Date(Date)) %>%
  mutate(wyWeek = week(as.Date((waterDay(Date)-1), origin = as.Date("2024-01-01")))) %>%
  filter(Date == max(Date)) %>%
  pull(FPT)

sjr <- cdec_query('VNS', '41', 'D', "2025-01-01", Sys.Date()) %>% #San Joaquin flows at Vernalis
  filter(!is.na(parameter_value)) %>% select(Date = 3, VNS = 5) %>%
  mutate(Date = as.Date(Date)) %>%
  mutate(wyWeek = week(as.Date((waterDay(Date)-1), origin = as.Date("2024-01-01")))) %>%
  filter(Date == max(Date)) %>%
  pull(VNS)

#species loss
stlhd_loss <- steelhead %>% #calculates the previous weeks weekly steelhead loss for model
  filter(wyWeek == (wy_week-1)) %>% 
  summarise(total_loss = sum(Loss, na.rm = TRUE)) %>% 
  pull(total_loss)

wr_loss <- wr %>% #calculates the previous weeks weekly wr loss for model
  filter(wyWeek == (wy_week-1)) %>%
  summarise(total_loss = sum(LOSS, na.rm = TRUE)) %>%
  pull(total_loss)
  
#temp data for MAL
mal <- cdec_query('MAL', '25', 'H', "2025-01-01", Sys.Date()) %>% #pulls in the temp data at MAL
  mutate(Date = as.Date(datetime)) %>%
  filter(!is.na(parameter_value)) %>% filter(Date == max(Date, na.rm = TRUE)) %>%
  summarize(temp = mean((parameter_value-32)*(5/9))) %>% 
  pull(temp) #isolates most recent water temp

#pulling in covariates and adding species to Tillotson inputs
data <- expand.grid(OMR = c(-2500, -3000, -3500, -4000, -4500, -5000), species = c('Steelhead', 'Winter-run')) %>%
  mutate(fpt = sac,
         vrn = sjr,
         export = exports)


predictionsList <- list() #list to store model predictions

for (i in 1:nrow(data)){ #loop through each row of tillotson inputs
  omr <- data$OMR[i] #pull OMR value in current row
  sac <- data$fpt[i] #pull Sac flows in current row
  sjr <- data$vrn[i] #pull San Joaquin flows in current row
  species <- data$species[i] #pull species in current row
  export <- data$export[i] #pull and combine exports in current row
  
  #sets certain parameters based on species
  if (species == 'Steelhead') {
    speciesName <- 'stlhd_loss.pw' #set name for parameter
    speciesLoss <- stlhd_loss  # Assuming 'stlhd_loss' is defined elsewhere
    model <- Stlhd_Simple_Combined[[2]]  # Use Steelhead model
  } else if (species == 'Winter-run') {
    speciesName <- 'winter.pw' #set name for parameter
    speciesLoss <- wr_loss  # Assuming 'wr_loss' is defined elsewhere
    model <- WR_Simple_Combined[[2]]  # Use Winter-run model
  } else {
    # Handle other cases (if any)
    speciesName <- NA
    speciesLoss <- NA
    model <- NULL
  }
  
  #builds the data.frame for model predictions
  NewData_Combined <- data.frame(wy_week = wy_week,
                                 mal_temp = mal, 
                                 precip = 1,
                                 OMR= omr,
                                 sac = sac,
                                 sjr = sjr,
                                 dcc = "closed",
                                 daily_exports = export,
                                 speciesName = speciesLoss)
  colnames(NewData_Combined)[9] <- speciesName #renames columns for model to find
  Predictions <- data.frame(predict(model, #making predictions and cleaning output up
                                    newdata = NewData_Combined, 
                                    what=c(.25,.5,.75))) %>%
    mutate(OMR = omr, Exports = export, 
           'San Joaquin Flow' = sjr, 'Sacramento Flow' = sac, Species = species) %>%
    select(8, 4, 5, 6, 7,qtl25 = 1, median = 2, qtl75 = 3)
  predictionsList[[i]] <- Predictions #storing predictions in list
}
all <- bind_rows(predictionsList) %>% #combining predictions into single table
  mutate('Median Daily Loss' = round(median/7,2)) %>%
  rename('25th Percentile Weekly Loss' = 'qtl25', #renaming columns
         'Median Weekly Loss' = 'median', 
         '75th Percentile Weekly Loss' = 'qtl75')

#graph
tillotson_graph <- bind_rows(predictionsList) %>% 
  filter(Species == 'Winter-run') %>%
  ggplot(aes(x = factor(OMR))) +
  geom_crossbar(aes(y = median, ymin = qtl25, ymax = qtl75),  fill = '#CC9900') +
  labs(title = 'Salvage and Loss Predictions at different OMRI values', x = 'OMRI', y = 'Predicted Weekly Loss', fill = 'Median Predicted Loss') +
  theme_bw() +
  theme(plot.margin = ggplot2::margin(0.5, 0.5, 0.1, 0.2, unit = 'cm'),
        title = element_text(size = 14),
        axis.title.y = element_text(margin = ggplot2::margin(r = 15), size = 15),
        axis.title.x = element_text(margin = ggplot2::margin(t = 10), size = 15),
        strip.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = 'bottom')
tillotson_graph

tillotson_graph_sh <- bind_rows(predictionsList) %>% 
  filter(Species == 'Steelhead') %>%
  ggplot(aes(x = factor(OMR))) +
  geom_crossbar(aes(y = median, ymin = qtl25, ymax = qtl75),  fill = '#CC9900') +
  labs(title = 'Salvage and Loss Predictions at different OMRI values', x = 'OMRI', y = 'Predicted Weekly Loss', fill = 'Median Predicted Loss') +
  theme_bw() +
  theme(plot.margin = ggplot2::margin(0.5, 0.5, 0.1, 0.2, unit = 'cm'),
        title = element_text(size = 14),
        axis.title.y = element_text(margin = ggplot2::margin(r = 15), size = 15),
        axis.title.x = element_text(margin = ggplot2::margin(t = 10), size = 15),
        strip.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = 'bottom')
tillotson_graph_sh
#write SH to csv file for input into assessment
all %>% filter(Species == 'Steelhead') %>% 
  select(-1) %>%
  write.csv(paste0('Tillotson/outputs/SHmodel_predictions_',Sys.Date(),'.csv'), row.names = FALSE)

#write WR to csv file for input into assessment
all %>% filter(Species == 'Winter-run') %>% 
  select(-1) %>%
  write.csv(paste0('Tillotson/outputs/WRmodel_predictions_',Sys.Date(),'.csv'), row.names = FALSE)


ggsave(tillotson_graph, file = paste0('Tillotson/outputs/WR_tillotson_',Sys.Date(),'.png'), width = 9, height = 6)
ggsave(tillotson_graph_sh, file = paste0('Tillotson/outputs/SH_tillotson_',Sys.Date(),'.png'), width = 9, height = 6)
