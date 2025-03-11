library(tidyverse)
library(readxl)
library(flextable)
library(officer)

today <- Sys.Date()
start <- format(today-7, '%m/%d')
end <- format(today, '%m/%d')

####################setting source folder and finding all filenames in the folder
source_folder <- 'Monitoring/' #where the current files for the summary script to process live
destination_folder <- 'Monitoring/Processed' #where past processed files live
file_names <- list.files(path = source_folder, 
                         pattern = "\\.xlsx$", full.names = TRUE)

####################Delta Sampling Files
DJFMP <- file_names[grep('datcall',file_names, ignore.case = TRUE)] #look for datcall file
if (length(DJFMP) > 0){
DJFMPsheets <- excel_sheets(DJFMP)
DJFMPsheets <- DJFMPsheets[DJFMPsheets %in% c('Beach Seines', 'Sacramento Trawls', #everthing but EDSM
                                              'Mossdale Trawls', 'Chipps Island Trawls')]

DJFMPlist <- list()

for(i in DJFMPsheets){
  temp <- read_excel(DJFMP, sheet = i) %>% 
    select('Date', 'Species', 'ADCLIP' = 'Mark', 'Run' = 'RaceByLength', 'Catch', 'WaterTemp', 'Turbidity') %>%
    mutate(Sample = i)
  DJFMPlist[[i]] <- temp
}

delta <- bind_rows(DJFMPlist) 

delta <- delta %>% filter(Species %in% c('CHN', 'RBT')) %>%
  select(1,2,5,4,3,8,6,7) %>%
  mutate(ADCLIP = if_else(ADCLIP == 'None', 'N', 'Y'),
         Run = gsub('n/p', NA, Run),
         Flow = NA,
         Species = as.character(Species))
} else {
  delta <- data.frame()
}

EDSM <- read_excel(DJFMP, sheet = 'EDSM Trawls') %>% 
  select('Date', 'Species' = 'OrganismCode', 'ADCLIP' = 'MarkCode', 
         'Run' = 'RaceByLength', 'Catch' = 'SumOfSumOfCount', 'WaterTemp' = 'Temp', 'Turbidity' = 'Turb') %>%
  mutate(Sample = 'EDSM') %>%
  select(1,2,5,4,3,8,6,7) %>%
  filter(Species %in% c('CHN', 'RBT')) %>%
  mutate(ADCLIP = if_else(ADCLIP == 'AdClipped', 'Y', 'N'),
         Run = gsub('n/p', NA, Run),
         Flow = NA,
         Species = as.character(Species))

delta <- bind_rows(delta, EDSM)

####################Mossdale Sampling Files
MSfile <- file_names[grep('weekly', file_names, ignore.case = TRUE)] #look for mossdale file
if (length(MSfile) > 0){
msCountdata <- read_excel(MSfile, sheet = 'Daily Count', skip = 9) %>%
  select(Date = 1, CHN = 3, RBT = 7) %>% 
  gather(key = 'Species', value = 'Catch', 2:3) %>%
  mutate(Run = NA, ADCLIP = if_else(Species == 'CHN', 'Y', 'N'))

msFLdata <- read_excel(MSfile, sheet = "CHN FL")
ladmodel <- read_excel("ReferenceCode/LADModel.xlsx", sheet = "Sheet2")

msdataLAD <-left_join(msFLdata,ladmodel, by = "Date") %>%  #adds run assignment to data
  mutate(Run = case_when(`Forklength (mm)` >= WRMinFL & `Forklength (mm)` <= WRMaxFL | 
                           `Forklength (mm)` >= WRaltMinFL & `Forklength (mm)` <= WRaltMaxFL ~ "Winter",
                         `Forklength (mm)` >= SRMinFL & `Forklength (mm)` <= SRMaxFL | 
                           `Forklength (mm)` >= SRaltMinFL & `Forklength (mm)` <= SRaltMaxFL ~ "Spring",
                         `Forklength (mm)` >= FRMinFL & `Forklength (mm)` <= FRMaxFL | 
                           `Forklength (mm)` >= FRaltMinFL & `Forklength (mm)` <= FRaltMaxFL ~ "Fall",
                         `Forklength (mm)` >= LFRMinFL & `Forklength (mm)` <= LFRMaxFL | 
                           `Forklength (mm)` >= LFRaltMinFL & `Forklength (mm)` <= LFRaltMaxFL ~ "LateFall",
                         TRUE ~ "Error")) %>% 
  group_by(Date, Run) %>%
  summarize(Catch = n()) %>%
  mutate(ADCLIP = 'N', Species = 'CHN') %>%
  select(1,5,3,2,4)

mossdale <- bind_rows(msCountdata, msdataLAD) %>%
  mutate(Sample = 'Mossdale', WaterTemp = NA, Turbidity = NA, Flow = NA)
} else {
  mossdale <- data.frame()
}

####################Knights Landing RST
KLfile <- file_names[grep('knights', file_names, ignore.case = TRUE)] #look for knights landing file
if (length(KLfile) > 0)  {
  KL <- read_excel(KLfile, skip = 8)
KLWild <- KL %>% 
  select(Date = 3, Fall = 19, Spring = 20, Winter = 21, LateFall = 22, WaterTemp = 15, Turbidity = 16, Flow = 14) %>%
  gather(key = 'Run', value = 'Catch', 2:5) %>% 
  mutate(ADCLIP = 'N', 
         Date = as.Date(Date),
         Species = 'CHN')
KLclipped <- KL %>%
  select(Date = 3, Fall = 25, Spring = 26, Winter = 27, LateFall = 28,WaterTemp = 15, Turbidity = 16, Flow = 14) %>%
  gather(key = 'Run', value = 'Catch', 2:5) %>% 
  mutate(ADCLIP = 'Y', 
         Date = as.Date(Date),
         Species = 'CHN')
KLSH <- KL %>%
  select(Date = 3, N = 29, Y = 30, WaterTemp = 15, Turbidity = 16, Flow = 14) %>%
  gather(key = 'ADCLIP', value = 'Catch', 2:3) %>%
  mutate(Run = NA, 
         Date = as.Date(Date),
         Species = 'RBT')
KLall <- bind_rows(KLWild, KLclipped, KLSH) %>%
  filter(Date >= "2024-08-01")

KLall <- KLall %>% 
  mutate(Sample = 'Knights Landing RST') %>%
  select(1,8,6,5,7,9,2,3,4)
} else {
  KLall <- data.frame()
}

####################Lower Sac RST
LSfile <- file_names[grep('Lower', file_names, ignore.case = TRUE)] #look for Lower sac RST file
if (length(LSfile) > 0)  {
  LS <- read_excel(LSfile, skip = 8)
LSWild <- LS %>% 
  select(Date = 3, Fall = 17, Spring = 18, Winter = 19, LateFall = 20, WaterTemp = 13, Turbidity = 14, Flow = 12) %>%
  gather(key = 'Run', value = 'Catch', 2:5) %>% 
  mutate(ADCLIP = 'N', 
         Date = as.Date(Date),
         Species = 'CHN')
LSclipped <- LS %>%
  select(Date = 3, Fall = 21, Spring = 22, Winter = 23, LateFall = 24, WaterTemp = 13, Turbidity = 14, Flow = 12) %>%
  gather(key = 'Run', value = 'Catch', 2:5) %>% 
  mutate(ADCLIP = 'Y', 
         Date = as.Date(Date),
         Species = 'CHN')
LSSH <-  LS %>%
  select(Date = 3, N = 25, Y = 26, WaterTemp = 13, Turbidity = 14, Flow = 12) %>%
  gather(key = 'ADCLIP', value = 'Catch', 2:3) %>%
  mutate(Run = NA, 
         Date = as.Date(Date),
         Species = 'RBT')
LSall <- bind_rows(LSWild, LSclipped, LSSH) %>%
  mutate(Sample = 'Lower Sacramento RST') %>%
  select(1,8,6,5,7,9,2,3,4)
} else {
  LSall <- data.frame()
}

####################Lower Feather River
LFRfile <- file_names[grep('LFR', file_names, ignore.case = TRUE)] #look for Lower Feather file
if (length(LFRfile) > 0)  {
  LFR <- read_excel(LFRfile, skip = 9)
  
  LFRWild <- LFR %>%
    select(Date = 3, Fall = 18, Spring = 19, Winter = 20, LateFall = 21, WaterTemp = 14, Turbidity = 15, Flow = 13) %>%
    gather(key = 'Run', value = 'Catch', 2:5) %>% 
    mutate(ADCLIP = 'N', 
           Date = as.Date(Date),
           Species = 'CHN')
  
  LFRclipped <- LFR %>%
    select(Date = 3, Fall = 23, Spring = 24, Winter = 25, LateFall = 26, WaterTemp = 14, Turbidity = 15, Flow = 13) %>%
    gather(key = 'Run', value = 'Catch', 2:5) %>% 
    mutate(ADCLIP = 'Y', 
           Date = as.Date(Date),
           Species = 'CHN')
  
  LFRSH <-  LFR %>%
    select(Date = 3, N = 27, Y = 28, WaterTemp = 14, Turbidity = 15, Flow = 13) %>%
    gather(key = 'ADCLIP', value = 'Catch', 2:3) %>%
    mutate(Run = NA, 
           Date = as.Date(Date),
           Species = 'RBT')
  
  LFRall <- bind_rows(LFRWild, LFRclipped, LFRSH) %>% 
    mutate(Sample = 'Lower Feather RST') %>%
    select(1, 8, 6, 5, 7, 9, 2, 3, 4) %>%
    mutate(WaterTemp = as.numeric(WaterTemp),
           Turbidity = as.numeric(Turbidity)) %>%
    filter(Date <= today)
} else {
  LFRall <- data.frame()
}

####################Upper feather RSTs
FRfile <- file_names[grep('feather', file_names, ignore.case = TRUE)] #look for upper Feather file
if (length(FRfile) > 0)  {
  FRherr <- read_excel(FRfile, skip = 5, sheet = 'Herringer') %>%
    select(Date = 3, Fall = 21, Spring = 22, LateFall = 24, 
           ADCLIP = 23, SH = 25, SHClip = 26,
           WaterTemp = 17, Turbidity = 18, Flow = 16) %>%
    gather(key = 'Run', value = 'Catch', 2:7) %>%
    mutate(ADCLIP = if_else(Run %in% c('ADCLIP', 'SHClip'), 'Y', 'N')) %>%
    mutate(Species = if_else(Run %in% c('SH', 'SHClip'), 'RBT', 'CHN')) %>%
    mutate(Run = if_else(Run %in% c('SH', 'SHClip', 'ADCLIP'), NA, Run)) %>%
    mutate(Sample = 'Feather River (Herringer)')
  
  FReye <- read_excel(FRfile, skip = 5, sheet = 'Eye Side Channel') %>%
    select(Date = 3, Fall = 14, Spring = 15, LateFall = 16, SH = 17, SHClip = 18,
           WaterTemp = 10, Turbidity = 11, Flow = 9) %>%
    gather(key = 'Run', value = 'Catch', 2:6) %>%
    mutate(ADCLIP = if_else(Run == 'SHClip', 'Y', 'N')) %>%
    mutate(Species = if_else(Run %in% c('SH', 'SHClip'), 'RBT', 'CHN')) %>%
    mutate(Run = if_else(Run %in% c('SH', 'SHClip'), NA, Run)) %>%
    mutate(Sample = 'Feather River (Eye Side)')
  FRall <- bind_rows(FReye, FRherr) %>%
    mutate(Date = as.Date(Date)) %>%
    select(1,8,6,5,7,9,2,3,4)
} else {
  FRall <- data.frame()
}

####################Yuba RSTs
YubaFile <- file_names[grep('yuba', file_names, ignore.case = TRUE)] #look for Yuba file
if (length(YubaFile) > 0)  {
  Yuba <- read_excel(YubaFile, skip = 5, sheet = 'Hallwood') %>%
    select(Date = 3, Fall = 28, Spring = 29, LateFall = 30, SH = 31, SHClip = 32,
           WaterTemp = 24, Turbidity = 25, Flow = 23) %>%
    gather(key = 'Run', value = 'Catch', 2:6) %>%
    mutate(ADCLIP = if_else(Run == 'SHClip', 'Y', 'N')) %>%
    mutate(Species = if_else(Run %in% c('SH', 'SHClip'), 'RBT', 'CHN')) %>%
    mutate(Run = if_else(Run %in% c('SH', 'SHClip'), NA, Run)) %>%
    mutate(Sample = 'Yuba') %>%
    select(1,8,6,5,7,9,2,3,4) %>%
    mutate(Date = as.Date(Date))
} else {
  Yuba <- data.frame()
}

####################Butte Creek RST
ButteFile <- file_names[grep('butte', file_names, ignore.case = TRUE)] #look for Butte File
if (length(ButteFile) > 0)  {
  Butte <- read_excel(ButteFile, skip = 8) %>%
    select(Date = 3, Spring = 14, SH = 15, 
           WaterTemp = 10, Turbidity = 11, Flow = 9) %>%
    gather(key = 'Run', value = 'Catch', 2:3) %>%
    mutate(ADCLIP = 'N') %>%
    mutate(Species = if_else(Run == 'SH', 'RBT', 'CHN')) %>%
    mutate(Run = if_else(Run =='SH', NA, Run)) %>%
    mutate(Sample = 'Butte') %>%
    select(1,8,6,5,7,9,2,3,4) %>%
    mutate(Date = as.Date(Date)) %>%
    mutate(WaterTemp = as.numeric(WaterTemp), Turbidity = as.numeric(Turbidity))
} else {
  Butte <- data.frame()
}

####################Tisdale RST
TDfile <- file_names[grep('Tisdale', file_names, ignore.case = TRUE)] #look for Tisdale file
if (length(TDfile) > 0)  {
  TD <- read_excel(TDfile, skip = 8) 
  TDWild <- TD %>% 
    select(Date = 3, Fall = 17, Spring = 18, Winter = 19, LateFall = 20, WaterTemp = 13, Turbidity = 14, Flow = 12) %>%
    gather(key = 'Run', value = 'Catch', 2:5) %>% 
    mutate(ADCLIP = 'N', 
           Date = as.Date(Date),
           Species = 'CHN',
           WaterTemp = as.numeric(WaterTemp))
  TDclipped <- TD %>%
    select(Date = 3, Fall = 22, Spring = 23, Winter = 24, LateFall = 25,WaterTemp = 13, Turbidity = 14, Flow = 12) %>%
    gather(key = 'Run', value = 'Catch', 2:5) %>% 
    mutate(ADCLIP = 'Y', 
           Date = as.Date(Date),
           Species = 'CHN',
           WaterTemp = as.numeric(WaterTemp))
  TDSH <- TD %>%
    select(Date = 3, N = 26, Y = 27, WaterTemp = 15, Turbidity = 16, Flow = 14) %>%
    gather(key = 'ADCLIP', value = 'Catch', 2:3) %>%
    mutate(Run = NA, 
           Date = as.Date(Date),
           Species = 'RBT',
           WaterTemp = as.numeric(WaterTemp))
  TDall <- bind_rows(TDWild, TDclipped, TDSH) 
  
  TDall <- TDall %>% 
    mutate(Sample = 'Tisdale RST') %>%
    select(1,8,6,5,7,9,2,3,4)
} else {
  TDall <- data.frame()
}
####################Bind all datasets together
dataframes <- list(delta, KLall, LSall, LFRall, mossdale, FRall, Yuba, Butte, TDall) #put all dataframes in a list
AllFish <- dplyr::bind_rows(dataframes)

####################create base table for joining with file names for some conditional stuff
baseTable <- data.frame(Location = c('Beach Seines', 'Sacramento Trawls', 'Mossdale Trawls', 
                                     'Chipps Island Trawls','Knights Landing RST', 'Lower Sacramento RST', 
                                     'Lower Feather RST', 'Feather River (Herringer)','Feather River (Eye Side)',
                                     'Yuba', 'Butte', 'Tisdale RST', 'EDSM'),
                        Type = c('Delta', 'Delta', 'Delta', 'Delta', 
                                 'RST', 'RST', 'RST', 'RST', 'RST', 'RST', 'RST', 'RST', 'Delta'),
                        file = c(if(length(DJFMP) == 0) NA else DJFMP, 
                                 if(length(DJFMP) == 0) NA else DJFMP, 
                                 if(length(DJFMP) == 0) NA else DJFMP, 
                                 if(length(DJFMP) == 0) NA else DJFMP, 
                                 if(length(KLfile) == 0) NA else KLfile, 
                                 if(length(LSfile) == 0) NA else LSfile, 
                                 if(length(LFRfile) == 0) NA else LFRfile, 
                                 if(length(FRfile) == 0) NA else FRfile, 
                                 if(length(FRfile) == 0) NA else FRfile, 
                                 if(length(YubaFile) == 0) NA else YubaFile, 
                                 if(length(ButteFile) == 0) NA else ButteFile, 
                                 if(length(TDfile) == 0) NA else TDfile,
                                 if(length(DJFMP) == 0) NA else DJFMP)) %>%
  mutate(char = nchar(file)) %>%
  mutate(FileExists = if_else(char > 0, 1, 0))

####################Summarize Data
Summary <- AllFish %>% 
  filter(Species != 'RBT') %>%
  filter(Date <= today & Date >= today - 6) %>%
  mutate(Run = if_else(ADCLIP == 'Y', 'Ad-clip', Run)) %>%
  group_by(Sample, Run) %>%
  summarize(Catch = sum(Catch, na.rm = TRUE)) %>% 
  pivot_wider(
    names_from = Run,
    values_from = Catch) %>%
  replace(is.na(.), 0) %>%
  right_join(baseTable, by = c('Sample' = 'Location')) %>%
  mutate(FileExists = if_else(is.na(FileExists), 0, FileExists)) %>%
  mutate(across(where(is.numeric), ~ if_else(FileExists == 1 & is.na(.), 0, .))) %>%
  mutate(Dates = paste0(start, '-', end))

####################putting data into a flextable
RSTtable <- Summary %>% 
  filter(Type == 'RST') %>%
  select(Location = 1,11,6,3,2,5,4) %>%
  mutate(across(everything(), ~ as.character(.))) %>%
  mutate(across(everything(), ~ replace_na(., "na"))) %>%
  flextable() %>% 
  align(j = c(1:7), align = 'center', part = 'all') %>%
  fontsize(size = 12, part = 'body') %>%
  fontsize(size = 12, part = 'header') %>%
  font(fontname = 'Calibri', part = 'all') %>%
  bold(bold = TRUE, part = 'header') %>%
  vline()%>%
  hline() %>%
  border_outer() %>%
  width(width = c(1.75, 1, 0.75, 0.75, 0.75, 0.75, 0.75)) %>%
  height_all(height = 0.20, part = 'body') %>%
  height(height = 0.75, part = 'header') %>%
  hrule(rule = 'atleast', part = 'body')
RSTtable

Deltatable <- Summary %>% 
  filter(Type == 'Delta') %>%
  select(Location = 1,11,6,3,2,5,4) %>%
  mutate(across(everything(), ~ as.character(.))) %>%
  mutate(across(everything(), ~ replace_na(., "na"))) %>%
  flextable() %>%
  align(j = c(1:7), align = 'center', part = 'all') %>%
  fontsize(size = 12, part = 'body') %>%
  fontsize(size = 12, part = 'header') %>%
  font(fontname = 'Calibri', part = 'all') %>%
  bold(bold = TRUE, part = 'header') %>%
  vline()%>%
  hline() %>%
  border_outer() %>%
  width(width = c(1.75, 1, 0.75, 0.75, 0.75, 0.75, 0.75)) %>%
  height_all(height = 0.20, part = 'body') %>%
  height(height = 0.75, part = 'header') %>%
  hrule(rule = 'atleast', part = 'body')
Deltatable

####################save relevant files
save_as_docx(RSTtable, Deltatable, path = paste0('FishTable/catch_',today,'.docx')) #saves the flextable to a word document
write.csv(Summary, paste0('FishTable/fishTable',today,'.csv'), row.names = FALSE) #saves the summary table to csv
write.csv(AllFish, paste0('FishTable/OldData/allfish_',today,'.csv'), row.names = FALSE) #saves the all fish table that appends data onto previous data

####################putting data into a flextable
new_file_path <- file.path(destination_folder, basename(file_names))
file.rename(file_names, new_file_path)

