#01_MossdaleCDFW_LAD_assignment.R----

# Nick Bertrand
# Start Date: Mon Apr 15 10:19:12 2024

#About----
#Project:Realtime operations  
#Purpose: This is to address the lack of LAD assignment in the data files provided by CDFW 
#when they are doing the mossdale trawl. 

#Libraries ----
library(tidyverse)
library(readr)
library(readxl)

# Set working directory ----
#set to location of root object to highest tier directory
getwd()
root <- "C:/Users/nbertrand/OneDrive - DOI/Desktop/Bertrand/GitHub/SAMTLossSummary"
setwd(root)
getwd()
#these root object use directories 
data_root<-file.path(root,"Data")
code_root <- file.path(root,"R_scripts")
#table_output_root <- file.path(root,"Table_Output")
#viz_output_root <- file.path(root,"Viz_Output")

#Import Data----
msdata <- read_excel(file.path(data_root,"MossdaleTrawl/Weekly Update_2024.xlsx"), sheet = "CHN FL")
View(msdata)

ladmodel <- read_excel(file.path(data_root,"MossdaleTrawl/LADModel.xlsx"), sheet = "Sheet2")
View(ladmodel)


#Assigment----
#The lad model values for all runs are joined to the the mossdale trawl date and length data 
msdata_join <-left_join(msdata,ladmodel, by = "Date") 
#view(msdata_join)

#the Mossdale observed FL is then compared to the ranges of each run
#LAD model assumes that there is no overlap in LAD
#any value outside of the run ranges should generate an "Error" value
assigned <- msdata_join %>% 
  mutate(Run = case_when(`Forklength (mm)` >= WRMinFL & `Forklength (mm)` <= WRMaxFL | 
                          `Forklength (mm)` >= WRaltMinFL & `Forklength (mm)` <= WRaltMaxFL ~ "WR",
                        `Forklength (mm)` >= SRMinFL & `Forklength (mm)` <= SRMaxFL | 
                          `Forklength (mm)` >= SRaltMinFL & `Forklength (mm)` <= SRaltMaxFL ~ "SR",
                        `Forklength (mm)` >= FRMinFL & `Forklength (mm)` <= FRMaxFL | 
                          `Forklength (mm)` >= FRaltMinFL & `Forklength (mm)` <= FRaltMaxFL ~ "FR",
                        `Forklength (mm)` >= LFRMinFL & `Forklength (mm)` <= LFRMaxFL | 
                          `Forklength (mm)` >= LFRaltMinFL & `Forklength (mm)` <= LFRaltMaxFL ~ "LFR",
                        TRUE ~ "Error"))
#view(assigned)

#Summarize Run Counts----
#Summarize the run assignments to count of the previous 7 days based on the date this script is run.

#Create start and end dates
startdate <- Sys.Date() - 7
startdate
enddate <- Sys.Date()
enddate 

#Summarize observations
summaryMSD <- assigned %>% filter(Date >= startdate & Date<= enddate) %>% 
  group_by(Run) %>% 
  count()

view(summaryMSD)


