# Install packages if you don't have them
# install.packages("readr")
# install.packages("dplyr")

library(readr)
library(dplyr)

# Define the input filename
# Note: You mentioned an Excel file, but the uploaded file is a CSV.
# This code assumes you are reading the CSV.
input_file <- "CodeFiles/rubias_assignment_SaMT-WY2025-20250724.2.csv"

# Read the data
all_data <- read_csv(input_file)

# Filter for CVP and write to CSV
cvp_data <- all_data %>%
  filter(Facility == "CVP")

write_csv(cvp_data, "CodeFiles/CVP_genetics.csv")

# Filter for SWP and write to CSV
swp_data <- all_data %>%
  filter(Facility == "SWP")

write_csv(swp_data, "CodeFiles/SWP_genetics.csv")

print("Files CVP_genetics.csv and SWP_genetics.csv have been created.")