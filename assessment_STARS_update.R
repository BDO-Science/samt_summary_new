# Load required libraries
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)

# Read CSV file
data <- read_csv("Data.WR.Mar.17.2025.10_59_45.csv")

# Convert Date column to Date format
data$Date <- as.Date(data$Date, format="%m/%d/%Y")

# Rename columns for easier reference
colnames(data) <- gsub(" ", "_", colnames(data))

# Select only the overall survival columns
data_overall <- data %>%
  select(Date, 
         Survival_Overall_Est, 
         Survival_Overall_LCL_80, 
         Survival_Overall_UCL_80) %>%
  rename(Estimate = Survival_Overall_Est, 
         LCL = Survival_Overall_LCL_80, 
         UCL = Survival_Overall_UCL_80)

# Convert survival values to numeric
data_overall <- data_overall %>%
  mutate(across(c(Estimate, LCL, UCL), as.numeric))

# Plot overall survival estimates with 80% credible intervals
ggplot(data_overall, aes(x = Date, y = Estimate)) +
  geom_line(color = "black", linewidth = 1) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), fill = "grey", alpha = 0.2) +
  labs(title = NULL,
       x = NULL,
       y = "Survival Estimate") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add a box around the plot
    text = element_text(size = 14, face = "bold")  # Make text bigger and bold
  )


# Assuming your dataframe is named `data` and has columns:
# 'Date', 'Estimate', 'LCL', and 'UCL'

# Convert Date column to Date format if it's not already
data$Date <- as.Date(data$Date)

# Filter for February to the end of March
filtered_data <- data %>%
  filter(Date >= as.Date("2025-02-01"))  # Replace YYYY with the actual year

# Compute the range (min/max) for survival estimates and credible intervals
range_estimate <- range(filtered_data$Survival_Overall_Est, na.rm = TRUE)
range_LCL <- range(filtered_data$Survival_Overall_LCL_80, na.rm = TRUE)
range_UCL <- range(filtered_data$Survival_Overall_UCL_80, na.rm = TRUE)

# Print results
cat("Survival Estimate Range:", range_estimate, "\n")
cat("80% LCL Range:", range_LCL, "\n")
cat("80% UCL Range:", range_UCL, "\n")



