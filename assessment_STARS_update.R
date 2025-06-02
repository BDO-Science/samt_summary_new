# Load required libraries
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)

# Get a list of all files that start with "Data.WR" for winter-run or "Data.LF" for late-fall run
file_list <- list.files(pattern = "^Data\\.WR.*\\.csv$", full.names = TRUE)

# Check if there are any files that match the pattern
if (length(file_list) > 0) {
  # Get the most recent file based on modification time
  most_recent_file <- file_list[which.max(file.info(file_list)$mtime)]
  
  # Read the most recent file into a data frame
  data <- read_csv(most_recent_file)
} else {
  # Handle the case where no files are found
  data <- NULL
  message("No files found matching the pattern.")
}

# Convert Date column to Date format
data$Date <- as.Date(data$Date, format="%m/%d/%Y")

# Rename columns for easier reference
colnames(data) <- gsub(" ", "_", colnames(data))

routing <- data |>
  select(contains("Routing"), Date) |>
  filter(month(Date) == 4)  # Filter for March (3)

# Pivot longer
routing_long <- routing |>
  pivot_longer(
    cols = -c(Date, starts_with("Delta")),  # Exclude "Date" and "Delta_*" columns
    names_to = c("Location", ".value"),
    names_pattern = "(?i)Routing_probability_([A-Za-z_]+)_(Est|LCL|UCL)"
  ) 

routing_long$Year <- year(routing_long$Date)

# Plot
ggplot(routing_long, aes(x = Date, y = Est, color = Location)) +
  geom_line() +  # Line for estimate
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Location), alpha = 0.2) +  # 80% Credible Interval
  facet_wrap(~ Year, scales = "free_x") +  # Facet by Year and Location
  theme_minimal() +
  labs(
    title = "Routing Probability Estimates by Year and Location",
    x = NULL,
    y = "Estimate",
    color = "Location",
    fill = "Location"
  ) +
  theme(legend.position = "bottom")

survival <- data |>
  select(contains("Survival"), Date) |>
  filter(month(Date) == 4)  # Filter for March (3)

# Pivot longer
survival_long <- survival |>
  pivot_longer(
    cols = -c(Date, starts_with("Delta")),  # Exclude "Date" and "Delta_*" columns
    names_to = c("Location", ".value"),
    names_pattern = "(?i)Survival_([A-Za-z_]+)_(Est|LCL|UCL)"
  ) 

survival_long$Year <- year(survival_long$Date)

# Plot
ggplot(survival_long, aes(x = Date, y = Est, color = Location)) +
  geom_line() +  # Line for estimate
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Location), alpha = 0.2) +  # 80% Credible Interval
  facet_wrap(~ Year, scales = "free_x") +  # Facet by Year
  theme_minimal() +
  labs(
    title = "Survival Probability Estimates by Year and Location",
    x = NULL,
    y = "Estimate",
    color = "Location",
    fill = "Location"
  ) +
  theme(legend.position = "bottom")

travel_time <- data |>
  select(contains("Median_Travel_Time"), Date) |>
  filter(month(Date) == 4)  # Filter for March (3)

# Pivot longer
travel_time_long <- travel_time |>
  pivot_longer(
    cols = -c(Date, starts_with("Delta")),  # Exclude "Date" and "Delta_*" columns
    names_to = c("Location", ".value"),
    names_pattern = "(?i)Median_Travel_Time_([A-Za-z_]+)_(Est|LCL|UCL)"
  ) 

travel_time_long$Year <- year(travel_time_long$Date)

# Plot
ggplot(travel_time_long, aes(x = Date, y = Est, color = Location)) +
  geom_line() +  # Line for estimate
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Location), alpha = 0.2) +  # 80% Credible Interval
  facet_wrap(~ Year, scales = "free_x") +  # Facet by Year
  theme_minimal() +
  labs(
    title = "Median Travel Time Estimates by Year and Location",
    x = NULL,
    y = "Estimate (days)",
    color = "Location",
    fill = "Location"
  ) +
  theme(legend.position = "bottom")


data_overall_survival <- data |>
  select(
    Date,
    any_of(c("Survival_Overall_Est", "Survival_Overall_LCL_80", "Survival_Overall_UCL_80")),
    any_of(c("Survival_All_Reaches_Est", "Survival_All_Reaches_LCL", "Survival_All_Reaches_UCL"))
  ) |>
  select(Date, everything()) |>  # Ensure Date is first
  # Conditionally rename columns based on their existence
  rename(
    Estimate = any_of(c("Survival_Overall_Est", "Survival_All_Reaches_Est")),
    LCL = any_of(c("Survival_Overall_LCL_80", "Survival_All_Reaches_LCL")),
    UCL = any_of(c("Survival_Overall_UCL_80", "Survival_All_Reaches_UCL"))
  ) |>
  mutate(across(c(Estimate, LCL, UCL), as.numeric)) |>
  filter(Date >= as.Date("2024-10-01"))


# Plot overall survival estimates with 80% credible intervals
STARS_plot <- ggplot(data_overall_survival, aes(x = Date, y = Estimate)) +
  geom_line(color = "black", linewidth = 1) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), fill = "grey", alpha = 0.2) +
  labs(title = NULL,
       x = NULL,
       y = "Survival Estimate") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add a box around the plot
    text = element_text(size = 14, face = "bold")  # Make text bigger and bold
  )+
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b")  # Adjust the breaks and labels
STARS_plot

# Assuming your dataframe is named `data` and has columns:
# 'Date', 'Estimate', 'LCL', and 'UCL'

# Convert Date column to Date format if it's not already
data$Date <- as.Date(data$Date)

survival_range <- data_overall_survival |>
  filter(Date >= as.Date("2025-02-01")) |>
  summarise(
    range_estimate = range(Estimate, na.rm = TRUE),
    range_LCL = range(LCL, na.rm = TRUE),
    range_UCL = range(UCL, na.rm = TRUE)
  )



