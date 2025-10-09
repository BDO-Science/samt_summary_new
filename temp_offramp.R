library(CDECRetrieve)
library(dplyr)
library(ggplot2)
library(lubridate)

#–– PARAMETERS
stations    <- c(MSD = "Mossdale", PPT = "Prisoners Point")
sensor      <- "25"             # water-temp sensor (°F)
dur_code    <- "D"              # daily
start_date  <- "2025-06-01"
threshold   <- 22.2             # °C

#–– 1) PULL & PREP
all_temps <- bind_rows(
  lapply(names(stations), function(code) {
    cdec_query(
      station    = code,
      sensor_num = sensor,
      dur_code   = dur_code,
      start_date = start_date
    ) %>%
      # convert datetime→Date, parameter_value→numeric Temp_F
      mutate(
        Date   = as_date(datetime),                  # <— here!
        Temp_F = as.numeric(parameter_value),
        Station= stations[code],
        Temp_C = (Temp_F - 32) * 5/9,
        exceed = Temp_C > threshold
      ) %>%
      select(Date, Station, Temp_F, Temp_C, exceed)
  })
)

#–– 2) LABELS
label_data <- all_temps %>%
  group_by(Station) %>%
  summarize(
    n_exceed = sum(exceed),
    x        = min(Date),
    y        = min(Temp_C) - 0.5,
    .groups  = "drop"
  ) %>%
  mutate(label = paste0("n = ", n_exceed))

#–– 3) PLOT
ggplot(all_temps, aes(Date, Temp_C, group = Station)) +
  # threshold
  geom_hline(yintercept = threshold, color = "red", linetype = "dashed", size = 0.8) +
  # line and colored points
  geom_line(size = 0.8) +
  geom_point(aes(color = exceed), size = 3) +
  # color scale
  scale_color_manual(
    values = c(`FALSE` = "forestgreen", `TRUE` = "red"),
    labels = c("≤ threshold", "> threshold"),
    name   = NULL
  ) +
  # facet
  facet_wrap(~ Station, ncol = 1) +
  # labels
  geom_text(
    data    = label_data,
    aes(x = x, y = y, label = label),
    hjust   = 0, 
    vjust   = 1,
    fontface = "bold",
    size     = 4
  ) +
  # x axis from June 1 → last date in data
  scale_x_date(
    date_labels = "%b %d",
    date_breaks = "1 week",
    limits     = c(as.Date("2025-06-01"), max(all_temps$Date))
  ) +
  # y axis 15–30 °C
  scale_y_continuous(
    limits = c(15, 30),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Water Temperature (°C)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    text            = element_text(face = "bold"),
    axis.title      = element_text(face = "bold", size = 14),
    axis.text       = element_text(face = "bold", size = 12),
    strip.text      = element_text(face = "bold", size = 14),
    legend.text     = element_text(face = "bold", size = 12),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.margin     = margin(5, 5, 5, 5)
  )