library(ggplot2)
library(dplyr)
library(tidyr)

data <- read.csv("exoplanets.csv")

# ==============================
# Data cleaning and selection
# ==============================

# Select relevant variables for analysis
data_clean <- data %>%
  select(
    koi_disposition,
    koi_pdisposition,
    koi_period,
    koi_time0bk,
    koi_impact,
    koi_duration,
    koi_depth,
    koi_prad,
    koi_teq,
    koi_insol,
    koi_model_snr,
    koi_steff,
    koi_slogg,
    koi_srad,
    ra,
    dec,
    koi_kepmag,
    koi_fpflag_nt,
    koi_fpflag_ss,
    koi_fpflag_co,
    koi_fpflag_ec
  )

# ==============================
# Outlier analysis using boxplots
# ==============================

# Select continuous variables of interest
selected_data <- data_clean %>%
  select(
    koi_prad,
    koi_period,
    koi_impact,
    koi_duration,
    koi_depth,
    koi_teq,
    koi_insol
  ) %>%
  na.omit()  # remove missing values for cleaner plots

# Convert data to long format for faceted plotting
selected_data_long <- selected_data %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  )

# ==============================
# Boxplot visualization
# ==============================

# Each panel shows the distribution of one variable
# Boxplots help identify outliers and spread of the data
ggplot(selected_data_long, aes(x = "", y = value)) +
  geom_boxplot(fill = "lightblue") +
  facet_wrap(~variable, scales = "free") +
  labs(
    title = "Boxplots of Selected Variables",
    subtitle = "Visualization of distributions and potential outliers",
    x = "",
    y = "Value"
  ) +
  theme_minimal()