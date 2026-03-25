library(ggplot2)
library(dplyr)

data <- read.csv("exoplanets.csv")
data

#clean
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

#6. Outliers


library(dplyr)
library(tidyr)
library(ggplot2)

selected_data <- data_clean %>%
  select(koi_prad, koi_period, koi_impact, koi_duration, koi_depth, koi_teq, koi_insol) %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "value")

ggplot(selected_data, aes(x = "", y = value)) +
  geom_boxplot(fill = "lightblue") +
  facet_wrap(~variable, scales = "free") +
  labs(
    title = "Boxplots of Selected Variables",
    x = "",
    y = "Value"
  ) +
  theme_minimal()

