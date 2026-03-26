library(tidyverse)
library(ggplot2)

getwd()
dataset <- read.csv("exoplanets.csv")

# ==============================
# Data selection and cleaning
# ==============================

# Select relevant variables for analysis
data_clean <- dataset %>%
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
# Distribution of planetary radius
# ==============================

# Remove missing values and apply log transformation
radius <- na.omit(data_clean$koi_prad)
log_radius <- log10(radius)

# Plot histogram with median line
ggplot(data.frame(log_radius), aes(x = log_radius)) +
  geom_histogram(bins = 40, fill = "lightblue", color = "white") +
  geom_vline(xintercept = median(log_radius), color = "red", linewidth = 0.8) +
  labs(
    title = "Distribution of Exoplanet Radii (Log Scale)",
    x = "log10(Planetary Radius in Earth Radii)",
    y = "Frequency"
  ) +
  theme_minimal()

# ==============================
# Distribution of orbital period
# ==============================

# Remove missing values and apply log transformation
period <- na.omit(data_clean$koi_period)
log_period <- log10(period)

# Plot histogram
ggplot(data.frame(log_period), aes(x = log_period)) +
  geom_histogram(bins = 40, fill = "lightblue", color = "white") +
  geom_vline(xintercept = median(log_period), color = "red", linewidth = 0.8) +
  labs(
    title = "Distribution of Orbital Period (Log Scale)",
    x = "log10(Orbital Period [days])",
    y = "Frequency"
  ) +
  theme_minimal()

# ==============================
# Distribution of transit duration
# ==============================

# Remove missing values and apply log transformation
duration <- na.omit(data_clean$koi_duration)
log_duration <- log10(duration)

# Plot histogram
ggplot(data.frame(log_duration), aes(x = log_duration)) +
  geom_histogram(bins = 40, fill = "lightblue", color = "white") +
  geom_vline(xintercept = median(log_duration), color = "red", linewidth = 0.8) +
  labs(
    title = "Distribution of Transit Duration (Log Scale)",
    x = "log10(Transit Duration [hours])",
    y = "Frequency"
  ) +
  theme_minimal()

# ==============================
# Distribution of transit depth
# ==============================

# Remove missing values and apply log transformation
depth <- na.omit(data_clean$koi_depth)
log_depth <- log10(depth)

# Plot histogram
ggplot(data.frame(log_depth), aes(x = log_depth)) +
  geom_histogram(bins = 40, fill = "lightblue", color = "white") +
  geom_vline(xintercept = median(log_depth), color = "red", linewidth = 0.8) +
  labs(
    title = "Distribution of Transit Depth (Log Scale)",
    x = "log10(Transit Depth [ppm])",
    y = "Frequency"
  ) +
  theme_minimal()

# ==============================
# Distribution of stellar temperature
# ==============================

# Remove missing values (no log transformation applied)
star_temp <- na.omit(data_clean$koi_steff)

# Plot histogram
ggplot(data.frame(star_temp), aes(x = star_temp)) +
  geom_histogram(bins = 40, fill = "lightblue", color = "white") +
  geom_vline(xintercept = median(star_temp), color = "red", linewidth = 0.8) +
  labs(
    title = "Distribution of Stellar Effective Temperature",
    x = "Stellar Effective Temperature [K]",
    y = "Frequency"
  ) +
  theme_minimal()

# ==============================
# Distribution of stellar radius
# ==============================

# Remove missing values and apply log transformation
srad <- na.omit(data_clean$koi_srad)
log_srad <- log10(srad)

# Plot histogram
ggplot(data.frame(log_srad), aes(x = log_srad)) +
  geom_histogram(bins = 40, fill = "lightblue", color = "white") +
  geom_vline(xintercept = median(log_srad), color = "red", linewidth = 0.8) +
  labs(
    title = "Distribution of Stellar Radius (Log Scale)",
    x = "log10(Stellar Radius [Solar radii])",
    y = "Frequency"
  ) +
  theme_minimal()