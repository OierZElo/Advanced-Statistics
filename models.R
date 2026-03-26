library(tidyverse)
library(ggplot2)
library(dplyr)

getwd()
data <- read.csv("exoplanets.csv")

# Select variables used in the project
df <- data %>%
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
# Simple Linear Regression
# ==============================

# Prepare data
df_lin_simple <- df %>%
  select(koi_prad, koi_period) %>%
  na.omit()

# ==============================
# Outlier removal
# ==============================

# Outliers for koi_prad
Q1_prad <- quantile(df_lin_simple$koi_prad, 0.25, na.rm = TRUE)
Q3_prad <- quantile(df_lin_simple$koi_prad, 0.75, na.rm = TRUE)
IQR_value_prad <- Q3_prad - Q1_prad
lower_prad <- Q1_prad - 1.5 * IQR_value_prad
upper_prad <- Q3_prad + 1.5 * IQR_value_prad

# Outliers for koi_period
Q1_period <- quantile(df_lin_simple$koi_period, 0.25, na.rm = TRUE)
Q3_period <- quantile(df_lin_simple$koi_period, 0.75, na.rm = TRUE)
IQR_value_period <- Q3_period - Q1_period
lower_period <- Q1_period - 1.5 * IQR_value_period
upper_period <- Q3_period + 1.5 * IQR_value_period

# Remove outliers
df_no_outliers <- df_lin_simple %>%
  filter(koi_prad >= lower_prad & koi_prad <= upper_prad) %>%
  filter(koi_period >= lower_period & koi_period <= upper_period)

# ==============================
# Scatter plot
# ==============================

ggplot(df_no_outliers, aes(x = log10(koi_period), y = koi_prad)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", col = "red") +
  labs(
    title = "Planet Radius vs Orbital Period",
    x = "Orbital Period (days) *log10*",
    y = "Planet Radius (Earth radii)"
  ) +
  theme_minimal()

# ==============================
# Linear regression model
# ==============================

# Fit simple linear regression model
model_lin_simple <- lm(koi_prad ~ koi_period, data = df_no_outliers)

# Model summary
summary(model_lin_simple)

# Base R diagnostic plots
plot(model_lin_simple)

# ==============================
# Diagnostic data
# ==============================

diag_df_simple <- data.frame(
  fitted = fitted(model_lin_simple),
  residuals = resid(model_lin_simple),
  std_residuals = rstandard(model_lin_simple),
  leverage = hatvalues(model_lin_simple)
)

# ==============================
# Diagnostic plots
# ==============================

# Residuals vs fitted values
ggplot(diag_df_simple, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(se = FALSE, method = "loess", color = "blue") +
  labs(
    title = "Residuals vs Fitted Values",
    x = "Fitted values",
    y = "Residuals"
  ) +
  theme_minimal()

# Normal Q-Q plot
ggplot(diag_df_simple, aes(sample = std_residuals)) +
  stat_qq(alpha = 0.3) +
  stat_qq_line(color = "red") +
  labs(
    title = "Normal Q-Q Plot",
    x = "Theoretical Quantiles",
    y = "Standardized Residuals"
  ) +
  theme_minimal()

# Scale-location plot
ggplot(diag_df_simple, aes(x = fitted, y = sqrt(abs(std_residuals)))) +
  geom_point(alpha = 0.3) +
  geom_smooth(se = FALSE, method = "loess", color = "blue") +
  labs(
    title = "Scale-Location Plot",
    x = "Fitted values",
    y = expression(sqrt("|Standardized residuals|"))
  ) +
  theme_minimal()

# Residuals vs leverage
ggplot(diag_df_simple, aes(x = leverage, y = std_residuals)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Residuals vs Leverage",
    x = "Leverage",
    y = "Standardized Residuals"
  ) +
  theme_minimal()

# ==============================
# Multiple Linear Regression
# ==============================

# Prepare data
df_lin_multi <- df %>%
  select(koi_teq, koi_insol, koi_steff, koi_period) %>%
  na.omit()

# ==============================
# Outlier removal function
# ==============================

remove_outliers <- function(data, column) {
  Q1 <- quantile(data[[column]], 0.25, na.rm = TRUE)
  Q3 <- quantile(data[[column]], 0.75, na.rm = TRUE)
  IQR_value <- Q3 - Q1
  
  lower <- Q1 - 1.5 * IQR_value
  upper <- Q3 + 1.5 * IQR_value
  
  data %>% filter(data[[column]] >= lower & data[[column]] <= upper)
}

# Remove outliers variable by variable
df_lin_multi <- remove_outliers(df_lin_multi, "koi_teq")
df_lin_multi <- remove_outliers(df_lin_multi, "koi_insol")
df_lin_multi <- remove_outliers(df_lin_multi, "koi_steff")
df_lin_multi <- remove_outliers(df_lin_multi, "koi_period")

# Fit multiple linear regression model
model_lin_multi <- lm(koi_teq ~ koi_insol + koi_steff + koi_period,
                      data = df_lin_multi)

# Model summary
summary(model_lin_multi)

# Predicted values
df_lin_multi$pred <- predict(model_lin_multi)

# ==============================
# Observed vs predicted plot
# ==============================

ggplot(df_lin_multi, aes(x = pred, y = koi_teq)) +
  geom_point(alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  labs(
    title = "Observed vs Predicted Equilibrium Temperature",
    x = "Predicted Temperature",
    y = "Observed Temperature"
  ) +
  theme_minimal()

# ==============================
# Diagnostic data
# ==============================

diag_df <- data.frame(
  fitted = fitted(model_lin_multi),
  residuals = resid(model_lin_multi),
  std_residuals = rstandard(model_lin_multi),
  leverage = hatvalues(model_lin_multi)
)

# ==============================
# Diagnostic plots
# ==============================

# Residuals vs fitted values
ggplot(diag_df, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(se = FALSE, method = "loess", color = "blue") +
  labs(
    title = "Residuals vs Fitted",
    x = "Fitted values",
    y = "Residuals"
  ) +
  theme_minimal()

# Normal Q-Q plot
ggplot(diag_df, aes(sample = std_residuals)) +
  stat_qq(alpha = 0.3) +
  stat_qq_line(color = "red") +
  labs(
    title = "Normal Q-Q Plot",
    x = "Theoretical Quantiles",
    y = "Standardized Residuals"
  ) +
  theme_minimal()

# Scale-location plot
ggplot(diag_df, aes(x = fitted, y = sqrt(abs(std_residuals)))) +
  geom_point(alpha = 0.3) +
  geom_smooth(se = FALSE, method = "loess", color = "blue") +
  labs(
    title = "Scale-Location Plot",
    x = "Fitted values",
    y = expression(sqrt("|Standardized residuals|"))
  ) +
  theme_minimal()

# Residuals vs leverage
ggplot(diag_df, aes(x = leverage, y = std_residuals)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Residuals vs Leverage",
    x = "Leverage",
    y = "Standardized Residuals"
  ) +
  theme_minimal()

# ==============================
# Simple Logistic Regression
# ==============================

# Prepare data
df_no_na <- df %>%
  select(koi_disposition, koi_prad) %>%
  na.omit()

# Sort values for a cleaner probability curve
df_no_na <- df_no_na[order(df_no_na$koi_prad), ]

# Create binary response variable
df_no_na$confirmed <- ifelse(df_no_na$koi_disposition == "CONFIRMED", 1, 0)

# Fit logistic regression model
model_log_simple <- glm(confirmed ~ koi_prad, data = df_no_na, family = binomial)

# Model summary and odds ratios
summary(model_log_simple)
exp(coef(model_log_simple))

# Predicted probabilities
df_no_na$pred_prob <- predict(model_log_simple, type = "response")

# ==============================
# Probability plot
# ==============================

ggplot(df_no_na, aes(x = koi_prad, y = confirmed)) +
  geom_jitter(height = 0.02, alpha = 0.2) +
  geom_line(aes(y = pred_prob), color = "red", linewidth = 1) +
  scale_x_log10() +
  labs(
    x = "Planet Radius (log scale)",
    y = "Probability of being Confirmed",
    title = "Simple Logistic Regression: koi_prad -> confirmed"
  ) +
  theme_minimal()

# ==============================
# Diagnostic data
# ==============================

diag_log <- data.frame(
  fitted = fitted(model_log_simple),
  residuals = residuals(model_log_simple, type = "deviance"),
  std_residuals = rstandard(model_log_simple),
  leverage = hatvalues(model_log_simple)
)

# ==============================
# Diagnostic plots
# ==============================

# Residuals vs fitted values
ggplot(diag_log, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(
    title = "Residuals vs Fitted Values",
    x = "Fitted probabilities",
    y = "Deviance residuals"
  ) +
  theme_minimal()

# Observed vs predicted probabilities
ggplot(df_no_na, aes(x = pred_prob, y = confirmed)) +
  geom_jitter(height = 0.02, alpha = 0.2) +
  labs(
    title = "Observed vs Predicted Probabilities",
    x = "Predicted probability",
    y = "Observed outcome"
  ) +
  theme_minimal()

# Residuals vs leverage
ggplot(diag_log, aes(x = leverage, y = std_residuals)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(se = FALSE, method = "loess", color = "blue") +
  labs(
    title = "Residuals vs Leverage",
    x = "Leverage",
    y = "Standardized residuals"
  ) +
  theme_minimal()

# ==============================
# Classification thresholds
# ==============================

df_no_na$pred_010 <- ifelse(df_no_na$pred_prob > 0.10, 1, 0)
df_no_na$pred_020 <- ifelse(df_no_na$pred_prob > 0.20, 1, 0)
df_no_na$pred_030 <- ifelse(df_no_na$pred_prob > 0.30, 1, 0)

# ==============================
# Confusion matrix function
# ==============================

plot_confusion <- function(pred, actual, title) {
  
  # Force a complete 2x2 confusion matrix
  conf_mat <- table(
    factor(pred, levels = c(0, 1)),
    factor(actual, levels = c(0, 1))
  )
  
  # Convert table to dataframe
  conf_df <- as.data.frame(conf_mat)
  
  # Add labels in the correct order
  conf_df$label <- c("TN", "FN", "FP", "TP")
  
  # Plot confusion matrix
  ggplot(conf_df, aes(x = Var2, y = Var1, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = paste(label, "\n", Freq)), size = 5) +
    labs(
      x = "Actual",
      y = "Predicted",
      title = title
    ) +
    theme_minimal()
}

# ==============================
# Confusion matrices
# ==============================

plot_confusion(df_no_na$pred_010, df_no_na$confirmed, "Confusion Matrix (threshold = 0.10)")
plot_confusion(df_no_na$pred_020, df_no_na$confirmed, "Confusion Matrix (threshold = 0.20)")
plot_confusion(df_no_na$pred_030, df_no_na$confirmed, "Confusion Matrix (threshold = 0.30)")

# ==============================
# Multiple Logistic Regression
# ==============================

# Prepare data (select variables and remove NA)
df_no_na_multi <- df %>%
  select(koi_disposition,
         koi_prad, koi_period, koi_insol,
         koi_depth, koi_model_snr, koi_steff) %>%
  na.omit()

# Create binary response variable
df_no_na_multi$confirmed <- ifelse(df_no_na_multi$koi_disposition == "CONFIRMED", 1, 0)

# ==============================
# Fit logistic regression model
# ==============================

model_log_multi <- glm(
  confirmed ~ koi_prad + koi_period + koi_insol +
    koi_depth + koi_model_snr + koi_steff,
  data = df_no_na_multi,
  family = binomial
)

# Model summary and odds ratios
summary(model_log_multi)
exp(coef(model_log_multi))

# ==============================
# Predicted probabilities
# ==============================

df_no_na_multi$pred_prob <- predict(model_log_multi, type = "response")

# ==============================
# Predicted probability curves
# ==============================

# Choose representative levels for orbital period and signal-to-noise ratio
period_levels <- quantile(df_no_na_multi$koi_period,
                          probs = c(0.25, 0.50, 0.75),
                          na.rm = TRUE)

snr_levels <- quantile(df_no_na_multi$koi_model_snr,
                       probs = c(0.25, 0.50, 0.75),
                       na.rm = TRUE)

# Create grid of values
plot_data_multi <- expand.grid(
  koi_prad = seq(
    quantile(df_no_na_multi$koi_prad, 0.01, na.rm = TRUE),
    quantile(df_no_na_multi$koi_prad, 0.99, na.rm = TRUE),
    length.out = 200
  ),
  koi_period = period_levels,
  koi_model_snr = snr_levels
)

# Fix remaining variables at median values
plot_data_multi$koi_insol <- median(df_no_na_multi$koi_insol, na.rm = TRUE)
plot_data_multi$koi_depth <- median(df_no_na_multi$koi_depth, na.rm = TRUE)
plot_data_multi$koi_steff <- median(df_no_na_multi$koi_steff, na.rm = TRUE)

# Predict probabilities
plot_data_multi$pred_prob <- predict(
  model_log_multi,
  newdata = plot_data_multi,
  type = "response"
)

# Create readable labels
plot_data_multi$period_group <- factor(
  plot_data_multi$koi_period,
  levels = period_levels,
  labels = c("Low period", "Median period", "High period")
)

plot_data_multi$snr_group <- factor(
  plot_data_multi$koi_model_snr,
  levels = snr_levels,
  labels = c("Low SNR", "Median SNR", "High SNR")
)

# Plot predicted probabilities
ggplot(plot_data_multi, aes(x = koi_prad, y = pred_prob, color = period_group)) +
  geom_line(linewidth = 1) +
  scale_x_log10() +
  facet_wrap(~ snr_group) +
  labs(
    title = "Predicted Probability of Confirmation",
    subtitle = "Effect of planet radius across orbital period and signal-to-noise levels",
    x = "Planet Radius (log scale)",
    y = "Predicted probability",
    color = "Orbital period"
  ) +
  theme_minimal()

# ==============================
# Diagnostic data
# ==============================

diag_log_multi <- data.frame(
  fitted = fitted(model_log_multi),
  residuals = residuals(model_log_multi, type = "deviance"),
  std_residuals = rstandard(model_log_multi),
  leverage = hatvalues(model_log_multi)
)

# ==============================
# Diagnostic plots
# ==============================

# Residuals vs fitted
ggplot(diag_log_multi, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(
    title = "Residuals vs Fitted Values",
    x = "Fitted probabilities",
    y = "Deviance residuals"
  ) +
  theme_minimal()

# Observed vs predicted
ggplot(df_no_na_multi, aes(x = pred_prob, y = confirmed)) +
  geom_jitter(height = 0.02, alpha = 0.2) +
  labs(
    title = "Observed vs Predicted Probabilities",
    x = "Predicted probability",
    y = "Observed outcome"
  ) +
  theme_minimal()

# Residuals vs leverage
ggplot(diag_log_multi, aes(x = leverage, y = std_residuals)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(se = FALSE, method = "loess", color = "blue") +
  labs(
    title = "Residuals vs Leverage",
    x = "Leverage",
    y = "Standardized residuals"
  ) +
  theme_minimal()

# ==============================
# Classification thresholds
# ==============================

# Create predictions for different thresholds
df_no_na_multi$pred_010 <- ifelse(df_no_na_multi$pred_prob > 0.10, 1, 0)
df_no_na_multi$pred_020 <- ifelse(df_no_na_multi$pred_prob > 0.20, 1, 0)
df_no_na_multi$pred_030 <- ifelse(df_no_na_multi$pred_prob > 0.30, 1, 0)

# ==============================
# Confusion matrices
# ==============================

plot_confusion(df_no_na_multi$pred_010, df_no_na_multi$confirmed, "Confusion Matrix (threshold = 0.10)")
plot_confusion(df_no_na_multi$pred_020, df_no_na_multi$confirmed, "Confusion Matrix (threshold = 0.20)")
plot_confusion(df_no_na_multi$pred_030, df_no_na_multi$confirmed, "Confusion Matrix (threshold = 0.30)")
