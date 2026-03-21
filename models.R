library(tidyverse)
library(ggplot2)
library(dplyr)

getwd()
data <- read.csv("exoplanets.csv")
data

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

df

#SIMPLE LINEAR REGRESSION
df_lin_simple <- df %>%
  select(koi_prad, koi_period) %>%
  na.omit
  #koi_prad outliers removal
Q1_prad <- quantile(df_lin_simple$koi_prad, 0.25, na.rm = TRUE)
Q3_prad <- quantile(df_lin_simple$koi_prad, 0.75, na.rm = TRUE)
IQR_value_prad <- Q3_prad - Q1_prad
lower_prad <- Q1_prad - 1.5 * IQR_value_prad
upper_prad <- Q3_prad + 1.5 * IQR_value_prad
  #koi_period outliers removal  
Q1_period <- quantile(df_lin_simple$koi_period, 0.25, na.rm = TRUE)
Q3_period <- quantile(df_lin_simple$koi_period, 0.75, na.rm = TRUE)
IQR_value_period <- Q3_period - Q1_period
lower_period <- Q1_period - 1.5 * IQR_value_period
upper_period <- Q3_period + 1.5 * IQR_value_period

df_no_outliers <- df_lin_simple %>%
  filter(koi_prad >= lower_prad & koi_prad <= upper_prad) %>%
  filter(koi_period >= lower_period & koi_period <= upper_period)

plot(df_no_outliers$koi_prad, df_no_outliers$koi_period)

model_lin_simple <- lm(koi_prad ~ koi_period, data = df_no_outliers)
summary(model_lin_simple)
abline(model_lin_simple, col = "red", lwd=2)

#MULTIPLE LINEAR REGRESSION
df_lin_multi <- df %>%
  select(koi_teq, koi_insol, koi_steff, koi_period) %>%
  na.omit()


remove_outliers <- function(data, column) {
  Q1 <- quantile(data[[column]], 0.25, na.rm = TRUE)
  Q3 <- quantile(data[[column]], 0.75, na.rm = TRUE)
  IQR_value <- Q3 - Q1
  
  lower <- Q1 - 1.5 * IQR_value
  upper <- Q3 + 1.5 * IQR_value
  
  data %>% filter(data[[column]] >= lower & data[[column]] <= upper)
}

df_lin_multi <- remove_outliers(df_lin_multi, "koi_teq")
df_lin_multi <- remove_outliers(df_lin_multi, "koi_insol")
df_lin_multi <- remove_outliers(df_lin_multi, "koi_steff")
df_lin_multi <- remove_outliers(df_lin_multi, "koi_period")

model_lin_multi <- lm(koi_teq ~ koi_insol + koi_steff + koi_period,
                      data = df_lin_multi)

summary(model_lin_multi)

# Predicciones
df_lin_multi$pred <- predict(model_lin_multi)

# Scatter plot
plot(df_lin_multi$pred, df_lin_multi$koi_teq,
     xlab = "Predicted Equilibrium Temperature",
     ylab = "Observed Equilibrium Temperature",
     main = "Observed vs Predicted Equilibrium Temperature")
abline(0,1, col="red", lwd=2)  # línea de identidad

#SIMPLE LOGISTIC REGRESSION
df_no_na <- df[!is.na(df$koi_prad), ]
df_no_na$confirmed <- ifelse(df_no_na$koi_disposition == "CONFIRMED", 1, 0)
model_log_simple <- glm(confirmed ~ koi_prad, data = df_no_na, family = binomial)
summary(model_log_simple)
exp(coef(model_log_simple))

# Predicciones del modelo
df_no_na$pred_prob <- predict(model_log_simple, type = "response")

# Scatter plot con línea de probabilidad
ggplot(df_no_na, aes(x = koi_prad, y = confirmed)) +
  geom_jitter(height = 0.02, alpha = 0.2) +
  geom_line(aes(y = pred_prob), color = "red", size = 1) +
  scale_x_log10() +
  labs(x = "Planet Radius (log10)",
       y = "Probability of being Confirmed",
       title = "Simple Logistic Regression: koi_prad → confirmed") +
  theme_minimal()

#MULTIPLE LOGISTIC REGRESSION

# Filtrar filas sin NA en los predictores y respuesta
df_no_na_multi <- df_no_na[!is.na(df_no_na$koi_prad) &
                             !is.na(df_no_na$koi_period) &
                             !is.na(df_no_na$koi_insol), ]

# Ajustar modelo multivariate
model_log_multi <- glm(confirmed ~ koi_prad + koi_period + koi_insol,
                       data = df_no_na_multi,
                       family = binomial)

# Ver summary
summary(model_log_multi)

# Odds ratios
exp(coef(model_log_multi))

# Crear dataframe para la visualización
plot_data <- data.frame(
  koi_prad = seq(quantile(df_no_na_multi$koi_prad, 0.01),
                 quantile(df_no_na_multi$koi_prad, 0.99),
                 length.out = 100),
  koi_period = median(df_no_na_multi$koi_period, na.rm = TRUE),
  koi_insol = median(df_no_na_multi$koi_insol, na.rm = TRUE)
)

# Predicciones del modelo multivariate
plot_data$pred_prob <- predict(model_log_multi, newdata = plot_data, type = "response")

# Gráfico
ggplot(plot_data, aes(x = koi_prad, y = pred_prob)) +
  geom_line(color = "red", size = 1) +
  scale_x_log10() +
  labs(
    x = "Planet Radius (log10 koi_prad)",
    y = "Predicted Probability of being Confirmed",
    title = "Multivariate Logistic Regression: Predicted Probability vs Planet Radius"
  ) +
  theme_minimal()
