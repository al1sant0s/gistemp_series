setwd("C:/Users/Kleber Santos/Dropbox/PesquisaAlissonKleber/scriptR")
source("gbparma.R")

# Pacotes -----------------------------------------------------------------
library(forecast)
library(ggplot2)
library(patchwork)
library(PTSR)
library(extraDistr)
library(tidyverse)
library(tsibble)

install.packages("pacman")
library("pacman")
p_load("tidyverse",
  "rio",
  "forecast",
  "tsibble",
  "lubridate",
  "fabletools",
  "ggtime",
  "gt",
  "rlang",
  "seastests",
  "tseries",
  "lmtest",
  "Kendall",
  "randtests",
  "nortest",
  "FinTS",
  "trend",
  "extraDistr",
  "katex"
)

# Funcoes para medidas de acuracia ----------------------------------------
mse <- function(y_true, y_pred) {
  mean((y_true - y_pred)^2)
}

mae <- function(y_true, y_pred) {
  mean(abs(y_true - y_pred))
}

smape <- function(y_true, y_pred) {
  numerator <- 2 * abs(y_true - y_pred)
  denominator <- abs(y_true) + abs(y_pred)
  mean(numerator / denominator) * 100
}

mape <- function(y_true, y_pred) {
  n <- length(y_true)
  mape <- (1/n) * sum(abs((y_true - y_pred) / y_true) * 100)
  return(mape)
}

# Importando os dados -----------------------------------------------------
avg_gtemp <- import("C:/Users/Kleber Santos/Dropbox/PesquisaAlissonKleber/data/temperatures.csv", HEADER = TRUE) |>
  select(1:13) |>
  mutate(across(.cols = -Year, .fns = as.numeric)) |>
  pivot_longer(cols = -Year, names_to = "Month", values_to = "Mean", values_drop_na = TRUE) |>
  mutate(Date = yearmonth(paste(Year, Month)), .before = Mean) |>
  filter(Date >= yearmonth("2010 Jan"), Date < yearmonth("2026 Apr")) |>
  as_tsibble(index = Date)

x <- avg_gtemp$Mean

y <- ts(data = x, start = c(2010,1), end = c(2026,3), frequency = 12)

# Descritivas -------------------------------------------------------------
summary(y)
moments::kurtosis(y)
moments::skewness(y)
sd(y)

# Plots da serie ----------------------------------------------------------
serie <- autoplot(y, ylab = "Flow")
sub_serie <- ggsubseriesplot(y, ylab = "Flow")

bacf <- acf(y, plot = F, lag.max = 24)
bacfdf <- with(bacf, data.frame(lag, acf))[-1,]

acf_serie <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  geom_hline(aes(yintercept = 2/sqrt(length(y))), linetype = 2, color = 'black') +
  geom_hline(aes(yintercept = -2/sqrt(length(y))), linetype = 2, color = 'black') +
  ylab("ACF") + xlab("Lag")

bpacf <- pacf(y, plot = F, lag.max = 24)
bpacfdf <- with(bpacf, data.frame(lag, acf))

pacf_serie <- ggplot(data = bpacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  geom_hline(aes(yintercept = 2/sqrt(length(y))), linetype = 2, color = 'black') +
  geom_hline(aes(yintercept = -2/sqrt(length(y))), linetype = 2, color = 'black') +
  ylab("PACF") + xlab("Lag")

(serie + sub_serie) / (acf_serie + pacf_serie)

# Serie da amostra --------------------------------------------------------
n <- length(y)
h1 <- 6
n1 <- n-h1

y1 <- ts(y[1:n1], start = c(2010,1), frequency = 12)


# Variaveis harmonicas ----------------------------------------------------
t <- 1:n1 # in sample
t_hat <- (n1+1):(n1+h1) # out of sample

# deterministic seazonality
C<-cos(2*pi*t/12) # in sample 
C_hat<-cos(2*pi*t_hat/12) # out of sample

S<-sin(2*pi*t/12) # in sample
S_hat<-sin(2*pi*t_hat/12) # out of sample

nina34 <- import("C:/Users/Kleber Santos/Dropbox/PesquisaAlissonKleber/data/nina34.anom.csv", HEADER = TRUE) |>
  mutate(Date = yearmonth(Date)) |>
  filter(Date >= yearmonth("2010-01"), Date < yearmonth("2026-04")) |>
  mutate(Mean = parse_number(Mean, locale = locale(decimal_mark = ","))) |>
  na.omit()

nina34_treino <- nina34$Mean[1:n1]
nina34_teste <- tail(nina34$Mean,h1)

# More than on covariate
mX<-cbind(S,C,t,nina34_treino) # in sample
mX_hat<-cbind(S_hat,C_hat,t_hat,nina34_teste) # out of sample

# Modelo proposto ---------------------------------------------------------
fit_gbparma <- gbparma(y = y1, ar = c(1), ma = c(1), h1 = h1, X = mX, 
                       X_hat = mX_hat, resid = 1, diag = 0)

# Portmanteau -------------------------------------------------------------

Box.test(fit_gbparma$resid1, lag = round(sqrt(n1)), type = "Ljung-Box", fitdf = 1)

# Analise de diagnostico --------------------------------------------------

bacf <- acf(fit_gbparma$resid3[-1], plot = F, lag.max = 24)
bacfdf <- with(bacf, data.frame(lag, acf))[-1,]

acf_res <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  geom_hline(aes(yintercept = 1.96/sqrt(n-1)), linetype = 2, color = 'black') +
  geom_hline(aes(yintercept = -1.96/sqrt(n-1)), linetype = 2, color = 'black') +
  ylab("ACF") + xlab("Lag")

bpacf <- pacf(fit_gbparma$resid3[-1], plot = F, lag.max = 24)
bpacfdf <- with(bpacf, data.frame(lag, acf))

pacf_res <- ggplot(data = bpacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  geom_hline(aes(yintercept = 1.96/sqrt(n-1)), linetype = 2, color = 'black') +
  geom_hline(aes(yintercept = -1.96/sqrt(n-1)), linetype = 2, color = 'black') +
  ylab("PACF") + xlab("Lag")

(acf_res + pacf_res)

# Plot do ajuste ----------------------------------------------------------
b <- seq.Date(as.Date("2010/1/1"), as.Date("2025/9/1"), "month")
x1 <- data.frame(Data = b, Observado = y1,
                 GBPARMA = fit_gbparma$fitted)

ajuste <- ggplot(x1, aes(x=Data)) +
  geom_line(aes(y = Observado, colour = "Observed data")) +
  geom_line(aes(y = GBPARMA, colour = "Fitted values")) +
  ggplot2::scale_colour_manual(
    values = c('Observed data' = "gray", 'Fitted values' = "black"),
  ) +
  ggplot2::labs(colour = NULL) +
  xlab("Time") +
  ylab("Flow"); ajuste

cor(fit_gbparma$fitted[-c(1:3)],y1[-c(1:3)])

# plots da precisao --------------------------------------------------------
summary(fit_gbparma$fitted_prec[-c(1:3)])

b <- data.frame(y = fit_gbparma$fitted_prec[-c(1:3)],
                x = 4:(length(fit_gbparma$fitted_prec)))

ggplot(b, aes(x = x, y = y)) +
  geom_line() +
  # geom_hline(yintercept = fit_bparma$prec, linetype = "dashed") +
  ylab(expression(hat(phi)[~t])) +
  xlab("t") 

# Teste da RV -------------------------------------------------------------
xx <- 2*(fit_gbparma$loglik - fit_bparma$loglik)
(1-pchisq(xx,1))
qchisq(0.95,1)

# Critério de selecao de modelos ------------------------------------------
fit_gbparma$aic;fit_gbparma$bic


# Demais modelos ----------------------------------------------------------
SARIMA <- auto.arima(y1, ic = "aic",allowdrift = F);SARIMA
prev_sarima <- predict(SARIMA,6)

# medidas de acuracia das previsões ---------------------------------------

## MAE
a1 <- mae(y_true = y[n1+1], y_pred = fit_gbparma$forecast[1])
a2 <- mae(y_true = y[(n1+1):(n1+2)], y_pred = fit_gbparma$forecast[1:2])
a3 <- mae(y_true = y[(n1+1):(n1+3)], y_pred = fit_gbparma$forecast[1:3])
a4 <- mae(y_true = y[(n1+1):(n1+4)], y_pred = fit_gbparma$forecast[1:4])
a5 <- mae(y_true = y[(n1+1):(n1+5)], y_pred = fit_gbparma$forecast[1:5])
a6 <- mae(y_true = y[(n1+1):(n1+6)], y_pred = fit_gbparma$forecast[1:6])

c1 <- mae(y_true = y[n1+1], y_pred = prev_sarima$pred[1])
c2 <- mae(y_true = y[(n1+1):(n1+2)], y_pred = prev_sarima$pred[1:2])
c3 <- mae(y_true = y[(n1+1):(n1+3)], y_pred = prev_sarima$pred[1:3])
c4 <- mae(y_true = y[(n1+1):(n1+4)], y_pred = prev_sarima$pred[1:4])
c5 <- mae(y_true = y[(n1+1):(n1+5)], y_pred = prev_sarima$pred[1:5])
c6 <- mae(y_true = y[(n1+1):(n1+6)], y_pred = prev_sarima$pred[1:6])

mae1 <- data.frame(
  'h=1' = c(a1,c1),
  'h=2' = c(a2,c2),
  'h=3' = c(a3,c3),
  'h=4' = c(a4,c4),
  'h=5' = c(a5,c5),
  'h=6' = c(a6,c6)
)

rownames(mae1) <- c("Generalized BPARMA", "SARIMA");mae1
round(mae1,4)

# xtable::xtable(mae1, digits = 4)

## MAPE
a1 <- mape(y_true = y[n1+1], y_pred = fit_gbparma$forecast[1])
a2 <- mape(y_true = y[(n1+1):(n1+2)], y_pred = fit_gbparma$forecast[1:2])
a3 <- mape(y_true = y[(n1+1):(n1+3)], y_pred = fit_gbparma$forecast[1:3])
a4 <- mape(y_true = y[(n1+1):(n1+4)], y_pred = fit_gbparma$forecast[1:4])
a5 <- mape(y_true = y[(n1+1):(n1+5)], y_pred = fit_gbparma$forecast[1:5])
a6 <- mape(y_true = y[(n1+1):(n1+6)], y_pred = fit_gbparma$forecast[1:6])

c1 <- mape(y_true = y[n1+1], y_pred = prev_sarima$pred[1])
c2 <- mape(y_true = y[(n1+1):(n1+2)], y_pred = prev_sarima$pred[1:2])
c3 <- mape(y_true = y[(n1+1):(n1+3)], y_pred = prev_sarima$pred[1:3])
c4 <- mape(y_true = y[(n1+1):(n1+4)], y_pred = prev_sarima$pred[1:4])
c5 <- mape(y_true = y[(n1+1):(n1+5)], y_pred = prev_sarima$pred[1:5])
c6 <- mape(y_true = y[(n1+1):(n1+6)], y_pred = prev_sarima$pred[1:6])

mape1 <- data.frame(
  'h=1' = c(a1,c1),
  'h=2' = c(a2,c2),
  'h=3' = c(a3,c3),
  'h=4' = c(a4,c4),
  'h=5' = c(a5,c5),
  'h=6' = c(a6,c6)
)

rownames(mape1) <- c("Generalized BPARMA", "SARIMA");mape1
round(mape1,4)

# xtable::xtable(mape1, digits = 4)

## sMAPE
a1 <- smape(y_true = y[n1+1], y_pred = fit_gbparma$forecast[1])
a2 <- smape(y_true = y[(n1+1):(n1+2)], y_pred = fit_gbparma$forecast[1:2])
a3 <- smape(y_true = y[(n1+1):(n1+3)], y_pred = fit_gbparma$forecast[1:3])
a4 <- smape(y_true = y[(n1+1):(n1+4)], y_pred = fit_gbparma$forecast[1:4])
a5 <- smape(y_true = y[(n1+1):(n1+5)], y_pred = fit_gbparma$forecast[1:5])
a6 <- smape(y_true = y[(n1+1):(n1+6)], y_pred = fit_gbparma$forecast[1:6])

c1 <- smape(y_true = y[n1+1], y_pred = prev_sarima$pred[1])
c2 <- smape(y_true = y[(n1+1):(n1+2)], y_pred = prev_sarima$pred[1:2])
c3 <- smape(y_true = y[(n1+1):(n1+3)], y_pred = prev_sarima$pred[1:3])
c4 <- smape(y_true = y[(n1+1):(n1+4)], y_pred = prev_sarima$pred[1:4])
c5 <- smape(y_true = y[(n1+1):(n1+5)], y_pred = prev_sarima$pred[1:5])
c6 <- smape(y_true = y[(n1+1):(n1+6)], y_pred = prev_sarima$pred[1:6])

smape1 <- data.frame(
  'h=1' = c(a1,c1),
  'h=2' = c(a2,c2),
  'h=3' = c(a3,c3),
  'h=4' = c(a4,c4),
  'h=5' = c(a5,c5),
  'h=6' = c(a6,c6)
)

rownames(smape1) <- c("Generalized BPARMA", "SARIMA");smape1
round(smape1,4)

# xtable::xtable(smape1, digits = 4)
