source("src/util.R")
source("src/Analise de reisduos.R")
source("src/modelo_est.R")

library(zoo)
library(quantmod)
library(imputeTS)
library(ggplot2)
library(dplyr)

# Carregando dados e graficos ---------------------------------------------

load("dados/NYSE.RData")
BegSample <- '2004-01-01'
EndSample <- '2010-12-31'
crises <- as.Date(c("2007-07-01", "2008-11-28")) # Nao usar segunda data

NYSE <- NYSE %>% fortify.zoo %>% as_tibble 

ggplot(NYSE, aes(x = Index, y = 100 * nyse)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Preço", title = "NYSE") +
  theme_minimal()

NYSE <- NYSE %>% 
  filter(Index >= as.Date(BegSample), Index <= as.Date(EndSample)) %>% 
  mutate(id = row_number(), nyse = 100 * nyse) %>%
  select(Index, id, everything())

yt <- NYSE$nyse %>% as.vector()
Varyt <- var(yt[1:50])

# Graficos ----------------------------------------------------------------

p1 <- ggplot(NYSE, aes(x = Index, y = NYSE)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Preço", title = "NYSE") +
  theme_minimal()

p2 <- ggplot(NYSE, aes(x = Index, y = nyse)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "NYSE") +
  theme_minimal() + 
  geom_vline(xintercept = as.Date(c("2007-07-01", "2008-11-28")), 
             colour = 'red', size = 1.5, linetype = "dashed")

p3 <- ggplot(NYSE, aes(x = Index, y = nyse)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "NYSE com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(crises[1],
                            as.Date(c("2008-09-15","2008-09-19", "2008-10-08"))), 
             colour = 'red', size = 1.5, linetype = "dashed")
p3
gridExtra::grid.arrange(p1, p2, p3, ncol = 1)

ggplot(NYSE, aes(x = Index, y = 100 * nyse)) +
  geom_line(size = 1L, colour = "#112446") + 
  labs(x = "Tempo", y = "Retorno", title = "NYSE") +
  theme_minimal() 
ggsave(r"{graficos\USA\usa_serie.png}", width = 6, height = 3.5)

acf(yt, plot = F) %>% autoplot() + ylim(c(-1,1)) + ggtitle("") +
  theme_minimal() 
ggsave(r"{graficos\USA\usa_fac_serie.png}", width = 6, height = 3.5)
pacf(yt, plot = F) %>% autoplot() + ylim(c(-1,1)) + ggtitle("") +
  theme_minimal() 
ggsave(r"{graficos\USA\usa_facp_serie.png}", width = 6, height = 3.5)

acf(yt^2, plot = F) %>% autoplot() + ylim(c(-1,1)) + ggtitle("") +
  theme_minimal() 
ggsave(r"{graficos\USA\usa_fac_quad.png}", width = 6, height = 3.5)
pacf(yt^2, plot = F) %>% autoplot() + ylim(c(-1,1)) + ggtitle("") +
  theme_minimal() 
ggsave(r"{graficos\USA\usa_facp_quad.png}", width = 6, height = 3.5)

# Modelo 00 AR(1)-GARCH(1,1) ----------------------------------------------

pars <- list(
  psi1 = log(1),
  psi2 = log(.15),
  psi3 = log(.84),
  ar = .2,
  deltaMedia = 0.0
)

alpha_order <- length(pars$psi2)
beta_order <- length(pars$psi3)
kmed <- length(pars$deltaMedia)
n <- length(yt) # Tamanho da serie

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
# Ordens e Parametros - FIM

# Estimando e residuos - INICIO
(opt0 <- estimando(llike_garch, pars))

media_cond_mod0 <- esp_cond_garch(
  data = yt,
  est = opt0,
  dummy1 = dummy1,
  alpha_order = alpha_order,
  beta_order = beta_order,
  kmed = kmed,
  n = n
)

var_cond_mod0 <- var_cond_garch(
  data = yt,
  est = opt0,
  dummy1 = dummy1,
  Varyt = Varyt,
  alpha_order = alpha_order,
  beta_order = beta_order,
  kmed = kmed,
  n = n
)

resid_pad_mod0 <- (yt - media_cond_mod0)/sqrt(var_cond_mod0)
resid_pad_mod0 <- resid_pad_mod0[-(1:50)]

resid_pad_data <- data.frame(resid_pad = resid_pad_mod0, 
                             time = seq_along(resid_pad_mod0))
resid_pad_data <- resid_pad_data[-1, ]

plot(resid_pad_mod0, type = 'l')

mean(resid_pad_data$resid_pad)
var(resid_pad_data$resid_pad)

# Estimando e residuos - FIM

# FAC e FACP - INICIO
acf(resid_pad_data$resid_pad, plot = F) %>% autoplot() + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot() + ylim(c(-1,1))

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot() + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot() + ylim(c(-1,1))
# FAC e FACP - FIM

# Poder preditivo - INICIO
poder_pred(yt, media_cond_mod0, var_cond_mod0)$rmse
cor(var_cond_mod0[-(1:50)], ((yt - media_cond_mod0)^2)[-(1:50)])^2
# Poder preditivo - FIM

# QQplot e Histograma - INICIO
ggplot(resid_pad_data, aes(sample = resid_pad)) + 
  stat_qq() + 
  geom_abline(slope = 1, intercept = 0) + 
  ylim(-6,6) + 
  scale_x_continuous(limits = c(-6, 6),  breaks = c(-6, -4, -2, 0, 2, 4, 6))

ggplot(resid_pad_data, aes(x = resid_pad)) + 
  geom_histogram(aes(y =..density..), fill = "#0c4c8a") +
  theme_minimal() +
  labs(x = "Residuos padronizados", y = 'Densidade') + 
  scale_x_continuous(limits = c(-6, 6),  breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  stat_function(fun = dnorm, args = list(0, 1), color = 'red')
# QQplot e Histograma - FIM

# TH - INICIO

Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box', lag = 30)
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box', lag = 30)

shapiro.test(resid_pad_data$resid_pad)
tseries::jarque.bera.test(resid_pad_data$resid_pad)
nortest::ad.test(resid_pad_data$resid_pad)

(dw <- sum(diff(yt - media_cond_mod0)^2)/sum((yt - media_cond_mod0)^2))
moments::kurtosis(resid_pad_mod0)
moments::skewness(resid_pad_mod0)

# TH - FIM

# Graficos de linha para esp_cond e var_cond - INICIO
data <- data.frame(
  yt = yt,
  one_step_predict = media_cond_mod0,
  var_cond = var_cond_mod0,
  time = 1:n
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

ggplot(data, aes(x = time, y = sqrt(var_cond))) +
  labs(x = "Tempo", y = "Desvio Condicional") + 
  geom_line(size = 1L, colour = "red") + 
  geom_line(aes(x = time, y = abs(yt)), colour = "blue")
ggsave(r"{graficos\USA\desvio_cond_modelo0.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM


# Modelo 01 ---------------------------------------------------------------

# Ordens e Parametros - INICIO
pars <- list(
  psi2 = log(.15),
  psi3 = log(.84),
  ar = .2,
  deltaMedia = 0.0,
  deltaVar = c(-3, -3, -1, -3, -3)
)

alpha_order <- length(pars$psi2)
beta_order <- length(pars$psi3)
kmed <- length(pars$deltaMedia)
kvar <- length(pars$deltaVar)
n <- length(yt) # Tamanho da serie

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 881, 1185, 1189, 1202),
                                 c(880, 1184, 1188, 1201, n)))
# Ordens e Parametros - FIM

# Estimando e residuos - INICIO

(opt1 <- estimando(llike_model_garch, pars))

media_cond_mod1 <- esp_cond_model(
  data = yt,
  est = opt1,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_cond_mod1 <- var_cond_model(
  data = yt,
  est = opt1,
  dummy1 = dummy1,
  dummy2 = dummy2,
  Varyt = Varyt,
  alpha_order = alpha_order,
  beta_order = beta_order,
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_incond_mod1 <- var_indcond(
  data = yt,
  est = opt1,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  kmed = kmed,
  kvar = kvar
)

resid_pad_mod1 <- (yt - media_cond_mod1)/sqrt(var_cond_mod1)
resid_pad_mod1 <- resid_pad_mod1[-(1:50)]

resid_pad_data <- data.frame(resid_pad = resid_pad_mod1, 
                             time = seq_along(resid_pad_mod1))
resid_pad_data <- resid_pad_data[-1, ]

plot(resid_pad_mod1, type = 'l')
plot(var_incond_mod1, type = 'l')

mean(resid_pad_data$resid_pad)
var(resid_pad_data$resid_pad)

# Estimando e residuos - FIM

# FAC e FACP - INICIO
acf(resid_pad_data$resid_pad, plot = F) %>% autoplot() + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot() + ylim(c(-1,1))

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot() + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot() + ylim(c(-1,1))
# FAC e FACP - FIM

poder_pred(yt, media_cond_mod1, var_cond_mod1)$rmse
cor(var_cond_mod1[-(1:50)], ((yt - media_cond_mod1)^2)[-(1:50)])^2

Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box', lag = 30)
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box', lag = 30)

(dw <- sum(diff(yt - media_cond_mod1)^2)/sum((yt - media_cond_mod1)^2))


# QQplot e Histograma - INICIO
ggplot(resid_pad_data, aes(sample = resid_pad)) + 
  stat_qq() + 
  geom_abline(slope = 1, intercept = 0) + 
  ylim(-6,6) + 
  scale_x_continuous(limits = c(-6, 6),  breaks = c(-6, -4, -2, 0, 2, 4, 6))

ggplot(resid_pad_data, aes(x = resid_pad)) + 
  geom_histogram(aes(y =..density..), fill = "#0c4c8a") +
  theme_minimal() +
  labs(x = "Residuos padronizados", y = 'Densidade') + 
  scale_x_continuous(limits = c(-6, 6),  breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  stat_function(fun = dnorm, args = list(0, 1), color = 'red')
# QQplot e Histograma - FIM

# TH - INICIO

shapiro.test(resid_pad_data$resid_pad)
tseries::jarque.bera.test(resid_pad_data$resid_pad)
nortest::ad.test(resid_pad_data$resid_pad)

moments::kurtosis(resid_pad_mod1)
moments::skewness(resid_pad_mod1)

# TH - FIM

# Graficos de linha para esp_cond e var_cond - INICIO
data <- data.frame(
  yt = yt,
  one_step_predict = media_cond_mod1,
  var_incond = var_incond_mod1,
  var_cond = var_cond_mod1,
  time = 1:n
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

ggplot(data, aes(x = time, y = sqrt(var_cond))) +
  labs(y = "Tempo", x = "Desvio Condicional") + 
  geom_line(size = 1L, colour = "red") + 
  geom_line(aes(x = time, y = abs(yt-one_step_predict)), colour = "blue", alpha = .5)
ggsave(r"{graficos\USA\desvio_cond_modelo1.png}", width = 20, height = 10)

ggplot(data, aes(x = time, y = sqrt(var_incond))) +
  labs(x = "Tempo", y = "Desvio Incondicional") + 
  geom_line(size = 1L, colour = "red") + 
  geom_line(aes(x = time, y = abs(yt-opt1$deltaMedia)), colour = "blue", alpha = .5) +
  xlim(1180, 1190)
ggsave(r"{graficos\USA\desvio_incond_modelo1.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM


# Modelo 02 ---------------------------------------------------------------

# Ordens e Parametros - INICIO
pars <- list(
  psi2 = log(.15),
  psi3 = log(.84),
  ar = .2,
  deltaMedia = 0.0,
  deltaVar = c(-3, -3, -3, -3)
)

alpha_order <- length(pars$psi2)
beta_order <- length(pars$psi3)
kmed <- length(pars$deltaMedia)
kvar <- length(pars$deltaVar)
n <- length(yt) # Tamanho da serie

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 881, 1185, 1202),
                                 c(880, 1184, 1201, n)))
# Ordens e Parametros - FIM

# Estimando e residuos - INICIO

(opt2 <- estimando(llike_model_garch, pars))

media_cond_mod2 <- esp_cond_model(
  data = yt,
  est = opt2,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_cond_mod2 <- var_cond_model(
  data = yt,
  est = opt2,
  dummy1 = dummy1,
  dummy2 = dummy2,
  Varyt = Varyt,
  alpha_order = alpha_order,
  beta_order = beta_order,
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_incond_mod2 <- var_indcond(
  data = yt,
  est = opt2,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  kmed = kmed,
  kvar = kvar
)

resid_pad_mod2 <- (yt - media_cond_mod2)/sqrt(var_cond_mod2)
resid_pad_mod2 <- resid_pad_mod2[-(1:50)]

resid_pad_data <- data.frame(resid_pad = resid_pad_mod2, 
                             time = seq_along(resid_pad_mod2))
resid_pad_data <- resid_pad_data[-1, ]

plot(resid_pad_mod2, type = 'l')
plot(var_incond_mod2, type = 'l')

mean(resid_pad_data$resid_pad)
var(resid_pad_data$resid_pad)

# Estimando e residuos - FIM

# FAC e FACP - INICIO
acf(resid_pad_data$resid_pad, plot = F) %>% autoplot() + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot() + ylim(c(-1,1))

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot() + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot() + ylim(c(-1,1))
# FAC e FACP - FIM

# Poder preditivo - INICIO
poder_pred(yt, media_cond_mod2, var_cond_mod2)$rmse
cor(var_cond_mod2[-(1:50)], ((yt - media_cond_mod2)^2)[-(1:50)])^2
# Poder preditivo - FIM

# QQplot e Histograma - INICIO
ggplot(resid_pad_data, aes(sample = resid_pad)) + 
  stat_qq() + 
  geom_abline(slope = 1, intercept = 0) + 
  ylim(-6,6) + 
  scale_x_continuous(limits = c(-6, 6),  breaks = c(-6, -4, -2, 0, 2, 4, 6))

ggplot(resid_pad_data, aes(x = resid_pad)) + 
  geom_histogram(aes(y =..density..), fill = "#0c4c8a") +
  theme_minimal() +
  labs(x = "Residuos padronizados", y = 'Densidade') + 
  scale_x_continuous(limits = c(-6, 6),  breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  stat_function(fun = dnorm, args = list(0, 1), color = 'red')
# QQplot e Histograma - FIM

# TH - INICIO

Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box', lag = 30)
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box', lag = 30)

shapiro.test(resid_pad_data$resid_pad)
tseries::jarque.bera.test(resid_pad_data$resid_pad)
nortest::ad.test(resid_pad_data$resid_pad)

(dw <- sum(diff(yt - media_cond_mod2)^2)/sum((yt - media_cond_mod2)^2))
moments::kurtosis(resid_pad_mod2)
moments::skewness(resid_pad_mod2)

# TH - FIM

# Graficos de linha para esp_cond e var_cond - INICIO
data <- data.frame(
  yt = yt,
  one_step_predict = media_cond_mod2,
  var_incond = var_incond_mod2,
  var_cond = var_cond_mod2,
  time = 1:n
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

ggplot(data, aes(x = time, y = sqrt(var_cond))) +
  labs(y = "Tempo", x = "Desvio Condicional") + 
  geom_line(size = 1L, colour = "red") + 
  geom_line(aes(x = time, y = abs(yt)), colour = "blue", alpha = .5)
ggsave(r"{graficos\USA\desvio_cond_modelo2.png}", width = 20, height = 10)

ggplot(data, aes(x = time, y = sqrt(var_incond))) +
  labs(x = "Tempo", y = "Desvio Incondicional") + 
  geom_line(size = 1L, colour = "red") + 
  geom_line(aes(x = time, y = abs(yt)), colour = "blue", alpha = .5)
ggsave(r"{graficos\USA\desvio_incond_modelo2.png}", width = 20, height = 10)


# Modelo 03 ---------------------------------------------------------------

# Ordens e Parametros - INICIO
pars <- list(
  psi2 = log(.15),
  psi3 = log(.84),
  ar = .2,
  deltaMedia = 0.0,
  deltaVar = c(-3, -3, -3, -3)
)

alpha_order <- length(pars$psi2)
beta_order <- length(pars$psi3)
kmed <- length(pars$deltaMedia)
kvar <- length(pars$deltaVar)
n <- length(yt) # Tamanho da serie
delta_ind <- 3
t_ast <- 1201
t_til <- 1375

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 881, 1185, 1375),
                                 c(880, 1184, 1201, n)))

(opt3 <- estimando(llike_suave, pars))

media_cond_mod3 <- esp_cond_sauve(
  data = yt,
  est = opt3,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = 3,
  t_ast = 1201,
  t_til = 1375,
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_cond_mod3 <- var_cond_sauve(
  data = yt,
  est = opt3,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  Varyt = Varyt,
  delta_ind = 3,
  t_ast = 1201,
  t_til = 1375,
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_incond_mod3 <- var_indcond_sauve(
  data = yt,
  est = opt3,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = 3,
  t_ast = 1201,
  t_til = 1375,
  kmed = kmed,
  kvar = kvar
)

resid_pad_mod3 <- (yt - media_cond_mod3)/sqrt(var_cond_mod3)
resid_pad_mod3 <- resid_pad_mod3[-(1:50)]

resid_pad_data <- data.frame(resid_pad = resid_pad_mod3, 
                             time = seq_along(resid_pad_mod3))
resid_pad_data <- resid_pad_data[-1, ]

plot(resid_pad_mod3, type = 'l')
plot(var_incond_mod3, type = 'l')

mean(resid_pad_data$resid_pad)
var(resid_pad_data$resid_pad)

# Estimando e residuos - FIM

# FAC e FACP - INICIO
acf(resid_pad_data$resid_pad, plot = F) %>% autoplot() + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot() + ylim(c(-1,1))

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot() + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot() + ylim(c(-1,1))
# FAC e FACP - FIM

poder_pred(yt, media_cond_mod3, var_cond_mod3)$rmse
cor(var_cond_mod3[-(1:50)], ((yt - media_cond_mod3)^2)[-(1:50)])^2

Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box', lag = 30)
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box', lag = 30)

(dw <- sum(diff(yt - media_cond_mod3)^2)/sum((yt - media_cond_mod3)^2))

# QQplot e Histograma - INICIO
ggplot(resid_pad_data, aes(sample = resid_pad)) + 
  stat_qq() + 
  geom_abline(slope = 1, intercept = 0) + 
  ylim(-6,6) + 
  scale_x_continuous(limits = c(-6, 6),  breaks = c(-6, -4, -2, 0, 2, 4, 6))

ggplot(resid_pad_data, aes(x = resid_pad)) + 
  geom_histogram(aes(y =..density..), fill = "#0c4c8a") +
  theme_minimal() +
  labs(x = "Residuos padronizados", y = 'Densidade') + 
  scale_x_continuous(limits = c(-6, 6),  breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  stat_function(fun = dnorm, args = list(0, 1), color = 'red')
# QQplot e Histograma - FIM

# TH - INICIO

shapiro.test(resid_pad_data$resid_pad)
tseries::jarque.bera.test(resid_pad_data$resid_pad)
nortest::ad.test(resid_pad_data$resid_pad)

moments::kurtosis(resid_pad_mod3)
moments::skewness(resid_pad_mod3)

# TH - FIM

# Graficos de linha para esp_cond e var_cond - INICIO
data <- data.frame(
  yt = yt,
  one_step_predict = media_cond_mod3,
  var_incond = var_incond_mod3,
  var_cond = var_cond_mod3,
  time = 1:n
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

ggplot(data, aes(x = time, y = sqrt(var_cond))) +
  labs(y = "Tempo", x = "Desvio Condicional") + 
  geom_line(size = 1L, colour = "red") + 
  geom_line(aes(x = time, y = abs(yt)), colour = "blue", alpha = .5)
ggsave(r"{graficos\USA\desvio_cond_modelo3.png}", width = 20, height = 10)

ggplot(data, aes(x = time, y = sqrt(var_incond))) +
  labs(x = "Tempo", y = "Desvio Incondicional") + 
  geom_line(size = 1L, colour = "red") + 
  geom_line(aes(x = time, y = abs(yt)), colour = "blue", alpha = .5)
ggsave(r"{graficos\USA\desvio_incond_modelo3.png}", width = 20, height = 10)


# Modelo 04 ---------------------------------------------------------------

# Ordens e Parametros - INICIO
pars <- list(
  psi2 = log(.15),
  psi3 = log(.84),
  ar = .2,
  deltaMedia = 0.0,
  deltaVar = c(-3, -3, -3, -3)
)

alpha_order <- length(pars$psi2)
beta_order <- length(pars$psi3)
kmed <- length(pars$deltaMedia)
kvar <- length(pars$deltaVar)
n <- length(yt) # Tamanho da serie
delta_ind <- c(2, 3)
t_ast <- c(1184, 1201)
t_til <- c(1188, 1375)

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 881, 1188, 1375),
                                 c(880, 1184, 1201, n)))

(opt4 <- estimando(llike_suave, pars))

media_cond_mod4 <- esp_cond_sauve(
  data = yt,
  est = opt4,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = c(2, 3),
  t_ast = c(1184, 1201),
  t_til = c(1188, 1375),
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_cond_mod4 <- var_cond_sauve(
  data = yt,
  est = opt4,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  Varyt = Varyt,
  delta_ind = c(2, 3),
  t_ast = c(1184, 1201),
  t_til = c(1188, 1375),
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_incond_mod4 <- var_indcond_sauve(
  data = yt,
  est = opt4,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = c(2, 3),
  t_ast = c(1184, 1201),
  t_til = c(1188, 1375),
  kmed = kmed,
  kvar = kvar
)

resid_pad_mod4 <- (yt - media_cond_mod4)/sqrt(var_cond_mod4)
resid_pad_mod4 <- resid_pad_mod4[-(1:50)]

resid_pad_data <- data.frame(resid_pad = resid_pad_mod4, 
                             time = seq_along(resid_pad_mod4))
resid_pad_data <- resid_pad_data[-1, ]

plot(resid_pad_mod4, type = 'l')
plot(var_incond_mod4, type = 'l')

mean(resid_pad_data$resid_pad)
var(resid_pad_data$resid_pad)

# Estimando e residuos - FIM

# FAC e FACP - INICIO
acf(resid_pad_data$resid_pad, plot = F) %>% autoplot() + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot() + ylim(c(-1,1))

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot() + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot() + ylim(c(-1,1))
# FAC e FACP - FIM

poder_pred(yt, media_cond_mod4, var_cond_mod4)$rmse
cor(var_cond_mod4[-(1:50)], ((yt - media_cond_mod4)^2)[-(1:50)])^2

Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box', lag = 30)
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box', lag = 30)

(dw <- sum(diff(yt - media_cond_mod4)^2)/sum((yt - media_cond_mod4)^2))

# QQplot e Histograma - INICIO
ggplot(resid_pad_data, aes(sample = resid_pad)) + 
  stat_qq() + 
  geom_abline(slope = 1, intercept = 0) + 
  ylim(-6,6) + 
  scale_x_continuous(limits = c(-6, 6),  breaks = c(-6, -4, -2, 0, 2, 4, 6))

ggplot(resid_pad_data, aes(x = resid_pad)) + 
  geom_histogram(aes(y =..density..), fill = "#0c4c8a") +
  theme_minimal() +
  labs(x = "Residuos padronizados", y = 'Densidade') + 
  scale_x_continuous(limits = c(-6, 6),  breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  stat_function(fun = dnorm, args = list(0, 1), color = 'red')
# QQplot e Histograma - FIM

# TH - INICIO

shapiro.test(resid_pad_data$resid_pad)
tseries::jarque.bera.test(resid_pad_data$resid_pad)
nortest::ad.test(resid_pad_data$resid_pad)

moments::kurtosis(resid_pad_mod4)
moments::skewness(resid_pad_mod4)

# TH - FIM

# Graficos de linha para esp_cond e var_cond - INICIO
data <- data.frame(
  yt = yt,
  one_step_predict = media_cond_mod4,
  var_incond = var_incond_mod4,
  var_cond = var_cond_mod4,
  time = 1:n
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

ggplot(data, aes(x = time, y = sqrt(var_cond))) +
  labs(y = "Tempo", x = "Desvio Condicional") + 
  geom_line(size = 1L, colour = "red") + 
  geom_line(aes(x = time, y = abs(yt)), colour = "blue", alpha = .5)
ggsave(r"{graficos\USA\desvio_cond_modelo4.png}", width = 20, height = 10)

ggplot(data, aes(x = time, y = sqrt(var_incond))) +
  labs(x = "Tempo", y = "Desvio Incondicional") + 
  geom_line(size = 1L, colour = "red") + 
  geom_line(aes(x = time, y = abs(yt-opt4$deltaMedia)), colour = "blue", alpha = .5)
ggsave(r"{graficos\USA\desvio_incond_modelo4.png}", width = 20, height = 10)


# Modelo 05 ---------------------------------------------------------------

# Ordens e Parametros - INICIO
pars <- list(
  psi2 = log(.15),
  psi3 = log(.84),
  ar = .2,
  deltaMedia = 0.0,
  deltaVar = c(-3, -3, -3, -3, -3)
)

alpha_order <- length(pars$psi2)
beta_order <- length(pars$psi3)
kmed <- length(pars$deltaMedia)
kvar <- length(pars$deltaVar)
n <- length(yt) # Tamanho da serie
delta_ind <- 2
t_ast <- 1184
t_til <- 1201

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 881, 1201, 1241, 1376),
                                 c(880, 1184, 1240, 1375, n)))

(opt5 <- estimando(llike_model_garch, pars))

media_cond_mod5 <- esp_cond_sauve(
  data = yt,
  est = opt5,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = 2,
  t_ast = 1184,
  t_til = 1201,
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_cond_mod5 <- var_cond_sauve(
  data = yt,
  est = opt5,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  Varyt = Varyt,
  delta_ind = 2,
  t_ast = 1184,
  t_til = 1201,
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_incond_mod5 <- var_indcond_sauve(
  data = yt,
  est = opt5,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = 2,
  t_ast = 1184,
  t_til = 1201,
  kmed = kmed,
  kvar = kvar
)

resid_pad_mod5 <- (yt - media_cond_mod5)/sqrt(var_cond_mod5)
resid_pad_mod5 <- resid_pad_mod5[-(1:50)]

resid_pad_data <- data.frame(resid_pad = resid_pad_mod5, 
                             time = seq_along(resid_pad_mod5))
resid_pad_data <- resid_pad_data[-1, ]

plot(resid_pad_mod5, type = 'l')
plot(var_incond_mod5, type = 'l')

mean(resid_pad_data$resid_pad)
var(resid_pad_data$resid_pad)

# Estimando e residuos - FIM

# FAC e FACP - INICIO
acf(resid_pad_data$resid_pad, plot = F) %>% autoplot() + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot() + ylim(c(-1,1))

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot() + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot() + ylim(c(-1,1))
# FAC e FACP - FIM

# Poder preditivo - INICIO
poder_pred(yt, media_cond_mod5, var_cond_mod5)$rmse
cor(var_cond_mod5[-(1:50)], ((yt - media_cond_mod5)^2)[-(1:50)])^2
# Poder preditivo - FIM

# QQplot e Histograma - INICIO
ggplot(resid_pad_data, aes(sample = resid_pad)) + 
  stat_qq() + 
  geom_abline(slope = 1, intercept = 0) + 
  ylim(-6,6) + 
  scale_x_continuous(limits = c(-6, 6),  breaks = c(-6, -4, -2, 0, 2, 4, 6))

ggplot(resid_pad_data, aes(x = resid_pad)) + 
  geom_histogram(aes(y =..density..), fill = "#0c4c8a") +
  theme_minimal() +
  labs(x = "Residuos padronizados", y = 'Densidade') + 
  scale_x_continuous(limits = c(-6, 6),  breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  stat_function(fun = dnorm, args = list(0, 1), color = 'red')
# QQplot e Histograma - FIM

# TH - INICIO

Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box', lag = 30)
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box', lag = 30)

shapiro.test(resid_pad_data$resid_pad)
tseries::jarque.bera.test(resid_pad_data$resid_pad)
nortest::ad.test(resid_pad_data$resid_pad)

(dw <- sum(diff(yt - media_cond_mod5)^2)/sum((yt - media_cond_mod5)^2))
moments::kurtosis(resid_pad_mod5)
moments::skewness(resid_pad_mod5)

# TH - FIM

# Graficos de linha para esp_cond e var_cond - INICIO
data <- data.frame(
  yt = yt,
  one_step_predict = media_cond_mod5,
  var_incond = var_incond_mod5,
  var_cond = var_cond_mod5,
  time = 1:n
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

ggplot(data, aes(x = time, y = sqrt(var_cond))) +
  labs(y = "Tempo", x = "Desvio Condicional") + 
  geom_line(size = 1L, colour = "red") + 
  geom_line(aes(x = time, y = abs(yt)), colour = "blue", alpha = .5)
ggsave(r"{graficos\USA\desvio_cond_modelo5.png}", width = 20, height = 10)

ggplot(data, aes(x = time, y = sqrt(var_incond))) +
  labs(x = "Tempo", y = "Desvio Incondicional") + 
  geom_line(size = 1L, colour = "red") + 
  geom_line(aes(x = time, y = abs(yt)), colour = "blue", alpha = .5)
ggsave(r"{graficos\USA\desvio_incond_modelo5.png}", width = 20, height = 10)

# Modelo 06 ---------------------------------------------------------------

# Ordens e Parametros - INICIO
pars <- list(
  psi2 = log(.15),
  psi3 = log(.84),
  ar = .2,
  deltaMedia = 0.0,
  deltaVar = c(-3, -3, -3, -3, -3)
)

alpha_order <- length(pars$psi2)
beta_order <- length(pars$psi3)
kmed <- length(pars$deltaMedia)
kvar <- length(pars$deltaVar)
n <- length(yt) # Tamanho da serie
delta_ind <- c(2, 3)
t_ast <- c(1184, 1201)
t_til <- c(1201, 1375)

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 881, 1201, 1375, 1496),
                                 c(880, 1184, 1201, 1495, n)))

(opt6 <- estimando(llike_suave, pars))

media_cond_mod6 <- esp_cond_sauve(
  data = yt,
  est = opt6,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = c(2, 3),
  t_ast = c(1184, 1201),
  t_til = c(1201, 1375),
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_cond_mod6 <- var_cond_sauve(
  data = yt,
  est = opt6,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  Varyt = Varyt,
  delta_ind = c(2, 3),
  t_ast = c(1184, 1201),
  t_til = c(1201, 1375),
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_incond_mod6 <- var_indcond_sauve(
  data = yt,
  est = opt6,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = c(2, 3),
  t_ast = c(1184, 1201),
  t_til = c(1201, 1375),
  kmed = kmed,
  kvar = kvar
)

resid_pad_mod6 <- (yt - media_cond_mod6)/sqrt(var_cond_mod6)
resid_pad_mod6 <- resid_pad_mod6[-(1:50)]

resid_pad_data <- data.frame(resid_pad = resid_pad_mod6, 
                             time = seq_along(resid_pad_mod6))
resid_pad_data <- resid_pad_data[-1, ]

plot(resid_pad_mod6, type = 'l')
plot(var_incond_mod6, type = 'l')

mean(resid_pad_data$resid_pad)
var(resid_pad_data$resid_pad)

# Estimando e residuos - FIM

# FAC e FACP - INICIO
acf(resid_pad_data$resid_pad, plot = F) %>% autoplot() + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot() + ylim(c(-1,1))

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot() + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot() + ylim(c(-1,1))
# FAC e FACP - FIM

poder_pred(yt, media_cond_mod6, var_cond_mod6)$rmse
cor(var_cond_mod6[-(1:50)], ((yt - media_cond_mod6)^2)[-(1:50)])^2

Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box', lag = 30)
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box', lag = 30)

(dw <- sum(diff(yt - media_cond_mod6)^2)/sum((yt - media_cond_mod6)^2))

# QQplot e Histograma - INICIO
ggplot(resid_pad_data, aes(sample = resid_pad)) + 
  stat_qq() + 
  geom_abline(slope = 1, intercept = 0) + 
  ylim(-6,6) + 
  scale_x_continuous(limits = c(-6, 6),  breaks = c(-6, -4, -2, 0, 2, 4, 6))

ggplot(resid_pad_data, aes(x = resid_pad)) + 
  geom_histogram(aes(y =..density..), fill = "#0c4c8a") +
  theme_minimal() +
  labs(x = "Residuos padronizados", y = 'Densidade') + 
  scale_x_continuous(limits = c(-6, 6),  breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  stat_function(fun = dnorm, args = list(0, 1), color = 'red')
# QQplot e Histograma - FIM

# TH - INICIO

shapiro.test(resid_pad_data$resid_pad)
tseries::jarque.bera.test(resid_pad_data$resid_pad)
nortest::ad.test(resid_pad_data$resid_pad)

moments::kurtosis(resid_pad_mod6)
moments::skewness(resid_pad_mod6)

# TH - FIM

# Graficos de linha para esp_cond e var_cond - INICIO
data <- data.frame(
  yt = yt,
  one_step_predict = media_cond_mod6,
  var_incond = var_incond_mod6,
  var_cond = var_cond_mod6,
  time = 1:n
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

ggplot(data, aes(x = time, y = sqrt(var_cond))) +
  labs(y = "Tempo", x = "Desvio Condicional") + 
  geom_line(size = 1L, colour = "red") + 
  geom_line(aes(x = time, y = abs(yt)), colour = "blue", alpha = .5)
ggsave(r"{graficos\USA\desvio_cond_modelo6.png}", width = 20, height = 10)

ggplot(data, aes(x = time, y = sqrt(var_incond))) +
  labs(x = "Tempo", y = "Desvio Incondicional") + 
  geom_line(size = 1L, colour = "red") + 
  geom_line(aes(x = time, y = abs(yt)), colour = "blue", alpha = .5) +
  geom_vline(xintercept = 1375)
ggsave(r"{graficos\USA\desvio_incond_modelo6.png}", width = 20, height = 10)

# Resultados --------------------------------------------------------------

medidas <- function(modelo, nome){
  modelo %>% select(llike, AIC, BIC) %>% mutate(Modelo = nome) %>% 
    select(Modelo, llike, AIC, BIC)
  
}

resultado <- rbind(
  medidas(opt1, "opt1"),
  medidas(opt2, "opt2"),
  medidas(opt3, "opt3"),
  medidas(opt4, "opt4"),
  medidas(opt5, "opt5"),
  medidas(opt6, "opt6")
)

resultado
resultado %>% arrange(AIC)
resultado %>% arrange(BIC)

## Teste LR

teste_lr(opt1, opt0)
teste_lr(opt2, opt0)

teste_lr(opt1, opt2)
teste_lr(opt6, opt4)

## Poder preditivo
poder_pred(yt, media_cond_mod0, var_cond_mod0)$rmse
poder_pred(yt, media_cond_mod1, var_cond_mod1)$rmse
poder_pred(yt, media_cond_mod2, var_cond_mod2)$rmse
poder_pred(yt, media_cond_mod3, var_cond_mod3)$rmse
poder_pred(yt, media_cond_mod4, var_cond_mod4)$rmse
poder_pred(yt, media_cond_mod5, var_cond_mod5)$rmse
poder_pred(yt, media_cond_mod6, var_cond_mod6)$rmse

## Cor
cor(var_cond_mod0[-(1:50)], ((yt - media_cond_mod0)^2)[-(1:50)])^2
cor(var_cond_mod1[-(1:50)], ((yt - media_cond_mod1)^2)[-(1:50)])^2
cor(var_cond_mod2[-(1:50)], ((yt - media_cond_mod2)^2)[-(1:50)])^2
cor(var_cond_mod3[-(1:50)], ((yt - media_cond_mod3)^2)[-(1:50)])^2
cor(var_cond_mod4[-(1:50)], ((yt - media_cond_mod4)^2)[-(1:50)])^2
cor(var_cond_mod5[-(1:50)], ((yt - media_cond_mod5)^2)[-(1:50)])^2
cor(var_cond_mod6[-(1:50)], ((yt - media_cond_mod6)^2)[-(1:50)])^2

