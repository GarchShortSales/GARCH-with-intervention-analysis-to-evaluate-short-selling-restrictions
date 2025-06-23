source("src/util.R")
source("src/Analise de reisduos.R")
source("src/modelo_est.R")
source("src/modelo7_dax.R")
source("src/graficos.R")

library(zoo)
library(quantmod)
library(ggplot2)
library(dplyr)
library(imputeTS)

# Carregando dados --------------------------------------------------------

load("dados/DAX.RData")
BegSample <- '2004-01-01'
EndSample <- '2009-12-31'
crises <- as.Date(c("2007-07-01", "2008-11-28")) # Nao usar segunda data

DAX <- DAX %>% fortify.zoo %>% as_tibble 

ggplot(DAX, aes(x = Index, y = 100 * dax)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "DAX") + 
  theme_minimal()

DAX <- DAX %>% 
  filter(Index >= as.Date(BegSample), Index <= as.Date(EndSample)) %>% 
  mutate(id = row_number(), dax = 100*dax) %>% 
  select(Index, id, everything())

yt <- DAX$dax %>% as.vector()
Varyt <- var(yt[1:50])

# Graficos ----------------------------------------------------------------

p1 <- ggplot(DAX, aes(x = Index, y = DAX)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Preço", title = "DAX") +
  theme_minimal()

p2 <- ggplot(DAX, aes(x = Index, y = dax)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "DAX") +
  theme_minimal() 

p3 <- ggplot(DAX, aes(x = Index, y = dax)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "DAX com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(crises[1], as.Date("2008-09-15"), as.Date("2008-09-20")), 
             colour = c('red', 'red', 'red'), size = 1.5,
             linetype = "dashed")
# p2
gridExtra::grid.arrange(p2, p3, ncol = 1)

ggplot(DAX, aes(x = Index, y = 100 * dax)) +
  geom_line(size = 1L, colour = "#112446") + 
  labs(x = "Tempo", y = "Retorno", title = "DAX") +
  theme_minimal()  + tema
ggsave(r"{graficos\DAX\ger_serie.png}", width = 20, height = 10)

acf(yt, plot = F) %>% autoplot(main='') + ylim(c(-1,1)) + ggtitle("") + tema
ggsave(r"{graficos\DAX\ger_fac_serie.png}", width = 10, height = 10)

pacf(yt, plot = F) %>% autoplot(main='') + ylim(c(-1,1)) + ggtitle("") + tema
ggsave(r"{graficos\DAX\ger_facp_serie.png}", width = 10, height = 10)

acf(yt^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1)) + ggtitle("") + tema
ggsave(r"{graficos\DAX\ger_fac_quad.png}", width = 10, height = 10)

pacf(yt^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1)) + ggtitle("") + tema
ggsave(r"{graficos\DAX\ger_facp_quad.png}", width = 10, height = 10) 

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
acf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_fac_modelo0_serie.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_facp_modelo0_serie.png}", width = 10, height = 10)

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_fac_modelo0_quad.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_facp_modelo0_quad.png}", width = 10, height = 10)
# FAC e FACP - FIM

poder_pred(yt, media_cond_mod0, var_cond_mod0)$rmse
cor(var_cond_mod0[-(1:50)], ((yt - media_cond_mod0)^2)[-(1:50)])^2

Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box', lag = 30)
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box', lag = 30)

(dw <- sum(diff(yt - media_cond_mod0)^2)/sum((yt - media_cond_mod0)^2))

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

moments::kurtosis(resid_pad_mod0)
moments::skewness(resid_pad_mod0)

# TH - FIM

# Graficos de linha para esp_cond e var_cond - INICIO
data <- data.frame(
  yt = yt,
  one_step_predict = media_cond_mod0,
  var_cond = var_cond_mod0,
  time = DAX$Index
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

grafico_var_cond(data)
#ggsave(r"{graficos\DAX\desvio_cond_modelo0.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM

# Modelo 01 ---------------------------------------------------------------

ggplot(DAX, aes(x = Index, y = dax)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "DAX com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(crises[1], as.Date("2008-09-20")), 
             colour = c('red', 'red'), size = 1.5,
             linetype = "dashed")

# Ordens e Parametros - INICIO
pars <- list(
  psi2 = log(.15),
  psi3 = log(.84),
  ar = .2,
  deltaMedia = 0.0,
  deltaVar = c(-3, -3, -3)
)

alpha_order <- length(pars$psi2)
beta_order <- length(pars$psi3)
kmed <- length(pars$deltaMedia)
kvar <- length(pars$deltaVar)
n <- length(yt) # Tamanho da serie

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 896, 1208),
                                 c(895, 1207, n)))
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
acf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
# FAC e FACP - FIM

# Poder preditivo - INICIO
poder_pred(yt, media_cond_mod1, var_cond_mod1)$rmse
cor(var_cond_mod1[-(1:50)], ((yt - media_cond_mod1)^2)[-(1:50)])^2
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

(dw <- sum(diff(yt - media_cond_mod1)^2)/sum((yt - media_cond_mod1)^2))
moments::kurtosis(resid_pad_mod1)
moments::skewness(resid_pad_mod1)

# TH - FIM

# Graficos de linha para esp_cond e var_cond - INICIO
data <- data.frame(
  yt = yt,
  one_step_predict = media_cond_mod1,
  var_cond = var_cond_mod1,
  var_incond = var_incond_mod1,
  med_incond = opt1$deltaMedia,
  time = DAX$Index
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

grafico_var_cond(data)
#ggsave(r"{graficos\DAX\desvio_cond_modelo1.png}", width = 20, height = 10)

grafico_var_incond(data)
#ggsave(r"{graficos\DAX\desvio_incond_modelo1.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM

# Modelo 02 ---------------------------------------------------------------

ggplot(DAX, aes(x = Index, y = dax)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "DAX com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(crises[1], as.Date("2008-09-15")), 
             colour = c('red', 'red'), size = 1.5,
             linetype = "dashed")

# Ordens e Parametros - INICIO
pars <- list(
  psi2 = log(.15),
  psi3 = log(.84),
  ar = .2,
  deltaMedia = 0.0,
  deltaVar = c(-3, -3, -3)
)

alpha_order <- length(pars$psi2)
beta_order <- length(pars$psi3)
kmed <- length(pars$deltaMedia)
kvar <- length(pars$deltaVar)
n <- length(yt) # Tamanho da serie

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 896, 1203),
                                 c(895, 1202, n)))
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
acf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
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
  med_incond = opt2$deltaMedia,
  time = DAX$Index
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

grafico_var_cond(data)
#ggsave(r"{graficos\DAX\desvio_cond_modelo2.png}", width = 20, height = 10)

grafico_var_incond(data)
#ggsave(r"{graficos\DAX\desvio_incond_modelo2.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM

# Modelo 03 ---------------------------------------------------------------

ggplot(DAX, aes(x = Index, y = dax)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "DAX com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(crises[1], 
                            as.Date(c("2008-09-15", "2008-09-20"))), 
             colour = c('red', 'red', 'red'), size = 1.5,
             linetype = "dashed")

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
dummy2 <- as.matrix(dummy_on_off(n, c(1, 896, 1203, 1208),
                                 c(895, 1202, 1207, n)))
# Ordens e Parametros - FIM

# Estimando e residuos - INICIO

(opt3 <- estimando(llike_model_garch, pars))

media_cond_mod3 <- esp_cond_model(
  data = yt,
  est = opt3,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_cond_mod3 <- var_cond_model(
  data = yt,
  est = opt3,
  dummy1 = dummy1,
  dummy2 = dummy2,
  Varyt = Varyt,
  alpha_order = alpha_order,
  beta_order = beta_order,
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_incond_mod3 <- var_indcond(
  data = yt,
  est = opt3,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
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
acf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_fac_modelo3_serie.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_facp_modelo3_serie.png}", width = 10, height = 10)

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_fac_modelo3_quad.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_facp_modelo3_quad.png}", width = 10, height = 10)
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
  med_incond = opt3$deltaMedia,
  time = DAX$Index
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

grafico_var_cond(data)
#ggsave(r"{graficos\DAX\desvio_cond_modelo3.png}", width = 20, height = 10)

grafico_var_incond(data)
#ggsave(r"{graficos\DAX\desvio_incond_modelo3.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM

# Modelo 04 ---------------------------------------------------------------

ggplot(DAX, aes(x = Index, y = dax)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "DAX com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(crises[1], 
                            as.Date(c("2008-09-15", "2008-09-20"))), 
             colour = c('red', 'grey', 'yellow'), size = 1.5,
             linetype = "dashed")

# Ordens e Parametros - INICIO
pars <- list(
  psi2 = log(.15),
  psi3 = log(.84),
  ar = .2,
  deltaMedia = 0.0,
  deltaVar = c(-3, -3, -3)
)

alpha_order <- length(pars$psi2)
beta_order <- length(pars$psi3)
kmed <- length(pars$deltaMedia)
kvar <- length(pars$deltaVar)
n <- length(yt) # Tamanho da serie
delta_ind <- 2
t_ast <- 1202
t_til <- 1208

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 896, 1208),
                                 c(895, 1202, n)))
# Ordens e Parametros - FIM

# Estimando e residuos - INICIO

(opt4 <- estimando(llike_suave, pars))


media_cond_mod4 <- esp_cond_sauve(
  data = yt,
  est = opt4,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = 2,
  t_ast = 1202,
  t_til = 1208,
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
  delta_ind = 2,
  Varyt = Varyt,
  t_ast = 1202,
  t_til = 1208,
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
  delta_ind = 2,
  t_ast = 1202,
  t_til = 1208,
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
acf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_fac_modelo4_serie.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_facp_modelo4_serie.png}", width = 10, height = 10)

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_fac_modelo4_quad.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_facp_modelo4_quad.png}", width = 10, height = 10)
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
  med_incond = opt4$deltaMedia,
  time = DAX$Index
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

grafico_var_cond(data)
#ggsave(r"{graficos\DAX\desvio_cond_modelo4.png}", width = 20, height = 10)

grafico_var_incond(data)
#ggsave(r"{graficos\DAX\desvio_incond_modelo4.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM

# Modelo 05 ---------------------------------------------------------------

ggplot(DAX, aes(x = Index, y = dax)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "DAX com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(crises[1], 
                            as.Date(c("2008-09-15", "2008-09-20", "2008-12-11", "2009-05-26"))), 
             colour = c('red', 'red', 'red', 'red', 'red'), size = 1.5,
             linetype = "dashed")

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
delta_ind <- 2
t_ast <- 1202
t_til <- 1208

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 896, 1208, 1266),
                                 c(895, 1202, 1265, n)))
# Ordens e Parametros - FIM

# Estimando e residuos - INICIO

(opt5 <- estimando(llike_suave, pars))

media_cond_mod5 <- esp_cond_sauve(
  data = yt,
  est = opt5,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = 2,
  t_ast = 1202,
  t_til = 1208,
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
  delta_ind = 2,
  Varyt = Varyt,
  t_ast = 1202,
  t_til = 1208,
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
  t_ast = 1202,
  t_til = 1208,
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
acf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_fac_modelo5_serie.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_facp_modelo5_serie.png}", width = 10, height = 10)

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_fac_modelo5_quad.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_facp_modelo5_quad.png}", width = 10, height = 10)
# FAC e FACP - FIM

poder_pred(yt, media_cond_mod5, var_cond_mod5)$rmse
cor(var_cond_mod5[-(1:50)], ((yt - media_cond_mod5)^2)[-(1:50)])^2

Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box', lag = 30)
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box', lag = 30)

(dw <- sum(diff(yt - media_cond_mod5)^2)/sum((yt - media_cond_mod5)^2))

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

moments::kurtosis(resid_pad_mod5)
moments::skewness(resid_pad_mod5)

# TH - FIM

# Graficos de linha para esp_cond e var_cond - INICIO
data <- data.frame(
  yt = yt,
  one_step_predict = media_cond_mod5,
  var_incond = var_incond_mod5,
  var_cond = var_cond_mod5,
  med_incond = opt5$deltaMedia,
  time = DAX$Index
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

grafico_var_cond(data)
#ggsave(r"{graficos\DAX\desvio_cond_modelo5.png}", width = 20, height = 10)

grafico_var_incond(data)
#ggsave(r"{graficos\DAX\desvio_incond_modelo5.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM

# Modelo 06 ---------------------------------------------------------------

ggplot(DAX, aes(x = Index, y = dax)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "DAX com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(crises[1], 
                            as.Date(c("2008-09-15", "2008-09-20", 
                                      "2008-12-11", "2009-05-26"))), 
             colour = c('red', 'grey', 'yellow', 'red', 'red'), size = 1.5,
             linetype = "dashed")

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
t_ast <- 1202
t_til <- 1208

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 896, 1208, 1266, 1376),
                                 c(895, 1202, 1265, 1375, n)))

# Ordens e Parametros - FIM

# Estimando e residuos - INICIO

(opt6 <- estimando(llike_suave, pars))

media_cond_mod6 <- esp_cond_sauve(
  data = yt,
  est = opt6,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = 2,
  t_ast = 1202,
  t_til = 1208,
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
  delta_ind = 2,
  Varyt = Varyt,
  t_ast = 1202,
  t_til = 1208,
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
  delta_ind = 2,
  t_ast = 1202,
  t_til = 1208,
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
acf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_fac_modelo6_serie.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_facp_modelo6_serie.png}", width = 10, height = 10)

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_fac_modelo6_quad.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\DAX\ger_facp_modelo6_quad.png}", width = 10, height = 10)
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
  med_incond = opt6$deltaMedia,
  time = DAX$Index
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

grafico_var_cond(data)
#ggsave(r"{graficos\DAX\desvio_cond_modelo6.png}", width = 20, height = 10)

grafico_var_incond(data)
#ggsave(r"{graficos\DAX\desvio_incond_modelo6.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM


# Modelo 07 ---------------------------------------------------------------

ggplot(DAX, aes(x = Index, y = dax)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "DAX com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(crises[1],
                            as.Date(c("2008-09-15", "2008-09-20",
                                      "2009-05-26"))), 
             colour = c('red', 'grey', 'yellow', 'red'), size = 1.5,
             linetype = "dashed")

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
t_ast <- c(1202, 1208)
t_til <- c(1208, 1376)

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 896, 1208, 1376),
                                 c(895, 1202, 1208, n)))

# Estimando e residuos - INICIO

(opt7 <- estimando(llike_suave, pars))

media_cond_mod7 <- esp_cond_sauve(
  data = yt,
  est = opt7,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = c(2, 3),
  t_ast = c(1202, 1208),
  t_til = c(1208, 1376),
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_cond_mod7 <- var_cond_sauve(
  data = yt,
  est = opt7,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = c(2, 3),
  Varyt = Varyt,
  t_ast = c(1202, 1208),
  t_til = c(1208, 1376),
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_incond_mod7 <- var_indcond_sauve(
  data = yt,
  est = opt7,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = c(2, 3),
  t_ast = c(1202, 1208),
  t_til = c(1208, 1376),
  kmed = kmed,
  kvar = kvar
)

resid_pad_mod7 <- (yt - media_cond_mod7)/sqrt(var_cond_mod7)
resid_pad_mod7 <- resid_pad_mod7[-(1:50)]

resid_pad_data <- data.frame(resid_pad = resid_pad_mod7, 
                             time = seq_along(resid_pad_mod7))
resid_pad_data <- resid_pad_data[-1, ]

plot(resid_pad_mod7, type = 'l')
plot(var_incond_mod7, type = 'l')

mean(resid_pad_data$resid_pad)
var(resid_pad_data$resid_pad)

# Estimando e residuos - FIM

# FAC e FACP - INICIO
acf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
# FAC e FACP - FIM

# Poder preditivo - INICIO
poder_pred(yt, media_cond_mod7, var_cond_mod7)$rmse
cor(var_cond_mod7[-(1:50)], ((yt - media_cond_mod7)^2)[-(1:50)])^2
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

(dw <- sum(diff(yt - media_cond_mod7)^2)/sum((yt - media_cond_mod7)^2))
moments::kurtosis(resid_pad_mod7)
moments::skewness(resid_pad_mod7)

# TH - FIM

# Graficos de linha para esp_cond e var_cond - INICIO
data <- data.frame(
  yt = yt,
  one_step_predict = media_cond_mod7,
  var_incond = var_incond_mod7,
  var_cond = var_cond_mod7,
  med_incond = opt7$deltaMedia,
  time = DAX$Index
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

grafico_var_cond(data)
#ggsave(r"{graficos\DAX\desvio_cond_modelo7.png}", width = 20, height = 10)

grafico_var_incond(data)
#ggsave(r"{graficos\DAX\desvio_incond_modelo7.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM

# Modelo 08 -------------------------------------------------------------

ggplot(DAX, aes(x = Index, y = dax)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "DAX com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(as.Date(c("2008-09-15", "2008-09-20",
                                      "2008-12-11", "2009-05-26"))), 
             colour = c('grey', 'yellow', 'red', 'red'), size = 1.5,
             linetype = "dashed")

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
delta_ind <- 1
t_ast <- 1202
t_til <- 1208

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1,  1208, 1266, 1376),
                                 c(1202, 1265, 1375, n)))
# Ordens e Parametros - FIM

# Estimando e residuos - INICIO

(opt8 <- estimando(llike_suave, pars))

media_cond_mod8 <- esp_cond_sauve(
  data = yt,
  est = opt8,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = 1,
  t_ast = 1202,
  t_til = 1208,
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_cond_mod8 <- var_cond_sauve(
  data = yt,
  est = opt8,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = 1,
  Varyt = Varyt,
  t_ast = 1202,
  t_til = 1208,
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_incond_mod8 <- var_indcond_sauve(
  data = yt,
  est = opt8,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = 1,
  t_ast = 1202,
  t_til = 1208,
  kmed = kmed,
  kvar = kvar
)

resid_pad_mod8 <- (yt - media_cond_mod8)/sqrt(var_cond_mod8)
resid_pad_mod8 <- resid_pad_mod8[-(1:50)]

resid_pad_data <- data.frame(resid_pad = resid_pad_mod8, 
                             time = seq_along(resid_pad_mod8))
resid_pad_data <- resid_pad_data[-1, ]

plot(resid_pad_mod8, type = 'l')
plot(var_incond_mod8, type = 'l')

mean(resid_pad_data$resid_pad)
var(resid_pad_data$resid_pad)

# Estimando e residuos - FIM

# FAC e FACP - INICIO
acf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
# FAC e FACP - FIM

# Poder preditivo - INICIO
poder_pred(yt, media_cond_mod8, var_cond_mod8)$rmse
cor(var_cond_mod8[-(1:50)], ((yt - media_cond_mod8)^2)[-(1:50)])^2
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

(dw <- sum(diff(yt - media_cond_mod8)^2)/sum((yt - media_cond_mod8)^2))
moments::kurtosis(resid_pad_mod8)
moments::skewness(resid_pad_mod8)

# TH - FIM

# Graficos de linha para esp_cond e var_cond - INICIO
data <- data.frame(
  yt = yt,
  one_step_predict = media_cond_mod8,
  var_incond = var_incond_mod8,
  var_cond = var_cond_mod8,
  med_incond = opt8$deltaMedia,
  time = DAX$Index
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

grafico_var_cond(data)
#ggsave(r"{graficos\DAX\desvio_cond_modelo8.png}", width = 20, height = 10)

grafico_var_incond(data)
#ggsave(r"{graficos\DAX\desvio_incond_modelo8.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM

# Modelo 09 ---------------------------------------------------------------

ggplot(DAX, aes(x = Index, y = dax)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "DAX com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(crises[1], 
                            as.Date(c("2008-09-15", "2008-09-20", 
                                      "2008-12-11", "2009-05-26"))), 
             colour = c('red', 'grey', 'yellow', 'red', 'red'), size = 1.5,
             linetype = "dashed")

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
delta_ind <- 2
t_ast <- 1202
t_til <- 1208

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 896, 1208, 1266, 1376),
                                 c(895, 1202, 1265, 1375, n)))
# Ordens e Parametros - FIM

# Estimando e residuos - INICIO

opt9 <- estimando_dax(llike_suave_dax, pars)

media_cond_mod9 <- opt9$media_cond
var_cond_mod9 <- opt9$dp_cond^2
var_incond_mod9 <- opt9$var_indcond

resid_pad_mod9 <- (yt - media_cond_mod9)/sqrt(var_cond_mod9)
resid_pad_mod9 <- resid_pad_mod9[-(1:50)]

resid_pad_data <- data.frame(resid_pad = resid_pad_mod9, 
                             time = seq_along(resid_pad_mod9))
resid_pad_data <- resid_pad_data[-1, ]

plot(resid_pad_mod9, type = 'l')
plot(var_incond_mod9, type = 'l')

mean(resid_pad_data$resid_pad)
var(resid_pad_data$resid_pad)

# Estimando e residuos - FIM

# FAC e FACP - INICIO
acf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1)) + tema
ggsave(r"{graficos\DAX\ger_fac_modelo9_serie.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1)) + tema
ggsave(r"{graficos\DAX\ger_facp_modelo9_serie.png}", width = 10, height = 10)

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1)) + tema
ggsave(r"{graficos\DAX\ger_fac_modelo9_quad.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1)) + tema
ggsave(r"{graficos\DAX\ger_facp_modelo9_quad.png}", width = 10, height = 10)
# FAC e FACP - FIM

poder_pred(yt, media_cond_mod9, var_cond_mod9)$rmse
cor(var_cond_mod9[-(1:50)], ((yt - media_cond_mod9)^2)[-(1:50)])^2

Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box', lag = 30, fitdf = 8)
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box', lag = 30, fitdf = 8)

(dw <- sum(diff(yt - media_cond_mod9)^2)/sum((yt - media_cond_mod9)^2))

# QQplot e Histograma - INICIO

grafico_qqplot(resid_pad_data)
grafico_hist(resid_pad_data)

juntando_hist_qq(resid_pad_data)
ggsave(r"{graficos\DAX\qqplot_hist_modelo9.png}", width = 20, height = 10)
# QQplot e Histograma - FIM

# TH - INICIO

shapiro.test(resid_pad_data$resid_pad)
tseries::jarque.bera.test(resid_pad_data$resid_pad)
nortest::ad.test(resid_pad_data$resid_pad)

moments::kurtosis(resid_pad_mod9)
moments::skewness(resid_pad_mod9)

# TH - FIM

# Graficos de linha para esp_cond e var_cond - INICIO
data <- data.frame(
  yt = yt,
  one_step_predict = media_cond_mod9,
  var_incond = var_incond_mod9,
  var_cond = var_cond_mod9,
  med_incond = opt9$data$deltaMedia,
  time = DAX$Index
)

data %>% 
  select(time, var_incond) %>% 
  write.csv(file=r"{dados\Volatilidade\dax_vol.csv}")

data %>% 
  select(time, var_incond) %>% 
  saveRDS(file=r"{dados\Volatilidade\dax_vol.rds}")

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

grafico_var_cond(data)
ggsave(r"{graficos\DAX\desvio_cond_modelo9.png}", width = 20, height = 10)

grafico_var_incond(data) + theme_bw() + 
  geom_vline(xintercept = c(crises[1], 
                          as.Date(c("2008-09-15",
                                    "2008-12-11", 
                                    "2009-05-26"))), 
           colour = c('darkorange', "yellow", 'green', 'green'), size = 1.5,
           linetype = "dashed") +
  geom_rect(data=data, 
            mapping=aes(xmin= as.Date("2008-09-20"), 
                        xmax=as.Date('2009-12-31'),
                        ymin=0, ymax=max(abs(yt-med_incond))), 
            color="grey", alpha=0.002) + tema +
  annotate(geom = "text",
           x = as.Date(c("2009-01-05", "2009-06-21")), y = 7.5, 
           label = c("11-12-2008", "26-05-2009"),
           color = "red", size = 10,
           angle = 90)
ggsave(r"{graficos\DAX\desvio_incond_modelo9.png}", width = 20, height = 10)

# Graficos de linha para esp_cond e var_cond - FIM

# Ajuste Fino - Modelo 06 -------------------------------------------------

# Ordens e Parametros - INICIO
pars <- list(
  psi2 = log(c(.05, .1)),
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
delta_ind <- 2
t_ast <- 1202
t_til <- 1208

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 896, 1208, 1266, 1376),
                                 c(895, 1202, 1265, 1375, n)))

# Ordens e Parametros - FIM

# Estimando e residuos - INICIO

(opt6_1 <- estimando(llike_suave, pars))


# Resultado ---------------------------------------------------------------

medidas <- function(modelo, nome){
  modelo %>% select(llike, AIC, BIC) %>% mutate(Modelo = nome) %>% 
    select(Modelo, llike, AIC, BIC)
}

resultado <- rbind(
  medidas(opt0, "opt0"),
  medidas(opt1, "opt1"),
  medidas(opt2, "opt2"),
  medidas(opt3, "opt3"),
  medidas(opt4, "opt4"),
  medidas(opt5, "opt5"),
  medidas(opt6, "opt6"),
  medidas(opt7, "opt7"),
  medidas(opt8, "opt8"),
  medidas(opt9$data, "opt9")

)

resultado
resultado %>% arrange(AIC)
resultado %>% arrange(BIC)

## Teste LR

teste_lr(opt1, opt0)
teste_lr(opt2, opt0)
teste_lr(opt3, opt0)
teste_lr(opt4, opt0)
teste_lr(opt5, opt0)
teste_lr(opt6, opt0)
teste_lr(opt7, opt0)
teste_lr(opt8$data, opt0)
teste_lr(opt9$data, opt0)

teste_lr(opt5, opt4)
teste_lr(opt10, opt5)
teste_lr(opt6, opt5)
teste_lr(opt10, opt8$data)

teste_lr(opt6, opt9$data)

## Poder preditivo
poder_pred(yt, media_cond_mod0, var_cond_mod0)$rmse
poder_pred(yt, media_cond_mod1, var_cond_mod1)$rmse
poder_pred(yt, media_cond_mod2, var_cond_mod2)$rmse
poder_pred(yt, media_cond_mod3, var_cond_mod3)$rmse
poder_pred(yt, media_cond_mod4, var_cond_mod4)$rmse
poder_pred(yt, media_cond_mod5, var_cond_mod5)$rmse
poder_pred(yt, media_cond_mod6, var_cond_mod6)$rmse
poder_pred(yt, media_cond_mod7, var_cond_mod7)$rmse
poder_pred(yt, media_cond_mod8, var_cond_mod8)$rmse
poder_pred(yt, media_cond_mod9, var_cond_mod9)$rmse

## Cor
cor(var_cond_mod0[-(1:50)], ((yt - media_cond_mod0)^2)[-(1:50)])^2
cor(var_cond_mod1[-(1:50)], ((yt - media_cond_mod1)^2)[-(1:50)])^2
cor(var_cond_mod2[-(1:50)], ((yt - media_cond_mod2)^2)[-(1:50)])^2
cor(var_cond_mod3[-(1:50)], ((yt - media_cond_mod3)^2)[-(1:50)])^2
cor(var_cond_mod4[-(1:50)], ((yt - media_cond_mod4)^2)[-(1:50)])^2
cor(var_cond_mod5[-(1:50)], ((yt - media_cond_mod5)^2)[-(1:50)])^2
cor(var_cond_mod6[-(1:50)], ((yt - media_cond_mod6)^2)[-(1:50)])^2
cor(var_cond_mod7[-(1:50)], ((yt - media_cond_mod7)^2)[-(1:50)])^2
cor(var_cond_mod8[-(1:50)], ((yt - media_cond_mod8)^2)[-(1:50)])^2
cor(var_cond_mod9[-(1:50)], ((yt - media_cond_mod9)^2)[-(1:50)])^2

