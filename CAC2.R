source("src/util.R")
source("src/Analise de reisduos.R")
source("src/modelo_est.R")
source("src/modelo8_france.R")
source("src/graficos.R")

library(zoo)
library(quantmod)
library(ggplot2)
library(dplyr)
library(imputeTS)

# Carregando dados --------------------------------------------------------

load("dados/CAC.RData")
BegSample <- '2004-01-01'
EndSample <- '2009-12-31'
crises <- as.Date(c("2007-07-01", "2008-11-28"))

CAC <- CAC %>% fortify.zoo %>% as_tibble

ggplot(CAC, aes(x = Index, y = 100 * cac)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "CAC") +
  theme_minimal()

CAC <- CAC %>%
  filter(Index >= as.Date(BegSample), Index <= as.Date(EndSample)) %>% 
  mutate(id = row_number(), cac = 100*cac) %>%
  select(Index, id, everything())

yt <- CAC$cac %>% as.vector()
Varyt <- var(yt[1:50])

# Graficos ----------------------------------------------------------------

p1 <- ggplot(CAC, aes(x = Index, y = CAC)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Preço", title = "CAC") +
  theme_minimal()

p2 <- ggplot(CAC, aes(x = Index, y = cac)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "CAC") +
  theme_minimal()

p3 <- ggplot(CAC, aes(x = Index, y = cac)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "CAC com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(crises, as.Date("2008-09-22")), 
             colour = c('red', "blue", "grey"), size = 1.5,
             linetype = "dashed")
p3 + xlim(c(as.Date("2008-04-01"), as.Date("2009-12-31")))
gridExtra::grid.arrange(p2, p3, ncol = 1)

ggplot(CAC, aes(x = Index, y = 100 * cac)) +
  geom_line(size = 1L, colour = "#112446") + 
  labs(x = "Tempo", y = "Retorno", title = "") +
  theme_minimal() + tema
ggsave(r"{graficos\CAC\fra_serie.png}", width = 20, height = 10)

acf(yt, plot = F) %>% autoplot(main='') + ylim(c(-1,1)) + ggtitle("") + tema
ggsave(r"{graficos\CAC\fra_fac_serie.png}", width = 10, height = 10)

pacf(yt, plot = F) %>% autoplot(main='') + ylim(c(-1,1)) + ggtitle("") + tema
ggsave(r"{graficos\CAC\fra_facp_serie.png}", width = 10, height = 10)

acf(yt^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1)) + ggtitle("") + tema
ggsave(r"{graficos\CAC\fra_fac_quad.png}", width = 10, height = 10)

pacf(yt^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1)) + ggtitle("") + tema
ggsave(r"{graficos\CAC\fra_facp_quad.png}", width = 10, height = 10)

# Modelo 00 AR(1)-GARCH(1,1) ----------------------------------------------

# Ordens e Parametros - INICIO
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
#ggsave(r"{graficos\CAC\fra_fac_modelo0_serie.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_facp_modelo0_serie.png}", width = 10, height = 10)

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_fac_modelo0_quad.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_facp_modelo0_quad.png}", width = 10, height = 10)
# FAC e FACP - FIM

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
  time = CAC$Index
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

ggplot(data, aes(x = time, y = sqrt(var_cond))) +
  labs(x = "Tempo", y = "Desvio Condicional") + 
  geom_line(size = 1L, colour = "red") + 
  geom_line(aes(x = time, y = abs(yt)), colour = "blue", alpha = .5)
#ggsave(r"{graficos\CAC\desvio_cond_modelo0.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM

# Modelo 01 ----------------------------------------------------------------

ggplot(CAC, aes(x = Index, y = cac)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "CAC com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(crises[1], as.Date("2008-09-22")), 
             colour = c('red', "red"), size = 1.5,
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
dummy2 <- as.matrix(dummy_on_off(n, c(1, 899, 1213),
                                 c(898, 1212, n)))
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
#ggsave(r"{graficos\CAC\fra_fac_modelo1_serie.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_facp_modelo1_serie.png}", width = 10, height = 10)

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_fac_modelo1_quad.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_facp_modelo1_quad.png}", width = 10, height = 10)
# FAC e FACP - FIM

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
  var_incond = var_incond_mod1,
  var_cond = var_cond_mod1,
  med_incond = opt1$deltaMedia,
  time = CAC$Index
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

grafico_var_cond(data)
#ggsave(r"{graficos\CAC\desvio_cond_modelo1.png}", width = 20, height = 10)

grafico_var_incond(data)
#ggsave(r"{graficos\CAC\desvio_incond_modelo1.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM

# Modelo 02 ---------------------------------------------------------------

ggplot(CAC, aes(x = Index, y = cac)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "CAC com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(crises[1], as.Date("2008-09-15"), 
                            as.Date("2008-09-22")), 
             colour = c('red', "red", "red"), size = 1.5,
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
dummy2 <- as.matrix(dummy_on_off(n, c(1, 899, 1207, 1212),
                                 c(898, 1206, 1211, n)))

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
#ggsave(r"{graficos\CAC\fra_fac_modelo2_serie.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_facp_modelo2_serie.png}", width = 10, height = 10)

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_fac_modelo2_quad.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_facp_modelo2_quad.png}", width = 10, height = 10)
# FAC e FACP - FIM

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

poder_pred(yt, media_cond_mod2, var_cond_mod2)$rmse
cor(var_cond_mod2[-(1:50)], ((yt - media_cond_mod2)^2)[-(1:50)])^2

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
  time = CAC$Index
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

grafico_var_cond(data)
#ggsave(r"{graficos\CAC\desvio_cond_modelo2.png}", width = 20, height = 10)

grafico_var_incond(data)
#ggsave(r"{graficos\CAC\desvio_incond_modelo2.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM

# Modelo 03 ---------------------------------------------------------------

ggplot(CAC, aes(x = Index, y = cac)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "CAC com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(crises[1], as.Date("2008-09-15"), 
                            as.Date("2008-09-22")), 
             colour = c('red', "grey", "yellow"), size = 1.5,
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
delta_ind = 2
t_ast = 1206
t_til = 1212

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 899, 1212, 1213),
                                 c(898, 1206, 1212, n)))

# Estimando e residuos - INICIO

(opt3 <- estimando(llike_suave, pars))

media_cond_mod3 <- esp_cond_sauve(
  data = yt,
  est = opt3,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = 2,
  t_ast = 1206,
  t_til = 1212,
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_cond_mod3 <- var_cond_sauve(
  data = yt,
  est = opt3,
  dummy1 = dummy1,
  dummy2 = dummy2,
  Varyt = Varyt,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = 2,
  t_ast = 1206,
  t_til = 1212,
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
  delta_ind = 2,
  t_ast = 1206,
  t_til = 1212,
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
#ggsave(r"{graficos\CAC\fra_fac_modelo3_serie.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_facp_modelo3_serie.png}", width = 10, height = 10)

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_fac_modelo3_quad.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_facp_modelo3_quad.png}", width = 10, height = 10)
# FAC e FACP - FIM

poder_pred(yt, media_cond_mod3, var_cond_mod3)$rmse
cor(var_cond_mod3[-(1:50)], ((yt - media_cond_mod3)^2)[-(1:50)])^2

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

poder_pred(yt, media_cond_mod3, var_cond_mod3)$rmse
cor(var_cond_mod3[-(1:50)], ((yt - media_cond_mod3)^2)[-(1:50)])^2

# TH - INICIO

Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box', lag = 30)
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box', lag = 30)

shapiro.test(resid_pad_data$resid_pad)
tseries::jarque.bera.test(resid_pad_data$resid_pad)
nortest::ad.test(resid_pad_data$resid_pad)

(dw <- sum(diff(yt - media_cond_mod3)^2)/sum((yt - media_cond_mod3)^2))
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
  time = CAC$Index
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

grafico_var_cond(data)
#ggsave(r"{graficos\CAC\desvio_cond_modelo3.png}", width = 20, height = 10)

grafico_var_incond(data)
#ggsave(r"{graficos\CAC\desvio_incond_modelo3.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM

# Modelo 04 ---------------------------------------------------------------

# Modelo Triangulo |>
ggplot(CAC, aes(x = Index, y = cac)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "CAC com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(crises[1], as.Date("2008-09-15"), 
                            as.Date("2008-09-22"), as.Date("2009-04-07")), 
             colour = c('red', "grey", "yellow", "yellow"), size = 1.5,
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
delta_ind = c(2, 3)
t_ast = c(1206, 1212)
t_til = c(1212, 1350)

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 899, 1212, 1350),
                                 c(898, 1206, 1212, n)))

# Estimando e residuos - INICIO

(opt4 <- estimando(llike_suave, pars))

media_cond_mod4 <- esp_cond_sauve(
  data = yt,
  est = opt4,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = c(2, 3),
  t_ast = c(1206, 1212),
  t_til = c(1212, 1350),
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_cond_mod4 <- var_cond_sauve(
  data = yt,
  est = opt4,
  dummy1 = dummy1,
  dummy2 = dummy2,
  Varyt = Varyt,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = c(2, 3),
  t_ast = c(1206, 1212),
  t_til = c(1212, 1350),
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
  t_ast = c(1206, 1212),
  t_til = c(1212, 1350),
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
#ggsave(r"{graficos\CAC\fra_fac_modelo4_serie.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_facp_modelo4_serie.png}", width = 10, height = 10)

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_fac_modelo4_quad.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_facp_modelo4_quad.png}", width = 10, height = 10)
# FAC e FACP - FIM

poder_pred(yt, media_cond_mod4, var_cond_mod4)$rmse
cor(var_cond_mod4[-(1:50)], ((yt - media_cond_mod4)^2)[-(1:50)])^2

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

poder_pred(yt, media_cond_mod4, var_cond_mod4)$rmse
cor(var_cond_mod4[-(1:50)], ((yt - media_cond_mod4)^2)[-(1:50)])^2

# TH - INICIO

Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box', lag = 30)
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box', lag = 30)

shapiro.test(resid_pad_data$resid_pad)
tseries::jarque.bera.test(resid_pad_data$resid_pad)
nortest::ad.test(resid_pad_data$resid_pad)

(dw <- sum(diff(yt - media_cond_mod4)^2)/sum((yt - media_cond_mod4)^2))
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
  time = CAC$Index
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

grafico_var_cond(data)
#ggsave(r"{graficos\CAC\desvio_cond_modelo4.png}", width = 20, height = 10)

grafico_var_incond(data)
#ggsave(r"{graficos\CAC\desvio_incond_modelo4.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM

# Modelo 05 ---------------------------------------------------------------

ggplot(CAC, aes(x = Index, y = cac)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "CAC com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(crises[1], as.Date("2008-09-15"), 
                            as.Date("2008-09-22"), 
                            as.Date("2008-12-18"), 
                            as.Date("2009-04-07")), 
             colour = c('red', "grey", "yellow", 'red', 'red'), size = 1.5,
             linetype = "dashed")

# Ordens e Parametros - INICIO
pars <- list(
  psi2 = log(.05),
  psi3 = log(.85),
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
t_ast <- 1207
t_til <- 1212

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 899, 1212, 1271, 1351),
                                 c(898, 1207, 1270, 1350, n)))

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
  t_ast = 1207,
  t_til = 1212,
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_cond_mod5 <- var_cond_sauve(
  data = yt,
  est = opt5,
  dummy1 = dummy1,
  dummy2 = dummy2,
  Varyt = Varyt,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = 2,
  t_ast = 1207,
  t_til = 1212,
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
  t_ast = 1207,
  t_til = 1212,
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
#ggsave(r"{graficos\CAC\fra_fac_modelo5_1_serie.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_facp_modelo5_serie.png}", width = 10, height = 10)

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_fac_modelo5_1_quad.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_facp_modelo5_1_quad.png}", width = 10, height = 10)
# FAC e FACP - FIM

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

poder_pred(yt, media_cond_mod5, var_cond_mod5)$rmse
cor(var_cond_mod5[-(1:50)], ((yt - media_cond_mod5)^2)[-(1:50)])^2

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
  med_incond = opt5$deltaMedia,
  time = CAC$Index
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

grafico_var_cond(data)
#ggsave(r"{graficos\CAC\desvio_cond_modelo5.png}", width = 20, height = 10)

grafico_var_incond(data)
#ggsave(r"{graficos\CAC\desvio_incond_modelo5.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM

# Modelo 06 ---------------------------------------------------------------

ggplot(CAC, aes(x = Index, y = cac)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "CAC com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(crises[1], as.Date("2008-09-15"), 
                            as.Date("2008-09-22"), 
                            as.Date("2008-12-18"), 
                            as.Date("2009-04-07")), 
             colour = c('red', "grey", "yellow", 'red', 'red'), size = 1.5,
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
t_ast <- 1207
t_til <- 1212

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 899, 1212, 1271),
                                 c(898, 1207, 1270, n)))

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
  t_ast = 1207,
  t_til = 1212,
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_cond_mod6 <- var_cond_sauve(
  data = yt,
  est = opt6,
  dummy1 = dummy1,
  dummy2 = dummy2,
  Varyt = Varyt,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = 2,
  t_ast = 1207,
  t_til = 1212,
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
  t_ast = 1207,
  t_til = 1212,
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
#ggsave(r"{graficos\CAC\fra_fac_modelo6_serie.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_facp_modelo6_serie.png}", width = 10, height = 10)

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_fac_modelo6_quad.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_facp_modelo6_quad.png}", width = 10, height = 10)
# FAC e FACP - FIM

poder_pred(yt, media_cond_mod6, var_cond_mod6)$rmse
cor(var_cond_mod6[-(1:50)], ((yt - media_cond_mod6)^2)[-(1:50)])^2

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

poder_pred(yt, media_cond_mod6, var_cond_mod6)$rmse
cor(var_cond_mod6[-(1:50)], ((yt - media_cond_mod6)^2)[-(1:50)])^2

# TH - INICIO

Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box', lag = 30)
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box', lag = 30)

shapiro.test(resid_pad_data$resid_pad)
tseries::jarque.bera.test(resid_pad_data$resid_pad)
nortest::ad.test(resid_pad_data$resid_pad)

(dw <- sum(diff(yt - media_cond_mod6)^2)/sum((yt - media_cond_mod6)^2))
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
  time = CAC$Index
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
#ggsave(r"{graficos\CAC\desvio_cond_modelo6.png}", width = 20, height = 10)

ggplot(data, aes(x = time, y = sqrt(var_incond))) +
  labs(x = "Tempo", y = "Desvio Incondicional") + 
  geom_line(size = 1L, colour = "red") + 
  geom_line(aes(x = time, y = abs(yt-med_incond)), colour = "blue", alpha = .5)
#ggsave(r"{graficos\CAC\desvio_incond_modelo6.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM

# Modelo 07 ---------------------------------------------------------------

ggplot(CAC, aes(x = Index, y = cac)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "CAC com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(crises[1], as.Date("2008-09-15"), 
                            as.Date("2008-09-22"), 
                            as.Date("2008-12-18"), 
                            as.Date("2009-06-19")), 
             colour = c('red', "grey", "yellow", 'red', 'red'), size = 1.5,
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
t_ast <- 1206
t_til <- 1212

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 899, 1212, 1271, 1401),
                                 c(898, 1207, 1270, 1400, n)))

# Estimando e residuos - INICIO

(opt7 <- estimando(llike_suave, pars))

media_cond_mod7 <- esp_cond_sauve(
  data = yt,
  est = opt7,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = 2,
  t_ast = 1206,
  t_til = 1212,
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_cond_mod7 <- var_cond_sauve(
  data = yt,
  est = opt7,
  dummy1 = dummy1,
  dummy2 = dummy2,
  Varyt = Varyt,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = 2,
  t_ast = 1206,
  t_til = 1212,
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
  delta_ind = 2,
  t_ast = 1206,
  t_til = 1212,
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
#ggsave(r"{graficos\CAC\fra_fac_modelo7_serie.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_facp_modelo7_serie.png}", width = 10, height = 10)

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_fac_modelo7_quad.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_facp_modelo7_quad.png}", width = 10, height = 10)
# FAC e FACP - FIM

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

poder_pred(yt, media_cond_mod7, var_cond_mod7)$rmse
cor(var_cond_mod7[-(1:50)], ((yt - media_cond_mod7)^2)[-(1:50)])^2

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
  time = CAC$Index
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

grafico_var_cond(data)
#ggsave(r"{graficos\CAC\desvio_cond_modelo7.png}", width = 20, height = 10)

grafico_var_incond(data)
#ggsave(r"{graficos\CAC\desvio_incond_modelo7.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM

# Ajuste Fino - Modelo 05 -------------------------------------------------

pars <- list(
  psi2 = log(c(.05, .05)),
  psi3 = log(c(.4, .4)),
  ar = .2,
  deltaMedia = 0.0,
  deltaVar = c(-3, -3, 3, -3, -3)
)

alpha_order <- length(pars$psi2)
beta_order <- length(pars$psi3)
kmed <- length(pars$deltaMedia)
kvar <- length(pars$deltaVar)
n <- length(yt) # Tamanho da serie
delta_ind <- 2
t_ast <- 1206
t_til <- 1212

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 899, 1212, 1271, 1351),
                                 c(898, 1206, 1270, 1350, n)))

# Estimando e residuos - INICIO

(opt5_1 <- estimando(llike_suave, pars))

media_cond_mod5_1 <- esp_cond_sauve(
  data = yt,
  est = opt5_1,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = 2,
  t_ast = 1206,
  t_til = 1212,
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_cond_mod5_1 <- var_cond_sauve(
  data = yt,
  est = opt5_1,
  dummy1 = dummy1,
  dummy2 = dummy2,
  Varyt = Varyt,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = 2,
  t_ast = 1206,
  t_til = 1212,
  kmed = kmed,
  kvar = kvar,
  n = n
)

var_incond_mod5_1 <- var_indcond_sauve(
  data = yt,
  est = opt5_1,
  dummy1 = dummy1,
  dummy2 = dummy2,
  alpha_order = alpha_order,
  beta_order = beta_order,
  delta_ind = 2,
  t_ast = 1206,
  t_til = 1212,
  kmed = kmed,
  kvar = kvar
)

resid_pad_mod5_1 <- (yt - media_cond_mod5_1)/sqrt(var_cond_mod5_1)
resid_pad_mod5_1 <- resid_pad_mod5_1[-(1:50)]

resid_pad_data <- data.frame(resid_pad = resid_pad_mod5_1, 
                             time = seq_along(resid_pad_mod5_1))
resid_pad_data <- resid_pad_data[-1, ]

plot(resid_pad_mod5_1, type = 'l')
plot(var_incond_mod5_1, type = 'l')

mean(resid_pad_data$resid_pad)
var(resid_pad_data$resid_pad)

# Estimando e residuos - FIM

# FAC e FACP - INICIO
acf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_fac_modelo5_1_serie.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_facp_modelo5_1_serie.png}", width = 10, height = 10)

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_fac_modelo5_1_quad.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1))
#ggsave(r"{graficos\CAC\fra_facp_modelo5_1_quad.png}", width = 10, height = 10)
# FAC e FACP - FIM

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

poder_pred(yt, media_cond_mod5_1, var_cond_mod5_1)$rmse
cor(var_cond_mod5_1[-(1:50)], ((yt - media_cond_mod5_1)^2)[-(1:50)])^2

# TH - INICIO

Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box', lag = 30)
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box', lag = 30)

shapiro.test(resid_pad_data$resid_pad)
tseries::jarque.bera.test(resid_pad_data$resid_pad)
nortest::ad.test(resid_pad_data$resid_pad)

(dw <- sum(diff(yt - media_cond_mod5_1)^2)/sum((yt - media_cond_mod5_1)^2))
moments::kurtosis(resid_pad_mod5_1)
moments::skewness(resid_pad_mod5_1)

# TH - FIM

# Graficos de linha para esp_cond e var_cond - INICIO
data <- data.frame(
  yt = yt,
  one_step_predict = media_cond_mod5_1,
  var_incond = var_incond_mod5_1,
  var_cond = var_cond_mod5_1,
  time = CAC$Index
)

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

ggplot(data, aes(x = time, y = sqrt(var_cond))) +
  labs(y = "Tempo", x = "Desvio Condicional") + 
  geom_line(size = 1L, colour = "red") + 
  geom_line(aes(x = time, y = abs(yt)), colour = "blue")

ggplot(data, aes(x = time, y = sqrt(var_incond))) +
  labs(x = "Tempo", y = "Desvio Incondicional") + 
  geom_line(size = 1L, colour = "red") + 
  geom_line(aes(x = time, y = abs(yt)), colour = "blue", alpha = .5)
# Graficos de linha para esp_cond e var_cond - FIM

# Modelo 08 ---------------------------------------------------------------

ggplot(CAC, aes(x = Index, y = cac)) +
  geom_line(size = 1L, colour = "#112446") +
  labs(x = "Tempo", y = "Retorno", title = "CAC com as Restrições") +
  theme_minimal() + 
  geom_vline(xintercept = c(crises[1], as.Date("2008-09-15"), 
                            as.Date("2008-09-22"), 
                            as.Date("2008-12-18"), 
                            as.Date("2009-06-19")), 
             colour = c('red', "grey", "yellow", 'red', 'red'), size = 1.5,
             linetype = "dashed")

# Ordens e Parametros - INICIO
pars <- list(
  psi2 = log(.15),
  psi3 = log(.84),
  ar = .2,
  deltaMedia = 0.0,
  deltaVar = c(-3, -3, -1, -1) # Chute 
)

alpha_order <- length(pars$psi2)
beta_order <- length(pars$psi3)
kmed <- length(pars$deltaMedia)
kvar <- length(pars$deltaVar)
n <- length(yt) # Tamanho da serie
delta_ind <- 2
t_ast <- 1207
t_til <- 1212

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 899, 1212, 1271, 1351),
                                 c(898, 1207, 1270, 1350, n)))

# Estimando e residuos - INICIO

opt8 <- estimando_france(llike_suave_france, pars)

media_cond_mod8 <- opt8$media_cond
var_cond_mod8 <- opt8$dp_cond^2
var_incond_mod8 <- opt8$var_indcond

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
acf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1)) + tema
ggsave(r"{graficos\CAC\fra_fac_modelo8_serie.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad, plot = F) %>% autoplot(main='') + ylim(c(-1,1))+ tema
ggsave(r"{graficos\CAC\fra_facp_modelo8_serie.png}", width = 10, height = 10)

acf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1)) + tema
ggsave(r"{graficos\CAC\fra_fac_modelo8_quad.png}", width = 10, height = 10)
pacf(resid_pad_data$resid_pad^2, plot = F) %>% autoplot(main='') + ylim(c(-1,1)) + tema
ggsave(r"{graficos\CAC\fra_facp_modelo8_quad.png}", width = 10, height = 10)
# FAC e FACP - FIM

# QQplot e Histograma - INICIO
grafico_qqplot(resid_pad_data)
grafico_hist(resid_pad_data)

juntando_hist_qq(resid_pad_data)
ggsave(r"{graficos\CAC\qqplot_hist_modelo8.png}", width = 20, height = 10)

# QQplot e Histograma - FIM

poder_pred(yt, media_cond_mod8, var_cond_mod8)$rmse
cor(var_cond_mod8[-(1:50)], ((yt - media_cond_mod8)^2)[-(1:50)])^2

# TH - INICIO

Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box', lag = 30, fitdf = 8)
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box', lag = 30, fitdf = 8)

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
  med_incond = opt8$data$deltaMedia,
  time = CAC$Index
)

data %>% 
  select(time, var_incond) %>% 
  write.csv(file=r"{dados\Volatilidade\cac_vol.csv}")

data %>% 
  select(time, var_incond) %>% 
  saveRDS(file=r"{dados\Volatilidade\cac_vol.rds}")
  

ggplot(data, aes(x = time, y = yt)) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') 

grafico_var_cond(data)
ggsave(r"{graficos\CAC\desvio_cond_modelo8.png}", width = 20, height = 10)

grafico_var_incond(data) + 
  theme_bw() + 
  geom_vline(xintercept = c(crises[1], as.Date("2008-09-15"),
                            as.Date("2008-12-18"), 
                            as.Date("2009-04-07")), 
             colour = c('darkorange', "yellow", 'green', 'green'), size = 1.5,
             linetype = "dashed") +
  geom_rect(data=data, 
            mapping=aes(xmin= as.Date("2008-09-22"), 
                        xmax=as.Date('2009-12-31'),
                        ymin=0, ymax=max(abs(yt-med_incond))), 
            color="grey", alpha=0.003) + tema +
  annotate(geom = "text",
           x = as.Date(c("2009-01-13", "2009-05-2")), y = 7.5, 
           label = c("18-12-2008", "07-04-2009"),
           color = "red", size = 10,
           angle = 90)
ggsave(r"{graficos\CAC\desvio_incond_modelo8.png}", width = 20, height = 10)
# Graficos de linha para esp_cond e var_cond - FIM

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
                   medidas(opt8$data, "opt8")
)
resultado
resultado %>% arrange(AIC)
resultado %>% arrange(BIC)

## Teste LR
teste_lr(opt2, opt1)
teste_lr(opt5, opt6)
teste_lr(opt5, opt8$data)

teste_lr(opt7, opt0)
teste_lr(opt8$data, opt0)



## Poder preditivo
poder_pred(yt, media_cond_mod1, var_cond_mod1)$rmse
poder_pred(yt, media_cond_mod2, var_cond_mod2)$rmse
poder_pred(yt, media_cond_mod3, var_cond_mod3)$rmse
poder_pred(yt, media_cond_mod4, var_cond_mod4)$rmse
poder_pred(yt, media_cond_mod5, var_cond_mod5)$rmse
poder_pred(yt, media_cond_mod6, var_cond_mod6)$rmse
poder_pred(yt, media_cond_mod7, var_cond_mod7)$rmse
poder_pred(yt, media_cond_mod8, var_cond_mod8)$rmse

## Cor
cor(var_cond_mod1[-(1:50)], ((yt - media_cond_mod1)^2)[-(1:50)])^2
cor(var_cond_mod2[-(1:50)], ((yt - media_cond_mod2)^2)[-(1:50)])^2
cor(var_cond_mod3[-(1:50)], ((yt - media_cond_mod3)^2)[-(1:50)])^2
cor(var_cond_mod4[-(1:50)], ((yt - media_cond_mod4)^2)[-(1:50)])^2
cor(var_cond_mod5[-(1:50)], ((yt - media_cond_mod5)^2)[-(1:50)])^2
cor(var_cond_mod6[-(1:50)], ((yt - media_cond_mod6)^2)[-(1:50)])^2
cor(var_cond_mod7[-(1:50)], ((yt - media_cond_mod7)^2)[-(1:50)])^2
cor(var_cond_mod8[-(1:50)], ((yt - media_cond_mod8)^2)[-(1:50)])^2