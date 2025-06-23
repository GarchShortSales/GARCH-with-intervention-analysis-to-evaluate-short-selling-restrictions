source("src/modelo_est.R")
source("src/util.R")
source("src/modelo_sim.R")
source("src/graficos.R")
source("src/Analise de reisduos.R")

# Informacoes da serie ----------------------------------------------------

# Informacoes - INICIO
n <- 2400

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))

dummy2 <- as.matrix(dummy_on_off(n, c(1, 1901),
                                 c(1900, n)))
kc <- 2
# Informacoes - FIM

# Definicao parametros verdadeiros - INICIO
pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = c(-5, -6.5))
# Definicao parametros verdadeiros - FIM

# Geracao de uma serie - INICIO
set.seed(87544)
data <- modelo(pars, dummy1, dummy2, n)
line(data, "time", "yt", "serie gerada") + 
  geom_vline(xintercept = c(1900),
             linetype="dashed", color = "blue", size = .9)
yt <-  data$yt
# Geracao de uma serie - FIM

# Otimizacao --------------------------------------------------------------

# Chute inicial - INICIO
pars_init <- list(psi2 = log(.13), psi3 = log(.85),
                  ar = .2, 
                  deltaMedia = 10,
                  deltaVar = c(-5, -6.5))
# Chute inicial - FIM

# Chamando otimizador e contando o tempo - INICIO 
opt <- estimando(llike_model_ar_dumvar, pars_init)
opt
unlist(pars)[names(opt)[1:9]]

# Analise de residuos -----------------------------------------------------

# Residuos pad - INICIO
par <- unlist(opt)
media_condyt <- esp_cond_completo(yt, unname(par), dummy1, dummy2, kc, n)
var_condyt <- var_cond_completo(yt, unname(par), dummy1, dummy2, kc, n)
resid_pad <- (yt - media_condyt)/sqrt(var_condyt)

resid_pad_data <- data.frame(resid_pad = resid_pad[51:n], 
                             time = seq_along(resid_pad)[51:n])

mean(resid_pad_data$resid_pad)
var(resid_pad_data$resid_pad)
# Residuos pad - FIM

# Graficos - INICIO
line(resid_pad_data, 'time', 'resid_pad', 'Residous padronizados')

## FAC e FACP  - INICIO
acf_plot(resid_pad_data$resid_pad, "FAC do residuo")
pacf_plot(resid_pad_data$resid_pad, "FACP do residuo")
## FAC e FACP - FIM

## FAC e FACP do quadrado - INICIO
acf_plot(resid_pad_data$resid_pad^2, "FAC do residuo")
pacf_plot(resid_pad_data$resid_pad^2, "FACP do quadrado")
## FAC e FACP - FIM

## QQplot e Histograma - INICIO
ggplot(resid_pad_data, aes(sample = resid_pad)) + 
  stat_qq() + 
  geom_abline(slope = 1, intercept = 0) + 
  tema + ylim(-6,6) + 
  scale_x_continuous(limits = c(-6, 6),  breaks = c(-6, -4, -2, 0, 2, 4, 6))

ggplot(resid_pad_data, aes(x = resid_pad)) + 
  geom_histogram(aes(y =..density..), fill = "#0c4c8a") +
  theme_minimal() + tema + 
  labs(x = "Residuos padronizados", y = 'Densidade') + 
  scale_x_continuous(limits = c(-6, 6),  breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  stat_function(fun = dnorm, args = list(0, 1), color = 'red')
## QQplot e Histograma - FIM

# TH - INICIO

## Teste de Ljung-Box - INICIO
Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box', lag = 30)
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box', lag = 30)
## Teste de Ljung-Box - FIM

## TH para normalidade - INICIO
shapiro.test(resid_pad_data$resid_pad)
tseries::jarque.bera.test(resid_pad_data$resid_pad)
nortest::ad.test(resid_pad_data$resid_pad)
## TH para normalidade - FIM

# TH - FIM

# Graficos de linha para esp_cond e var_cond - INICIO
data <- data.frame(yt = yt, one_step_predict = media_condyt, 
                   var_cond = var_condyt, time = 1:n)

ggplot(data) +
  aes(x = time, y = yt) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  tema

ggplot(data, aes(x = time, y = sqrt(var_cond))) +
  labs(y = "Variancia Condicional", x = "Tempo") + 
  geom_line(size = 1L, colour = 'red') +
  geom_line(aes(y = abs(yt - one_step_predict)), colour = "#2F97EC",
            alpha = .5) + 
  tema + xlim(1700, 2050) + 
  geom_vline(xintercept = c(1850, 1900, 1950, 2000, n),
             linetype="dashed", color = "blue", size = .9)
# Graficos de linha para esp_cond e var_cond - FIM

