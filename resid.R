source("src/modelo_est.R")
source("src/util.R")
source("src/modelo_sim.R")
source("src/graficos.R")
source("src/Analise de reisduos.R")

# Exemplo 1 ---------------------------------------------------------------

# Informacoes dos dados - INICIO
n <- 2500
dummy1 <- as.matrix(dummy_step(n, 1, "Media"));
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1500), c(1499, n)))
kc <- 2
# Informacoes dos dados - FIM

# Definicao parametros verdadeiros - INICIO
pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = c(1, 3))
# Definicao parametros verdadeiros - FIM

# Geracao de uma serie - INICIO
data <- modelo(pars, dummy1, dummy2, n)
line(data, "time", "yt", "serie gerada") + 
  geom_vline(xintercept = 1500, linetype="dashed", color = "blue", size = .9)
yt <-  data$yt
# Geracao de uma serie - FIM

# Otimizacao - INICIO 

# Chute inicial - INICIO
pars_init <- list(psi2 = log(.13), psi3 = log(.85),
                  ar = 0, 
                  deltaMedia = 10,
                  deltaVar = c(1, 3))
# Chute inicial - FIM

# Chamando otimizador e contando o tempo - INICIO 
opt <- estimando(llike_model_ar_dumvar, pars_init)
opt
unlist(pars)[names(opt)]
# Chamando otimizador e contando o tempo - FIM

# Analise de residuos - INICIO

# Residuos pad - INICIO
par <- unlist(opt)
media_condyt <- esp_cond_completo(yt, unname(par), dummy1, dummy2, kc, n)
var_condyt <- var_cond_completo(yt, unname(par), dummy1, dummy2, kc, n)
resid_pad <- (yt - media_condyt)/sqrt(var_condyt)

resid_pad_data <- data.frame(resid_pad = resid_pad, time = seq_along(resid_pad))
resid_pad_data <- resid_pad_data[-1, ]

mean(resid_pad_data$resid_pad)
var(resid_pad_data$resid_pad)
# Residuos pad - FIM

# Erro pad - INICIO

# Chute inicial - INICIO
pars_verd <- list(alpha = .13, beta = .85,
                  ar = 0, 
                  deltaMedia = 10,
                  deltaVar = c(1, 3))
# Chute inicial - FIM

pars_verd <- unlist(pars_verd)
media_condyt2 <- esp_cond_completo(yt, unname(pars_verd), dummy1, dummy2, kc, n)
var_condyt2 <- var_cond_completo(yt, unname(pars_verd), dummy1, dummy2, kc, n)
resid_pad2 <- (yt - media_condyt2)/sqrt(var_condyt2)

erro_pad_data <- data.frame(erro_pad = resid_pad2, time = seq_along(resid_pad))
erro_pad_data <- erro_pad_data[-1, ]

mean(erro_pad_data$erro_pad)
var(erro_pad_data$erro_pad)
# Erro pad - FIM

# Graficos - INICIO
line(resid_pad_data, 'time', 'resid_pad', 'Residous padronizados')
line(erro_pad_data, 'time', 'erro_pad', 'Erro padronizado')

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
Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box')
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box')
## Teste de Ljung-Box - FIM

## TH para normalidade - INICIO
shapiro.test(resid_pad_data$resid_pad)
tseries::jarque.bera.test(resid_pad_data$resid_pad)
nortest::ad.test(resid_pad_data$resid_pad)
## TH para normalidade - FIM

# TH - FIM

# Graficos de linha para esp_cond e var_cond - INICIO
data <- data.frame(yt = yt, one_step_predict = media_condyt, 
                   var_cond = var_condyt, time = 1:2500)

ggplot(data) +
  aes(x = time, y = yt) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  tema

ggplot(data) +
  aes(x = time, y = var_cond) +
  labs(y = "Tempo", x = "Variancia Condicional") + 
  geom_line(size = 1L, colour = "#0c4c8a") +
  tema
# Graficos de linha para esp_cond e var_cond - FIM

# Exemplo 2 ---------------------------------------------------------------

# Informacoes dos dados - INICIO
n <- 2500
dummy1 <- as.matrix(dummy_step(n, 1, "Media"));
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1000), c(999, n)))
kc <- 2
# Informacoes dos dados - FIM

# Definicao parametros verdadeiros - INICIO
pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = c(1, 3))
# Definicao parametros verdadeiros - FIM

# Geracao de uma serie - INICIO
data <- modelo(pars, dummy1, dummy2, n)
line(data, "time", "yt", "serie gerada") + 
  geom_vline(xintercept = 1000, linetype="dashed", color = "blue", size = .9)
yt <-  data$yt
# Geracao de uma serie - FIM

# Otimizacao - INICIO 

# Chute inicial - INICIO
pars_init <- list(psi2 = log(.13), psi3 = log(.85),
                  ar = 0, 
                  deltaMedia = 10,
                  deltaVar = c(3, 1))
# Chute inicial - FIM

# Chamando otimizador e contando o tempo - INICIO 
opt <- estimando(llike_model_ar_dumvar, pars_init)
opt
unlist(pars)[names(opt)]
# Chamando otimizador e contando o tempo - FIM

# Analise de residuos - INICIO

# Residuos pad - INICIO
par <- unlist(opt)
media_condyt <- esp_cond_completo(yt, unname(par), dummy1, dummy2, kc, n)
var_condyt <- var_cond_completo(yt, unname(par), dummy1, dummy2, kc, n)
resid_pad <- (yt - media_condyt)/sqrt(var_condyt)

resid_pad_data <- data.frame(resid_pad = resid_pad, time = seq_along(resid_pad))
resid_pad_data <- resid_pad_data[-1, ]

mean(resid_pad_data$resid_pad)
var(resid_pad_data$resid_pad)
# Residuos pad - FIM

# Erro pad - INICIO

# Chute inicial - INICIO
pars_verd <- list(alpha = .13, beta = .85,
                  ar = 0, 
                  deltaMedia = 10,
                  deltaVar = c(1, 3))
# Chute inicial - FIM

pars_verd <- unlist(pars_verd)
media_condyt2 <- esp_cond_completo(yt, unname(pars_verd), dummy1, dummy2, kc, n)
var_condyt2 <- var_cond_completo(yt, unname(pars_verd), dummy1, dummy2, kc, n)
resid_pad2 <- (yt - media_condyt2)/sqrt(var_condyt2)

erro_pad_data <- data.frame(erro_pad = resid_pad2, time = seq_along(resid_pad))
erro_pad_data <- erro_pad_data[-1, ]

mean(erro_pad_data$erro_pad)
var(erro_pad_data$erro_pad)
# Erro pad - FIM

# Graficos - INICIO
line(resid_pad_data, 'time', 'resid_pad', 'Residous padronizados')
line(erro_pad_data, 'time', 'erro_pad', 'Erro padronizado')

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
Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box')
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box')
## Teste de Ljung-Box - FIM

## TH para normalidade - INICIO
shapiro.test(resid_pad_data$resid_pad)
tseries::jarque.bera.test(resid_pad_data$resid_pad)
nortest::ad.test(resid_pad_data$resid_pad)
## TH para normalidade - FIM

# TH - FIM

# Graficos de linha para esp_cond e var_cond - INICIO
data <- data.frame(yt = yt, one_step_predict = media_condyt, 
                   var_cond = var_condyt, time = 1:2500)

ggplot(data) +
  aes(x = time, y = yt) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  tema

ggplot(data) +
  aes(x = time, y = var_cond) +
  labs(y = "Tempo", x = "Variancia Condicional") + 
  geom_line(size = 1L, colour = "#0c4c8a") +
  tema
# Graficos de linha para esp_cond e var_cond - FIM

# Exemplo 3 ---------------------------------------------------------------

# Informacoes dos dados - INICIO
n <- 2500
dummy1 <- as.matrix(dummy_step(n, 1, "Media"));
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1400, 2000), c(1399, 1999, n)))
kc <- 3
# Informacoes dos dados - FIM

# Definicao parametros verdadeiros - INICIO
pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = c(1, 3.5, 2.5))
# Definicao parametros verdadeiros - FIM

# Geracao de uma serie - INICIO
data <- modelo(pars, dummy1, dummy2, n)

line(data, "time", "yt", "serie gerada") + 
  geom_vline(xintercept = c(1400, 2000), linetype="dashed", color = "blue", size = .9)
yt <-  data$yt
# Geracao de uma serie - FIM

# Otimizacao - INICIO 

# Chute inicial - INICIO
pars_init <- list(psi2 = log(.13), psi3 = log(.85),
                  ar = 0, 
                  deltaMedia = 10,
                  deltaVar = c(1, 3.5, 2.5))
# Chute inicial - FIM

# Chamando otimizador e contando o tempo - INICIO 
opt <- estimando(llike_model_ar_dumvar, pars_init)
opt
unlist(pars)[names(opt)]
# Chamando otimizador e contando o tempo - FIM

# Analise de residuos - INICIO

# Residuos pad - INICIO
par <- unlist(opt)
media_condyt <- esp_cond_completo(yt, unname(par), dummy1, dummy2, kc, n)
var_condyt <- var_cond_completo(yt, unname(par), dummy1, dummy2, kc, n)
resid_pad <- (yt - media_condyt)/sqrt(var_condyt)

resid_pad_data <- data.frame(resid_pad = resid_pad, time = seq_along(resid_pad))
resid_pad_data <- resid_pad_data[-1, ]

mean(resid_pad_data$resid_pad)
var(resid_pad_data$resid_pad)
# Residuos pad - FIM

# Erro pad - INICIO

pars_verd <- list(alpha = .13, beta = .85,
                  ar = 0, 
                  deltaMedia = 10,
                  deltaVar = c(1, 3.5, 2.5))

pars_verd <- unlist(pars_verd)
media_condyt2 <- esp_cond_completo(yt, unname(pars_verd), dummy1, dummy2, kc, n)
var_condyt2 <- var_cond_completo(yt, unname(pars_verd), dummy1, dummy2, kc, n)
resid_pad2 <- (yt - media_condyt2)/sqrt(var_condyt2)

erro_pad_data <- data.frame(erro_pad = resid_pad2, time = seq_along(resid_pad))
erro_pad_data <- erro_pad_data[-1, ]

mean(erro_pad_data$erro_pad)
var(erro_pad_data$erro_pad)
# Erro pad - FIM

# Graficos - INICIO
line(resid_pad_data, 'time', 'resid_pad', 'Residous padronizados')
line(erro_pad_data, 'time', 'erro_pad', 'Erro padronizado')

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
Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box')
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box')
## Teste de Ljung-Box - FIM

## TH para normalidade - INICIO
shapiro.test(resid_pad_data$resid_pad)
tseries::jarque.bera.test(resid_pad_data$resid_pad)
nortest::ad.test(resid_pad_data$resid_pad)
## TH para normalidade - FIM

# TH - FIM

# Graficos de linha para esp_cond e var_cond - INICIO
data <- data.frame(yt = yt, one_step_predict = media_condyt, 
                   var_cond = var_condyt, time = 1:2500)

ggplot(data) +
  aes(x = time, y = yt) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  tema

ggplot(data) +
  aes(x = time, y = var_cond) +
  labs(y = "Tempo", x = "Variancia Condicional") + 
  geom_line(size = 1L, colour = "#0c4c8a") +
  tema
# Graficos de linha para esp_cond e var_cond - FIM

# Exemplo 4 ---------------------------------------------------------------

# Informacoes dos dados - INICIO
n <- 2500
dummy1 <- as.matrix(dummy_on_off(n, c(1, 1400, 2000), c(1399, 1999, n), "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1400, 2000), c(1399, 1999, n)))
kvar <- 3
kmed <- 3
# Informacoes dos dados - FIM

# Definicao parametros verdadeiros - INICIO
pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = c(1, 7, 3), 
             deltaVar = c(1, 3.5, 2.5))
# Definicao parametros verdadeiros - FIM

# Geracao de uma serie - INICIO
data <- modelo(pars, dummy1, dummy2, n)

line(data, "time", "yt", "serie gerada") + 
  geom_vline(xintercept = c(1400, 2000), linetype="dashed", color = "blue", size = .9)
yt <-  data$yt
Varyt <- var(yt[1:50])
# Geracao de uma serie - FIM

# Otimizacao - INICIO 

# Chute inicial - INICIO
pars_init <- list(psi2 = log(.13), psi3 = log(.85),
                  ar = 0, 
                  deltaMedia = c(1, 7, 3),
                  deltaVar = c(1, 3.5, 2.5))
# Chute inicial - FIM

# Chamando otimizador e contando o tempo - INICIO 
opt <-
  estimando_cpp_geral(
    pars_init = pars_init,
    data = yt,
    dummyMedia = dummy1,
    dummyVar = dummy2,
    kmed = kmed,
    kvar = kvar, 
    n = n,
    Varyt = Varyt
  )
opt
unlist(pars)[names(opt)]
# Chamando otimizador e contando o tempo - FIM

# Analise de residuos - INICIO

# Residuos pad - INICIO
par <- unlist(opt)

media_condyt <- esp_cond_geral(yt, unname(par), dummy1, dummy2, kmed, kvar, n)
var_condyt <- var_cond_geral(yt, unname(par), dummy1, dummy2, kmed, kvar, n)
resid_pad <- (yt - media_condyt)/sqrt(var_condyt)

resid_pad_data <- data.frame(resid_pad = resid_pad, time = seq_along(resid_pad))

mean(resid_pad_data$resid_pad)
var(resid_pad_data$resid_pad)
# Residuos pad - FIM

# Erro pad - INICIO
pars_verd <- list(alpha = .13, beta = .85,
                  ar = 0.2, 
                  deltaMedia = c(1, 7, 3), 
                  deltaVar = c(1, 3.5, 2.5))

pars_verd <- unlist(pars_verd)
media_condyt2 <- esp_cond_geral(yt, unname(pars_verd), dummy1, dummy2, kmed, kvar, n)
var_condyt2 <- var_cond_geral(yt, unname(pars_verd), dummy1, dummy2, kmed, kvar, n)
resid_pad2 <- (yt - media_condyt2)/sqrt(var_condyt2)

erro_pad_data <- data.frame(erro_pad = resid_pad2, time = seq_along(resid_pad))

mean(erro_pad_data$erro_pad)
var(erro_pad_data$erro_pad)
# Erro pad - FIM

# Graficos - INICIO
line(resid_pad_data, 'time', 'resid_pad', 'Residous padronizados')
line(erro_pad_data, 'time', 'erro_pad', 'Erro padronizado')

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
Box.test(resid_pad_data$resid_pad, type = 'Ljung-Box')
Box.test(resid_pad_data$resid_pad^2, type = 'Ljung-Box')
## Teste de Ljung-Box - FIM

## TH para normalidade - INICIO
shapiro.test(resid_pad_data$resid_pad)
tseries::jarque.bera.test(resid_pad_data$resid_pad)
nortest::ad.test(resid_pad_data$resid_pad)
## TH para normalidade - FIM
# TH - FIM

# Graficos de linha para esp_cond e var_cond - INICIO
data <- data.frame(yt = yt, one_step_predict = media_condyt, 
                   var_cond = var_condyt, time = 1:2500)

ggplot(data) +
  aes(x = time, y = yt) +
  geom_line(size = 1L, colour = "#0c4c8a") +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = 'Tempo') +
  geom_line(aes(y = one_step_predict), size = 1L, colour = "red") +
  tema

ggplot(data) +
  aes(x = time, y = var_cond) +
  labs(y = "Tempo", x = "Variancia Condicional") + 
  geom_line(size = 1L, colour = "#0c4c8a") +
  tema
# Graficos de linha para esp_cond e var_cond - FIM
