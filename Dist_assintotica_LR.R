source("src/modelo_sim.R")
source("src/util.R")
source("src/modelo_est.R")
source("src/graficos.R")

# Exemplo 1 t = 50 ----------------------------------------------------------------

# Definindo algumas estruturas das series - INICIO
n <- 2000

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 50), c(49, n)))
dummy2_red <- as.matrix(dummy_step(n, 1))

kc <- 2; kc_red <- 1

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = log(5))

## Chute inicial - INICIO

pars_init_red <- list(psi2 = log(.13), psi3 = log(.85),
                      ar = .2, 
                      deltaMedia = 10, 
                      deltaVar = log(5))

pars_init_comp <- list(psi2 = log(.13), psi3 = log(.85),
                       ar = .2, 
                       deltaMedia = 10, 
                       deltaVar = c(log(5), log(5)))

## Chute inicial - FIM

# Inciando MC - INICIO

## Definindo estruturas - INICIO
M <- 500
modelo_completo <- data.frame(alpha_com = rep(NA, M), beta_com = rep(NA, M), 
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M),
                              deltaVar2_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M),
                              cont_com = rep(NA, M))

modelo_reduzido <- data.frame(alpha = rep(NA, M), beta = rep(NA, M), 
                              ar = rep(NA, M),
                              deltaMedia = rep(NA, M),
                              deltaVar1 = rep(NA, M),
                              ite = rep(NA, M), llike_red = rep(NA, M),
                              cont_rep = rep(NA, M))
## Definindo estruturas - FIM

set.seed(95423)
inicio <- Sys.time()
for (cont in 1:M){
  
  data <- modelo(pars, dummy1, dummy2_red, n)
  yt <- data$yt
  
  modelo_completo[cont, ] <- c(estimando(llike_model_ar_dumvar, pars_init_comp), cont)
  modelo_reduzido[cont, ] <- c(estimando(llike_model_ar_dumvar_red, pars_init_red), cont)
  
  print(cont)
  
}
fim <- Sys.time()
fim - inicio
# Inciando MC - FIM

## Grafico e medidas descritivas das estimativas - INICIO
modelo_completo_ord <- modelo_completo %>% arrange(-ite_com)
modelo_completo_ret <- modelo_completo_ord[floor(.01*M):M, ]

media_est <- apply(modelo_completo_ret, 2, mean); media_est
sd_est <- apply(modelo_completo_ret, 2, sd); sd_est
cor_est <- cor(modelo_completo_ret[, 1:6]); cor_est
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar")]

est_pad <- modelo_completo_ret %>% select(alpha_com:deltaVar2_com) %>% pad()

## QQplot e histograma - INICIO
q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n))
h1 <- map(names(est_pad), ~histo_MC(est_pad, .x, M, n))
## QQplot e histograma - FIM

## Grafico e medidas descritivas das estimativas - FIM

# Estudando convergencia - INICIO
est_data <- est_data %>% select(alpha_com:deltaVar2_com)
media_est <- est_data %>% apply(2, mean); media_est
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar")]

sd_est <- est_data %>% apply(2, sd)
cv <- sd_est/media_est; cv
cor(est_data[, 1:6])

# Acessando distribuicao assintotica - INICIO
llike_modelos_retirado <- llike_retirado(modelo_completo, modelo_reduzido)

apply(llike_modelos_retiradoe[c('rejeita_h0_1', 
                                'rejeita_h0_5', 
                                'rejeita_h0_10')],
      2, mean) 


# Teste da razao - FIM

## Histogramas e QQplots - INICIO
QQplot_LR(llike_modelos_retirado, razao_llike_est, M, n)

histo_LR(llike_modelos_retirado, razao_llike_est, M, n)
## Histogramas e QQplots - FIM

# Acessando distribuicao assintotica - FIM

# Exemplo 2 t = 1950 ----------------------------------------------------------------

# Definindo algumas estruturas das series - INICIO
n <- 2000

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1950), c(1949, n)))
dummy2_red <- as.matrix(dummy_step(n, 1))

kc <- 2; kc_red <- 1

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = log(5))

## Chute inicial - INICIO

pars_init_red <- list(psi2 = log(.13), psi3 = log(.85),
                      ar = .2, 
                      deltaMedia = 10, 
                      deltaVar = log(5))

pars_init_comp <- list(psi2 = log(.13), psi3 = log(.85),
                       ar = .2, 
                       deltaMedia = 10, 
                       deltaVar = c(log(5), log(5)))

## Chute inicial - FIM

# Inciando MC - INICIO

## Definindo estruturas - INICIO
M <- 500
modelo_completo <- data.frame(alpha_com = rep(NA, M), beta_com = rep(NA, M), 
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M),
                              deltaVar2_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M))

modelo_red <- data.frame(alpha = rep(NA, M), beta = rep(NA, M), 
                         ar = rep(NA, M),
                         deltaMedia = rep(NA, M),
                         deltaVar1 = rep(NA, M),
                         ite = rep(NA, M), llike_red = rep(NA, M))
## Definindo estruturas - FIM

set.seed(8000)
inicio <- Sys.time()
for (cont in 1:M){
  
  yt <- modelo(pars, dummy1, dummy2_red, n)$yt
  
  modelo_completo[cont, ] <- estimando(llike_model_ar_dumvar, pars_init_comp)
  modelo_red[cont, ] <- estimando(llike_model_ar_dumvar_red, pars_init_red)
  
  print(cont)
  
}
fim <- Sys.time()
fim - inicio
# Inciando MC - FIM

# Tirando as primeiras estimativas - INICIO
est_data <- modelo_completo %>% arrange(-ite_com) 
est_data <- est_data[floor(.2*M):M, ]
# Tirando as primeiras estimativas- FIM

# Estudando convergencia - INICIO
est_data <- est_data %>% select(alpha_com:deltaVar2_com)
media_est <- est_data %>% apply(2, mean); media_est
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar")]

sd_est <- est_data %>% apply(2, sd)
cv <- sd_est/media_est; cv
cor(est_data[, 1:6])

# Acessando distribuicao assintotica - INICIO
valor_critico <- qchisq(.05, 1, lower.tail = FALSE)

llike_modelos <- cbind(modelo_completo, modelo_red) %>% 
  select(llike_com, llike, ite_com, everything(), -ite) %>% 
  mutate(
    comparando_llike = llike_com > llike,
    razao_llike_est = 2 * (llike_com - llike),
    rejeita_h0 = razao_llike_est > valor_critico
  )

llike_modelos <- llike_modelos %>% arrange(-ite_com) 
llike_modelos_retirado <- llike_modelos[floor(.2*M):M, ]
## Guardando as log-likes e medidas para dos modelos reduzido e completo - FIM 

sum(llike_modelos_retirado$rejeita_h0)/length(llike_modelos_retirado$rejeita_h0) 
# Teste da razao - FIM

# Graficos - INICIO
ggplot(llike_modelos_retirado, aes(sample = razao_llike_est)) + 
  stat_qq(distribution = stats::qchisq, dparams = 1) + 
  stat_qq_line(distribution = stats::qchisq, dparams = 1) +
  tema +
  labs(x = 'Quantil Teorico', y = 'Quantil Amostral', 
       title = glue("QQplot com M = {M} e n = {n}"))

ggplot(llike_modelos_retirado, aes(x = razao_llike_est))+ 
  geom_histogram(aes(y = ..density..), fill = "#0c4c8a") +
  tema + 
  labs(y = '', x = '',
       title = glue("QQplot com M = {M} e n = {n}")) + 
  theme(axis.title.y = element_text(angle=0, size = 15, vjust = .6)) +
  scale_x_continuous(limits = c(-1, 10),  breaks = seq(0, 10, 2)) +
  stat_function(fun = dchisq, args = list(1), color = 'red')
# Acessando distribuicao assintotica - FIM
# Exemplo 3 t = 1000 ---------------------------------------------------------------

# Definindo algumas estruturas das series - INICIO
n <- 2000

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1000), c(999, n)))
dummy2_red <- as.matrix(dummy_step(n, 1))

kc <- 2; kc_red <- 1

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = log(5))

## Chute inicial - INICIO

pars_init_red <- list(psi2 = log(.13), psi3 = log(.85),
                      ar = .2, 
                      deltaMedia = 10, 
                      deltaVar = log(5))

pars_init_comp <- list(psi2 = log(.13), psi3 = log(.85),
                       ar = .2, 
                       deltaMedia = 10, 
                       deltaVar = c(log(5), log(5)))

## Chute inicial - FIM

# Inciando MC - INICIO

## Definindo estruturas - INICIO
M <- 500
modelo_completo <- data.frame(alpha_com = rep(NA, M), beta_com = rep(NA, M), 
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M),
                              deltaVar2_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M))

modelo_red <- data.frame(alpha = rep(NA, M), beta = rep(NA, M), 
                         ar = rep(NA, M),
                         deltaMedia = rep(NA, M),
                         deltaVar1 = rep(NA, M),
                         ite = rep(NA, M), llike_red = rep(NA, M))
## Definindo estruturas - FIM

set.seed(8000)
inicio <- Sys.time()
for (cont in 1:M){
  
  yt <- modelo(pars, dummy1, dummy2_red, n)$yt
  
  modelo_completo[cont, ] <- estimando(llike_model_ar_dumvar, pars_init_comp)
  modelo_red[cont, ] <- estimando(llike_model_ar_dumvar_red, pars_init_red)
  
  print(cont)
  
}
fim <- Sys.time()
fim - inicio
# Inciando MC - FIM

# Tirando as primeiras estimativas - INICIO
est_data <- modelo_completo %>% arrange(ite_com) 
est_data <- est_data[floor(.2*M):M, ]
# Tirando as primeiras estimativas- FIM

# Estudando convergencia - INICIO
est_data <- est_data %>% select(alpha_com:deltaVar2_com)
media_est <- est_data %>% apply(2, mean); media_est
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar")]

sd_est <- est_data %>% apply(2, sd)
cv <- sd_est/media_est; cv
cor(est_data[, 1:6])

# Acessando distribuicao assintotica - INICIO
valor_critico <- qchisq(.05, 1, lower.tail = FALSE)

llike_modelos <- cbind(modelo_completo, modelo_red) %>% 
  select(llike_com, llike, ite_com, everything(), -ite) %>% 
  mutate(
    comparando_llike = llike_com > llike,
    razao_llike_est = 2 * (llike_com - llike),
    rejeita_h0 = razao_llike_est > valor_critico
  )

llike_modelos <- llike_modelos %>% arrange(-ite_com) 
llike_modelos_retirado <- llike_modelos[floor(.2*M):M, ]
## Guardando as log-likes e medidas para dos modelos reduzido e completo - FIM 

sum(llike_modelos_retirado$rejeita_h0)/length(llike_modelos_retirado$rejeita_h0) 
# Teste da razao - FIM

# Graficos - INICIO
ggplot(llike_modelos_retirado, aes(sample = razao_llike_est)) + 
  stat_qq(distribution = stats::qchisq, dparams = 1) + 
  stat_qq_line(distribution = stats::qchisq, dparams = 1) +
  tema +
  labs(x = 'Quantil Teorico', y = 'Quantil Amostral', 
       title = glue("QQplot com M = {M} e n = {n}"))

ggplot(llike_modelos_retirado, aes(x = razao_llike_est))+ 
  geom_histogram(aes(y = ..density..), fill = "#0c4c8a") +
  tema + 
  labs(y = '', x = '',
       title = glue("QQplot com M = {M} e n = {n}")) + 
  theme(axis.title.y = element_text(angle=0, size = 15, vjust = .6)) +
  scale_x_continuous(limits = c(-1, 10),  breaks = seq(0, 10, 2)) +
  stat_function(fun = dchisq, args = list(1), color = 'red')
# Acessando distribuicao assintotica - FIM
# Exemplo 4 t = 300 ---------------------------------------------------------------

# Definindo algumas estruturas das series - INICIO
n <- 2000

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 300), c(299, n)))
dummy2_red <- as.matrix(dummy_step(n, 1))

kc <- 2; kc_red <- 1

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = log(5))

## Chute inicial - INICIO

pars_init_red <- list(psi2 = log(.13), psi3 = log(.85),
                      ar = .2, 
                      deltaMedia = 10, 
                      deltaVar = log(5))

pars_init_comp <- list(psi2 = log(.13), psi3 = log(.85),
                       ar = .2, 
                       deltaMedia = 10, 
                       deltaVar = c(log(5), log(5)))

## Chute inicial - FIM

# Inciando MC - INICIO

## Definindo estruturas - INICIO
M <- 1500
modelo_completo <- data.frame(alpha_com = rep(NA, M), beta_com = rep(NA, M), 
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M),
                              deltaVar2_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M))

modelo_red <- data.frame(alpha = rep(NA, M), beta = rep(NA, M), 
                         ar = rep(NA, M),
                         deltaMedia = rep(NA, M),
                         deltaVar1 = rep(NA, M),
                         ite = rep(NA, M), llike_red = rep(NA, M))
## Definindo estruturas - FIM

set.seed(851966)
inicio <- Sys.time()
for (cont in 1:M){
  
  yt <- modelo(pars, dummy1, dummy2_red, n)$yt
  
  modelo_completo[cont, ] <- estimando(llike_model_ar_dumvar, pars_init_comp)
  modelo_red[cont, ] <- estimando(llike_model_ar_dumvar_red, pars_init_red)
  
  print(cont)
  
}
fim <- Sys.time()
fim - inicio
# Inciando MC - FIM

# Tirando as primeiras estimativas - INICIO
est_data <- modelo_completo %>% arrange(ite_com) 
est_data <- est_data[floor(.2*M):M, ]
# Tirando as primeiras estimativas- FIM

# Estudando convergencia - INICIO
est_data <- est_data %>% select(alpha_com:deltaVar2_com)
media_est <- est_data %>% apply(2, mean); media_est
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar")]

sd_est <- est_data %>% apply(2, sd)
cv <- sd_est/media_est; cv
cor(est_data[, 1:6])

# Acessando distribuicao assintotica - INICIO
valor_critico <- qchisq(.05, 1, lower.tail = FALSE)

llike_modelos <- cbind(modelo_completo, modelo_red) %>% 
  select(llike_com, llike, ite_com, everything(), -ite) %>% 
  mutate(
    comparando_llike = llike_com > llike,
    razao_llike_est = 2 * (llike_com - llike),
    rejeita_h0 = razao_llike_est > valor_critico
  )

llike_modelos <- llike_modelos %>% arrange(-ite_com) 
llike_modelos_retirado <- llike_modelos[floor(.2*M):M, ]
## Guardando as log-likes e medidas para dos modelos reduzido e completo - FIM 

sum(llike_modelos_retirado$rejeita_h0)/length(llike_modelos_retirado$rejeita_h0) 
# Teste da razao - FIM

# Graficos - INICIO
ggplot(llike_modelos_retirado, aes(sample = razao_llike_est)) + 
  stat_qq(distribution = stats::qchisq, dparams = 1) + 
  stat_qq_line(distribution = stats::qchisq, dparams = 1) +
  tema +
  labs(x = 'Quantil Teorico', y = 'Quantil Amostral', 
       title = glue("QQplot com M = {M} e n = {n}"))

ggplot(llike_modelos_retirado, aes(x = razao_llike_est))+ 
  geom_histogram(aes(y = ..density..), fill = "#0c4c8a") +
  tema + 
  labs(y = '', x = '',
       title = glue("QQplot com M = {M} e n = {n}")) + 
  theme(axis.title.y = element_text(angle=0, size = 15, vjust = .6)) +
  scale_x_continuous(limits = c(-1, 10),  breaks = seq(0, 10, 2)) +
  stat_function(fun = dchisq, args = list(1), color = 'red')
# Acessando distribuicao assintotica - FIM

# Exemplo 5 t = 1500 ---------------------------------------------------------------

# Definindo algumas estruturas das series - INICIO
n <- 2000

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1500), c(1499, n)))
dummy2_red <- as.matrix(dummy_step(n, 1))

kc <- 2; kc_red <- 1

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = log(5))

## Chute inicial - INICIO

pars_init_red <- list(psi2 = log(.13), psi3 = log(.85),
                      ar = .2, 
                      deltaMedia = 10, 
                      deltaVar = log(5))

pars_init_comp <- list(psi2 = log(.13), psi3 = log(.85),
                       ar = .2, 
                       deltaMedia = 10, 
                       deltaVar = c(log(5), log(5)))

## Chute inicial - FIM

# Inciando MC - INICIO

## Definindo estruturas - INICIO
M <- 500
modelo_completo <- data.frame(alpha_com = rep(NA, M), beta_com = rep(NA, M), 
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M),
                              deltaVar2_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M))

modelo_red <- data.frame(alpha = rep(NA, M), beta = rep(NA, M), 
                         ar = rep(NA, M),
                         deltaMedia = rep(NA, M),
                         deltaVar1 = rep(NA, M),
                         ite = rep(NA, M), llike_red = rep(NA, M))
## Definindo estruturas - FIM

set.seed(79600)
inicio <- Sys.time()
for (cont in 1:M){
  
  yt <- modelo(pars, dummy1, dummy2_red, n)$yt
  
  modelo_completo[cont, ] <- estimando(llike_model_ar_dumvar, pars_init_comp)
  modelo_red[cont, ] <- estimando(llike_model_ar_dumvar_red, pars_init_red)
  
  print(cont)
  
}
fim <- Sys.time()
fim - inicio
# Inciando MC - FIM

# Tirando as primeiras estimativas - INICIO
est_data <- modelo_completo %>% arrange(ite_com) 
est_data <- est_data[floor(.2*M):M, ]
# Tirando as primeiras estimativas- FIM

# Estudando convergencia - INICIO
est_data <- est_data %>% select(alpha_com:deltaVar2_com)
media_est <- est_data %>% apply(2, mean); media_est
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar")]

sd_est <- est_data %>% apply(2, sd)
cv <- sd_est/media_est; cv
cor(est_data[, 1:6])

# Acessando distribuicao assintotica - INICIO
valor_critico <- qchisq(.05, 1, lower.tail = FALSE)

llike_modelos <- cbind(modelo_completo, modelo_red) %>% 
  select(llike_com, llike, ite_com, everything(), -ite) %>% 
  mutate(
    comparando_llike = llike_com > llike,
    razao_llike_est = 2 * (llike_com - llike),
    rejeita_h0 = razao_llike_est > valor_critico
  )

llike_modelos <- llike_modelos %>% arrange(-ite_com) 
llike_modelos_retirado <- llike_modelos[floor(.2*M):M, ]
## Guardando as log-likes e medidas para dos modelos reduzido e completo - FIM 

sum(llike_modelos_retirado$rejeita_h0)/length(llike_modelos_retirado$rejeita_h0) 
# Teste da razao - FIM

# Graficos - INICIO
ggplot(llike_modelos_retirado, aes(sample = razao_llike_est)) + 
  stat_qq(distribution = stats::qchisq, dparams = 1) + 
  stat_qq_line(distribution = stats::qchisq, dparams = 1) +
  tema +
  labs(x = 'Quantil Teorico', y = 'Quantil Amostral', 
       title = glue("QQplot com M = {M} e n = {n}"))

ggplot(llike_modelos_retirado, aes(x = razao_llike_est))+ 
  geom_histogram(aes(y = ..density..), fill = "#0c4c8a") +
  tema + 
  labs(y = '', x = '',
       title = glue("QQplot com M = {M} e n = {n}")) + 
  theme(axis.title.y = element_text(angle=0, size = 15, vjust = .6)) +
  scale_x_continuous(limits = c(-1, 10),  breaks = seq(0, 10, 2)) +
  stat_function(fun = dchisq, args = list(1), color = 'red')
# Acessando distribuicao assintotica - FIM


# Exemplo 5 t = 1500, n = 5000 ---------------------------------------------------------------

# Definindo algumas estruturas das series - INICIO
n <- 5000

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1500), c(1499, n)))
dummy2_red <- as.matrix(dummy_step(n, 1))

kc <- 2; kc_red <- 1

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = log(5))

## Chute inicial - INICIO

pars_init_red <- list(psi2 = log(.13), psi3 = log(.85),
                      ar = .2, 
                      deltaMedia = 10, 
                      deltaVar = log(5))

pars_init_comp <- list(psi2 = log(.13), psi3 = log(.85),
                       ar = .2, 
                       deltaMedia = 10, 
                       deltaVar = c(log(5), log(5)))

## Chute inicial - FIM

# Inciando MC - INICIO

## Definindo estruturas - INICIO
M <- 500
modelo_completo <- data.frame(alpha_com = rep(NA, M), beta_com = rep(NA, M), 
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M),
                              deltaVar2_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M))

modelo_red <- data.frame(alpha = rep(NA, M), beta = rep(NA, M), 
                         ar = rep(NA, M),
                         deltaMedia = rep(NA, M),
                         deltaVar1 = rep(NA, M),
                         ite = rep(NA, M), llike_red = rep(NA, M))
## Definindo estruturas - FIM

set.seed(79600)
inicio <- Sys.time()
for (cont in 1:M){
  
  yt <- modelo(pars, dummy1, dummy2_red, n)$yt
  
  modelo_completo[cont, ] <- estimando(llike_model_ar_dumvar, pars_init_comp)
  modelo_red[cont, ] <- estimando(llike_model_ar_dumvar_red, pars_init_red)
  
  print(cont)
  
}
fim <- Sys.time()
fim - inicio
# Inciando MC - FIM

# Tirando as primeiras estimativas - INICIO
est_data <- modelo_completo %>% arrange(ite_com) 
est_data <- est_data[floor(.2*M):M, ]
# Tirando as primeiras estimativas- FIM

# Estudando convergencia - INICIO
est_data <- est_data %>% select(alpha_com:deltaVar2_com)
media_est <- est_data %>% apply(2, mean); media_est
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar")]

sd_est <- est_data %>% apply(2, sd)
cv <- sd_est/media_est; cv
cor(est_data[, 1:6])

# Acessando distribuicao assintotica - INICIO
valor_critico <- qchisq(.05, 1, lower.tail = FALSE)

llike_modelos <- cbind(modelo_completo, modelo_red) %>% 
  select(llike_com, llike, ite_com, everything(), -ite) %>% 
  mutate(
    comparando_llike = llike_com > llike,
    razao_llike_est = 2 * (llike_com - llike),
    rejeita_h0 = razao_llike_est > valor_critico
  )

llike_modelos <- llike_modelos %>% arrange(-ite_com) 
llike_modelos_retirado <- llike_modelos[floor(.2*M):M, ]
## Guardando as log-likes e medidas para dos modelos reduzido e completo - FIM 

sum(llike_modelos_retirado$rejeita_h0)/length(llike_modelos_retirado$rejeita_h0) 
# Teste da razao - FIM

# Graficos - INICIO
ggplot(llike_modelos_retirado, aes(sample = razao_llike_est)) + 
  stat_qq(distribution = stats::qchisq, dparams = 1) + 
  stat_qq_line(distribution = stats::qchisq, dparams = 1) +
  tema +
  labs(x = 'Quantil Teorico', y = 'Quantil Amostral', 
       title = glue("QQplot com M = {M} e n = {n}"))

ggplot(llike_modelos_retirado, aes(x = razao_llike_est))+ 
  geom_histogram(aes(y = ..density..), fill = "#0c4c8a") +
  tema + 
  labs(y = '', x = '',
       title = glue("QQplot com M = {M} e n = {n}")) + 
  theme(axis.title.y = element_text(angle=0, size = 15, vjust = .6)) +
  scale_x_continuous(limits = c(-1, 10),  breaks = seq(0, 10, 2)) +
  stat_function(fun = dchisq, args = list(1), color = 'red')
# Acessando distribuicao assintotica - FIM
