source("src/modelo_sim.R")
source("src/util.R")
source("src/modelo_est.R")
source("src/graficos.R")

# M de Monte Carlo - INICIO
M <- 200
# M de Monte Carlo - FIM

# Definindo algumas estruturas das series - INICIO
kc <- 2
kc_red <- 1

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = -5)

## Chute inicial - INICIO
pars_init_completo <- list(psi1 = log(1), psi2 = log(.13), psi3 = log(.85),
                           ar = .2, 
                           deltaMedia = 10, 
                           deltaVar = c(-5, -5))

pars_init_reduzido <- list(psi1 = log(1), psi2 = log(.13), psi3 = log(.85),
                           ar = .2, 
                           deltaMedia = 10, 
                           deltaVar = -5)
## Chute inicial - FIM

# Definindo algumas estruturas das series - INICIO

# Caso 1,1 ---------------------------------------------------------------

## Definindo estruturas - INICIO
modelo_completo <- data.frame(omega_com = rep(NA, M), alpha_com = rep(NA, M), 
                              beta_com = rep(NA, M), 
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M),
                              deltaVar2_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M))

modelo_reduzido <- data.frame(omega = rep(NA, M), alpha = rep(NA, M), 
                              beta = rep(NA, M), 
                              ar = rep(NA, M), 
                              deltaMedia = rep(NA, M), 
                              deltaVar1 = rep(NA, M),
                              ite = rep(NA, M), llike_red = rep(NA, M))
## Definindo estruturas - FIM

n <- 3000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))

dummy2 <- as.matrix(dummy_step(n, c(1, 50)))
dummy2_red <- as.matrix(dummy_step(n, 1))

inicio <- Sys.time()
for (cont in 1:M){
  
  yt <- modelo(pars, dummy1, dummy2_red, n)$yt
  
  modelo_reduzido[cont, ] <- estimando2(llike_model_ar_dumvar_red_ant, 
                                        pars_init_reduzido)
  modelo_completo[cont, ] <- estimando2(llike_model_ar_dumvar_ant, 
                                        pars_init_completo)
  
  print(cont)
}
fim <- Sys.time()

# Grafico e medidas descritivas das estimativas - INICIO
modelo_completo_ord <- modelo_completo %>% arrange(-ite_com)
plot(modelo_completo_ord$ite_com, ylab = 'iterecoes', type = 'l')

modelo_completo_ret <- modelo_completo_ord[floor(.1*M):M, ]

## Medidas descritivas modelo completo - INICIO
media_est <- apply(modelo_completo_ret, 2, mean); media_est
var_est <- apply(modelo_completo_ret, 2, sd); var_est
cor_est <- cor(modelo_completo_ret[, 1:6]); cor_est
unlist(pars)[c('omega', "alpha", "beta", "ar", "deltaMedia", "deltaVar")]
## Medidas descritivas modelo completo - FIM

est_pad <- modelo_completo_ret %>% select(alpha_com:deltaVar2_com) %>% pad()

## QQplot e histograma - INICIO
q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

## Grafico e medidas descritivas das estimativas - FIM

## Medidas descritivas modelo reduzido - INICIO
modelo_reduzido_ord <- modelo_reduzido %>% arrange(-ite)
plot(modelo_reduzido_ord$ite, ylab = 'iterecoes', type = 'l')

modelo_reduzido_ret <- modelo_reduzido_ord[floor(.1*M):M, ]

media_est <- apply(modelo_reduzido_ret, 2, mean); media_est
var_est <- apply(modelo_reduzido_ret, 2, sd); var_est
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar")]
## Medidas descritivas modelo reduzido - FIM

## QQplot e histograma - INICIO
est_pad <- modelo_reduzido_ret %>% select(alpha:deltaVar1) %>% pad()

q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

# Grafico e medidas descritivas das estimativas - FIM

# Teste da razao - INICIOs
llike_modelos_retirado <- llike_retirado(modelo_completo, modelo_reduzido, 
                                         taxa_rem = .2)

freq_rel <- apply(llike_modelos_retirado[c('rejeita_h0_1', 
                                           'rejeita_h0_5', 
                                           'rejeita_h0_10')], 2, mean); freq_rel

## Histogramas e QQplots - INICIO
QQplot_LR(llike_modelos_retirado, razao_llike_est, M, n)

histo_LR(llike_modelos_retirado, razao_llike_est, M, n)
## Histogramas e QQplots - FIM

# Teste da razao - FIM

# Caso 1,2 ---------------------------------------------------------------

## Definindo estruturas - INICIO
modelo_completo <- data.frame(omega_com = rep(NA, M), alpha_com = rep(NA, M), 
                              beta_com = rep(NA, M), 
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M),
                              deltaVar2_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M))

modelo_reduzido <- data.frame(omega = rep(NA, M), alpha = rep(NA, M), 
                              beta = rep(NA, M), 
                              ar = rep(NA, M), 
                              deltaMedia = rep(NA, M), 
                              deltaVar1 = rep(NA, M),
                              ite = rep(NA, M), llike_red = rep(NA, M))
## Definindo estruturas - FIM

n <- 4000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))

dummy2 <- as.matrix(dummy_step(n, c(1, 50)))
dummy2_red <- as.matrix(dummy_step(n, 1))

inicio <- Sys.time()
for (cont in 1:M){
  
  yt <- modelo(pars, dummy1, dummy2_red, n)$yt
  
  modelo_reduzido[cont, ] <- estimando2(llike_model_ar_dumvar_red_ant, 
                                        pars_init_reduzido)
  modelo_completo[cont, ] <- estimando2(llike_model_ar_dumvar_ant, 
                                        pars_init_completo)
  
  print(cont)
}
fim <- Sys.time()

# Grafico e medidas descritivas das estimativas - INICIO
modelo_completo_ord <- modelo_completo %>% arrange(-ite_com)
plot(modelo_completo_ord$ite_com, ylab = 'iterecoes', type = 'l')

modelo_completo_ret <- modelo_completo_ord[floor(.1*M):M, ]

## Medidas descritivas modelo completo - INICIO
media_est <- apply(modelo_completo_ret, 2, mean); media_est
var_est <- apply(modelo_completo_ret, 2, sd); var_est
cor_est <- cor(modelo_completo_ret[, 1:6]); cor_est
unlist(pars)[c('omega', "alpha", "beta", "ar", "deltaMedia", "deltaVar")]
## Medidas descritivas modelo completo - FIM

est_pad <- modelo_completo_ret %>% select(alpha_com:deltaVar2_com) %>% pad()

## QQplot e histograma - INICIO
q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

## Grafico e medidas descritivas das estimativas - FIM

## Medidas descritivas modelo reduzido - INICIO
modelo_reduzido_ord <- modelo_reduzido %>% arrange(-ite)
plot(modelo_reduzido_ord$ite, ylab = 'iterecoes', type = 'l')

modelo_reduzido_ret <- modelo_reduzido_ord[floor(.1*M):M, ]

media_est <- apply(modelo_reduzido_ret, 2, mean); media_est
var_est <- apply(modelo_reduzido_ret, 2, sd); var_est
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar")]
## Medidas descritivas modelo reduzido - FIM

## QQplot e histograma - INICIO
est_pad <- modelo_reduzido_ret %>% select(alpha:deltaVar1) %>% pad()

q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

# Grafico e medidas descritivas das estimativas - FIM

# Teste da razao - INICIOs
llike_modelos_retirado <- llike_retirado(modelo_completo, modelo_reduzido, taxa_rem = .1)

freq_rel <- apply(llike_modelos_retirado[c('rejeita_h0_1', 
                                           'rejeita_h0_5', 
                                           'rejeita_h0_10')], 2, mean); freq_rel

## Histogramas e QQplots - INICIO
QQplot_LR(llike_modelos_retirado, razao_llike_est, M, n)

histo_LR(llike_modelos_retirado, razao_llike_est, M, n)
## Histogramas e QQplots - FIM

# Teste da razao - FIM

# Caso 1,3 ---------------------------------------------------------------

## Definindo estruturas - INICIO
modelo_completo <- data.frame(omega_com = rep(NA, M), alpha_com = rep(NA, M), 
                              beta_com = rep(NA, M), 
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M),
                              deltaVar2_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M))

modelo_reduzido <- data.frame(omega = rep(NA, M), alpha = rep(NA, M), 
                              beta = rep(NA, M), 
                              ar = rep(NA, M), 
                              deltaMedia = rep(NA, M), 
                              deltaVar1 = rep(NA, M),
                              ite = rep(NA, M), llike_red = rep(NA, M))
## Definindo estruturas - FIM

n <- 5000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))

dummy2 <- as.matrix(dummy_step(n, c(1, 50)))
dummy2_red <- as.matrix(dummy_step(n, 1))

inicio <- Sys.time()
for (cont in 1:M){
  
  yt <- modelo(pars, dummy1, dummy2_red, n)$yt
  
  modelo_reduzido[cont, ] <- estimando2(llike_model_ar_dumvar_red_ant, 
                                        pars_init_reduzido)
  modelo_completo[cont, ] <- estimando2(llike_model_ar_dumvar_ant, 
                                        pars_init_completo)
  
  print(cont)
}
fim <- Sys.time()

# Grafico e medidas descritivas das estimativas - INICIO
modelo_completo_ord <- modelo_completo %>% arrange(-ite_com)
plot(modelo_completo_ord$ite_com, ylab = 'iterecoes', type = 'l')

modelo_completo_ret <- modelo_completo_ord[floor(.1*M):M, ]

## Medidas descritivas modelo completo - INICIO
media_est <- apply(modelo_completo_ret, 2, mean); media_est
var_est <- apply(modelo_completo_ret, 2, sd); var_est
cor_est <- cor(modelo_completo_ret[, 1:6]); cor_est
unlist(pars)[c('omega', "alpha", "beta", "ar", "deltaMedia", "deltaVar")]
## Medidas descritivas modelo completo - FIM

est_pad <- modelo_completo_ret %>% select(alpha_com:deltaVar2_com) %>% pad()

## QQplot e histograma - INICIO
q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

## Grafico e medidas descritivas das estimativas - FIM

## Medidas descritivas modelo reduzido - INICIO
modelo_reduzido_ord <- modelo_reduzido %>% arrange(-ite)
plot(modelo_reduzido_ord$ite, ylab = 'iterecoes', type = 'l')

modelo_reduzido_ret <- modelo_reduzido_ord[floor(.1*M):M, ]

media_est <- apply(modelo_reduzido_ret, 2, mean); media_est
var_est <- apply(modelo_reduzido_ret, 2, sd); var_est
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar")]
## Medidas descritivas modelo reduzido - FIM

## QQplot e histograma - INICIO
est_pad <- modelo_reduzido_ret %>% select(alpha:deltaVar1) %>% pad()

q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

# Grafico e medidas descritivas das estimativas - FIM

# Teste da razao - INICIOs
llike_modelos_retirado <- llike_retirado(modelo_completo, modelo_reduzido, taxa_rem = .1)

freq_rel <- apply(llike_modelos_retirado[c('rejeita_h0_1', 
                                           'rejeita_h0_5', 
                                           'rejeita_h0_10')], 2, mean); freq_rel

## Histogramas e QQplots - INICIO
QQplot_LR(llike_modelos_retirado, razao_llike_est, M, n)

histo_LR(llike_modelos_retirado, razao_llike_est, M, n)
## Histogramas e QQplots - FIM

# Teste da razao - FIM

# Caso 2,1 ---------------------------------------------------------------

## Definindo estruturas - INICIO
modelo_completo <- data.frame(omega_com = rep(NA, M), alpha_com = rep(NA, M), 
                              beta_com = rep(NA, M), 
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M),
                              deltaVar2_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M))

modelo_reduzido <- data.frame(omega = rep(NA, M), alpha = rep(NA, M), 
                              beta = rep(NA, M), 
                              ar = rep(NA, M), 
                              deltaMedia = rep(NA, M), 
                              deltaVar1 = rep(NA, M),
                              ite = rep(NA, M), llike_red = rep(NA, M))
## Definindo estruturas - FIM

# Monte Carlo - INICIO
n <- 3000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))

dummy2 <- as.matrix(dummy_step(n, c(1, 1500)))
dummy2_red <- as.matrix(dummy_step(n, 1))

inicio <- Sys.time()
for (cont in 1:M){
  
  yt <- modelo(pars, dummy1, dummy2_red, n)$yt
  
  modelo_reduzido[cont, ] <- estimando2(llike_model_ar_dumvar_red_ant, 
                                        pars_init_reduzido)
  modelo_completo[cont, ] <- estimando2(llike_model_ar_dumvar_ant, 
                                        pars_init_completo)
  
  print(cont)
}
fim <- Sys.time()

# Grafico e medidas descritivas das estimativas - INICIO
modelo_completo_ord <- modelo_completo %>% arrange(-ite_com)
plot(modelo_completo_ord$ite_com, ylab = 'iterecoes', type = 'l')

modelo_completo_ret <- modelo_completo_ord[floor(.1*M):M, ]

## Medidas descritivas modelo completo - INICIO
media_est <- apply(modelo_completo_ret, 2, mean); media_est
var_est <- apply(modelo_completo_ret, 2, sd); var_est
cor_est <- cor(modelo_completo_ret[, 1:6]); cor_est
unlist(pars)[c('omega', "alpha", "beta", "ar", "deltaMedia", "deltaVar")]
## Medidas descritivas modelo completo - FIM

est_pad <- modelo_completo_ret %>% select(alpha_com:deltaVar2_com) %>% pad()

## QQplot e histograma - INICIO
q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

## Grafico e medidas descritivas das estimativas - FIM

## Medidas descritivas modelo reduzido - INICIO
modelo_reduzido_ord <- modelo_reduzido %>% arrange(-ite)
plot(modelo_reduzido_ord$ite, ylab = 'iterecoes', type = 'l')

modelo_reduzido_ret <- modelo_reduzido_ord[floor(.1*M):M, ]

media_est <- apply(modelo_reduzido_ret, 2, mean); media_est
var_est <- apply(modelo_reduzido_ret, 2, sd); var_est
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar")]
## Medidas descritivas modelo reduzido - FIM

## QQplot e histograma - INICIO
est_pad <- modelo_reduzido_ret %>% select(alpha:deltaVar1) %>% pad()

q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

# Grafico e medidas descritivas das estimativas - FIM

# Teste da razao - INICIOs
llike_modelos_retirado <- llike_retirado(modelo_completo, modelo_reduzido, taxa_rem = .1)

freq_rel <- apply(llike_modelos_retirado[c('rejeita_h0_1', 
                                           'rejeita_h0_5', 
                                           'rejeita_h0_10')], 2, mean); freq_rel

## Histogramas e QQplots - INICIO
QQplot_LR(llike_modelos_retirado, razao_llike_est, M, n)

histo_LR(llike_modelos_retirado, razao_llike_est, M, n)
## Histogramas e QQplots - FIM

# Teste da razao - FIM

# Caso 2,2 ---------------------------------------------------------------


## Definindo estruturas - INICIO
modelo_completo <- data.frame(omega_com = rep(NA, M), alpha_com = rep(NA, M), 
                              beta_com = rep(NA, M), 
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M),
                              deltaVar2_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M))

modelo_reduzido <- data.frame(omega = rep(NA, M), alpha = rep(NA, M), 
                              beta = rep(NA, M), 
                              ar = rep(NA, M), 
                              deltaMedia = rep(NA, M), 
                              deltaVar1 = rep(NA, M),
                              ite = rep(NA, M), llike_red = rep(NA, M))
## Definindo estruturas - FIM

# Monte Carlo - INICIO
n <- 4000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))

dummy2 <- as.matrix(dummy_step(n, c(1, 1500)))
dummy2_red <- as.matrix(dummy_step(n, 1))

inicio <- Sys.time()
for (cont in 1:M){
  
  yt <- modelo(pars, dummy1, dummy2_red, n)$yt
  
  modelo_reduzido[cont, ] <- estimando2(llike_model_ar_dumvar_red_ant, 
                                        pars_init_reduzido)
  modelo_completo[cont, ] <- estimando2(llike_model_ar_dumvar_ant, 
                                        pars_init_completo)
  
  print(cont)
}
fim <- Sys.time()

# Grafico e medidas descritivas das estimativas - INICIO
modelo_completo_ord <- modelo_completo %>% arrange(-ite_com)
plot(modelo_completo_ord$ite_com, ylab = 'iterecoes', type = 'l')

modelo_completo_ret <- modelo_completo_ord[floor(.1*M):M, ]

## Medidas descritivas modelo completo - INICIO
media_est <- apply(modelo_completo_ret, 2, mean); media_est
var_est <- apply(modelo_completo_ret, 2, sd); var_est
cor_est <- cor(modelo_completo_ret[, 1:6]); cor_est
unlist(pars)[c('omega', "alpha", "beta", "ar", "deltaMedia", "deltaVar")]
## Medidas descritivas modelo completo - FIM

est_pad <- modelo_completo_ret %>% select(alpha_com:deltaVar2_com) %>% pad()

## QQplot e histograma - INICIO
q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

## Grafico e medidas descritivas das estimativas - FIM

## Medidas descritivas modelo reduzido - INICIO
modelo_reduzido_ord <- modelo_reduzido %>% arrange(-ite)
plot(modelo_reduzido_ord$ite, ylab = 'iterecoes', type = 'l')

modelo_reduzido_ret <- modelo_reduzido_ord[floor(.1*M):M, ]

media_est <- apply(modelo_reduzido_ret, 2, mean); media_est
var_est <- apply(modelo_reduzido_ret, 2, sd); var_est
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar")]
## Medidas descritivas modelo reduzido - FIM

## QQplot e histograma - INICIO
est_pad <- modelo_reduzido_ret %>% select(alpha:deltaVar1) %>% pad()

q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

# Grafico e medidas descritivas das estimativas - FIM

# Teste da razao - INICIOs
llike_modelos_retirado <- llike_retirado(modelo_completo, modelo_reduzido, taxa_rem = .1)

freq_rel <- apply(llike_modelos_retirado[c('rejeita_h0_1', 
                                           'rejeita_h0_5', 
                                           'rejeita_h0_10')], 2, mean); freq_rel

## Histogramas e QQplots - INICIO
QQplot_LR(llike_modelos_retirado, razao_llike_est, M, n)

histo_LR(llike_modelos_retirado, razao_llike_est, M, n)
## Histogramas e QQplots - FIM

# Teste da razao - FIM

# Caso 2,3 ---------------------------------------------------------------

## Definindo estruturas - INICIO
modelo_completo <- data.frame(omega_com = rep(NA, M), alpha_com = rep(NA, M), 
                              beta_com = rep(NA, M), 
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M),
                              deltaVar2_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M))

modelo_reduzido <- data.frame(omega = rep(NA, M), alpha = rep(NA, M), 
                              beta = rep(NA, M), 
                              ar = rep(NA, M), 
                              deltaMedia = rep(NA, M), 
                              deltaVar1 = rep(NA, M),
                              ite = rep(NA, M), llike_red = rep(NA, M))
## Definindo estruturas - FIM

# Monte Carlo - INICIO
n <- 5000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))

dummy2 <- as.matrix(dummy_step(n, c(1, 1500)))
dummy2_red <- as.matrix(dummy_step(n, 1))

inicio <- Sys.time()
for (cont in 1:M){
  
  yt <- modelo(pars, dummy1, dummy2_red, n)$yt
  
  modelo_reduzido[cont, ] <- estimando2(llike_model_ar_dumvar_red_ant, 
                                        pars_init_reduzido)
  modelo_completo[cont, ] <- estimando2(llike_model_ar_dumvar_ant, 
                                        pars_init_completo)
  
  print(cont)
}
fim <- Sys.time()

# Grafico e medidas descritivas das estimativas - INICIO
modelo_completo_ord <- modelo_completo %>% arrange(-ite_com)
plot(modelo_completo_ord$ite_com, ylab = 'iterecoes', type = 'l')

modelo_completo_ret <- modelo_completo_ord[floor(.1*M):M, ]

## Medidas descritivas modelo completo - INICIO
media_est <- apply(modelo_completo_ret, 2, mean); media_est
var_est <- apply(modelo_completo_ret, 2, sd); var_est
cor_est <- cor(modelo_completo_ret[, 1:6]); cor_est
unlist(pars)[c('omega', "alpha", "beta", "ar", "deltaMedia", "deltaVar")]
## Medidas descritivas modelo completo - FIM

est_pad <- modelo_completo_ret %>% select(alpha_com:deltaVar2_com) %>% pad()

## QQplot e histograma - INICIO
q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

## Grafico e medidas descritivas das estimativas - FIM

## Medidas descritivas modelo reduzido - INICIO
modelo_reduzido_ord <- modelo_reduzido %>% arrange(-ite)
plot(modelo_reduzido_ord$ite, ylab = 'iterecoes', type = 'l')

modelo_reduzido_ret <- modelo_reduzido_ord[floor(.1*M):M, ]

media_est <- apply(modelo_reduzido_ret, 2, mean); media_est
var_est <- apply(modelo_reduzido_ret, 2, sd); var_est
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar")]
## Medidas descritivas modelo reduzido - FIM

## QQplot e histograma - INICIO
est_pad <- modelo_reduzido_ret %>% select(alpha:deltaVar1) %>% pad()

q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

# Grafico e medidas descritivas das estimativas - FIM

# Teste da razao - INICIOs
llike_modelos_retirado <- llike_retirado(modelo_completo, modelo_reduzido, taxa_rem = .1)

freq_rel <- apply(llike_modelos_retirado[c('rejeita_h0_1', 
                                           'rejeita_h0_5', 
                                           'rejeita_h0_10')], 2, mean); freq_rel

## Histogramas e QQplots - INICIO
QQplot_LR(llike_modelos_retirado, razao_llike_est, M, n)

histo_LR(llike_modelos_retirado, razao_llike_est, M, n)
## Histogramas e QQplots - FIM

# Teste da razao - FIM

# Caso 3,1 ----------------------------------------------------------------

## Definindo estruturas - INICIO
modelo_completo <- data.frame(omega_com = rep(NA, M), alpha_com = rep(NA, M), 
                              beta_com = rep(NA, M), 
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M),
                              deltaVar2_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M))

modelo_reduzido <- data.frame(omega = rep(NA, M), alpha = rep(NA, M), 
                              beta = rep(NA, M), 
                              ar = rep(NA, M), 
                              deltaMedia = rep(NA, M), 
                              deltaVar1 = rep(NA, M),
                              ite = rep(NA, M), llike_red = rep(NA, M))
## Definindo estruturas - FIM

inicio <- Sys.time()
n <- 3000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))

dummy2 <- as.matrix(dummy_step(n, c(1, 2950)))
dummy2_red <- as.matrix(dummy_step(n, 1))

inicio <- Sys.time()
for (cont in 1:M){
  
  yt <- modelo(pars, dummy1, dummy2_red, n)$yt
  
  modelo_reduzido[cont, ] <- estimando2(llike_model_ar_dumvar_red_ant, 
                                        pars_init_reduzido)
  modelo_completo[cont, ] <- estimando2(llike_model_ar_dumvar_ant, 
                                        pars_init_completo)
  
  print(cont)
}
fim <- Sys.time()

# Grafico e medidas descritivas das estimativas - INICIO
modelo_completo_ord <- modelo_completo %>% arrange(-ite_com)
plot(modelo_completo_ord$ite_com, ylab = 'iterecoes', type = 'l')

modelo_completo_ret <- modelo_completo_ord[floor(.1*M):M, ]

## Medidas descritivas modelo completo - INICIO
media_est <- apply(modelo_completo_ret, 2, mean); media_est
var_est <- apply(modelo_completo_ret, 2, sd); var_est
cor_est <- cor(modelo_completo_ret[, 1:6]); cor_est
unlist(pars)[c('omega', "alpha", "beta", "ar", "deltaMedia", "deltaVar")]
## Medidas descritivas modelo completo - FIM

est_pad <- modelo_completo_ret %>% select(alpha_com:deltaVar2_com) %>% pad()

## QQplot e histograma - INICIO
q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

## Grafico e medidas descritivas das estimativas - FIM

## Medidas descritivas modelo reduzido - INICIO
modelo_reduzido_ord <- modelo_reduzido %>% arrange(-ite)
plot(modelo_reduzido_ord$ite, ylab = 'iterecoes', type = 'l')

modelo_reduzido_ret <- modelo_reduzido_ord[floor(.1*M):M, ]

media_est <- apply(modelo_reduzido_ret, 2, mean); media_est
var_est <- apply(modelo_reduzido_ret, 2, sd); var_est
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar")]
## Medidas descritivas modelo reduzido - FIM

## QQplot e histograma - INICIO
est_pad <- modelo_reduzido_ret %>% select(alpha:deltaVar1) %>% pad()

q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

# Grafico e medidas descritivas das estimativas - FIM

# Teste da razao - INICIOs
llike_modelos_retirado <- llike_retirado(modelo_completo, modelo_reduzido, taxa_rem = .1)

freq_rel <- apply(llike_modelos_retirado[c('rejeita_h0_1', 
                                           'rejeita_h0_5', 
                                           'rejeita_h0_10')], 2, mean); freq_rel

## Histogramas e QQplots - INICIO
QQplot_LR(llike_modelos_retirado, razao_llike_est, M, n)

histo_LR(llike_modelos_retirado, razao_llike_est, M, n)
## Histogramas e QQplots - FIM

# Teste da razao - FIM

# Caso 3,2 ---------------------------------------------------------------

## Definindo estruturas - INICIO
modelo_completo <- data.frame(omega_com = rep(NA, M), alpha_com = rep(NA, M), 
                              beta_com = rep(NA, M), 
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M),
                              deltaVar2_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M))

modelo_reduzido <- data.frame(omega = rep(NA, M), alpha = rep(NA, M), 
                              beta = rep(NA, M), 
                              ar = rep(NA, M), 
                              deltaMedia = rep(NA, M), 
                              deltaVar1 = rep(NA, M),
                              ite = rep(NA, M), llike_red = rep(NA, M))
## Definindo estruturas - FIM

inicio <- Sys.time()
n <- 4000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))

dummy2 <- as.matrix(dummy_step(n, c(1, 2950)))
dummy2_red <- as.matrix(dummy_step(n, 1))

inicio <- Sys.time()
for (cont in 1:M){
  
  yt <- modelo(pars, dummy1, dummy2_red, n)$yt
  
  modelo_reduzido[cont, ] <- estimando2(llike_model_ar_dumvar_red_ant, 
                                        pars_init_reduzido)
  modelo_completo[cont, ] <- estimando2(llike_model_ar_dumvar_ant, 
                                        pars_init_completo)
  
  print(cont)
}
fim <- Sys.time()

# Grafico e medidas descritivas das estimativas - INICIO
modelo_completo_ord <- modelo_completo %>% arrange(-ite_com)
plot(modelo_completo_ord$ite_com, ylab = 'iterecoes', type = 'l')

modelo_completo_ret <- modelo_completo_ord[floor(.1*M):M, ]

## Medidas descritivas modelo completo - INICIO
media_est <- apply(modelo_completo_ret, 2, mean); media_est
var_est <- apply(modelo_completo_ret, 2, sd); var_est
cor_est <- cor(modelo_completo_ret[, 1:6]); cor_est
unlist(pars)[c('omega', "alpha", "beta", "ar", "deltaMedia", "deltaVar")]
## Medidas descritivas modelo completo - FIM

est_pad <- modelo_completo_ret %>% select(alpha_com:deltaVar2_com) %>% pad()

## QQplot e histograma - INICIO
q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

## Grafico e medidas descritivas das estimativas - FIM

## Medidas descritivas modelo reduzido - INICIO
modelo_reduzido_ord <- modelo_reduzido %>% arrange(-ite)
plot(modelo_reduzido_ord$ite, ylab = 'iterecoes', type = 'l')

modelo_reduzido_ret <- modelo_reduzido_ord[floor(.1*M):M, ]

media_est <- apply(modelo_reduzido_ret, 2, mean); media_est
var_est <- apply(modelo_reduzido_ret, 2, sd); var_est
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar")]
## Medidas descritivas modelo reduzido - FIM

## QQplot e histograma - INICIO
est_pad <- modelo_reduzido_ret %>% select(alpha:deltaVar1) %>% pad()

q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

# Grafico e medidas descritivas das estimativas - FIM

# Teste da razao - INICIOs
llike_modelos_retirado <- llike_retirado(modelo_completo, modelo_reduzido, taxa_rem = .1)

freq_rel <- apply(llike_modelos_retirado[c('rejeita_h0_1', 
                                           'rejeita_h0_5', 
                                           'rejeita_h0_10')], 2, mean); freq_rel

## Histogramas e QQplots - INICIO
QQplot_LR(llike_modelos_retirado, razao_llike_est, M, n)

histo_LR(llike_modelos_retirado, razao_llike_est, M, n)
## Histogramas e QQplots - FIM

# Teste da razao - FIM

# Caso 3,3 --------------------------------------------------------------

## Definindo estruturas - INICIO
modelo_completo <- data.frame(omega_com = rep(NA, M), alpha_com = rep(NA, M), 
                              beta_com = rep(NA, M), 
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M),
                              deltaVar2_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M))

modelo_reduzido <- data.frame(omega = rep(NA, M), alpha = rep(NA, M), 
                              beta = rep(NA, M), 
                              ar = rep(NA, M), 
                              deltaMedia = rep(NA, M), 
                              deltaVar1 = rep(NA, M),
                              ite = rep(NA, M), llike_red = rep(NA, M))
## Definindo estruturas - FIM

inicio <- Sys.time()
n <- 5000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))

dummy2 <- as.matrix(dummy_step(n, c(1, 2950)))
dummy2_red <- as.matrix(dummy_step(n, 1))

inicio <- Sys.time()
for (cont in 1:M){
  
  yt <- modelo(pars, dummy1, dummy2_red, n)$yt
  
  modelo_reduzido[cont, ] <- estimando2(llike_model_ar_dumvar_red_ant, 
                                        pars_init_reduzido)
  modelo_completo[cont, ] <- estimando2(llike_model_ar_dumvar_ant, 
                                        pars_init_completo)
  
  print(cont)
}
fim <- Sys.time()

# Grafico e medidas descritivas das estimativas - INICIO
modelo_completo_ord <- modelo_completo %>% arrange(-ite_com)
plot(modelo_completo_ord$ite_com, ylab = 'iterecoes', type = 'l')

modelo_completo_ret <- modelo_completo_ord[floor(.1*M):M, ]

## Medidas descritivas modelo completo - INICIO
media_est <- apply(modelo_completo_ret, 2, mean); media_est
var_est <- apply(modelo_completo_ret, 2, sd); var_est
cor_est <- cor(modelo_completo_ret[, 1:6]); cor_est
unlist(pars)[c('omega', "alpha", "beta", "ar", "deltaMedia", "deltaVar")]
## Medidas descritivas modelo completo - FIM

est_pad <- modelo_completo_ret %>% select(alpha_com:deltaVar2_com) %>% pad()

## QQplot e histograma - INICIO
q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

## Grafico e medidas descritivas das estimativas - FIM

## Medidas descritivas modelo reduzido - INICIO
modelo_reduzido_ord <- modelo_reduzido %>% arrange(-ite)
plot(modelo_reduzido_ord$ite, ylab = 'iterecoes', type = 'l')

modelo_reduzido_ret <- modelo_reduzido_ord[floor(.1*M):M, ]

media_est <- apply(modelo_reduzido_ret, 2, mean); media_est
var_est <- apply(modelo_reduzido_ret, 2, sd); var_est
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar")]
## Medidas descritivas modelo reduzido - FIM

## QQplot e histograma - INICIO
est_pad <- modelo_reduzido_ret %>% select(alpha:deltaVar1) %>% pad()

q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

# Grafico e medidas descritivas das estimativas - FIM

# Teste da razao - INICIOs
llike_modelos_retirado <- llike_retirado(modelo_completo, modelo_reduzido, taxa_rem = .1)

freq_rel <- apply(llike_modelos_retirado[c('rejeita_h0_1', 
                                           'rejeita_h0_5', 
                                           'rejeita_h0_10')], 2, mean); freq_rel

## Histogramas e QQplots - INICIO
QQplot_LR(llike_modelos_retirado, razao_llike_est, M, n)

histo_LR(llike_modelos_retirado, razao_llike_est, M, n)
## Histogramas e QQplots - FIM

# Teste da razao - FIM