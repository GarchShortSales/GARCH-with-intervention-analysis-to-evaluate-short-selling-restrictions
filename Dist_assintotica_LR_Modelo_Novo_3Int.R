source("src/modelo_sim.R")
source("src/util.R")
source("src/modelo_est.R")
source("src/graficos.R")

# M de Monte Carlo

M <- 50

# Cenario 1 ---------------------------------------------------------------

# Definindo algumas estruturas das series - INICIO
kc <- 3
kc_red <- 1

pars <- list(ar = -.1,
             omega = 1, alpha = .09, beta = .89,
             deltaMedia = .1,
             deltaVar = -3.5)

## Chute inicial - INICIO
pars_init_completo <- list(psi2 = log(.09), psi3 = log(.89),
                           ar = -.1,
                           deltaMedia = .1,
                           deltaVar = c(-3.5, -3.5, -3.5))

pars_init_reduzido <- list(psi2 = log(.09), psi3 = log(.89),
                           ar = -.1,
                           deltaMedia = .1,
                           deltaVar = -3.5)
## Chute inicial - FIM

# Definindo algumas estruturas das series - INICIO

## Definindo estruturas - INICIO
modelo_completo <- data.frame(alpha_com = rep(NA, M), beta_com = rep(NA, M),
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M),
                              deltaVar2_com = rep(NA, M),
                              deltaVar3_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M),
                              AIC_com = rep(NA, M),
                              BIC_com = rep(NA, M))

modelo_reduzido <- data.frame(alpha = rep(NA, M), beta = rep(NA, M),
                              ar = rep(NA, M),
                              deltaMedia = rep(NA, M),
                              deltaVar1 = rep(NA, M),
                              ite = rep(NA, M), llike_red = rep(NA, M), 
                              AIC = rep(NA, M),
                              BIC = rep(NA, M))

## Definindo estruturas - FIM

set.seed(89540)

n <- 2000

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))

dummy2 <- as.matrix(dummy_on_off(n, c(1, 1301, 1701),
                                 c(1300, 1700, n)))
dummy2_red <- as.matrix(dummy_step(n, 1))

inicio <- Sys.time()
for (cont in 1:M){

  yt <- modelo(pars, dummy1, as.matrix(dummy_step(n, 1)), n)$yt

  modelo_reduzido[cont, ] <- estimando(llike_model_ar_dumvar_red, pars_init_reduzido)
  modelo_completo[cont, ] <- estimando(llike_model_ar_dumvar, pars_init_completo)

  print(cont)
}
fim <- Sys.time()

### Grafico e medidas descritivas das estimativas do modelo completo

# Grafico e medidas descritivas das estimativas - INICIO
modelo_completo_ord <- modelo_completo %>% arrange(-ite_com)
modelo_completo_ret <- modelo_completo_ord[floor(.01 * M):M, ]

media_est <- apply(modelo_completo_ret, 2, mean); media_est
var_est <- apply(modelo_completo_ret, 2, sd);  var_est
cor_est <- cor(modelo_completo_ret[, 1:6]);  cor_est

## Grafico de linha - INICIO
plot(modelo_completo_ord$ite_com, ylab = "iterecoes", type = "l")
## Grafico de linha - FIM

est_pad <- modelo_completo %>% select(alpha_com:deltaVar3_com) %>% pad()

## QQplot e histograma - INICIO
q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

# Grafico e medidas descritivas das estimativas - FIM


### Grafico e medidas descritivas das estimativas do modelo reduzido

# Grafico e medidas descritivas das estimativas - INICIO
modelo_reduzido_ord <- modelo_reduzido %>% arrange(-ite)
modelo_reduzido_ret <- modelo_reduzido_ord[floor(.01 * M):M, ]

media_est <- apply(modelo_reduzido_ret, 2, mean);media_est
var_est <- apply(modelo_reduzido_ret, 2, sd); var_est

## Grafico de linha - INICIO
plot(modelo_reduzido_ord$ite, ylab = 'iterecoes', type = 'l')
## Grafico de linha - FIM

est_pad <- modelo_reduzido_ret %>% select(alpha:deltaVar1) %>% pad()

## QQplot e histograma - INICIO
q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

# Grafico e medidas descritivas das estimativas - FIM


### Teste LR acessando dist assintotica

# Teste da razao - INICIO
llike_modelos_retirado <- llike_retirado(modelo_completo, 
                                         modelo_reduzido, 
                                         df = kc - kc_red)

freq_rel <- apply(llike_modelos_retirado[c('rejeita_h0_1', 
                                           'rejeita_h0_5', 
                                           'rejeita_h0_10')], 2, mean); freq_rel

apply(llike_modelos_retirado['AIC_comparar'], 2, mean)
apply(llike_modelos_retirado['BIC_comparar'], 2, mean)

## Histogramas e QQplots - INICIO
QQplot_LR(llike_modelos_retirado, razao_llike_est, M, n, df = kc - kc_red)

histo_LR(llike_modelos_retirado, razao_llike_est, M, n, df = kc - kc_red)
## Histogramas e QQplots - FIM

# Teste da razao - FIM

# Cenario 2 ---------------------------------------------------------------

# Definindo algumas estruturas das series - INICIO
kc <- 3
kc_red <- 2

n <- 2000

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))

dummy2 <- as.matrix(dummy_on_off(n, c(1, 1301, 1701),
                                 c(1300, 1700, n)))
dummy2_red <- as.matrix(dummy_on_off(n, c(1, 1301), 
                                     c(1300, n)))

pars <- list(ar = -.1, 
             omega = 1, alpha = .09, beta = .89, 
             deltaMedia = .1, 
             deltaVar = c(-3.5, -2.6))

## Chute inicial - INICIO
pars_init_completo <- list(psi2 = log(.09), psi3 = log(.89),
                           ar = -.1, 
                           deltaMedia = .1, 
                           deltaVar = c(-3.5, -2.6, -2.6))

pars_init_reduzido <- list(psi2 = log(.09), psi3 = log(.89),
                           ar = -.1, 
                           deltaMedia = .1, 
                           deltaVar = c(-3.5, -2.6))
## Chute inicial - FIM

# Definindo algumas estruturas das series - INICIO

## Definindo estruturas - INICIO
modelo_completo <- data.frame(alpha_com = rep(NA, M), beta_com = rep(NA, M),
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M), 
                              deltaVar2_com = rep(NA, M),
                              deltaVar3_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M))

modelo_reduzido <- data.frame(alpha = rep(NA, M), beta = rep(NA, M),
                              ar = rep(NA, M),
                              deltaMedia = rep(NA, M),
                              deltaVar1 = rep(NA, M),
                              deltaVar2 = rep(NA, M),
                              ite = rep(NA, M), llike_red = rep(NA, M))
## Definindo estruturas - FIM

set.seed(97560)

inicio <- Sys.time()
for (cont in 1:M){
  
  yt <- modelo(pars, dummy1, dummy2_red, n)$yt
  
  modelo_reduzido[cont, ] <- estimando(llike_model_ar_dumvar_red, pars_init_reduzido)
  modelo_completo[cont, ] <- estimando(llike_model_ar_dumvar, pars_init_completo)
  
  print(cont)
}
fim <- Sys.time()


### Grafico e medidas descritivas das estimativas do modelo completo

# Grafico e medidas descritivas das estimativas - INICIO
modelo_completo_ord <- modelo_completo %>% arrange(-ite_com)
modelo_completo_ret <- modelo_completo_ord[floor(.01*M):M, ]

media_est <- apply(modelo_completo_ret, 2, mean); media_est
var_est <- apply(modelo_completo_ret, 2, sd); var_est
cor_est <- cor(modelo_completo_ret[, 1:6]); cor_est

## Grafico de linha - INICIO
plot(modelo_completo_ord$ite_com, ylab = 'iterecoes', type = 'l')
## Grafico de linha - FIM

est_pad <- modelo_completo %>% select(alpha_com:deltaVar3_com) %>% pad()

## QQplot e histograma - INICIO
q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

# Grafico e medidas descritivas das estimativas - FIM

### Grafico e medidas descritivas das estimativas do modelo reduzido

# Grafico e medidas descritivas das estimativas - INICIO
modelo_reduzido_ord <- modelo_reduzido %>% arrange(-ite)
modelo_reduzido_ret <- modelo_reduzido_ord[floor(.01 * M):M, ]

media_est <- apply(modelo_reduzido_ret, 2, mean); media_est
var_est <- apply(modelo_reduzido_ret, 2, sd); var_est
print("Verdadeiros parametros, ")
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar")]

## Grafico de linha - INICIO
plot(modelo_reduzido_ord$ite, ylab = 'iterecoes', type = 'l')
## Grafico de linha - FIM

est_pad <- modelo_reduzido_ret %>% select(alpha:deltaVar1) %>% pad()

## QQplot e histograma - INICIO
q1 <- map(names(est_pad), ~QQplot_MC(est_pad, .x, M, n)); q1
## QQplot e histograma - FIM

# Grafico e medidas descritivas das estimativas - FIM


### Teste LR acessando dist assintotica

# Teste da razao - INICIO

## Guardando as log-likes e medidas para dos modelos reduzido e completo - INICIO
llike_modelos_retirado <- llike_retirado(modelo_completo, modelo_reduzido, df = 1)

freq_rel <- apply(llike_modelos_retirado[c("rejeita_h0_1",
                                           "rejeita_h0_5",
                                           "rejeita_h0_10")],
                  2, mean)

cat("Frequencia Relativa de rejeição de H0 com alpha 1%, 5% e 10%:\n",
    freq_rel)

## Histogramas e QQplots - INICIO
QQplot_LR(llike_modelos_retirado, razao_llike_est, M, n, df = kc - kc_red)

histo_LR(llike_modelos_retirado, razao_llike_est, M, n, df = kc - kc_red)
## Histogramas e QQplots - FIM

# Teste da razao - FIM