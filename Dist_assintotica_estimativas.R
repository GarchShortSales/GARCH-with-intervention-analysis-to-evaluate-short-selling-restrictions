source("src/modelo_sim.R")
source("src/util.R")
source("src/modelo_est.R")
source("src/graficos.R")

require(purrr)
require(dplyr)

# Exemplo 1 ---------------------------------------------------------------

# Definindo algumas estruturas das series - INICIO
n <- 1000
dummy1 <- as.matrix(dummy_on_off(n, c(1, 500), c(499, n), "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 500), c(499, n)))
kc <- 2

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = c(10, 5),  
             deltaVar = c(log(5), 5))

## Chute inicial - INICIO
pars_init<- list(psi2 = log(.13), psi3 = log(.85),
                           ar = .2, 
                           deltaMedia = 10, 
                           deltaVar = c(log(5), 5))
## Chute inicial - FIM

# Inciando MC - INICIO

## Definindo estruturas - INICIO
M <- 10
data <- matrix(nrow = n, ncol = M)
modelo_completo <- data.frame(alpha_com = rep(NA, M), beta_com = rep(NA, M), 
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M), deltaVar2_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M))
## Definindo estruturas - FIM

set.seed(5500559)
inicio <- Sys.time()
for (cont in 1:M){
  
  yt <- modelo(pars, dummy1, dummy2, n)$yt
  modelo_completo[cont, ] <- estimando(llike_model_ar_dumvar, pars_init)
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
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar1", "deltaVar2")]

sd_est <- est_data %>% apply(2, sd)
cv <- sd_est/media_est; cv

# Padronizando - INICIO
MC_pad <- est_data %>% pad()
apply(MC_pad, 2, mean)
apply(MC_pad, 2, sd)
# Padronizando - FIM

# QQplot e histograma - INICIO
q1 <- map(names(MC_pad), ~QQplot_MC(MC_pad, .x, M, n))
h1 <- map(names(MC_pad), ~histo_MC(MC_pad, .x, M, n))
# QQplot e histograma - FIM

map(MC_pad, ~shapiro.test(.x))
map(MC_pad, ~tseries::jarque.bera.test(.x))
# Estudando convergencia - FIM

# Exemplo 2 ---------------------------------------------------------------

# Definindo algumas estruturas das series - INICIO
n <- 2000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1900), c(1899, n)))
kc <- 2

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = c(log(5), 5))

## Chute inicial - INICIO
pars_init<- list(psi2 = log(.13), psi3 = log(.85),
                 ar = .2, 
                 deltaMedia = 10, 
                 deltaVar = c(log(5), 5))
## Chute inicial - FIM

# Inciando MC - INICIO

## Definindo estruturas - INICIO
M <- 1000
data <- matrix(nrow = n, ncol = M)
modelo_completo <- data.frame(alpha_com = rep(NA, M), beta_com = rep(NA, M), 
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M), deltaVar2_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M))
## Definindo estruturas - FIM

set.seed(1000)
inicio <- Sys.time()
for (cont in 1:M){
  
  yt <- modelo(pars, dummy1, dummy2, n)$yt
  modelo_completo[cont, ] <- estimando(llike_model_ar_dumvar, pars_init)
  print(cont)
  
}
fim <- Sys.time()
fim - inicio
# Inciando MC - FIM

# Tirando as primeiras estimativas - INICIO
est_data <- modelo_completo %>% arrange(-ite_com) 
est_data <- est_data[floor(.15*M):M, ]
# Tirando as primeiras estimativas- FIM

# Estudando convergencia - INICIO
est_data <- est_data %>% select(alpha_com:deltaVar2_com)
media_est <- est_data %>% apply(2, mean); media_est
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar1", "deltaVar2")]

sd_est <- est_data %>% apply(2, sd)
cv <- sd_est/media_est; cv

# Padronizando - INICIO
MC_pad <- est_data %>% pad()
apply(MC_pad, 2, mean)
apply(MC_pad, 2, sd)
# Padronizando - FIM

# QQplot e histograma - INICIO
q1 <- map(names(MC_pad), ~QQplot_MC(MC_pad, .x, M, n))
h1 <- map(names(MC_pad), ~histo_MC(MC_pad, .x, M, n))
# QQplot e histograma - FIM

map(MC_pad, ~shapiro.test(.x))
map(MC_pad, ~tseries::jarque.bera.test(.x))
# Estudando convergencia - FIM
# Exemplo 3 ---------------------------------------------------------------

# Definindo algumas estruturas das series - INICIO
n <- 2000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 100), c(99, n)))
kc <- 2

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = c(log(5), 5))

## Chute inicial - INICIO
pars_init<- list(psi2 = log(.13), psi3 = log(.85),
                 ar = .2, 
                 deltaMedia = 10, 
                 deltaVar = c(log(5), 5))
## Chute inicial - FIM

# Inciando MC - INICIO

## Definindo estruturas - INICIO
M <- 1000
data <- matrix(nrow = n, ncol = M)
modelo_completo <- data.frame(alpha_com = rep(NA, M), beta_com = rep(NA, M), 
                              ar_com = rep(NA, M),
                              deltaMedia_com = rep(NA, M),
                              deltaVar1_com = rep(NA, M), deltaVar2_com = rep(NA, M),
                              ite_com = rep(NA, M), llike_com = rep(NA, M))
## Definindo estruturas - FIM

set.seed(1000)
inicio <- Sys.time()
for (cont in 1:M){
  
  yt <- modelo(pars, dummy1, dummy2, n)$yt
  modelo_completo[cont, ] <- estimando(llike_model_ar_dumvar, pars_init)
  print(cont)
  
}
fim <- Sys.time()
fim - inicio
# Inciando MC - FIM

# Tirando as primeiras estimativas - INICIO
est_data <- modelo_completo %>% arrange(-ite_com) 
est_data <- est_data[floor(.15*M):M, ]
# Tirando as primeiras estimativas- FIM

# Estudando convergencia - INICIO
est_data <- est_data %>% select(alpha_com:deltaVar2_com)
media_est <- est_data %>% apply(2, mean); media_est
unlist(pars)[c("alpha", "beta", "ar", "deltaMedia", "deltaVar1", "deltaVar2")]

sd_est <- est_data %>% apply(2, sd)
cv <- sd_est/media_est; cv

# Padronizando - INICIO
MC_pad <- est_data %>% pad()
apply(MC_pad, 2, mean)
apply(MC_pad, 2, sd)
# Padronizando - FIM

# QQplot e histograma - INICIO
q1 <- map(names(MC_pad), ~QQplot_MC(MC_pad, .x, M, n))
h1 <- map(names(MC_pad), ~histo_MC(MC_pad, .x, M, n))
# QQplot e histograma - FIM

map(MC_pad, ~shapiro.test(.x))
map(MC_pad, ~tseries::jarque.bera.test(.x))
# Estudando convergencia - FIM