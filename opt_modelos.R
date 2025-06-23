source("src/modelo_sim.R")
source("src/util.R")
source("src/modelo_est.R")
source("src/graficos.R")

# Exemplo 1 ---------------------------------------------------------------

# Definindo algumas estruturas das series - INICIO
n <- 2000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1500), c(1499, n)))
kc <- 2

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = c(1, 5))
# Definindo algumas estruturas das series - FIM

## Chute inicial - INICIO
pars_init<- list(psi2 = log(.13), psi3 = log(.6),
                 ar = .2, 
                 deltaMedia = 10, 
                 deltaVar = c(1, 5))
## Chute inicial - FIM

# pars <- unlist(pars_init, use.names = F)
# Chute inicial - INICIO
dados <- modelo(pars, dummy1, dummy2, n)
yt <- dados$yt
line(dados, "time", "yt", "Serie simulado") + 
  geom_vline(xintercept = 1500, 
             linetype="dashed", color = "blue", size = .9)
# Chute inicial - FIM

# Estimando - INICIO
inicio <- Sys.time()
opt <- estimando(llike_model_ar_dumvar, pars_init); opt
Sys.time() - inicio
# Estimando - FIM

# Calculando Verossimilhança - INICIO
llike_normal_model_ar_dumvar(unlist(opt, use.names = F))
opt
# Calculando Verossimilhança - FIM

# Exemplo 2 ---------------------------------------------------------------

# Definindo algumas estruturas das series - INICIO
n <- 2000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 500, 700), c(499, 699, n)))
kc <- 3

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = c(4, 5, 3))
# Definindo algumas estruturas das series - FIM

# Chute inicial - INICIO
pars_init<- list(psi2 = log(.13), psi3 = log(.85),
                 ar = .2, 
                 deltaMedia = 10, 
                 deltaVar = c(4, 5, 3))
# Chute inicial - FIM

# Gerando dados - INICIO
dados <- modelo(pars, dummy1, dummy2, n)
yt <- dados$yt
line(dados, "time", "yt", "Serie simulado") +
  geom_vline(xintercept = c(500, 700), 
             linetype="dashed", color = "blue", size = .9)
# Gerando dados - FIM

# Estimando - INICIO
inicio <- Sys.time()
opt <- estimando(llike_model_ar_dumvar, pars_init); opt
fim <- Sys.time()
fim - inicio

unlist(pars)[names(opt)[1:7]]
# Estimando - FIM

# Calculando Verossimilhança - INICIO
llike_normal_model_ar_dumvar(unlist(opt, use.names = F))
opt
# Calculando Verossimilhança - FIM

# Exemplo 3 ---------------------------------------------------------------

# Definindo algumas estruturas das series - INICIO
n <- 2000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 500, 900), c(499, 899, n)))
kc <- 3

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = c(.5, .3, .1))
# Definindo algumas estruturas das series - FIM

## Chute inicial - INICIO
pars_init<- list(psi2 = log(.13), psi3 = log(.85),
                 ar = .2, 
                 deltaMedia = 10, 
                 deltaVar = c(.5, .3, .1))
## Chute inicial - FIM

# Chute inicial - INICIO
dados <- modelo(pars, dummy1, dummy2, n)
yt <- dados$yt
line(dados, "time", "yt", "Serie simulado") +
  geom_vline(xintercept = c(500, 900), 
             linetype="dashed", color = "blue", size = .9)
# Chute inicial - FIM

# Estimando - INICIO
inicio <- Sys.time()
opt <- estimando(llike_model_ar_dumvar, pars_init); opt
fim <- Sys.time()
fim - inicio

unlist(pars)[names(opt)[1:7]]
# Estimando - FIM

# Calculando Verossimilhança - INICIO
llike_normal_model_ar_dumvar(unlist(opt, use.names = F))
opt
# Calculando Verossimilhança - FIM
# Exemplo 4 ---------------------------------------------------------------

# Definindo algumas estruturas das series - INICIO
n <- 2000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_step(n, 1))
kc <- 2

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = -c(10, 10))
# Definindo algumas estruturas das series - FIM

## Chute inicial - INICIO
pars_init<- list(psi2 = log(.25), psi3 = log(.7),
                 ar = .2, 
                 deltaMedia = 10, 
                 deltaVar = -c(8, 7))
## Chute inicial - FIM

# Chute inicial - INICIO
set.seed(5050)
dados <- modelo(pars, dummy1, dummy2, n)
yt <- dados$yt
line(dados, "time", "yt", "Serie simulado") +
  geom_vline(xintercept = 50, 
             linetype="dashed", color = "blue", size = .9)
# Chute inicial - FIM

media <- 10 + .2*(yt[1:(n-1)] - 10)
acf((yt[2:n] - media)^2)
pacf((yt[2:n] - media)^2)

# Estimando - INICIO
inicio <- Sys.time()
opt <- estimando(llike_model_ar_dumvar, pars_init_completo)
fim <- Sys.time()
fim - inicio

opt
unlist(pars)[names(opt)[1:7]]
# Estimando - FIM

# Calculando Verossimilhança - INICIO
llike_normal_model_ar_dumvar(unlist(opt, use.names = F))
opt
# Calculando Verossimilhança - FIM