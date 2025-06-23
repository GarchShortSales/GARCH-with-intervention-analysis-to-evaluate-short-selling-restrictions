source("src/modelo_sim.R")
source("src/util.R")
source("src/modelo_est.R")
source("src/graficos.R")

# Definindo algumas estruturas das series - INICIO
kc <- 2

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = -6)

## Chute inicial - INICIO
pars_init_completo <- list(psi2 = log(.13), psi3 = log(.85),
                           ar = .2, 
                           deltaMedia = 10, 
                           deltaVar = -c(6, 6))

# Chute inicial - INICIO
n <- 6000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))

dummy2 <- as.matrix(dummy_on_off(n, c(1, 3000), c(2999, n)))
dummy2_red <- as.matrix(dummy_step(n, 1))

set.seed(87891212)
dados <- modelo(pars, dummy1, dummy2_red, n)
yt <- dados$yt
line(dados, "time", "yt", "Serie simulado")
# Chute inicial - FIM
acf(yt)
pacf(yt)

media <- 10 + .2*(yt[1:(n-1)] - 10)
plot((yt[2:n] - media)^2, type = 'l')
acf((yt[2:n] - media)^2)
pacf((yt[2:n] - media)^2)

# Estimando - INICIO
inicio <- Sys.time()
opt <- estimando(llike_model_ar_dumvar, pars_init_completo)
fim <- Sys.time()
fim - inicio

opt
unlist(pars)
# Estimando - FIM