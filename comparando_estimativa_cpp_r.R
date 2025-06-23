source("src/modelo_sim.R")
source("src/util.R")
source("src/modelo_est.R")

# Cenario - complicado ----------------------------------------------------

n <- 1500
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1451),
                                 c(1450, n)))
kc <- 2

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = -c(1, 5))
# Definindo algumas estruturas das series - FIM

## Chute inicial - INICIO
pars_init<- list(psi2 = log(.13), psi3 = log(.85),
                 ar = .2, 
                 deltaMedia = 10, 
                 deltaVar = -c(1, 5))
## Chute inicial - FIM

# pars <- unlist(pars_init, use.names = F)
# Chute inicial - INICIO
set.seed(1000)

inicio <- Sys.time()
dados <- modelo(pars, dummy1, dummy2, n)
fim <- Sys.time()

yt <- dados$yt
Varyt <- var(yt[1:50])

inicio1 <- Sys.time()
opt <- estimando(llike_model_ar_dumvar, pars_init)
fim1 <- Sys.time()

inicio2 <- Sys.time()
opt2 <- estimando_cpp(pars_init, yt, dummy1, dummy2, kc, n, Varyt)
fim2 <- Sys.time()

fim - inicio
fim1 - inicio1
fim2 - inicio2

# Cenario 1 + complicado ----------------------------------------------------

n <- 5000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1500, 1901, 3001, 4951),
                                 c(1499, 1900, 3000, 4950, n)))
kc <- 5

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = -c(1, 5, 8, 3, 4))
# Definindo algumas estruturas das series - FIM

## Chute inicial - INICIO
pars_init<- list(psi2 = log(.13), psi3 = log(.85),
                 ar = .2, 
                 deltaMedia = 10, 
                 deltaVar = -c(1, 5, 8, 3, 4))
## Chute inicial - FIM

# pars <- unlist(pars_init, use.names = F)
# Chute inicial - INICIO
dados <- modelo(pars, dummy1, dummy2, n)
yt <- dados$yt

inicio2 <- Sys.time()
opt2 <- estimando_cpp(pars_init, yt, dummy1, dummy2, kc); opt2
fim2 <- Sys.time()

inicio1 <- Sys.time()
opt <- estimando(llike_model_ar_dumvar, pars_init); opt
fim1 <- Sys.time()

fim1 - inicio1
fim2 - inicio2


# Cenario 2 + complicado ----------------------------------------------------

n <- 1500
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1000, 1201, 1301, 1351),
                                 c(999, 1200, 1300, 1350, n)))
kc <- 5

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = -c(1, 5, 8, 3, 4))
# Definindo algumas estruturas das series - FIM

## Chute inicial - INICIO
pars_init<- list(psi2 = log(.13), psi3 = log(.85),
                 ar = .2, 
                 deltaMedia = 10, 
                 deltaVar = -c(1, 5, 8, 3, 4))
## Chute inicial - FIM

# pars <- unlist(pars_init, use.names = F)
# Chute inicial - INICIO
dados <- modelo(pars, dummy1, dummy2, n)
yt <- dados$yt

inicio2 <- Sys.time()
opt2 <- estimando_cpp(pars_init, yt, dummy1, dummy2, kc); opt2
fim2 <- Sys.time()

inicio1 <- Sys.time()
opt <- estimando(llike_model_ar_dumvar, pars_init); opt
fim1 <- Sys.time()

fim1 - inicio1
fim2 - inicio2



# Cenario com 2 int na media ----------------------------------------------

n <- 2500
dummy1 <- as.matrix(dummy_on_off(n, c(1, 1401),
                                 c(1400, n), 'Media'))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1401),
                                 c(1400, n)))
kvar <- 2 
kmed <- 2

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = c(4, 5), 
             deltaVar = -c(1, 5))
# Definindo algumas estruturas das series - FIM

## Chute inicial - INICIO
pars_init <- list(psi2 = log(.13), psi3 = log(.85),
                  ar = .2, 
                  deltaMedia = c(4, 5), 
                  deltaVar = -c(1, 5))
## Chute inicial - FIM

# pars <- unlist(pars_init, use.names = F)
# Chute inicial - INICIO
set.seed(8517)

dados <- modelo(pars, dummy1, dummy2, n)

yt <- dados$yt
Varyt <- var(yt[1:50])

inicio1 <- Sys.time()
opt <- estimando(llike_model_geral, pars_init)
fim1 <- Sys.time()

inicio2 <- Sys.time()
opt2 <- estimando_cpp_geral(pars_init, yt, dummy1, dummy2, 
                            kmed, kvar, n, Varyt)
fim2 <- Sys.time()

fim1 - inicio1
fim2 - inicio2

opt
opt2

# Cenario com 2 int na media ----------------------------------------------

n <- 1500
dummy1 <- as.matrix(dummy_on_off(n, c(1, 1301),
                                 c(1300, n), 'Media'))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1001),
                                 c(1000, n)))
kvar <- 2 
kmed <- 2

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = c(1, 5), 
             deltaVar = -c(1, 5))
# Definindo algumas estruturas das series - FIM

## Chute inicial - INICIO
pars_init <- list(psi2 = log(.13), psi3 = log(.85),
                  ar = .2, 
                  deltaMedia = c(1, 5), 
                  deltaVar = -c(1, 5))
## Chute inicial - FIM

# pars <- unlist(pars_init, use.names = F)
# Chute inicial - INICIO
set.seed(8888)

dados <- modelo(pars, dummy1, dummy2, n)

yt <- dados$yt
Varyt <- var(yt[1:50])

inicio1 <- Sys.time()
opt <- estimando(llike_model_geral, pars_init)
fim1 <- Sys.time()

inicio2 <- Sys.time()
opt2 <- estimando_cpp_geral(pars_init, yt, dummy1, dummy2, 
                            kvar, kmed, n, Varyt)
fim2 <- Sys.time()

fim1 - inicio1
fim2 - inicio2

opt
opt2

# Cenario com 2 int na media ----------------------------------------------

n <- 2500
dummy1 <- as.matrix(dummy_on_off(n, c(1, 1401),
                                 c(1400, n), 'Media'))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1401),
                                 c(1400, n)))
kvar <- 2 
kmed <- 2

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = c(3, 5), 
             deltaVar = -c(5, 1))
# Definindo algumas estruturas das series - FIM

## Chute inicial - INICIO
pars_init <- list(psi2 = log(.13), psi3 = log(.85),
                  ar = .2, 
                  deltaMedia = c(1, 5), 
                  deltaVar = -c(5, 1))
## Chute inicial - FIM

# pars <- unlist(pars_init, use.names = F)
# Chute inicial - INICIO
set.seed(78417)

dados <- modelo(pars, dummy1, dummy2, n)

yt <- dados$yt
Varyt <- var(yt[1:50])

inicio1 <- Sys.time()
opt <- estimando(llike_model_geral, pars_init)
fim1 <- Sys.time()

inicio2 <- Sys.time()
opt2 <- estimando_cpp_geral(pars_init, yt, dummy1, dummy2, 
                            kmed, kvar, n, Varyt)
fim2 <- Sys.time()

fim1 - inicio1
fim2 - inicio2

opt
opt2

# Cenario com 3 int na media e 2 var ----------------------------------------------

n <- 1500
dummy1 <- as.matrix(dummy_on_off(n, c(1, 1301, 1401),
                                 c(1300, 1400, n), 'Media'))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1001),
                                 c(1000, n)))
kvar <- 2 
kmed <- 3

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = c(1, 5, 3), 
             deltaVar = -c(1, 5))
# Definindo algumas estruturas das series - FIM

## Chute inicial - INICIO
pars_init <- list(psi2 = log(.13), psi3 = log(.85),
                  ar = .2, 
                  deltaMedia = c(1, 5, 3), 
                  deltaVar = -c(1, 5))
## Chute inicial - FIM

# pars <- unlist(pars_init, use.names = F)
# Chute inicial - INICIO
set.seed(451206)

dados <- modelo(pars, dummy1, dummy2, n)

yt <- dados$yt
Varyt <- var(yt[1:50])

inicio1 <- Sys.time()
opt <- estimando(llike_model_geral, pars_init)
fim1 <- Sys.time()

inicio2 <- Sys.time()
opt2 <- estimando_cpp_geral(pars_init, yt, dummy1, dummy2, 
                            kmed, kvar, n, Varyt)
fim2 <- Sys.time()

fim1 - inicio1
fim2 - inicio2

opt
opt2


# Cenario suave -----------------------------------------------------------
n <- 2000
kvar <- 2

t_ast <- 1900
delta_t <- 25
t_til <- t_ast + delta_t

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, 
                                 c(1, t_til - ceiling(delta_t/2)), 
                                 c(t_ast + floor(delta_t/2) - 1, n)))

pars <- list(ar = .2,
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 1, 
             deltaVar = c(-5, -4))

pars_init <- list(psi2 = log(.13), psi3 = log(.85),
                  ar = .2,  
                  deltaMedia = 1, 
                  deltaVar = c(-5, -4))

set.seed(100)
dados <- modelo_cpp_suave(pars, dummy1, dummy2, t_ast, t_til, n)
plot(dados$int2, type = 'l')
plot(dados$int2)
plot(dados$int2, xlim = c(1895, 1910))

set.seed(100)
dados2 <- modelo_smooth(pars, dummy1, dummy2, t_ast, t_til, n)

yt <- dados$yt
Varyt <- var(yt[1:50])

inicio1 <- Sys.time()
opt1 <- estimando_cpp_geral_suave(pars_init, yt, dummy1, dummy2,
                                  t_ast, t_til,
                                  1, kvar, n, Varyt)
fim1 <- Sys.time()

inicio2 <- Sys.time()
opt2 <- estimando(llike_suave, pars_init)
fim2 <- Sys.time()

fim1 - inicio1
fim2 - inicio2

# Cenario com suave -------------------------------------------------------


