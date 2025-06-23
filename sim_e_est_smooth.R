source("src/modelo_sim.R")
source("src/util.R")
source("src/modelo_est.R")
source("src/graficos.R")
source('src/modelo_est.R')

# Cenario 1 ---------------------------------------------------------------

n <- 2500
kvar <- 2
t_ast <- 1980
t_til <- 2000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, t_til), c(t_ast, n)))

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = -c(5, 4))
# Definindo algumas estruturas das series - FIM

# Gerando serie - INICIO
set.seed(100000)
dados <- modelo_smooth(pars, dummy1, dummy2, t_ast, t_til, n)
yt <- dados$yt
Varyt <- var(yt[1:50])

plot(dados$time, dados$int2, type = 'l')

line(dados, "time", "yt", "Serie simulado") + 
  geom_vline(xintercept = c(t_ast, t_til), 
             linetype = "dashed", color = "blue", size = .9) +
  xlim(1950, 2050)
# Gerando serie - FIM

## Chute inicial - INICIO
pars_init<- list(psi2 = log(.13), psi3 = log(.6),
                 ar = .2, 
                 deltaMedia = 10, 
                 deltaVar = -c(5, 1))
## Chute inicial - FIM

# pars <- unlist(pars_init, use.names = F)

# Estimando - INICIO
inicio <- Sys.time()
opt <- estimando(llike_suave, pars_init); opt
Sys.time() - inicio

inicio <- Sys.time()
opt2 <- estimando_cpp_geral_suave(pars_init, yt,
                                 dummy1, dummy2,
                                 t_ast, t_til,
                                 1, kvar, n, Varyt); opt2
Sys.time() - inicio
# Estimando - FIM

# Cenario 2 ---------------------------------------------------------------

n <- 2500
kvar <- 2
t_ast <- 1499
t_til <- 1550
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, t_til), c(t_ast, n)))

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = -c(5, 1))
# Definindo algumas estruturas das series - FIM

# Gerando serie - INICIO
dados <- modelo_smooth(pars, dummy1, dummy2, t_ast, t_til, n)
yt <- dados$yt
Varyt <- var(yt[1:50])

plot(dados$time, dados$int2, type = 'l')

line(dados, "time", "yt", "Serie simulado") + 
  geom_vline(xintercept = c(t_ast, t_til), 
             linetype = "dashed", color = "blue", size = .9)
# Gerando serie - FIM

## Chute inicial - INICIO
pars_init<- list(psi2 = log(.13), psi3 = log(.6),
                 ar = .2, 
                 deltaMedia = 10, 
                 deltaVar = -c(5, 1))
## Chute inicial - FIM

# pars <- unlist(pars_init, use.names = F)

# Estimando - INICIO
inicio <- Sys.time()
opt <- estimando(llike_suave, pars_init); opt
Sys.time() - inicio

inicio <- Sys.time()
opt2 <- estimando_cpp_geral_suave(pars_init, yt,
                                  dummy1, dummy2,
                                  t_ast, t_til,
                                  1, kvar, n, Varyt); opt2
Sys.time() - inicio
# Estimando - FIM

# Cenario 3 ---------------------------------------------------------------

n <- 2500
kvar <- 2
t_ast <- 1499
t_til <- 2000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, t_til), c(t_ast, n)))

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = -c(1, 5))
# Definindo algumas estruturas das series - FIM

# Gerando serie - INICIO
dados <- modelo_smooth(pars, dummy1, dummy2, t_ast, t_til, n)
yt <- dados$yt
Varyt <- var(yt[1:50])

plot(dados$time, dados$int2, type = 'l')

line(dados, "time", "yt", "Serie simulado") + 
  geom_vline(xintercept = c(t_ast, t_til), 
             linetype = "dashed", color = "blue", size = .9)
# Gerando serie - FIM

## Chute inicial - INICIO
pars_init<- list(psi2 = log(.13), psi3 = log(.6),
                 ar = .2, 
                 deltaMedia = 10, 
                 deltaVar = -c(5, 1))
## Chute inicial - FIM

# Estimando - INICIO
inicio <- Sys.time()
opt <- estimando(llike_suave, pars_init); opt
Sys.time() - inicio

inicio <- Sys.time()
opt2 <- estimando_cpp_geral_suave(pars_init, yt,
                                  dummy1, dummy2,
                                  t_ast, t_til,
                                  1, kvar, n, Varyt); opt2
Sys.time() - inicio
# Estimando - FIM
# Cenario 4 ---------------------------------------------------------------

n <- 2500
kvar <- 2
t_ast <- 1499
t_til <- 1550
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, t_til), c(t_ast, n)))

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = -c(1, 5))
# Definindo algumas estruturas das series - FIM

# Gerando serie - INICIO
dados <- modelo_smooth(pars, dummy1, dummy2, t_ast, t_til, n)
yt <- dados$yt
Varyt <- var(yt[1:50])

plot(dados$time, dados$int2, type = 'l')

line(dados, "time", "yt", "Serie simulado") + 
  geom_vline(xintercept = c(t_ast, t_til), 
             linetype = "dashed", color = "blue", size = .9)
# Gerando serie - FIM

## Chute inicial - INICIO
pars_init<- list(psi2 = log(.13), psi3 = log(.6),
                 ar = .2, 
                 deltaMedia = 10, 
                 deltaVar = -c(5, 1))
## Chute inicial - FIM

# Estimando - INICIO
inicio <- Sys.time()
opt <- estimando(llike_suave, pars_init); opt
Sys.time() - inicio

inicio <- Sys.time()
opt2 <- estimando_cpp_geral_suave(pars_init, yt,
                                  dummy1, dummy2,
                                  t_ast, t_til,
                                  1, kvar, n, Varyt); opt2
Sys.time() - inicio
# Estimando - FIM



# Cenario 5 ---------------------------------------------------------------

n <- 2240
kvar <- 4
delta_t <- 20
t_ast <- c(1951, 2030)
t_til <- t_ast + delta_t
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1551, 1970, 2050), 
                                 c(1550, 1950, 2030, n)))

# teste <- dummy2 %>% as.data.frame()
# teste['nrow'] <- 1:n
# 
# ggplot(teste, aes(x = nrow, y = Dummy_Var1)) +
#   geom_step(size = 0.5, colour = "red") +
#   geom_step(aes(x = nrow, y = Dummy_Var2), size = 0.5, colour = "#112446") + 
#   geom_step(aes(x = nrow, y = Dummy_Var3), size = 0.5, colour = "blue") + 
#   geom_step(aes(x = nrow, y = Dummy_Var4), size = 0.5, colour = "green") 

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = -c(5, 4, 2, 3))
# Definindo algumas estruturas das series - FIM

# Gerando serie - INICIO
set.seed(87171)
dados <- modelo_smooth2(pars, dummy1, dummy2, t_ast, t_til, n)
yt <- dados$yt
Varyt <- var(yt[1:50])

set.seed(87171)
dados2 <- modelo_cpp_suave2(pars, dummy1, dummy2, t_ast, t_til, n)
max(abs(dados2$yt - dados$yt))

plot(dados$int2, type = 'l')
plot(dados$int2)

plot(dados2$int2, type = 'l')
plot(dados2$int2)

line(dados, "time", "yt", "Serie simulado") + 
  geom_vline(xintercept = c(t_ast, t_til), 
             linetype = "dashed", color = "blue", size = .9) + 
  xlim(1900, 2100)

# Gerando serie - FIM

## Chute inicial - INICIO
pars_init<- list(psi2 = log(.13), psi3 = log(.6),
                 ar = .2, 
                 deltaMedia = 10, 
                 deltaVar = -c(5, 4, 2, 3))
## Chute inicial - FIM

# Estimando - INICIO
inicio <- Sys.time()
opt <- estimando(llike_suave2, pars_init); opt
Sys.time() - inicio

inicio <- Sys.time()
opt2 <- estimando_cpp_geral_suave2(pars_init, yt,
                                  dummy1, dummy2,
                                  t_ast, t_til,
                                  1, kvar, n, Varyt); opt2
Sys.time() - inicio



# Cenario 6 ---------------------------------------------------------------

t1 <- 1550
t2 <- 1950
t3 <- t2 + 20
t4 <- t3 + 50
t5 <- t4 + 20
n <- t5 + 200

t_ast <- c(t2, t4)
t_til <- c(t3, t5)

dummy1 <- as.matrix(dummy_step(n, 1), "Media")
dummy2_red <- as.matrix(dummy_step(n, 1))

dummy2 <- as.matrix(dummy_on_off(n, 
                                 c(1, t1 + 1, t3, t5),
                                 c(t1, t2, t4, n)))
a <- 1.5
pars <- list(
  ar = .2,
  omega = 1,
  alpha = .13,
  beta = .85,
  deltaMedia = 1,
  deltaVar = c(-5, -5 + a, -5 + 3 * a, -5 + 2 * a)
)

pars_init <- list(
  psi2 = log(.13),
  psi3 = log(.85),
  ar = .2,
  deltaMedia = 1,
  deltaVar = c(-5, -5 + a, -5 + 3 * a, -5 + 2 * a)
)

set.seed(1000)
dados <- modelo_smooth2(pars, dummy1, dummy2, t_ast, t_til, n)

set.seed(1000)
dados2 <- modelo_cpp_suave2(pars, dummy1, as.matrix(dummy2), 
                           c(t2, t4), c(t3, t5), n)

plot(dados$int2, type = 'l')
plot(dados$int2)

plot(dados2$int2, type = 'l')
plot(dados2$int2)

max(abs(dados$yt - dados2$yt))

line(dados, "time", "yt", "Serie simulado") + 
  geom_vline(xintercept = c(t2, t4, t3, t5), 
             linetype = "dashed", color = "blue", size = .9) + 
  xlim(t2 - 100, t5 + 200)

# Benchmark ---------------------------------------------------------------
require(dplyr)
library(rbenchmark)

n <- 2000
kvar <- 2
t_ast <- 1499
t_til <- t_ast + 10
dummy1_suave <- as.matrix(dummy_step(n, 1, "Media"))
dummy2_suave <- as.matrix(dummy_on_off(n, c(1, t_til), c(t_ast, n)))

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = -c(5, 4))

pars_init <- list(psi2 = log(.13), psi3 = log(.85),
                  ar = .2, 
                  deltaMedia = 1, 
                  deltaVar = -c(5, 4))
# Definindo algumas estruturas das series - FIM

# Gerando serie - INICIO
set.seed(84561)
dados <- modelo_smooth(pars, dummy1_suave, dummy2_suave, t_ast, t_til, n)
yt <- dados$yt
Varyt <- var(yt[1:50])

line(dados, "time", "yt", "Serie simulado") + 
  geom_vline(xintercept = c(t_ast, t_til), 
             linetype = "dashed", color = "blue", size = .5)
# Gerando serie - FIM

# Estimando modelo intervencao abrupta - INICIO
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, t_til - 5 + 1), 
                                 c(t_ast + 5, n)))

opt <- estimando_cpp_geral(pars_init, yt, dummy1, dummy2, 
                           1, kvar, n, Varyt)
# Estimando modelo intervencao abrupta - FIM

# Estimando modelo intervencao suave - INICIO
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, t_til - 5 + 1), 
                                 c(t_ast + 5, n)))

opt2 <- estimando_cpp_geral_suave(pars_init, yt, dummy1, dummy2, 
                                  t_ast, t_til,
                                  1, kvar, n, Varyt)
# Estimando modelo intervencao sauve - FIM

# Benchmark - INICIO
bench::mark(
  t1 = {
    set.seed(5000)
    dados <- modelo_smooth(pars, dummy1_suave, dummy2_suave, t_ast, t_til, n)
    yt <- dados$yt
    estimando_cpp_geral(pars_init, yt, dummy1, dummy2, 
                        1, kvar, n, Varyt)
  },
  t2 = {
    set.seed(5000)
    dados <- modelo_smooth(pars, dummy1_suave, dummy2_suave, t_ast, t_til, n)
    yt <- dados$yt
    estimando_cpp_geral_suave(pars_init, yt, dummy1, dummy2, 
                              t_ast, t_til,
                              1, kvar, n, Varyt)
  },
  check = FALSE
)
# Benchmark - FIM



