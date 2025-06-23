library(Rcpp)
library(RcppArmadillo)
library(rbenchmark)
library(microbenchmark)

source('src/modelo_sim.R')
source('src/util.R')

# Comparando GARCH --------------------------------------------------------

pars <- c(.1, .2, .3)

inicio <- Sys.time()
set.seed(1000)
g1 <- garch11(pars, 10000)
fim <- Sys.time()

pars2 <- list(omega = .1, alpha = .2, beta = .3)
inicio2 <- Sys.time()
set.seed(1000)
g2 <- Garch(pars2, 10000)
fim2 <- Sys.time()
set.seed(1000)
inicio5 <- Sys.time()
pars5 <- list(omega = .1, alpha = .2, beta = .3)
g5 <- garch(pars5, 10000)
fim5 <- Sys.time()

fim - inicio
fim2 - inicio2
fim5 - inicio5

microbenchmark(
  t1 = Garch(pars2, 10000),
  t2 = garch11(pars, 10000),
  t3 = garch(pars5, 10000)
)

benchmark(
  t1 = Garch(pars2, 10000),
  t2 = garch11(pars, 10000),
  t3 = garch(pars5, 10000)
)
# Comparando AR-GARCH -----------------------------------------------------

set.seed(1000)
inicio3 <- Sys.time()
pars3 <- c(.1, .1, .2, .3)
g3 <- ar_garch11(pars3, 1000)
fim3 <- Sys.time()

set.seed(1000)
inicio4 <- Sys.time()
pars4 <- list(ar = .1, omega = .1, alpha = .2, beta = .3)
g4 <- AR_Garch(pars4, 1000)
fim4 <- Sys.time()

set.seed(1000)
inicio5 <- Sys.time()
pars4 <- list(ar = .1, omega = .1, alpha = .2, beta = .3)
g5 <- ar_garch(pars4, 1000)
fim5 <- Sys.time()

fim3 - inicio3
fim4 - inicio4
fim5 - inicio5

microbenchmark(
  t1 = ar_garch11(pars3, 10000),
  t2 = AR_Garch(pars4, 10000),
  t3 = ar_garch(pars4, 10000)
)

benchmark(
  t1 = ar_garch11(pars3, 10000),
  t2 = AR_Garch(pars4, 10000),
  t3 = ar_garch(pars4, 10000)
)

# Comparando modelo final --------------------------------------------------------

n <- 2000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1500), c(1499, n)))

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = c(-3, -5))

set.seed(12345679)
dados1 <- modelo(pars, dummy1, dummy2, n)

set.seed(12345679)
dados2 <- modelo_cpp(pars, dummy1, dummy2, n)

plot(dados1$yt, type = 'l')
plot(dados2$yt, type = 'l')

microbenchmark(
  t1 = modelo(pars, dummy1, dummy2, n),
  t2 = modelo_cpp(pars, dummy1, dummy2, n),
)

benchmark(
  t1 = modelo(pars, dummy1, dummy2, n),
  t2 = modelo_cpp(pars, dummy1, dummy2, n),
)

# Comparando modelo final --------------------------------------------------------

n <- 2000
dummy1 <- as.matrix(dummy_on_off(n, c(1, 1500), c(1499, n), "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 1500), c(1499, n)))

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = c(5, 10), 
             deltaVar = c(-3, -5))

set.seed(8741)
dados1 <- modelo(pars, dummy1, dummy2, n)

set.seed(8741)
dados2 <- modelo_cpp(pars, dummy1, dummy2, n)

plot(dados1$yt, type = 'l')
plot(dados2$yt, type = 'l')

microbenchmark(
  t1 = modelo(pars, dummy1, dummy2, n),
  t2 = modelo_cpp(pars, dummy1, dummy2, n),
)

benchmark(
  t1 = modelo(pars, dummy1, dummy2, n),
  t2 = modelo_cpp(pars, dummy1, dummy2, n),
)
