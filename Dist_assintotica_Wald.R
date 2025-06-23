source("src/modelo_sim.R")
source("src/util.R")
source("src/modelo_est.R")
source("src/graficos.R")

require(numDeriv)

# Exemplo 1 ---------------------------------------------------------------

# Definindo algumas estruturas das series - INICIO
n <- 1000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 500), c(499, n)))
dummy2_red <- as.matrix(dummy_step(n, 1))
kc <- 2
M <- 500

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = log(5))

## Chute inicial - INICIO
pars_init <- list(psi2 = log(.13), psi3 = log(.85),
                 ar = .2, 
                 deltaMedia = 10, 
                 deltaVar = c(log(5), log(5)))
## Chute inicial - FIM

# Monte Carlo - INICIO

## Definindo estruturas - INICIO
data <- matrix(nrow = n, ncol = M)
modelo_completo <- data.frame(alpha= rep(NA, M), beta = rep(NA, M), 
                              ar = rep(NA, M),
                              deltaMedia = rep(NA, M), 
                              deltaVar1 = rep(NA, M), deltaVar2 = rep(NA, M),
                              ite = rep(NA, M), llike = rep(NA, M), 
                              var_dummy = rep(NA, M))
## Definindo estruturas - FIM

set.seed(10003333)
cont <- 1
inicio <- Sys.time()
while (cont <= M){
  
  yt <- modelo(pars, dummy1, dummy2_red, n)$yt
  est <- estimando(llike_model_ar_dumvar, pars_init)
  
  FI <- solve(-1*hessian(llike_normal_model_ar_dumvar, 
                         unlist(est[1:6], use.names = FALSE)))
  
  var_dummy <- FI[6,6]
  
  if (var_dummy > 0){
    data[, cont] <- yt
    modelo_completo[cont, ] <- c(est, var_dummy)
    print(cont)
    cont <- cont + 1
  }
  
}
fim <- Sys.time()
fim - inicio

## Acessando distribuicao assintotica da estatistica de teste - INICIO 

valor_critico <- qnorm(.975)

### Retirando estimativas ruins - INICIO
modelo_comp_retirado <- modelo_completo %>% arrange(-ite) 
modelo_comp_retirado <- modelo_comp_retirado[floor(.15*M):M, ]
### Retirando estimativas ruins - FIM

modelo_comp_retirado <- modelo_comp_retirado %>% 
  mutate(test_wald_dummy = (deltaVar)/sqrt(var_dummy),
         rejeito_h0 =  abs(test_wald_dummy) >=  valor_critico)

sum(modelo_comp_retirado$rejeito_h0)/(.85*M)

pad(modelo_comp_retirado) %>% QQplot_MC( "test_wald_dummy", M, n)
pad(modelo_comp_retirado) %>% histo_MC("test_wald_dummy", M, n)

## Acessando distribuicao assintotica da estatistica de teste - FIM

# Monte Carlo - FIM