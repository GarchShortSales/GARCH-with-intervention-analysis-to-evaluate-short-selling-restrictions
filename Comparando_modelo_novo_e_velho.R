source("src/modelo_sim.R")
source("src/util.R")
source("src/modelo_est.R")
source("src/graficos.R")

# Confirgurando os cenarios - INICIO
n <- 2000
dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 700), c(699, n)))
dummy2_velho <- as.matrix(dummy_step(n, 700))
kc <- 2

pars_novo <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = c(-10, -10))

pars_velho <- list(ar = 0.2, 
                  omega = exp(-10), alpha = .13, beta = .85, 
                  deltaMedia = 10, 
                  deltaVar = 0)
# Confirgurando os cenarios - FIM

# Gerando serie - INICIO
set.seed(25088)
yt1 <- modelo(pars_novo, dummy1, dummy2, n)$yt # Modelo Novo

set.seed(25088)
yt2 <- modelo(pars_velho, dummy1, dummy2_velho, n)$yt # Modelo Velho
# Gerando serie - FIM

# Verificando - INICIO
yt1[1]
yt1[126]
yt1[1005]
yt1[1953]

yt2[1]
yt2[126]
yt2[1005]
yt2[1953]

max(abs(yt1 - yt2))
# Verificando - FIM

# Gerando outro serie - INICIO
yt <- modelo(pars_novo, dummy1, dummy2, n)$yt
# Gerando outro serie - FIM 

# Chutes iniciais - INICIO
pars_init_novo <- list(psi2 = log(.13), psi3 = log(.85), 
                       ar = 0.2, 
                       deltaMedia = 10, 
                       deltaVar = c(-10, -10))

pars_init_velho <- list(psi1 = -10, psi2 = log(.13), psi3 = log(.85), 
                        ar = 0.2, 
                        deltaMedia = 10, 
                        deltaVar = 0)
# Chutes iniciais - FIM

optvelho <- estimando2(llike_model_ar_dumvar_ant, pars_init_velho); optvelho
optnovo <- estimando(llike_model_ar_dumvar, pars_init_novo); optnovo

log(optvelho$omega)
optnovo$deltaVar1

optvelho$deltaVar
optnovo$deltaVar2 - optnovo$deltaVar1

# Comparando tempo - INICIO
M <- 500
modelo_novo <- data.frame(alpha = rep(NA, M), beta = rep(NA, M),
                          ar = rep(NA, M),
                          deltaMedia = rep(NA, M),
                          deltaVar1 = rep(NA, M),
                          deltaVar2 = rep(NA, M),
                          ite = rep(NA, M), llike = rep(NA, M))

modelo_velho <- data.frame(omega = rep(NA, M),
                          alpha = rep(NA, M), beta = rep(NA, M),
                          ar = rep(NA, M),
                          deltaMedia= rep(NA, M),
                          deltaVar1 = rep(NA, M),
                          ite = rep(NA, M), llike = rep(NA, M))

set.seed(10000)
inicio_velho <- Sys.time()
for (i in 1:M){
  yt <- modelo(pars_novo, dummy1, dummy2, n)$yt
  modelo_velho[i, ] <- estimando2(llike_model_ar_dumvar_ant, pars_init_velho);
}
fim_velho <- Sys.time()

set.seed(10000)
inicio_novo <- Sys.time()
for (i in 1:M){
  yt <- modelo(pars_novo, dummy1, dummy2, n)$yt
  modelo_novo[i, ] <- estimando(llike_model_ar_dumvar, pars_init_novo)
}
fim_novo <- Sys.time()
# Comparando tempo - FIM

# Medidas descritivas - INICO

media_est_novo <- apply(modelo_novo, 2, mean); media_est_novo
media_est_velho <- apply(modelo_velho, 2, mean); media_est_velho

sd_est_novo <- apply(modelo_novo, 2, sd); sd_est_novo
sd_est_velho <- apply(modelo_velho, 2, sd); sd_est_velho

tempo_exec_novo <- fim_novo - inicio_novo; tempo_exec_novo
tempo_exec_velho <- fim_velho - inicio_velho; tempo_exec_velho

mean(modelo_novo$deltaVar2 - modelo_novo$deltaVar1)
mean(modelo_velho$deltaVar1)
# Medidas descritivas - FIM

