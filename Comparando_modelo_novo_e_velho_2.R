source("src/modelo_sim.R")
source("src/util.R")
source("src/modelo_est.R")
source("src/graficos.R")

# Confirgurando os cenarios - INICIO
n <- 2000

dummy1 <- as.matrix(dummy_step(n, 1, "Media"))
dummy2 <- as.matrix(dummy_on_off(n, c(1, 700), 
                                 c(699, n)))
dummy2_velho <- as.matrix(dummy_step(n, 700))
kc <- 2

pars <- list(ar = 0.2, 
             omega = 1, alpha = .13, beta = .85, 
             deltaMedia = 10, 
             deltaVar = c(-5, -5))
# Confirgurando os cenarios - FIM

# Chutes iniciais - INICIO
deltaV <- c(-5, -8) # Facilitar o chute no modelo novo e velho

pars_init_novo <- list(psi2 = log(.13), psi3 = log(.85), 
                       ar = 0.2, 
                       deltaMedia = 10, 
                       deltaVar = deltaV)

pars_init_velho <- list(psi1 = deltaV[1], 
                        psi2 = log(.13), psi3 = log(.85), 
                        ar = 0.2, 
                        deltaMedia = 10, 
                        deltaVar = deltaV[2] - deltaV[1])
# Chutes iniciais - FIM

# Gerando serie - INICIO
set.seed(74152367)
yt <- modelo(pars, dummy1, dummy2, n)$yt
# Gerando serie - FIM

# Otimizando - INICIO
optvelho <- estimando2(llike_model_ar_dumvar_ant, pars_init_velho); optvelho
optnovo <- estimando(llike_model_ar_dumvar, pars_init_novo); optnovo
# Otimizando - FIM

# Comparando parametrizacao - INICIO
log(optvelho$omega)
optnovo$deltaVar1

optvelho$deltaVar
optnovo$deltaVar2 - optnovo$deltaVar1
# Comparando parametrizacao - FIM