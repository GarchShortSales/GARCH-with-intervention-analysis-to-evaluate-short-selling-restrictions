plot(modelo_completo$llike_com,
     modelo_completo$ite_com)
     
plot((modelo_completo$deltaVar2_com - log(5))^2, 
     modelo_completo$ite_com)

cor.test(modelo_completo$llike_com,
         modelo_completo$ite_com, 
         method = "kendall")

cor.test(modelo_completo$beta_com,
         modelo_completo$ite_com, 
         method = "spearman")

cor.test((modelo_completo$alpha_com-.13)^2, 
         modelo_completo$ite_com, 
         method = "kendall")

S <- 0*(modelo_completo$deltaVar2_com-log(5))^2 + 
  0*(modelo_completo$deltaVar1_com-log(5))^2 + 
  0*(modelo_completo$deltaMedia_com-10)^2 + 
  0*(modelo_completo$ar_com-.2)^2 + 
    (modelo_completo$beta_com-.85)^2 + 
  (modelo_completo$alpha_com-.13)^2 


plot(S, modelo_completo$ite_com)
plot(S^(1/2), modelo_completo$ite_com)

cor.test(S, modelo_completo$ite_com, 
         method = "spearman")

cor.test(S, modelo_completo$ite_com, 
         method = "spearman")

cor.test(modelo_completo$beta_com, 
         modelo_completo$ite_com,
         method = 'spearman')

cor(modelo_completo)
