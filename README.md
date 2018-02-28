# algoritmo-EM
## Construindo um algoritmo EM
#######################  Algoritmo EM  ###############################


setwd("C:\\Users\\Computador\\Documents\\Heidi\\Dissertação-Mestrado\\programas no R")  # mudar o diretorio
x=allelesub = load("alleleSub.rda")  # importar o banco de dados

  A = alleleAsub  ## conjunto de dados dos alelos A
  B = alleleBsub  ## conjunto de dados dos alelos B
  
  m=1180              ## número de pessoas observadas
  n= 1000             ## Número de posiÃ§Ãµes no genoma
  
  pis0= rep(1/3, 3)     ## Dando a mesma proporção de elementos em cada grupo
  coefs= cbind(c(0,0,0), c(1,0,-1))   ### Chute inicial dos coeficientes da regressão
  dp=c(1,1,1)         ### chute inicial dos desvios padrões
  
  for (i in 1:3){ 
    M = log2(A[,555])-log2(B[,555])
    S = (log2(A[,555])+log2(B[,555]))/2
    data=data.frame(S, M)  ## Montando nosso banco de dados
  
}

fitEM = function(M, S, pis0=rep(1/3, 3),
                 coefs=cbind(c(-1.37,-0.32,1.84), c(0.62, 0.034,-0.62)),
                 dp=c(1,1,1),
                 maxit=100){
  erro = 1
  it=1
  while (erro > 1e-6 & it<=maxit){
print(it)            ### ver o número de interações
    
    ### Passo E ####
    fs=matrix(data = NA, nrow = nrow(data), ncol = 3)    ## Matriz das densidades
    fs[,1]=dnorm(M, coefs[1,1]+S*coefs[1,2], dp[1])
    fs[,2]=dnorm(M, coefs[2,1]+S*coefs[2,2], dp[2])
    fs[,3]=dnorm(M, coefs[3,1]+S*coefs[3,2], dp[3])

    W=sweep(fs, 2, pis0, "*")             ## multiplicação entre fs e pis0
    #print(head(W))
    
    W=W/rowSums(W)                       ## Calculando os pesos medios de cada grupo
   print(head(W))
    plot(S,M, col=apply(W, 1, which.max), main = it)
    pis0= c(mean(W[,1]), mean(W[,2]), mean(W[,3]))     ## Montar um vetor com os novos pesos medios de cada grupo
        
    ### Passo M ###
        
    #fit1 = lm(M~S, data = data, weights = W[,1])   ## regressÃ£o usando os pesos do primeiro grupo
    #fit2 = lm(M~S, data = data, weights = W[,2])   ## regressÃ£o usando os pesos do segundo grupo
    #fit3 = lm(M~S, data = data, weights = W[,3])   ## regressÃ£o usando os pesos do terceiro grupo
    #abline(fit1)
    #abline(fit2, col=2)
    #abline(fit3, col=3)
    #newcoefs=t(cbind(fit1$coefficients, fit2$coefficients, fit3$coefficients))  ## vetor com os novos coeficientes
    
    beta1=c((sum(S*M*W[,1])-(sum(S*M*(W[,1]^2))/sum(W[,1])))/(sum(S^2*W[,1])-(sum(S^2*W[,1])/sum(W[,1]))), 0,
            (sum(S*M*W[,3])-(sum(S*M*(W[,3]^2))/sum(W[,3])))/(sum(S^2*W[,3])-(sum(S^2*W[,3])/sum(W[,3]))))
    
    alfa1=c((sum(W[,1]*M)-(beta1[1]*sum(W[,1]*S)))/sum(W[,1]),
             sum(W[,2]*M)/sum(W[,2]),
             (sum(W[,3]*M)-(beta1[3]*sum(W[,3]*S)))/sum(W[,3]))
    
    newcoefs=cbind(alfa1, beta1)
    
    erro = max(abs(coefs-newcoefs))
    
    coefs = newcoefs                ### substituindo o coefs pelos newcoefs calculados pelas novas regressÃµes
    
    rm(newcoefs)
    
    dp= c(sqrt(sum(((M-alfa1[1]-(beta1[1]*S))^2))/n),
          sqrt(sum(((M-alfa1[2]-(beta1[2]*S))^2))/n),
          sqrt(sum(((M-alfa1[3]-(beta1[3]*S))^2))/n) )       ## novos desvios padrões
    
    it=it+1
    
  } ##FECHA WHILE    
  return(list(coefs=coefs,W=head(W),pis0=pis0, it))      
}
   
res=fitEM(M,S)
