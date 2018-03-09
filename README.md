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
                 maxit=1000){
  llvelho=-Inf
  it=1
  mudanca=1
  while (mudanca > 1e-6 & it<=maxit){
    print(it)            ### ver o nÃºmero de interaÃ§Ãµes
    
    ### Passo E ####
    fs=matrix(data = NA, nrow = nrow(data), ncol = 3)    ## Matriz das densidades
    fs[,1]=dnorm(M, coefs[1,1]+S*coefs[1,2], dp[1])
    fs[,2]=dnorm(M, coefs[2,1]+S*coefs[2,2], dp[2])
    fs[,3]=dnorm(M, coefs[3,1]+S*coefs[3,2], dp[3])
    
    W=sweep(fs, 2, pis0, "*")             ## multiplicaÃ§Ã£o entre fs e pis0
    #print(head(W))
    
    W=W/rowSums(W)                       ## Calculando os pesos medios de cada grupo
    #print(head(W))
    #plot(S,M, col=apply(W, 1, which.max), main = it)
    pis0= c(mean(W[,1]), mean(W[,2]), mean(W[,3]))     ## Montar um vetor com os novos pesos medios de cada grupo
    
    ### Passo M ###
    
    #fit1 = lm(M~S, data = data, weights = W[,1])   ## regressÃ£o usando os pesos do primeiro grupo
    #fit2 = lm(M~S, data = data, weights = W[,2])   ## regressÃ£o usando os pesos do segundo grupo
    #fit3 = lm(M~S, data = data, weights = W[,3])   ## regressÃ£o usando os pesos do terceiro grupo
    
    #newcoefs=t(cbind(fit1$coefficients, fit2$coefficients, fit3$coefficients))  
    
    beta1=c((sum(S*M*W[,1])-(sum(S*M*(W[,1]^2))/sum(W[,1])))/(sum(S^2*W[,1])-(sum(S^2*W[,1]^2)/sum(W[,1]))),
            0,
            (sum(S*M*W[,3])-(sum(S*M*(W[,3]^2))/sum(W[,3])))/(sum(S^2*W[,3])-(sum(S^2*W[,3]^2)/sum(W[,3]))))
    
    alfa1=c((sum(W[,1]*M)-(beta1[1]*sum(W[,1]*S)))/sum(W[,1]),
            sum(W[,2]*M)/sum(W[,2]),
            (sum(W[,3]*M)-(beta1[3]*sum(W[,3]*S)))/sum(W[,3]))
    
    coefs=cbind(alfa1, beta1)         ## vetor com os novos coeficientes
    
    Rq=(M-(cbind(1,S)%*%t(coefs)))^2     ## calculando o residuo ao quadrado
    dp=sqrt(sum(W*Rq)/sum(W))            ## novos desvios padrões iguais para os 3 grupos
    dp=rep(dp,3)                  
    
    
    #abline(alfa1[1],beta1[1])
    #abline(alfa1[2],beta1[2],col=2)
    #abline(alfa1[3],beta1[3],col=3)
    
    fs[,1]=dnorm(M, coefs[1,1]+S*coefs[1,2], dp[1])
    fs[,2]=dnorm(M, coefs[2,1]+S*coefs[2,2], dp[2])
    fs[,3]=dnorm(M, coefs[3,1]+S*coefs[3,2], dp[3])
    
    log_verossimilhanca=sum(W * sweep(log(fs), 2, log(pis0), "+"))   ## calculando a log-verossimilhança
    
    mudanca = log_verossimilhanca - llvelho
    
    message("antigo:", llvelho)
    
    message("novo", log_verossimilhanca)
    
    llvelho = log_verossimilhanca
    
    if(mudanca<0) stop("verossimilhança diminuiu")
    
    print(log_verossimilhanca)
    
    it=it+1
    
    print(mudanca>1e-6)
    
  } ##FECHA WHILE    
  return(list(coefs=coefs,W=head(W),pis0=pis0, it))      
}


res=fitEM(M,S)
res

