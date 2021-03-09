rm(list = ls())

library(mclust)
library(devtools)
library(e1071)


Alizadeh    = read.csv("/Users/srahman/Documents/Clustering/alizadeh-2000-v2_database1.csv", header = T, row.names = 1)

D  = t(Alizadeh) 

source("RJMclust.R")

##### R-J Mclust-implementation ######

system.time({ 
#Z1 = scale(D, center = T, scale =  F)
#Z1 = scale(log(D), center = T, scale =  F)
Z1 = D
#boxplot(Z1[,1:100])
p = ncol(Z1)
n = nrow(Z1)

GG  = Z1%*%t(Z1)/p
gg = GG 
gg_wodiag       = gg - diag(diag(gg))
gg_wdiag        = cbind(gg_wodiag, diag(gg))
GG_new          = cbind(gg_wodiag + diag(diag(colMeans(gg_wodiag))), diag(gg))

############ Applying the Hierachical and EM Algorithm Raftery's Package #############
MclustGG = Mclust(GG_new, modelNames = "VVI")  # On GG 
summary(MclustGG)
})
groupG = MclustGG$class




source("bicRJ.R")

####### R-J exact implementation ####### 
system.time({
C_max = 10; bic_eval = NULL ; Z = Sigma = mu = list()
#C = 4
#init         =  Mclust(GG_new, modelNames = "VVI", G = C) 

  for(C in 2:C_max)
  { 
    system.time({
      init         =  Mclust(GG_new, modelNames = "VVI", G = C)})
    
    
    iter.max     =  50     
    d            =  n+1
    mu           =  t(init$parameters$mean)
    Sigma        =  init$parameters$variance$sigma
    N            =  init$n
    prob         =  init$parameters$pr
    #bic.val      =  bic.G(C, GG_new , mu , Sigma, N, prob, iter.max)
    bic.val      =  bic.G(C, GG_new , mu , Sigma, N, prob, iter.max)
    bic_eval     =  c(bic_eval, bic.val$bic.value)
    Z[[C-1]]     =  bic.val$z  
  }
  

bic_Pom   = bic_eval

Clust.no     = which(bic_eval == max(bic_eval)) + 1 ;
Clust.bic    = max(bic_eval)
Class.Matrix = Z[[which(bic_eval == max(bic_eval))]]
Class.Matrix = Z[[3]]
Class        = colSums(Class.Matrix)
Class.Sigma  = Gcov(GG_new, Clust.no, z = Class.Matrix, N)
Class.member = list(); Class.sigma = Sigma = Class.pro = list(); 
for(kk in 1:Clust.no)
{
  Class.member[[kk]]  = which(Class.Matrix[,kk]==1)
  #Class.mean[[kk]]    = colMeans(y[Class.member[[kk]], ]) 
  #Class.sigma[[kk]]   = var(y[Class.member[[kk]], ])
  #Class.pro[[kk]]     = length(Class.member[[kk]])/N  
}

})


