rm(list=ls(all=TRUE)) 

library(mclust)
library(RJcluster)

ami_hs1 = ami_bic = NULL 
K_hs1  = K_bic = NULL 

set.seed(44)
Seeds = sample(0:10000, 100, replace = FALSE)

for (ii in 1:100)
{
  # New Penalty
  set.seed(Seeds[ii])  
  n      = c(20,20,200,200)         # Unequal Cluster size settings
  p      = 220                    # first 4 being informative and remaining ones are non-informative 
  C      = 4                      # initializing every individual as their own clusters 
  sigma1 = 1                      # noise level in informative variables
  sigma2 = 1                      # noise level in uninformative variables
  ## sigma = 1 ( SNR : high signal)
  ## sigma = 2 ( SNR:  low  signal) 
  group = c(rep(1, n[1]), rep(2, n[2]), rep(3, n[3]), rep(4, n[4]))
  # Set1:  sigma = 1 ( SNR : high level) , and n  = c(20,20,20,20) 
  # Set2:  sigma = 2 ( SNR : low level) ,  and n  = c(20,20,20,20) 
  # Set1:  sigma = 1 ( SNR : high level) , and n  = c(20,20,200,200) (unbalanced)
  
  X    = matrix(rnorm(sum(n)*p,0, sigma2), nrow = sum(n), ncol = p, byrow = TRUE)
  
  #Cluster 1: N(2.5, sigma)(1-10), N(1.5, sigma)(11-20) 
  
  X[1:n[1],1:10]                                       =   rnorm(n[1]*10, 2.5, sigma1)
  X[1:n[1],(1 + 10):(10 + 10)]                         =   rnorm(n[1]*10, 1.5, sigma1)
  
  #Cluster 2: N(0, sigma)(1-10), N(1.5, sigma) (11-20)
  X[(n[1] + 1):(n[1] + n[2]),1:10]                         =   rnorm(n[2]*10, 0, sigma1)
  X[(n[1] + 1):(n[1] + n[2]),(1 + 10):(10 + 10)]           =   rnorm(n[2]*10, 1.5, sigma1)
  
  #Cluster 3: N(0, sigma)(1-10), N(-1.5,sigma)(11-20)
  X[(n[1] + n[2] + 1):(n[1] + n[2] + n[3]),1:10]               =   rnorm(n[3]*10, 0, sigma1)
  X[(n[1] + n[2] + 1):(n[1] + n[2] + n[3]),(1 + 10):(10 + 10)] =   rnorm(n[3]*10, -1.5, sigma1)
  
  #Cluster 4: N(-2.5,sigma)(1-10), N(-1.5, sigma)(11-20)
  X[(n[1] + n[2] + n[3] + 1):(n[1] + n[2] + n[3] + n[4]),1:10]                  =   rnorm(n[4]*10, -2.5, sigma1)
  X[(n[1] + n[2] + n[3] + 1):(n[1] + n[2] + n[3] + n[4]),(1 + 10):(10 + 10)]    =   rnorm(n[4]*10, -1.5, sigma1)
  
  
  N  = nrow(X)
  p  = ncol(X)
  GG          =  tcrossprod(X, X)/p
  gg_wodiag   =  GG - diag(diag(GG))
  GG_new      =  cbind(gg_wodiag + diag(colSums(gg_wodiag)/(N - 1)), diag(GG))
  Gmclust     =  Mclust(GG_new, modelNames = "VVI", verbose = F)
  K           =  Gmclust$G
  ami_bic     = c(ami_bic, Mutual_Information(Gmclust$classification, group)$ami)
  K_bic       = c(K_bic, K)
  
  Gclust      =  Mclust(GG_new, modelNames = "VVI", G = 1, verbose = F)
  M1          =  Gclust$parameters$mean  #N by 1 matrix
  RJMean1     =  RJ_mean(1, Gclust$class, GG)
  W1 = NULL 
  for (kk in 2:10)
  {
    Gclust      =  Mclust(GG_new, G = kk, modelNames = "VVI", verbose = F)
    RJMean      =  RJ_mean(kk, Gclust$class, GG, RJMean1)
    W1          =  c(W1, kmeans(RJMean, centers = RJMean1, iter.max = 1000, algorithm = "Lloyd")$tot.withinss)
    RJMean1     =  RJMean
  }

  W             =  W1 + (2:10)*(N + 1)/(p)
  K_hat         =  which.min(W)
  GG_W          =  Mclust(GG_new, modelNames = "VVI", G = K_hat , verbose = F)
  ami_hs1       =  c(ami_hs1, Mutual_Information(GG_W$classification, group)$ami)
  K_hs1         =  c(K_hs1, K_hat)
  
  
}


median(K_hs1)
sd(K_hs1)
median(K_bic)
sd(K_bic)

median(ami_hs1)
sd(ami_hs1)
median(ami_bic)
sd(ami_bic)




save(K_hs1, K_bic, ami_hs1, ami_bic, file = "RJ-unbalanced.RDATA")


cl        =  clustvarsel(X1, search = 'headlong', direction = 'forward', parallel = T ,
                         emModels1 = "E", emModels2 = mclust.options("VVI"), itermax = 50, allow.EEE = F)
raft.ami  =  Mutual_Information(cl$classification, group)$ami
raft.AMI  =  c(raft.AMI, raft.ami)
K.raft    =  c(K.raft, cl$K)


