library(mclust)
library(RJcluster)
library(clustvarsel)

### parallel
library(parallel)
library(doParallel)
library(foreach)
library(MASS)


set.seed(44)
Seeds = sample(0:10000, 100, replace = FALSE)

for (ii in 1:100)
{
  # New Penalty
  set.seed(Seeds[ii])  
  n      = c(20,20,20,20)         # Unequal Cluster size settings
  p      = 220                   # first 4 being informative and remaining ones are non-informative 
  C      = 4                      # initializing every individual as their own clusters 
  sigma1 = 1                      # noise level in informative variables
  sigma2 = 2                      # noise level in uninformative variables
  ## sigma = 1 ( SNR : high signal)
  ## sigma = 2 ( SNR:  low  signal) 
  group = c(rep(1, n[1]), rep(2, n[2]), rep(3, n[3]), rep(4, n[4]))
  # Set1:  sigma = 1 ( SNR : high level) , and n  = c(20,20,20,20) 
  # Set2:  sigma = 2 ( SNR : low level) ,  and n  = c(20,20,20,20) 
  # Set1:  sigma = 1 ( SNR : high level) , and n  = c(20,20,200,200) (unbalanced)
  
  X    = matrix(rnorm(sum(n)*p,0, sigma2), nrow = sum(n), ncol = p, byrow = TRUE)
  
  #Cluster 1: N(2.5, sigma)(1-10), N(1.5, sigma)(11-20) 
  info = 10
  
  X[1:n[1],1:info]                                =   rnorm(n[1]*info, 2.5, sigma1)
  X[1:n[1],(1 + info):(info + info)]              =   rnorm(n[1]*info, 1.5, sigma1)
  
  #Cluster 2: N(0, sigma)(1-10), N(1.5, sigma) (11-20)
  X[(n[1] + 1):(n[1] + n[2]),1:info]                    =   rnorm(n[2]*info, 0, sigma1)
  X[(n[1] + 1):(n[1] + n[2]),(1 + info):(info + info)]  =   rnorm(n[2]*info, 1.5, sigma1)
  
  #Cluster 3: N(0, sigma)(1-10), N(-1.5,sigma)(11-20)
  X[(n[1] + n[2] + 1):(n[1] + n[2] + n[3]),1:info]               =   rnorm(n[3]*info, 0, sigma1)
  X[(n[1] + n[2] + 1):(n[1] + n[2] + n[3]),(1 + info):(info + info)] =   rnorm(n[3]*info, -1.5, sigma1)
  
  #Cluster 4: N(-2.5,sigma)(1-10), N(-1.5, sigma)(11-20)
  X[(n[1] + n[2] + n[3] + 1):(n[1] + n[2] + n[3] + n[4]),1:info]                   =   rnorm(n[4]*info, -2.5, sigma1)
  X[(n[1] + n[2] + n[3] + 1):(n[1] + n[2] + n[3] + n[4]),(1 + info):(info + info)] =   rnorm(n[4]*info, -1.5, sigma1)
  
  
  #km.out   =  kmeans(X, centers = 4, iter.max = 100)
  #cl       =  km.out$cluster
  #km.ami   =  Mutual_Information(cl, group)$ami
  #km.AMI   =  c(km.AMI, km.ami)
  
  
  #km.out   =  KMeansSparseCluster.permute(X, K = 4, nperms = 25, wbounds = seq(1.1,10, by = 0.5),
  #silent = TRUE, nvals = 10, centers=NULL)
  #skmeans  =  KMeansSparseCluster(X, K = 4, wbounds = km.out$bestw)
  #cl       =  skmeans[[1]]$Cs
  #skm.ami  =  Mutual_Information(cl, group)$ami
  #skm.AMI  =  c(skm.AMI, skm.ami)
  
  #cl        =  SSC(X, k = 4)$cluster
  #ssc.ami   =  Mutual_Information(cl, group)$ami
  #ssc.AMI   =  c(ssc.AMI, ssc.ami)
  
  # hd.out   = hddc(X, K = 1:10, model = "ALL", init = "kmeans", threshold = 0.2, criterion = "bic", itermax = 1000, eps = 0.001, algo = "EM")
  # cl       = hd.out$class 
  # hd.ami   = Mutual_Information(cl, group)$ami
  # hd.AMI   = c(hd.AMI, hd.ami)
  # K.hd     = c(K.hd, hd.out$K)
  
  
  raft.ami  =  Mutual_Information(cl$classification, group)$ami
  raft.AMI  =  c(raft.AMI, raft.ami)
  K.raft    =  c(K.raft, cl$K)
  
  
  
}

save(K.raft, raft.AMI, file = "clustvarsel-HL.RData")

