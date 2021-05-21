library(RJcluster)
library(cluster)
library(sparcl)

### parallel
library(foreach)
library(parallel)
library(doParallel)
library(MASS)
numCores = detectCores()
#register a parallel backend using one of the packages that begin with do
registerDoParallel(numCores)


##### function 1

sim_data = function(n,p,sigma1, sigma2, seed){
  
  set.seed(seed)
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
  
  return(list(X = X, group = group))
  
}


###### function 2

sparcl_fx = function(X, k){
  
  km.out  = KMeansSparseCluster.permute(X, K = k, nperms = 25, wbounds = seq(1.1,10, by = 0.5),
                                        silent = TRUE, nvals = 10, centers=NULL)
  
  skmeans = KMeansSparseCluster(X, K = k, wbounds = km.out$bestw)
  cl      = skmeans[[1]]$Cs
  return(cl)
} 


##### function 3 

gap_fx = function(X){
  
  gap    = clusGap(X, FUN = pam, K.max = 10, B = nrow(X), d.power = 2)
  K      = maxSE(f = gap$Tab[, "gap"], SE.f = gap$Tab[, "SE.sim"])
  class  = sparcl_fx(X, K)
  return(list(class = class, K = K))
  
}


##### main function
set.seed(44)
Seeds = sample(0:10000, 100, replace = FALSE)


res = foreach(i = 1:100, .combine = rbind) %dopar% {
  
  sim   = sim_data(c(20,20,20,20), 220, 1,1, Seeds[i])
  X     = sim$X 
  group = sim$group
  
  cl    = gap_fx(X)
  c(Mutual_Information(cl$class, group)$ami , cl$K)
  
}


save(res, file = "gap_high_low.RData")