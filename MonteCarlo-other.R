K.hd = NULL 
km.AMI = skm.AMI = ssc.AMI = hd.AMI = NULL 

set.seed(44)
Seeds = sample(0:10000, 100, replace = FALSE)

for (ii in 1:100)
{
  # New Penalty
  set.seed(Seeds[ii])  
  n      = c(20,20,20,20)         # Unequal Cluster size settings
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
  
  
  
  km.out   = kmeans(X, centers = 4, iter.max = 100)
  cl       = km.out$cluster
  km.ami   = Mutual_Information(cl, group)$ami
  km.AMI   = c(km.AMI, km.ami)
  
  km.out   = KMeansSparseCluster.permute(X, K = 4, nperms = 25, wbounds = seq(1.1,10, by = 0.5),
                                        silent = FALSE, nvals = 10, centers=NULL)
  skmeans  = KMeansSparseCluster(X, K = 4, wbounds = km.out$bestw)
  cl       = skmeans[[1]]$Cs
  skm.ami  = Mutual_Information(cl, group)$ami
  skm.AMI   = c(skm.AMI, skm.ami)
  
  cl        = SSC(X, k = 4)$cluster
  ssc.ami   = Mutual_Information(cl, group)$ami
  ssc.AMI   = c(ssc.AMI, ssc.ami)
  
  hd.out   = hddc(X, K = 1:10, model = "ALL", init = "kmeans", threshold = 0.2, criterion = "bic", itermax = 1000, eps = 0.001, algo = "EM")
  cl       = hd.out$class 
  hd.ami   = Mutual_Information(cl, group)$ami
  hd.AMI   = c(hd.AMI, hd.ami)
  K.hd     = c(K.hd, hd.out$K)
  
  library(clustvarsel)
  
}

boxplot(km.AMI, skm.AMI, ssc.AMI, hd.AMI, names = c("Kmeans++", "sKmeans", "SSC", "HDDC"))
median(km.AMI)
sd(km.AMI)
median(skm.AMI)
sd(skm.AMI)
median(ssc.AMI)
sd(ssc.AMI)
median(na.omit(hd.AMI))
sd(na.omit(hd.AMI))



save(km.AMI, skm.AMI, ssc.AMI, hd.AMI, file = "sKmeans-HL.RDATA")


cl        =  clustvarsel(X1, search = 'headlong', direction = 'forward', parallel = T ,
                  emModels1 = "E", emModels2 = mclust.options("VVI"), itermax = 50, allow.EEE = F)
raft.ami  =  Mutual_Information(cl$classification, group)$ami
raft.AMI  =  c(raft.AMI, raft.ami)
K.raft    =  c(K.raft, cl$K)

N  = nrow(X)
p  = ncol(X)
GG          =  tcrossprod(X, X)/p
gg_wodiag   =  GG - diag(diag(GG))
GG_new      =  cbind(gg_wodiag + diag(colSums(gg_wodiag)/(N - 1)), diag(GG))
Gclust      =  Mclust(GG_new, modelNames = "VVI", G = 1, verbose = F)
M1          =  Gclust$parameters$mean  #N by 1 matrix
RJMean1     =  RJ_mean(1, Gclust$class, GG)
W = W1 = NULL 
#MM  = matrix(0, nrow = N, ncol = N+1)
for (kk in 2:10)
{
  Gclust      =  Mclust(GG_new, G = kk, modelNames = "VVI", verbose = F)
  RJMean      =  RJ_mean(kk, Gclust$class, GG, RJMean1)
  W1          =  c(W1, kmeans(RJMean, centers = RJMean1, iter.max = 1000, algorithm = "Lloyd")$tot.withinss)
  RJMean1     =  RJMean
}

W = NULL 

W        = W1 + (2:10)*(N + 1)/(p)
K_hat    = which.min(W)
GG_W     = Mclust(GG_new, modelNames = "VVI", G = K_hat , verbose = F)
#table(GG_W$classification, group)
#f_rez(GG_W$classification, group)$ami
ami_hs1  = c(ami_hs1, Mutual_Information(GG_W$classification, group)$ami)
K_hs1    = c(K_hs1, which.min(W))

install.packages('prclust')
library(prclust)

GCV(data,lambda1 = 2,lambda2 = 3,tau = 1,sigma = 0.25, B =100)
system.time({PR = PRclust(X, lambda1 = 1, lambda2 = 3, tau = 0.5)})


PR$group
PR$mu








