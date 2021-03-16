#Guo Simulation Setting 1 & 3  (Results without prescreening) 
#True number of groups = 2 

rm(list=ls(all=TRUE)) 
library(mclust)
library(RJcluster)



set.seed(44)
Seeds = sample(0:10000, 100, replace = FALSE)

AMI = C_hat = C_n = C_10 = C_20 = K0 = NULL
AMI_n = AMI_10 = AMI_20 = AMI_logp = NULL 

for (run in 1: 100)
{
  
set.seed(Seeds[run])
######### Bayesian Clustering Algorithm ############ 
#n    = c(20,20,50,50)         # Equal cluster size settings
n     = c(20,20,200,200)         # Unequal Cluster size settings
p     = 220                    # first 4 being informative and remaining ones are non-informative 
C     = 4                      # initializing every individual as their own clusters 
sigma = 1                      # noise level 
                               ## sigma = 1 ( SNR : high signal)
                               ## sigma = 2 ( SNR:  low  signal) 
group = c(rep(1,n[1]), rep(2,n[2]), rep(3,n[3]), rep(4,n[4]))
# Set1:  sigma = 1 ( SNR : high level) , and n      = c(20,20,20,20) 
# Set2:  sigma = 2 ( SNR : low level) ,  and n      = c(20,20,20,20) 
# Set1:  sigma = 1 ( SNR : high level) , and n      = c(20,20,200,200) (unbalanced)

X    = matrix(rnorm(sum(n)*p,0, 1), nrow = sum(n), ncol = p, byrow = TRUE)

#Cluster 1: N(2.5, sigma)(1-10), N(1.5, sigma)(11-20) 

X[1:n[1],1:10]                                       =   rnorm(n[1]*10, 2.5, sigma)
X[1:n[1],(1+10):(10+10)]                             =   rnorm(n[1]*10, 1.5, sigma)

#Cluster 2: N(0, sigma)(1-10), N(1.5, sigma) (11-20)
X[(n[1]+1):(n[1]+n[2]),1:10]                         =   rnorm(n[2]*10, 0, sigma)
X[(n[1]+1):(n[1]+n[2]),(1+10):(10+10)]               =   rnorm(n[2]*10, 1.5, sigma)

#Cluster 3: N(0, sigma)(1-10), N(-1.5,sigma)(11-20)
X[(n[1]+n[2]+1):(n[1]+n[2]+n[3]),1:10]               =   rnorm(n[3]*10, 0, sigma)
X[(n[1]+n[2]+1):(n[1]+n[2]+n[3]),(1+10):(10+10)]     =   rnorm(n[3]*10, -1.5, sigma)

#Cluster 4: N(-2.5,sigma)(1-10), N(-1.5, sigma)(11-20)
X[(n[1]+n[2]+n[3]+1):(n[1]+n[2]+n[3]+n[4]),1:10]               =   rnorm(n[4]*10, -2.5, sigma)
X[(n[1]+n[2]+n[3]+1):(n[1]+n[2]+n[3]+n[4]),(1+10):(10+10)]     =   rnorm(n[4]*10, -1.5, sigma)


#Z           =  scale(X, center = T, scale = T)

GG          =  tcrossprod(X, X)/p
N           =  sum(n)
gg_wodiag   =  GG - diag(diag(GG))
GG_new      =  cbind(gg_wodiag + diag(colSums(gg_wodiag)/(N-1)), diag(GG))
Gclust      =  Mclust(GG_new, verbose = F, modelNames = "VVI")
C_hat       =  c(C_hat, Gclust$G)
ami         =  f_rez(Gclust$classification, group)$ami
AMI         =  c(AMI,ami)

bic = NULL 
for(kk in 1:20){
  
  Gclust  = Mclust(GG_new, modelNames = "VVI", G = kk, verbose = F)
  loglik  = Gclust$loglik
  #print(loglik)
  nparams = nMclustParams(modelName = "VVI", d = ncol(GG_new) , G = kk)
  bic     = c(bic, 2 * loglik - nparams * log(p))
  
}
#plot(bic, type = "l", main = "High signal", xlab = "Clusters", ylab = "2.loglik - nparams.log(p)")
#abline(v = 4, col = "red", lwd = 2)
khat = which.max(bic)
K0  = c(K0, khat)
Gclust      =  Mclust(GG_new, G = khat, verbose = F, modelNames = "VVI")
ami         =  f_rez(Gclust$classification, group)$ami
AMI_logp    =  c(AMI_logp,ami)



boxplot(C_hat, K0, ylim = c(1, 10))
boxplot(AMI, AMI_logp)

G          =  tcrossprod(X, X)
L          =  sqrt(solve(diag(diag(G))))%*%G%*%sqrt(solve(diag(diag(G))))
SD         =  eigen(L)
Spec       =  SD$vectors
Smclust_n  = Mclust(Spec)
Smclust_10 = Mclust(Spec[,1:10])
Smclust_20 = Mclust(Spec[,1:20])

C_n        =  c(C_n, Smclust_n$G)
C_10       =  c(C_10, Smclust_10$G)
C_20       =  c(C_20, Smclust_20$G)

ami        =  f_rez(Smclust_n$classification, group)$ami
AMI_n      =  c(AMI_n,ami)
ami        =  f_rez(Smclust_10$classification, group)$ami
AMI_10     =  c(AMI_10,ami)
ami        =  f_rez(Smclust_20$classification, group)$ami
AMI_20     =  c(AMI_20,ami)


}



df  = list(AMI, AMI_10, AMI_20, AMI_n, C_hat, C_10, C_20, C_n)
save(df, file = "Guo_RJ_Set2_bal.RData")



boxplot(C_hat, C_n, C_20, C_10, names = c('RJ', 'Spec_N', 'Spec_20', 'Spec_10'), 
        main = "True number of clusters = 4")
abline(h = 4, lwd = 1, col = "blue", lty = "dashed")

boxplot(AMI, rep(0,100), rep(0,100), rep(0,100), names = c('RJ', 'Spec_N', 'Spec_20', 'Spec_10'), ylab = "Adjusted Mutual Index")

mean(C_hat)
sd(C_hat)


