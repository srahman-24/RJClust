library(mclust)
library(RJcluster)

# New Penalty
n     = c(20,20,20,20)         # Unequal Cluster size settings
p     = 220                    # first 4 being informative and remaining ones are non-informative 
C     = 4                      # initializing every individual as their own clusters 
sigma1 = 1                     # noise level in informative variables
sigma2 = 1                     # noise level in uninformative variables
## sigma = 1 ( SNR : high signal)
## sigma = 2 ( SNR:  low  signal) 
group = c(rep(1,n[1]), rep(2,n[2]), rep(3,n[3]), rep(4,n[4]))
# Set1:  sigma = 1 ( SNR : high level) , and n      = c(20,20,20,20) 
# Set2:  sigma = 2 ( SNR : low level) ,  and n      = c(20,20,20,20) 
# Set1:  sigma = 1 ( SNR : high level) , and n      = c(20,20,200,200) (unbalanced)

X    = matrix(rnorm(sum(n)*p,0, sigma2), nrow = sum(n), ncol = p, byrow = TRUE)

#Cluster 1: N(2.5, sigma)(1-10), N(1.5, sigma)(11-20) 

X[1:n[1],1:10]                                       =   rnorm(n[1]*10, 2.5, sigma1)
X[1:n[1],(1+10):(10+10)]                             =   rnorm(n[1]*10, 1.5, sigma1)

#Cluster 2: N(0, sigma)(1-10), N(1.5, sigma) (11-20)
X[(n[1]+1):(n[1]+n[2]),1:10]                         =   rnorm(n[2]*10, 0, sigma1)
X[(n[1]+1):(n[1]+n[2]),(1+10):(10+10)]               =   rnorm(n[2]*10, 1.5, sigma1)

#Cluster 3: N(0, sigma)(1-10), N(-1.5,sigma)(11-20)
X[(n[1]+n[2]+1):(n[1]+n[2]+n[3]),1:10]               =   rnorm(n[3]*10, 0, sigma1)
X[(n[1]+n[2]+1):(n[1]+n[2]+n[3]),(1+10):(10+10)]     =   rnorm(n[3]*10, -1.5, sigma1)

#Cluster 4: N(-2.5,sigma)(1-10), N(-1.5, sigma)(11-20)
X[(n[1]+n[2]+n[3]+1):(n[1]+n[2]+n[3]+n[4]),1:10]               =   rnorm(n[4]*10, -2.5, sigma1)
X[(n[1]+n[2]+n[3]+1):(n[1]+n[2]+n[3]+n[4]),(1+10):(10+10)]     =   rnorm(n[4]*10, -1.5, sigma1)


#Z           =  scale(X, center = T, scale = T)

GG          =  tcrossprod(X, X)/p
N           =  sum(n)
gg_wodiag   =  GG - diag(diag(GG))
GG_new      =  cbind(gg_wodiag + diag(colSums(gg_wodiag)/(N-1)), diag(GG))

#Clust_GAP   = clusGap(GG_new, FUN= pam, K.max = 20, B = 100, d.power = 2)
#plot(Clust_GAP)


BIC_GG      =  Mclust(GG_new, modelNames = "VVI")
table(BIC_GG$classification, group)
f_rez(BIC_GG$classification, group)
plot(BIC_GG$BIC)


Gclust      =  Mclust(GG_new, modelNames = "VVI", G = 1, verbose = F)
M1          =  Gclust$parameters$mean  #N by 1 matrix
RJMean1     =  RJ_mean(1, Gclust$class, GG)
W = W1 = NULL 
#MM  = matrix(0, nrow = N, ncol = N+1)
for (kk in 2:10)
{
  Gclust      =  Mclust(GG_new, modelNames = "VVI", G = kk, verbose = F)
  RJMean      =  RJ_mean(kk, Gclust$class, GG)
  #Mean        =  Gclust$parameters$mean
  #W           =  c(W, kmeans(t(Mean), centers = t(M1), iter.max = 1000, algorithm = "Lloyd")$tot.withinss)
  W1          =  c(W1, kmeans(RJMean, centers = RJMean1, iter.max = 1000, algorithm = "Lloyd")$tot.withinss)
  #M1          =  Mean
  RJMean1     =  RJMean
}

#plot(W, ylab = "|mu_(k+1) - mu_(k)|^2", xlab = "k", lwd = 2, pch = 2, col = "blue", main = "Data")
plot(W1, ylab = "|mu_(k+1) - mu_(k)|^2", xlab = "k", lwd = 2, pch = 2, col = "blue", main = "Data")

#plot(W1, ylab = "|mu_(k+1) - mu_(k)|^2", xlab = "k", lwd = 2, pch = 2, col = "blue", main = "Data")
abline(h = 100, lwd = 2, col = "red")
abline(h = 1, lwd = 2, col = "brown")





Gclust      =  Mclust(GG_new, modelNames = "VVI", G = 1, verbose = F)
M1          =  Gclust$parameters$mean  #N by 1 matrix 
W = ss = penalty = NULL
LL  = ss0
for (kk in 2:20)
{
Gclust      =  Mclust(GG_new, modelNames = "VVI", G = kk, verbose = F)
Mean        =  Gclust$parameters$mean   # N by kk matrix 

ss1         =  kmeans(GG_new, centers = t(Mean), iter.max = 1000, nstart = 1)$tot.withinss
ss          =  c(ss, ss0 - ss1)
W           =  c(W, kmeans(t(Mean), centers = t(M1), iter.max = 1000, algorithm = "Lloyd")$tot.withinss)
M1          =  Mean
ss0         =  ss1
LL          = c(LL, ss1)
penalty    =  c(penalty, (LL[ii] + log(N)*kk*(N+1)))

}

#*kk*(N+1)/(p)^(0.75)

plot(W, ylab = "|mu_(k+1) - mu_(k)|^2", xlab = "k", lwd = 2, pch = 2, col = "blue", main = "High signal (unbalanced) , K_true = 4")
abline(h = N/p, lwd = 2, col = "green")

plot(ss, pch = 2, lwd = 2, col = "blue", ylab = "W_(k+1)  -  W_(k)", xlab = "k")
#abline(h = N/p, lwd = 2, col = "red")

plot(LL)
plot(penalty)

Gclust      =  Mclust(GG_new, modelNames = "VVI", G = 4, verbose = F)
table(Gclust$classification, group)

#########################

set.seed(1000)
p = ncol(X)
N = nrow(X)
GG          =  tcrossprod(X, X)/p
gg_wodiag   =  GG - diag(diag(GG))
GG_new      =  cbind(gg_wodiag + diag(colSums(gg_wodiag)/(N-1)), diag(GG))
BIC_GG      =  Mclust(GG_new, modelNames = "VVI")
table(BIC_GG$classification, group)
f_rez(BIC_GG$classification, group)
plot(BIC_GG$BIC)

#library(cluster)
Clust_GAP = clusGap(GG_new, FUN= pam, K.max = 10, B = 100, d.power = 2)
plot(Clust_GAP)


Gclust      =  Mclust(GG_new, G = 1, verbose = F)
M1          =  Gclust$parameters$mean  #N by 1 matrix 
W = NULL 
MM  = matrix(0, nrow = N, ncol = N+1)
for (kk in 2:20)
{
  Gclust      =  Mclust(GG_new,  G = kk, verbose = F)
  Mean        =  Gclust$parameters$mean   # N by kk matrix
  for(ii in 1:N){
    MM[ii,]        =  t(Mean[ ,Gclust$class[ii]])
                }
  W           =  c(W, kmeans(MM, centers = t(M1), iter.max = 1000, algorithm = "Lloyd")$tot.withinss)
  M1          =  Mean
}

plot(W, ylab = "|mu_(k+1) - mu_(k)|^2", xlab = "k", lwd = 2, pch = 2, col = "blue", main = "Data")
abline(h = 1, lwd = 2, col = "red")
abline(h = 0.1, lwd = 2, col = "brown")


if(max(W) <= 1)
{
plot(W, ylab = "|mu_(k+1) - mu_(k)|^2", xlab = "k", lwd = 2, pch = 2, col = "blue", main = "Data", ylim = c(0,1))
abline(h = 0.1, lwd = 2, col = "red")
}
set.seed(1000)
Gclust      =  Mclust(GG_new, modelNames = "VVI", G = 3)
table(Gclust$classification, group)
f_rez(Gclust$classification, group)







