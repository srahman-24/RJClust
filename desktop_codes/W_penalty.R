rm(list=ls(all=TRUE)) 


library(mclust)
library(RJcluster)
library(latex2exp)

K_bic = K_aic = K_hs1 = K_hs2 = NULL 
ami_bic = ami_aic = ami_hs1 = ami_hs2 = NULL 

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


}


mu1 = as.matrix(c(rep(2.5, 10), rep(1.5, 10), rep(0, 200)))
mu2 = as.matrix(c(rep(0, 10), rep(1.5, 10), rep(0, 200)))
mu3 = as.matrix(c(rep(0, 10), rep(-1.5, 10), rep(0, 200)))
mu4 = as.matrix(c(rep(-2.5, 10), rep(-1.5, 10), rep(0, 200)))

MU = rbind(t(mu1), t(mu2), t(mu3), t(mu4))
dist(MU)

MM = tcrossprod(MU, MU)/p
dist(MM)
#N           =  sum(n)

#X = scale(X, center = T, scale = T)

#X = D
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



library(ggplot2)

pic = ggplot(data = NULL, aes(x = 1:9, y = W1, size = 3)) + geom_point(shape = 2, col = "blue", size = 3, alpha = 0.9) + 
  xlab('Number of Clusters (K)') + ylab(TeX('$\\widehat{V}_K$')) + 
  ggtitle("Hockey stick criterion")  + theme_minimal()
pic + scale_x_continuous(breaks = seq(1, 9, 1)) + theme(
  plot.title   = element_text(color = "black", size = 20, face="bold"),
  axis.title.x = element_text(color = "black", size = 18, face="bold"),
  axis.title.y = element_text(color = "black", size = 25, face="bold"), 
  axis.text.x  = element_text(face  = "bold",  size = 12),
  axis.text.y  = element_text(face  = "bold",  size = 12), 
  axis.line    = element_line(colour = "darkblue", size = 1, linetype = "solid"),
  #panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank()
   )  

pic = ggplot(data = NULL, aes(x = 1:9, y = W, size = 3)) + geom_point(shape = 2, col = "blue", size = 3, alpha = 0.9) + 
  xlab("Number of Clusters (K)") + ylab(TeX('$\\widehat{H}_K$')) + 
  ggtitle("Hockey Stick Penalty with l(P) = 1")  + theme_minimal()
pic
pic + scale_x_continuous(breaks = seq(1, 9, 1)) + theme(
  plot.title   = element_text(color = "black", size = 20, face="bold"),
  axis.title.x = element_text(color = "black", size = 18, face="bold"),
  axis.title.y = element_text(color = "black", size = 25, face="bold"), 
  axis.text.x  = element_text(face  = "bold",  size = 12),
  axis.text.y  = element_text(face  = "bold",  size = 12), 
  axis.line    = element_line(colour = "darkblue", size = 1, linetype = "solid"),
  #panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank()
)  





plot(W, ylab = "|mu_(k+1) - mu_(k)|^2", xlab = "K", lwd = 2, pch = 2, col = "blue", main = "Data", type = "l")


W       = W1 + (2:10)*(N + 1)*log(log(p))/p
K_hat   = which.min(W)
GG_W    = Mclust(GG_new, modelNames = "VVI", G = K_hat , verbose = F)
#table(GG_W$classification, group)
#f_rez(GG_W$classification, group)$ami
ami_hs2  = c(ami_hs2, Mutual_Information(GG_W$classification, group)$ami)
K_hs2    = c(K_hs2, which.min(W))


pic = ggplot(data = NULL, aes(x = 1:9, y = W, size = 3)) + geom_point(shape = 2, col = "blue", size = 3, alpha = 0.9) + 
  xlab("Number of Clusters (K)") + ylab(TeX('$V(\\widehat{\\Lambda}_{K+1} | \\widehat{\\Lambda}_{K}) + (K+1)N\\frac{loglog(P)}{P}$')) + 
  ggtitle(TeX('Hockey stick penalty $l(P) = loglog(P)$'))  + theme_minimal()
pic
pic + scale_x_continuous(breaks = seq(1, 9, 1)) + theme(
  plot.title   = element_text(color = "black", size = 14, face="bold"),
  axis.title.x = element_text(color = "black", size = 12, face="bold"),
  axis.title.y = element_text(color = "black", size = 14, face="bold"), 
  axis.text.x  = element_text(face  = "bold",  size = 12),
  axis.text.y  = element_text(face  = "bold",  size = 12), 
  axis.line    = element_line(colour = "darkblue", size = 1, linetype = "solid"),
  #panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank()
)  



bic = NULL ; aic = NULL ;
for (kk in 1:20) {
  
  Gclust  = Mclust(GG_new, modelNames = "VVI", G = kk, verbose = F)
  loglik  = Gclust$loglik
  nparams = nMclustParams(modelName = "VVI", d = ncol(GG_new) , G = kk)
  bic     = c(bic, 2 * loglik - nparams * log(N))
  aic     = c(aic, 2*loglik - 2 * nparams)
  
}

#plot(bic, type = "l", main = "High signal (BIC penalty)", xlab = "Clusters", ylab = "2.loglik - nparams.log(N)")
K_bic = c(K_bic, which.max(bic))
K_aic = c(K_aic, which.max(aic))

K_hat = which.max(bic)
GG_bic  = Mclust(GG_new, modelNames = "VVI", G = K_hat , verbose = F)
ami_bic = c(ami_bic, Mutual_Information(GG_bic$classification, group)$ami)
K_hat = which.max(aic)
GG_aic  = Mclust(GG_new, modelNames = "VVI", G = K_hat , verbose = F)
ami_aic = c(ami_aic, Mutual_Information(GG_aic$classification, group)$ami)


}

save(K_bic, K_aic, K_hs1, K_hs2, ami_bic, ami_aic, ami_hs1, ami_hs2, file = "HL(unbalanced).RData")



boxplot(K_bic, K_aic, K_hs1, K_hs2, names = c("BIC", "AIC", "HS", "HS-log"), main = "High signal - Low noise (unbalanced) ", 
        col = c("orange", "yellow", "gray", "cyan"), ylab = "Estimated number of clusters")
boxplot(ami_bic, ami_aic, ami_hs1,ami_hs2, names = c("BIC", "AIC", "HS", "HS-log"), main = "High signal - Low noise (unbalanced) ", 
        col = c("orange", "yellow", "gray", "cyan"), ylab = "Adjusted Mutual Index")


median(K_hs)
mean(K_hs)
sd(K_hs)

length(K_hs)


#plot(W, ylab = "|mu_(k+1) - mu_(k)|^2", xlab = "k", lwd = 2, pch = 2, col = "blue", main = "Data")
plot(W, ylab = "|mu_(k+1) - mu_(k)|^2", xlab = "K", lwd = 2, pch = 2, col = "blue", main = "Data", type = "l")



#Clust_GAP   = clusGap(GG_new, FUN= pam, K.max = 20, B = 100, d.power = 2)
#plot(Clust_GAP)

#####################################################################################

#BIC_GG      =  Mclust(GG_new, modelNames = "VVI")
#table(BIC_GG$classification, group)
#f_rez(BIC_GG$classification, group)
#plot(BIC_GG$BIC)

######################## BIC log(p) ######################




################################################################################

Gclust      =  Mclust(GG_new, modelNames = "VVI", G = 1, verbose = F)
M1          =  Gclust$parameters$mean  #N by 1 matrix
RJMean1     =  RJ_mean(1, Gclust$class, GG)
W = W1 = NULL 
#MM  = matrix(0, nrow = N, ncol = N+1)
for (kk in 2:10)
{
  Gclust      =  Mclust(GG_new, modelNames = "VVI", G = kk, verbose = F)
  RJMean      =  RJ_mean(kk, Gclust$class, GG)
  print(dim(RJMean))
  #Mean        =  Gclust$parameters$mean
  #W           =  c(W, kmeans(t(Mean), centers = t(M1), iter.max = 1000, algorithm = "Lloyd")$tot.withinss)
  W1          =  c(W1, kmeans(RJMean, centers = RJMean1, iter.max = 1000, algorithm = "Lloyd")$tot.withinss)
  #M1          =  Mean
  RJMean1     =  RJMean
}


#W = W1 + (2:10)*(N+1)/sqrt(p)
W = W1 + (2:10)*(N+1)/p
#W = W1 + (2:10)^2/p
#W = W1 + (2:10/p)

which.min(W)

#plot(W, ylab = "|mu_(k+1) - mu_(k)|^2", xlab = "k", lwd = 2, pch = 2, col = "blue", main = "Data")
plot(W, ylab = "|mu_(k+1) - mu_(k)|^2", xlab = "k", lwd = 2, pch = 2, col = "blue", main = "Data", type = "l")

#plot(W1, ylab = "|mu_(k+1) - mu_(k)|^2", xlab = "k", lwd = 2, pch = 2, col = "blue", main = "Data")
abline(h = 100, lwd = 2, col = "red")
abline(h = 1, lwd = 2, col = "brown")

BIC_GG      =  Mclust(GG_new, modelNames = "VVI")
table(BIC_GG$classification, group)
f_rez(BIC_GG$classification, group)
plot(BIC_GG$BIC)




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







