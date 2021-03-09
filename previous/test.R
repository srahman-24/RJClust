Z = D 
boxplot(Z[,1:100])
Z = Z-apply(Z,2, median)
boxplot(Z[,1:100])
Z = D
T1 = apply(Z,2,IQR)
Z = (Z-apply(Z,2, median))/T1
boxplot(Z[,1:100])
Z = D
Z = scale(Z, center = T, scale = F)
boxplot(Z[,1:100])

Z = log(D)
boxplot(Z[,1:100])
Z = Z-apply(Z,2, median)
boxplot(Z[,1:100])
Z = log(D)
T1 = apply(Z,2,IQR)
Z = (Z-apply(Z,2, median))/T1
boxplot(Z[,1:100])

Z = log(D)
Z = scale(Z, center = T, scale = T)
boxplot(Z[,1:100])
p = ncol(Z)
n = nrow(Z)

Z = sqrt(D)
boxplot(Z[,1:100])
Z = Z-apply(Z,2, median)
boxplot(Z[,1:100])

Z = Z1 = D
Z = log(D) 
boxplot(Z[,1:100])
Z2 = Z-apply(Z,2, median)
for(ii in 1:ncol(Z))
{
  Z1[,ii] = Z2[,ii]/mad(Z2[,ii])
}
boxplot(Z[,1:100])

source("RJMclust.R")
RJM =RJMclust(Z,group)
RJM$ami
table(group, RJM$groupG)

source("RJclust.R")
C.max = 5; iter.max = 100
RJ = RJclust(Z, group, C.max, iter.max)
plot(RJ$bic_eval)
RJ$ami 
RJ$n.clust 
table(group, RJ$groupG)


p           =    ncol(Z)
n           =    nrow(Z)      
GG          =    Z%*%t(Z)/p
gg          =    GG 
gg_wodiag   =    gg - diag(diag(gg))
gg_wdiag    =    cbind(gg_wodiag, diag(gg))
GG_new      =    cbind(gg_wodiag + diag(colSums(gg_wodiag)/(n-1)), diag(gg))
plot(mclustBIC(GG_new), ylab = "BIC", xlab = "Clusters", ylim = c(4000,12000))
lines(RJ$bic_eval, pch = 5, type = "b")
legend("topright", c("RJ-full"), pch = c(5))


#cbind(group, groupG)


