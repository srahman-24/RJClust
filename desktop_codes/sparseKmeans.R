### K-mean sparse clustering 

install.packages("sparcl")
library(sparcl)
# choose tuning parameter
system.time({
km.perm <- KMeansSparseCluster.permute(X,K=4,wbounds=seq(3,7,len=15),nperms=10)
print(km.perm)
plot(km.perm)
# run sparse k-means
km.out = KMeansSparseCluster(X, K=4, wbounds=km.perm$bestw)
plot(km.out)
print(km.out)
Mutual_Information(cl, group)$ami
})