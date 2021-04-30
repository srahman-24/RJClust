### Install T4cluster 
### library(T4cluster)
### 1. SSC = Sparse Subspace Clustering 
### 2. sparseKmeans 
### 3. spectralclustering 
### 4. Spectrum 
### 5. HDDC 
### 6. clustvarsel 

library(RJcluster)


km_fx  = function(X, K){
  
     km.out = kmeans(X, centers = K, iter.max = 1000)
     cl     = km.out$cluster
     return(cl)
}

Mutual_Information(cl, group)$ami


sparcl_fx = function(X, k){
library(sparcl)
km.out  = KMeansSparseCluster.permute(X, K = k, nperms = 25, wbounds = seq(1.1,10, by = 0.5),
                            silent = TRUE, nvals = 10, centers=NULL)

skmeans = KMeansSparseCluster(X, K = k, wbounds = km.out$bestw)
cl      = skmeans[[1]]$Cs
return(cl)
}

Mutual_Information(cl, group)$ami


SSC_fx = function(X, K){
library(T4cluster)
    cl = SSC(X, k = K)$cluster
    return(cl)
}

Mutual_Information(cl, group)$ami


#system.time({cl = MSM(X, k = 4)$cluster})
#Mutual_Information(cl, group)$ami

 
Spec_fx = function(X)
{ 
  library(Spectrum)
  res = Spectrum(t(X))
  return(res)
}

Mutual_Information(res$assignments, group)$ami


sec_fx = function(X, K){
         library(kernlab)
         spec = specc(X1, centers = K, kernel = "vanilladot")
         ## return ??
}



HDDC_fx = function(X){

library(HDclassif)

  hd.out = hddc(X, K = 1:10, model = "ALL", init = "kmeans", threshold = 0.2, criterion = "bic", itermax = 1000, eps = 0.001, algo = "EM")
  return(hd.out)

}

Mutual_Information(cl$class, group)$ami



Raftery_fx = function(X){
  library(clustvarsel)
  cl        =  clustvarsel(X, search = 'headlong', direction = 'forward', parallel = T ,
                           emModels1 = "E", emModels2 = mclust.options("VVI"), itermax = 40, allow.EEE = F)
  return(cl)
}

Mutual_Information(cl$class, group)$ami






