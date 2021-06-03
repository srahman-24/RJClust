RJ_hockey_stick = function(X){
  
  N           = nrow(X)
  p           = ncol(X)
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
}