RJ_bic = function(X){
  
  N           = nrow(X)
  p           = ncol(X)
  GG          =  tcrossprod(X, X)/p
  gg_wodiag   =  GG - diag(diag(GG))
  GG_new      =  cbind(gg_wodiag + diag(colSums(gg_wodiag)/(N - 1)), diag(GG))
  Gclust      =  Mclust(GG_new, modelNames = "VVI", verbose = F)
  return(list(class = Gclust$classification, K = Gclust$G))
  
}
