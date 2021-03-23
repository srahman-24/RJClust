default_RJ = function(X, penalty){

source("RJ_mean.R")  
  
p           =  ncol(X)
N           =  nrow(X)
GG          =  tcrossprod(X, X)/p
gg_wodiag   =  GG - diag(diag(GG))
GG_new      =  cbind(gg_wodiag + diag(colSums(gg_wodiag)/(N-1)), diag(GG))

bic = NULL ; aic = NULL ;
for(kk in 1:20){
  
  Gclust  = Mclust(GG_new, modelNames = "VVI", G = kk, verbose = F)
  loglik  = Gclust$loglik
  nparams = nMclustParams(modelName = "VVI", d = ncol(GG_new) , G = kk)
  if(penalty == "bic")
  {bic     = c(bic, 2 * loglik - nparams * log(p))}
  if(penalty == "aic")
  {aic     = c(aic, loglik - 2*nparams)}
  
}


if(penalty == "bic")
{
K = which.max(bic)
}
if(penalty == "aic")
{
K = which.max(aic)
}

results  = Mclust(GG_new, modelNames = "VVI", G = K, verbose = F)
class    = results$classification
prob     = results$parameters$pro
z        = results$z      
Mu       = RJ_mean(K, class, GG)


return(list(K = K, class = class, bic = bic, aic = aic, mean = Mu, prob = prob, z = z ))

}


