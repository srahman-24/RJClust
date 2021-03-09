### author: Shahina Rahman 


RJMclust = function(Z, group){
  
  source("AMI.R")  
  p = ncol(Z)
  n = nrow(Z)
  
  GG  = Z%*%t(Z)/p
  gg  = GG 
  gg_wodiag       = gg - diag(diag(gg))
  gg_wdiag        = cbind(gg_wodiag, diag(gg))
  GG_new          = cbind(gg_wodiag + diag(colSums(gg_wodiag)/(n-1)), diag(gg))
  
  
  system.time({
  MclustGG = Mclust(GG_new, G = 1:15, modelNames = "VVI")  
  groupG   = MclustGG$class
  })
  
  
  ami = f_rez(group, groupG)$ami
  
  
  return(list(groupG = groupG, n.clust = MclustGG$G, bic = MclustGG$bic, BIC = MclustGG$BIC, mean = MclustGG$parameters$mean, pro = MclustGG$parameters$parameters$pr, Sigma = MclustGG$parameters$variance$sigma, ami = ami))
}



