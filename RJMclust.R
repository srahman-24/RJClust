### author: Shahina Rahman 


RJMclust = function(Z){

source("AMI.R")  
p = ncol(Z1)
n = nrow(Z1)

GG  = Z1%*%t(Z1)/p
gg = GG 
gg_wodiag       = gg - diag(diag(gg))
gg_wdiag        = cbind(gg_wodiag, diag(gg))
GG_new          = cbind(gg_wodiag + diag(colSums(gg_wodiag)/(n-1)), diag(gg))

############ Applying the Hierachical and EM Algorithm Raftery's Package #############
MclustGG = Mclust(GG_new, modelNames = "VVI")  # On GG 
groupG   = MclustGG$class


ami = f_rez(group, groupG)$ami


return(list(groupG = groupG, n.clust = MclustGG$G, bic = MclustGG$bic, mean = MclustGG$parameters$mean, pro = MclustGG$parameters$parameters$pr, Sigma = MclustGG$parameters$variance$sigma, ami = ami))
}
