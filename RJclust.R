RJclust = function(Z1){

source("bicRJ.R")
source("AMI.R")  
  p           = ncol(Z1)
  n           = nrow(Z1)
  GG          = Z1%*%t(Z1)/p
  gg          = GG 
  gg_wodiag   = gg - diag(diag(gg))
  gg_wdiag    = cbind(gg_wodiag, diag(gg))
  GG_new      = cbind(gg_wodiag + diag(colSums(gg_wodiag)/(n-1)), diag(gg))

  C_max = 6; bic_eval = NULL ; Lat  = Sigma = mu = list(); iter.max =  100  
  d = n+1;  N  =  n
  
  for(C in 1:C_max)
  { 
    
    init                =  Mclust(GG_new, modelNames = "VVI", G = C)
    if(C == 1){bic_eval =  c(bic_eval, init$bic) }
    else
    {  
    mu           =  t(init$parameters$mean)
    Sigma        =  init$parameters$variance$sigma
    prob         =  init$parameters$pr
    bic.val      =  bic.G(C, GG_new, mu , Sigma, N, prob, iter.max)
    bic_eval     =  c(bic_eval, bic.val$bic.value)
    Lat[[C]]     =  bic.val$z
    }
  }

  n.clust      = which(bic_eval == max(bic_eval))  ;
  Clust.bic    = max(bic_eval)
  Class.Matrix = Lat[[which(bic_eval == max(bic_eval))]]
  Class        = colSums(Class.Matrix)
  Class.Sigma  = Gcov(GG_new, Clust.no, z = Class.Matrix, N)
  Class.member = list(); Class.sigma = Sigma = Class.pro = Class.mean = list(); 
  groupG = c(1:n)
  for(kk in 1:n.clust)
  {
    Class.member[[kk]]  = which(Class.Matrix[,kk]==1)
    groupG[Class.member[[kk]]] = kk
    Class.mean[[kk]]    = colMeans(GG_new[Class.member[[kk]], ]) 
    #Class.sigma[[kk]]   = var(y[Class.member[[kk]], ])
    Class.pro[[kk]]     = length(Class.member[[kk]])/N  
  }
  
  
  ami = f_rez(group, groupG)$ami
  
  return(c(bic = Clust.bic, bic.vector = bic_eval, n.clust = n.clust, groupG = groupG, Sigma = Class.Sigma, mean = Class.mean, pro = Class.pro, ami = ami))
  
  }
