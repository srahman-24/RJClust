RJclust = function(Z, group, C.max, iter.max){

source("bicRJ.R")
source("AMI.R")  
    p           =    ncol(Z)
    n           =    nrow(Z)      
    GG          =    Z%*%t(Z)/p
    gg          =    GG 
    gg_wodiag   =    gg - diag(diag(gg))
    gg_wdiag    =    cbind(gg_wodiag, diag(gg))
    GG_new      =    cbind(gg_wodiag + diag(colSums(gg_wodiag)/(n-1)), diag(gg))

  bic_eval = NULL ; Lat  = Sigma = mu = list();  
  d = (n+1);  N  =  n
  
  bic_eval = c(bic_eval, Mclust(GG_new, G = 1)$bic)
  
  for(C in 2:C.max)
  { 
    init         =  Mclust(GG_new, G = C)
    #init         =  Mclust(GG_new, modelNames = "VVI", G = C)
    mu           =  t(init$parameters$mean)
    Sigma        =  init$parameters$variance$sigma
    prob         =  init$parameters$pr
    bic.val      =  bic.G(C, GG_new, mu , Sigma, N, prob, iter.max)
    bic_eval     =  c(bic_eval, bic.val$bic.value)
    Lat[[C]]     =  bic.val$z
  }
  
  #plot(bic_eval)
  
  n.clust        =    which(bic_eval == max(bic_eval))  ;
  Clust.bic      =    max(bic_eval)
  Class.Matrix   =    Lat[[n.clust]]
  Class          =    colSums(Class.Matrix)
  Class.Sigma    =    Gcov(GG_new, n.clust, z = Class.Matrix, N)
  Class.member   =    list();     Class.sigma = Sigma = Class.pro = Class.mean = list(); 
  groupG         =    c(1:n)
  
  
  for(kk in 1:n.clust)
  {
    Class.member[[kk]]         = which(Class.Matrix[,kk]==1)
    groupG[Class.member[[kk]]] = kk
    Class.mean[[kk]]           = colMeans(GG_new[Class.member[[kk]], ]) 
    #Class.sigma[[kk]]         = var(y[Class.member[[kk]], ])
    Class.pro[[kk]]            = length(Class.member[[kk]])/N  
  }
  
  
  ami   =   f_rez(group, groupG)$ami
  
  
return(list(bic = Clust.bic, bic_eval = bic_eval, n.clust = n.clust, groupG = groupG, Sigma = Class.Sigma, mean = Class.mean, pro = Class.pro, ami = ami))
 
   
}


# Mclust(GG_new, modelNames = "VVI")$BIC 
# Class.Sigma[[1]][1:5, 1:5]
# det(Class.Sigma[[1]])
# table(group,groupG)


