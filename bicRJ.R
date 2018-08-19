##EM Algorithm for G-Matrix Cluster Estimation 
library(MASS)
library(mvtnorm)
library(mclust)



source("Gcov.R")

##EM function

weights_multi     = function(x, prob, mu, Sigma, ...)
{
  ww        = NULL 
  for(kk in 1:C)
    ww        =  c(ww, prob[kk]*dmvnorm(x, mean = mu[kk,], Sigma[,,kk], log = FALSE))
  return(ww/sum(ww))
}

emstep = function(y, C, mu, Sigma, prob, iter.max, N, ...)
{
  z = array(0, dim = c(N,C))  
  loglik_new  = sum(dmvnorm(GG_new, mu[1,], Sigma[,,1], log = TRUE))  
  for(ss in 1:iter.max)
  {
    #Estep:
    for(ii in 1:N)
    {
      ## print(weights(y[ii], p, mu, sigma))
      z[ii,] = t(rmultinom(1, 1, weights_multi(GG_new[ii,], prob, mu, Sigma)))
    }
    
    GCOV          =  Gcov(GG_new,C,z,N)
    
    ## Mstep Estimation:
    for(jj in 1:C)
    {
      if(sum(z[,jj]<=1))
      {break}
      else
      {
        y1           = y[which(z[,jj]==1), ]
        mu[jj,]      = colMeans(y1)
        Sigma[,,jj]  = GCOV[[jj]]                     #Calculated by G_cov function. 
        p[jj]        = nrow(y1)/nrow(y)
      }
    }
    
    loglik.diff     = abs(loglik_new - loglik.G(GG_new, z, mu, Sigma)) 
    if(loglik.diff< 1e-6) 
    {break}
    else  
    {
      loglik_new  = loglik(GG_new, z, mu, Sigma)
    }
    #print(Sigma)
  } 
  return(list(mu = mu, Sigma = Sigma, prob = prob, z = z))
}




loglik.G = function(y, z, mu, Sigma,...)
{
  loglik = 0 
  for(kk in 1:C)
  {
    loglik = loglik + sum(dmvnorm(y[which(z[,kk]==1),], mean= mu[kk,], Sigma[,,kk],log = TRUE))
  }
  return(loglik)
}


nparams.G = function(C,d)
{
  nparams = C * d                  # mu
  nparams = nparams + (C - 1)      # prob
  nparams = nparams + 4*C + 3*C*(C-1)/2 + C*(C-1)*(C-2)/3 # Sigma:  3 * C + 1 + (C + 1) * C * (C - 1)/3 
  #nparams = nparams + C^3
  return(nparams)
}



bic.G = function(C,y, mu, Sigma, N, prob, iter.max,...)
{
  
  mu    = emstep(y, C, mu, Sigma, prob, iter.max, N)$mu
  Sigma = emstep(y, C, mu, Sigma, prob, iter.max, N)$Sigma
  z     = emstep(y, C, mu, Sigma, prob, iter.max, N)$z
  if(any(colSums(z)<=1)) return(-1e6)
  return(list(bic.value = 2 * loglik.G(y, z, mu, Sigma) - nparams.G(C,ncol(GG)) * log(N), z = z, Sigma = Sigma, mu = mu))
  
}









