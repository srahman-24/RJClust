##EM Algorithm for G-Matrix Cluster Estimation 
library(MASS)
library(mvtnorm)
library(Rcpp)


#source("Gcov.R")
Rcpp::sourceCpp('Gcov.cpp')

##EM function

weights_multi     = function(x, prob, mu, Sigma,C, ...)
{
    ww        = NULL 
  for(kk in 1:C)
    ww        = c(ww, prob[kk]*dmvnorm(x, mean = mu[kk,], Sigma[,,kk], log = FALSE))
  return(ww/sum(ww))
}


emstep = function(y, C, mu, Sigma, prob, iter.max, N, ...)
{
  
  z = array(0, dim = c(N,C))  
  loglik_new  = sum(dmvnorm(y, mu[1,], Sigma[,,1], log = TRUE))  # multivariate normal log likelihood with 1 cluster.
  for(ss in 1:iter.max)
  {
    #Estep:
    for(ii in 1:N)
    {
      ## print(weights(y[ii], p, mu, sigma))
      z[ii,] = t(rmultinom(1, 1, weights_multi(y[ii,], prob, mu, Sigma,C)))
      #z1[ii,]= which[max(weights_multi(y[ii,], prob, mu, Sigma,C))]
    }
    
    #GCOV          =  Gcov(y,C,z,N)
    #print(GCOV$gamma[[2]][1])
    GCOV           =  GcovCPP(y, C, z, N)
    
    ## Mstep Estimation:
    for(jj in 1:C)
    {
      if(sum(z[,jj]<=1))
      {break}
      else
      {
        y1           = y[which(z[,jj]==1), ]
        mu[jj,]      = colMeans(y1)
        Sigma[,,jj]  = GCOV                     #Calculated by G_cov function. 
        prob[jj]     = nrow(y1)/nrow(y)
      }
    }
    
    loglik.diff     = abs(loglik_new - loglik.G(y, prob, mu, Sigma,C,N)) 
    if(loglik.diff< 1e-6) 
    {break}
    else  
    {
      loglik_new  = loglik.G(y, prob, mu, Sigma, C,N)
    }
    #print(Sigma)
  } 
  return(list(mu = mu, Sigma = Sigma, prob = prob, z = z))
  
}


loglik.G = function(y, prob, mu, Sigma,C, N,...)
{ 
  ll = NULL 
  for(ii in 1:N)
  {
    ww = NULL 
    for(kk in 1:C)
    {ww        =  c(ww, prob[kk]*dmvnorm(y[ii,], mean = mu[kk,], Sigma[,,kk], log = FALSE))}
    
    ll = c(ll,log(sum(ww)))
  }
  return(sum(ll))
}



# loglik.G = function(y, z, mu, Sigma,C,..)
# {
#   loglik = 0 
#   for(kk in 1:C)
#   {
#     loglik = loglik + sum(dmvnorm(y[which(z[,kk]==1),], mean= mu[kk,], Sigma[,,kk],log = TRUE))
#   }
#   return(loglik)
# }


nparams.G = function(C,d)
{
  nparams = C * d                  # mu
  nparams = nparams + (C - 1)      # prob
  nparams = nparams + 4*C + 3*C*(C-1)/2 + C*(C-1)*(C-2)/3 + (C + C + C*(C-1)/2) # Sigma:  3 * C + 1 + (C + 1) * C * (C - 1)/3 
  #nparams = nparams + C^3
  return(nparams)
}



bic.G = function(C,y, mu, Sigma, N, prob, iter.max,...)
{
  mu    = emstep(y, C, mu, Sigma, prob, iter.max, N)$mu
  Sigma = emstep(y, C, mu, Sigma, prob, iter.max, N)$Sigma
  z     = emstep(y, C, mu, Sigma, prob, iter.max, N)$z
  if(any(colSums(z)<=1)) return(-1e6)
  return(list(bic.value = 2 * loglik.G(y, prob, mu, Sigma,C, N) - (nparams.G(C,nrow(y)) * log(N)), z = z, Sigma = Sigma, mu = mu))
  
}









