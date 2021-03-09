##EM Algorithm for G-Matrix Cluster Estimation 
library(MASS)
library(mvtnorm)
library(Mclust)


### Data generate
d               =   4 ; N = 100
x               =   rbinom(N, 1, 0.50)
y               =   mvrnorm(length(x),  c(5,5,2,2), 1*diag(d))
y[which(x==0),] =   mvrnorm(length(which(x==0)), c(-1.5,0,0,0), 0.5*diag(d))
plot(density(y))

## Parameters: p,mu,sigma 
##Initial values: Mclust 
C = 2;               iter.max = 50
p = rep(1/C,C); 
mu = rbind(c(0,0,0,0),c(-1,1,-1,1)); Sigma = list(); 
Sigma[[1]] = Sigma[[2]] = diag(d)


##EM function

emstep = function(y, C, mu, Sigma, p, iter.max, N)
{
  z = array(0, dim = c(N,C))  
  loglik_new  = sum(dmvnorm(y, mu[1,], Sigma[[1]], log = TRUE))  
  for(ss in 1:iter.max)
  {
    #Estep:
    for(ii in 1:length(x))
    {
      #print(weights(y[ii], p, mu, sigma))
      z[ii,] = t(rmultinom(1, (C-1), prob = weights_multi(y[ii,], p, mu, Sigma)))
    }
    
    ## Mstep Estimation:
    for(jj in 1:C)
    {
      y1           = y[which(z[,jj]==1),]
      mu[jj,]      = colMeans(y1)
      Sigma[[jj]]  = var(y1)
      p[jj]        = nrow(y1)/nrow(y)
    }
    
    loglik.diff = abs(loglik_new - loglik(y, z, mu, Sigma)) 
    if(loglik.diff< 1e-6) 
    {break}
    else  
    {loglik_new  = loglik(y, z, mu, Sigma)}
    #print(Sigma)
  } 
  return(list(mu = mu, Sigma = Sigma, prob = p, z = z))
}


weights_multi     = function(x, p, mu, Sigma, ...)
{
  ww        = NULL 
  for(kk in 1:C)
    ww        =  c(ww, p[kk]*dmvnorm(x, mean = mu[kk,], Sigma[[kk]], log = FALSE))
  return(ww/sum(ww))
}

loglik = function(y, z, mu, Sigma)
{
  loglik = 0 
  for(kk in 1:C)
  {
    loglik = loglik + sum(dmvnorm(y[which(z[,kk]==1),], mean= mu[kk,], Sigma[[kk]],log = TRUE))
  }
  return(loglik)
}


bic = function(C,y, mu, Sigma, N, prob)
{
  
  mu    = emstep(y, C, mu, Sigma, prob)$mu
  Sigma = emstep(y, C, mu, Sigma, prob)$Sigma
  z     = emstep(y, C, mu, Sigma, prob)$z
  return(2 * loglik(y, z, mu, Sigma) - nparams(C,d) * log(N))
  
}


nparams = function(C,d)
{
  nparams = C * d                  # mu
  nparams = nparams + (C - 1)      # prob
  nparams = nparams + C*d*(d+1)/2  # Sigma
  return(nparams)
}


for(C in 1:C_max)
{ 
  #Initialization   
  mu     =  Mclust()$mu                    
  Sigma  =  Mclust()$Sigma
  prob   =  Mclust()$prob 
  bic_eval = c(bic_eval, bic(C, y, mu, Sigma, N, prob))
  
}




#Result: Compare
which(x==1) 
which(z[,1]==0)
