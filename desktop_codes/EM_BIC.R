### EM Algorithm with BIC
### EM Algorithm for multivariate case dim = 4, C = 2
install.packages("MASS")
install.packages("mvtnorm")

library(MASS)
library(mvtnorm)
library(mclust)

### Data generate
d                   =   4 ; N = 100
x                   =   t(rmultinom(N, 1, c(0.25,0.25,0.25,0.25)))
y                   =   mvrnorm(N,  c(0,0,0,0), 1*diag(d))
y[which(x[,1]==1),] =   mvrnorm(length(which(x[,1]==1)), c(-2.5,0,0,0),  0.5*diag(d))
y[which(x[,2]==1),] =   mvrnorm(length(which(x[,2]==1)), c(-2.5,0,0,-2.5), 0.5*diag(d))
y[which(x[,3]==1),] =   mvrnorm(length(which(x[,3]==1)), c(2.5,0,0,2.5),   0.5*diag(d))
y                   =   scale(y, center = T, scale =  T)
plot(density(y[,1]))
Mclust(y, modelNames = "EEI")$G

## Parameters: p,mu,sigma 


C_max = 8; bic_eval = NULL ; Z = list()

for(C in 2:C_max)
  { 
    init         =  Mclust(y, G = C) 
    iter.max     =  30
    mu           =  t(init$parameters$mean)
    Sigma        =  init$parameters$variance$sigma
    N            =  init$n
    prob         =  init$parameters$pr
    bic.val      =  bic.G(C, y , mu , Sigma, N, prob, iter.max)
    bic_eval     =  c(bic_eval, bic.val$bic.value)
    Z[[C-1]]     =  bic.val$z  
  }

Clust.no     = which(bic_eval == max(bic_eval)) + 1 ;
Class.Matrix = Z[[which(bic_eval == max(bic_eval))]]
Class        = colSums(Class.Matrix)
Class.member = Class.sigma = Class.mean = Class.pro = list(); 
for(kk in 1:Clust.no)
{
  Class.member[[kk]]  = which(Class.Matrix[,kk]==1)
  Class.mean[[kk]]    = colMeans(y[Class.member[[kk]], ]) 
  Class.sigma[[kk]]   = var(y[Class.member[[kk]], ])
  Class.pro[[kk]]     = length(Class.member[[kk]])/N  
}

plot( c(2:C_max), bic_eval, type = "l")
BIC = mclustBIC(y)
plot(BIC)
   


bic.G = function(C,y, mu, Sigma, N, prob, iter.max,...)
{
  
  mu    = emstep(y, C, mu, Sigma, prob, iter.max, N,d)$mu
  Sigma = emstep(y, C, mu, Sigma, prob, iter.max, N,d)$Sigma
  z     = emstep(y, C, mu, Sigma, prob, iter.max, N,d)$z
  if(any(colSums(z)<=1)) return(-1e16)
  return(list(bic.value = 2 * loglik.G(y, z, mu, Sigma) - nparams.G(C,d) * log(N), z = z))
  
}

##EM function

emstep = function(y, C, mu, Sigma, prob, iter.max, N, d, ...)
{
z = array(0, dim = c(N,C))  
loglik_new  = sum(dmvnorm(y, mu[1,], Sigma[,,1], log = TRUE))  
for(ss in 1:iter.max)
{
  #Estep:
  for(ii in 1:N)
  {
    #print(weights(y[ii], p, mu, sigma))
    z[ii,] = t(rmultinom(1, 1, weights_multi(y[ii,], prob, mu, Sigma)))
  }
  
  ## Mstep Estimation:
  for(jj in 1:C)
  {
    if(sum(z[,jj]<=1))
    {break}
    else
    {
    y1           = y[which(z[,jj]==1), ]
    mu[jj,]      = colMeans(y1)
    Sigma[,,jj]  = var(y1)
    p[jj]        = nrow(y1)/nrow(y)
    }
  }
  
  loglik.diff = abs(loglik_new - loglik(y, z, mu, Sigma)) 
  if(loglik.diff< 1e-6) 
  {break}
  else  
  {
    loglik_new  = loglik(y, z, mu, Sigma)
  }
  #print(Sigma)
} 
  return(list(mu = mu, Sigma = Sigma, prob = prob, z = z))
}


weights_multi     = function(x, prob, mu, Sigma, ...)
{
  ww        = NULL 
  for(kk in 1:C)
    ww        =  c(ww, prob[kk]*dmvnorm(x, mean = mu[kk,], Sigma[,,kk], log = FALSE))
  return(ww/sum(ww))
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
    nparams = C * d                        # mu
    nparams = nparams + (C - 1)            # prob
    nparams = nparams + C*d*(d+1)/2        # Sigma
  return(nparams)
}


for(C in 1:10)
{
  print(nparams.G(C,d)*log(N))
}


  
#Result: Compare
which(x==1) 
which(z[,1]==0)
