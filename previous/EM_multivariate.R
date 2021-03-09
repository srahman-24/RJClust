### EM Algorithm for multivariate case dim = 4, C = 2
install.packages("MASS")
install.packages("mvtnorm")

library(MASS)
library(mvtnorm)
### Data generate
d               =   4
x               =   rbinom(100, 1, 0.50)
y               =   mvrnorm(length(x),  c(5,5,2,2), 1*diag(d))
y[which(x==0),] =   mvrnorm(length(which(x==0)), c(-1.5,0,0,0), 0.5*diag(d))
plot(density(y))

## Parameters: p,mu,sigma 
##Initial values:  
C = 2;               iter.max = 50
p = rep(1/C,C); 
mu = rbind(c(0,0,0,0),c(-1,1,-1,1)); Sigma = list(); z = array(0, dim = c(length(x),C))
Sigma[[1]] = Sigma[[2]] = diag(d)
##Estep: Generating Latent Variables.

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
  
  print(mu)
  #print(Sigma)
}


weights_multi     = function(x, p, mu, Sigma, ...)
{
    ww        = NULL 
  for(kk in 1:C)
    ww        =  c(ww, p[kk]*dmvnorm(x, mean = mu[kk,], Sigma[[kk]], log = FALSE))
  return(ww/sum(ww))
}



#Result: Compare
which(x==1) 
which(z[,1]==0)
