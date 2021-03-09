### EM Algorithm for univariate case 

### Data generate 
x =  rbinom(100,1,0.75)
y  = rnorm(length(x),3,1)
y[which(x==0)] = rnorm(length(which(x==0)),-3,0.5)
plot(density(y))

## Parameters: p,mu,sigma 
##Initial values:  
C = 2; iter.max = 50
p = rep(0.5,C); mu = c(-2,2); sigma = rep(1,C); z = array(0, dim = c(length(x),C))

##Estep: Generating latent Variables.
for(ss in 1:iter.max)
{

for(ii in 1:length(x))
{
  #print(weights(y[ii], p, mu, sigma))
  z[ii,] = t(rmultinom(1, (C-1), prob = weights(y[ii], p, mu, sigma)))
}

## Mstep Estimation:
for(jj in 1:C)
{
  y1        = y[which(z[,jj]==1)]
  mu[jj]    = mean(y1)
  sigma[jj] = sd(y1)
  p[jj]     = length(y1)/length(y)
}

  print(mu)
  print(sigma)
  
}

weights     = function(x, p, mu, sigma)
{
  ww = NULL 
  for(kk in 1:C)
  ww        =  c(ww, p[kk]*dnorm(x, mean = mu[kk], sigma[kk]))
  return(ww/sum(ww))
}

#Result: Compare
length(intersect(which(x==1),which(z[,1]==0)))/length(which(x==1)) 


Loglik  = function()
{
   
}


bic   =  function()
{
  
}  


