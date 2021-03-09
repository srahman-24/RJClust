install.packages("rjmcmc")
install.packages("mclust")
install.packages("Rmixmod")
install.packages("EMCluster")

rm(list = ls()) 


library(lpSolve)
library(MASS)
library(Matrix)
library(Rcpp)
library(EMCluster)
library(madness)
library(rjmcmc)
library(mclust)
library(Rmixmod)

################## Simple Example ###############

#for(ii in 1:seed)
#{  
set.seed(20183)
######### True model (G Matrix)  ################
C  = 4
#n  = rep(20, C)   # Equal cluster size are simulated 
n  = c(20,20,200,200)
##################################################################
mu         = matrix(0,nrow = C, ncol = C)
sigma      = matrix(1,nrow = C, ncol = C)
for(ii in 1:(C-1))
{
  mu[ii,ii]     = 0.1   #homo-off diagonals
  sigma[ii,ii]  = 0.05
  for(jj in (ii+1):C)
  {
    mu[ii,jj]    = -0.1 #hetero-off diagonals
    sigma[ii,jj] =  0.05   
  }
}
mu[C,C] = 0.1; sigma[C,C] = 0.07

mu0  = sigma0 = rep(0,C)


########### Be careful about this one, it changes with C ############
mu0     = c(1, 1, 1, 1)       #diagonals 
sigma0  = c(0.1,0.1, 0.1,0.1) 



##################### Simulation of G matrix ###################
gg = matrix(0, nrow = sum(n), ncol = sum(n))
#diagonals
diag(gg)[1:n[1]] = rnorm(n[1], mu0[1], sigma0[1])
for(ii in 2: C)
{
  diag(gg)[(1+sum(n[1:(ii-1)])):sum(n[1:ii])] = rnorm(n[ii], mu0[ii], sigma0[ii]) 
}

#off diagonals


#homo Cluster 1
for(ii in 1:(n[1]-1))
{
  for(jj in (ii+1):n[1])
     gg[ii,jj] = rnorm(1, mu[1,1], sigma[1,1])
}
#hetero Cluster 2
for(kk in 2:C)
{
       gg[(1:n[1]),(sum(n[1:(kk-1)])+1):sum(n[1:kk])] = rnorm(1, mu[1,kk], sigma[1,kk])
}



for(ll in 2:C)
{
  for(kk in ll:C)
  {
   if(ll == kk)    #homo-off diagonals 
    {  
     for(ii in (sum(n[1:ll-1]) + 1): (sum(n[1:ll]) - 1))
      {   
        for(jj in (ii+1): sum(n[1:ll]))
        {
              print(c(ii,jj,ll,kk));
              gg[ii,jj] = rnorm(1, mu[ll,ll], sigma[ll,ll])
        }      
      }
    }
  else       #hetero-off diagonals 
  {
    print(c(ll,kk));
    gg[(sum(n[1:(ll-1)]) + 1): sum(n[1:ll]), (sum(n[1:(kk-1)]) + 1):  sum(n[1:kk])] = rnorm(n[ll]*n[kk], mu[ll,kk], sigma[ll,kk])
  }
    
  }
  
}



# all the diagonals are same
#g[1:n1, (n1+1):(n1+n2)] = rnorm(n1*n2, mu12, sigma12)
#g[1:n1, (n1+n2+1):(n1+n2+n3)] = rnorm(n1*n3, mu13, sigma13)
#g[(n1+1):(n1+n2), (n1+n2+1):(n1+n2+n3)] = rnorm(n2*n3, mu23, sigma23)
#for(ii in 1:(n-1))
#{
#  for(jj in (ii + 1):n1)
#  {
#    g[ii,jj]              = rnorm(1, mu11, sigma11)
#    g[n1+ii,n1+jj]        = rnorm(1, mu22, sigma22)
#    #g[n1+n2+ii, n1+n2+jj] = rnorm(1, mu33, sigma33)
#  }
#}
################ For symmetricity #####################
for(ii in 2:sum(n))
{
  for(jj in 1:(ii-1))
    gg[ii,jj] = gg[jj,ii]
}
##################################################
#as.vector(gg)
 

   


#mclustgg = mclustBIC(gg)
Mclustgg = Mclust(gg)
#summary(mclustgg)
summary(Mclustgg)
#matrix(Mclustgg$class, nrow = C, ncol = n[1], byrow = T)
#duplicated(which(Mclustgg$class== 1),1:n[1])
#Mclustgg$class


gg_wodiag        = gg - diag(diag(gg))
Mclustgg_wodiag  = Mclust(gg_wodiag)
summary(Mclustgg_wodiag)

gg_wdiag        = cbind(gg_wodiag, diag(gg))
Mclustgg_wdiag  = Mclust(gg_wdiag)
summary(Mclustgg_wdiag)

gg_new          = cbind(gg_wodiag + diag(diag(colMeans(gg_wodiag))), diag(gg))
Mclustgg_new  = Mclust(gg_new)
summary(Mclustgg_new)
BIC = mclustBIC(gg_new)
plot(BIC)


# mclustgg = mclustBIC(as.vector(gg))
# Mclustgg = Mclust(as.vector(gg))
# summary(mclustgg)
# summary(Mclustgg)
# matrix(Mclustgg$class, nrow = sum(n), ncol = sum(n))
# duplicated(which(Mclustgg$class== 1),1:n[1])
# Mclustgg$class







