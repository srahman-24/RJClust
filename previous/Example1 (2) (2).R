#Example 1 
#True number of groups = 2 
library(MCMCpack)

rm(list=ls(all=TRUE)) 
############################################################## Bayesian Clustering Algorithm ######################################################## 
n     = 10       # n > p settings
p     = 50       # first 4 being informative and remaining ones are non-informative 
C     = 2        # initializing every individual as their own clusters 
sigma = 1  
X     = matrix(0, nrow = n, ncol = p, byrow = TRUE)
for(ii in 1:n)
{
  for(jj in 1:p)
            X[ii,jj] = rnorm(1,0,1)
}
#Cluster 1 N(2.5, sigma^2)(1-10), N(1.5, sigma^2) 
#Cluster 2 N(0, sigma^2)(1-10), N(1.5, sigma^2) 
for(ii in 1:(n/2))
{
  for(jj in 1:2)
  {
    X[ii,jj]          =   rnorm(1, 2.5, sigma)
    X[ii, (jj+2)]     =   rnorm(1, 1.5, sigma)
    X[(ii+5), (jj+2)] =   rnorm(1, 2.5, sigma)
  }
}
n1 = n2 = 5

################################################################# Initialize Cluster ID #############################################################
clusterID            = matrix(-2, nrow = n , ncol = 1)   # n vector containing cluster id of each subject 
clusterNo            = matrix( 0, nrow = C , ncol = 1)   # gives number of subjects in each clusters 
clusterMember        = matrix( 0, nrow = C , ncol = n)   # gives the rows are clusters, columns subjects   
maxClusterID         = matrix( 0, nrow = n , ncol = 1)   # n vector containing the best cluster configuration


####  standardizing the columns of YY 
# for (ii in 1:p)
# {
#        X[,ii] = (X[,ii] - mean(X[,ii]))/sd(X[,ii])
# }
####  Form of G  and Z ## 
zz  = X*X
GG  = X%*%t(X)/p

################################################################# Initialize Gbar and Gbar_off (empirical diagonal and off-diagonal means)########################
Gbar                =   sum(diag(GG))/n
Gbar_off            =   sum(GG - diag(diag(GG)))/(n*(n-1))



C_start       = 2     # number of clusters. 
clusterID     = c(rep(1,3), rep(2,7))                  #cluster id 
mu_diag       = array(0, C_start)
sigma_diag    = array(1, C_start)
mu            = matrix(0, nrow = C_start, ncol = C_start)
sigma_start   = matrix(1, nrow = C_start, ncol = C_start)



################################################################# True G mu and G sigma #########################################################################
truemu_diag        = rep(0,2)
truesigma_diag     = rep(0,2)
true_mu            = matrix(0, nrow = C_start, ncol = C_start)
true_sigma         = matrix(1, nrow = C_start, ncol = C_start)
truemu_diag[1]     = sum(rep(2.5^2,2), rep(1.5^2,2), rep(0^2,(50-4)))/p+sum(rep(sigma^2,50))/p                                                                                         
truemu_diag[2]     = sum(rep(0^2,2), rep(1.5^2,2), rep(0^2,(50-4)))/p+sum(rep(sigma^2,50))/p
truesigma_diag[1]  = sqrt(2*sum(rep(1 + (2*2.5^2),2),rep(1 + (2*1.5^2),2), rep(1 + (2*0^2),(50-4)))/p^2) 
truesigma_diag[2]  = sqrt(2*sum(rep(1 + (2*0^2),2),rep(1 + (2*1.5^2),2), rep(1 + (2*0^2),(50-4)))/p^2) 
true_mu[1,1]       = sum(rep(2.5^2,2),rep(1.5^2,2),rep(0^2,(50-4)))/p 
true_mu[2,2]       = sum(rep(0^2,2),rep(1.5^2,2),rep(0^2,(50-4)))/p 
true_sigma[1,1]    = sqrt(truemu_diag[1]/2) 
true_sigma[2,2]    = sqrt(truemu_diag[2]/2)
true_mu[1,2]       = sum(rep(2.5*0,2),rep(1.5*1.5,2),rep(0*0,(50-4)))/p
true_sigma[1,2]    = sqrt(sum(rep(1 + (2.5^2) + (0^2),2),rep(1 + (1.5^2) + (1.5^2),2), rep(1 + (0^2) + (0^2),(50-4)))/p^2)
true_clusterID     = c(rep(1,5),rep(2,5))

################################################################### Likelihood function ######################################################################### 
###### When we consider the cluster ID ###############
uu = 1 
for(ii in 1:n)
{
          uu = uu*dnorm(diag(GG)[ii],mu_diag[clusterID[ii]],sigma_diag[clusterID[ii]]) 
}
ww = 1
for(ii in 1:(n-1))
{
  for(jj in (ii+1):n)
  {
        if(clusterID[ii] == clusterID[jj])
           ww = ww*dnorm(GG[ii,jj], mu[clusterID[ii],clusterID[ii]], sigma_start[clusterID[ii], clusterID[ii]]) 
        else
           ww = ww*dnorm(GG[ii,jj], mu[clusterID[ii],clusterID[jj]], sigma_start[clusterID[ii], clusterID[jj]]) 
  }
}
Likelihood = uu*ww


####################################################################  True_Likelihood  ###########################################################################
uu = 1 
for(ii in 1:n)
{
  uu = uu*dnorm(diag(GG)[ii],truemu_diag[true_clusterID[ii]],truesigma_diag[true_clusterID[ii]]) 
}
ww   = 1
for(ii in 1:(n-1))
{
  for(jj in (ii+1):n)
  {
    if(clusterID[ii] == clusterID[jj]) 
                   ww = ww*dnorm(GG[ii,jj],true_mu[true_clusterID[ii],true_clusterID[ii]],true_sigma[true_clusterID[ii],true_clusterID[ii]]) 
    else
                   ww = ww*dnorm(GG[ii,jj],true_mu[true_clusterID[ii],true_clusterID[jj]],true_sigma[true_clusterID[ii],true_clusterID[jj]]) 
  }
}
Likelihood_true       =  uu*ww


######################## Prior on mu with penalty ######################## 
mu_number            =  C*(C+3)/2
gamma_mat            =  as.matrix(rdirichlet(n, rep(0.5,C)))
mu_diag              =  rep(0,C)
mu_offdiag           =  matrix(rep(0,C^2),nrow = C, ncol = C)
sigma                =  matrix(rep(1,C^2),nrow = C, ncol = C)


mu_penalty    =  function(mu_diag, mu_offdiag, delta)
{
      pen  = NULL 
for(ii in 1:(C-1))
{
  for(jj in (ii+1):C)
  {
    pen    =  rbind(pen, 1*(abs(mu_diag[ii] - mu_diag[jj])>delta))              
    pen    =  rbind(pen, 1*(abs(mu_offdiag[ii,ii] - mu_offdiag[jj,jj])>delta))   
    pen    =  rbind(pen, 1*(abs(mu_offdiag[ii,ii] - mu_offdiag[ii,jj])>delta))  
    pen    =  rbind(pen, 1*(abs(mu_offdiag[jj,jj] - mu_offdiag[ii,jj])>delta))   
  for(kk in 1:C)
  {
    if((kk!=ii) && (kk!=jj))  
    { 
      pen  = rbind(pen, 1*(abs(mu_offdiag[ii,kk] - mu_offdiag[jj,kk]) > delta)) }                            
    }
  }
}
  if(sum(pen)>=1) 
    return(1)
  else     
    return(0) 
}







############################## Likelihood when the cluster ID are continuous  ######################### 


Likelihood = function(GG, mu_diag, mu_offdiag, sigma, gamma_mat)
{
LL  = 0
### Likelihood for diagonal
for(ii in 1:n)
  for(kk in 1:C)
  {
    LL = LL + gamma_mat[ii,kk]*dnorm(GG[ii,ii], mu_diag[kk], sigma[kk,kk],log=T)
    print(dnorm(GG[ii,ii], mu_diag[kk], sigma[kk,kk],log=T))
  }
### Likelihood for off-diagonals
for(ii in 1:(n-1)){
  for(jj in (ii+1):n) {
      for(kk1 in 1:C) {
        for(kk2 in 1:C) {
          if(kk1 == kk2)
          LL = LL + gamma_mat[ii,kk1]*gamma_mat[jj,kk2]*dnorm(GG[ii,jj], mu_offdiag[kk1,kk1], sigma[kk1,kk1]/sqrt(2),log = T)
          else
          LL = LL + gamma_mat[ii,kk1]*gamma_mat[jj,kk2]*dnorm(GG[ii,jj], mu_offdiag[kk1,kk2], sigma[kk1,kk2],log = T)  
                        }
                      }
                      }
}
return(LL)
}



################### Prior on cluster (Booth prior) : prior on Cluster probability ##################################


gamma_gibbs          =  gamma_mat
mu_diag_gibbs        =  rep(0,C)
mu_offdiag_gibbs     =  matrix(rep(0,C^2),nrow = C, ncol = C)
sigma_gibbs          =  matrix(rep(1,C^2),nrow = C, ncol = C)
################### MCMC steps #######################################

#### Gibbs step for mu 



for(kk1 in 1:C)
  { 
  aa = 0
  mu_diag_gibbs[kk1]   =   rnorm(1, sum(diag(GG)*(gamma[,kk1])^2)/sum(gamma[,kk1])^2, sqrt(sigma[kk1,kk1]^2/sum(gamma[,kk1])^2))
  aa                   =   sum((diag(GG) - mu_diag[kk1])^2*(gamma[,kk1])^2)/2 
  for(kk2 in 1:C)
  { 
    cc = dd = 0
    for(ii in 1:(n-1))
    { 
      for(jj in (ii+1):n)
      {
                   cc  =   cc +   GG[ii,jj]*gamma_gibbs[ii,kk1]*gamma_gibbs[jj,kk2]
                   dd  =   dd +   gamma_gibbs[ii,kk1]*gamma_gibbs[jj,kk2]
      }          
    }
                mu_offdiag_gibbs[kk1,kk2]  =   rnorm(1,  cc, sigma[kk1,kk2]/dd )
    }
  }





#Check the condition  
delta = 0.01
if(mu_penalty(mu_diag_gibbs, mu_offdiag_gibbs, delta) == 1) 
{ 
  print("yes")
  mu_diag    = mu_diag_gibbs;
  mu_offdiag = mu_offdiag_gibbs;
}else   print("no")

#### Gibbs step for sigma #####
MM1  = matrix(0, nrow = n, ncol = n)
MM2  = matrix(1, nrow = n, ncol = n)
MM2[lower.tri(MM2)] = 0
##### for the diagonals ######
for(kk1 in 1:C)
{
    for(ii in 1:n)
    {
    for(jj in ii:n)
    {
      if(jj==ii)
              MM1[ii,jj] = 0.5*(GG[ii,jj] - mu_diag[kk1])^2
      else
              MM1[ii,jj] = (GG[ii,jj] - mu_offdiag[kk1,kk1])^2
    }
    }
  sigma[kk1,kk1]   = sqrt(1/rgamma(1, shape = 0.5*t(gamma_mat[,kk1])%*%MM2%*%gamma_mat[,kk1] + 8, rate = t(gamma_mat[,kk1])%*%MM1%*%gamma_mat[,kk1] + (2.5/p)))
  for(kk2 in 1:C)     ##### for the offdiagonals of sigma ######                 
  {
    MM2  = matrix(1, nrow = n, ncol = n)
    MM2[lower.tri(MM2)] = 0
    MM1  = matrix(0, nrow = n, ncol = n)
    for(ii in 1:(n-1))
    {
      for(jj in (ii+1):n)
        {
         MM1[ii,jj] = 0.5*(GG[ii,jj] - mu_offdiag[kk1,kk2])^2
        }
    }
    sigma[kk1,kk2] = sqrt(1/rgamma(1, shape = 0.5*t(gamma_mat[,kk1])%*%MM2%*%gamma_mat[,kk2] + 8, rate = 0.5*(t(gamma_mat[,kk1])%*%MM1%*%gamma_mat[,kk2] + (2.5/p)))) 
  }
}





#### Gibbs step for gamma probability 





