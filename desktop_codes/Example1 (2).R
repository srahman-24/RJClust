#Example 1 
#True number of groups = 2 

rm(list=ls(all=TRUE)) 
##### Bayesian Clustering Algorithm ######## 
n     = 10       # n > p settings
p     = 50       # first 4 being informative and remaining ones are non-informative 
C     = n        # initializing every individual as their own clusters 
sigma = 1  
X    = matrix(0,nrow = n, ncol = p, byrow = TRUE)
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

##### Initialize Cluster ID ##########
clusterID     = matrix(-2, nrow = n , ncol = 1) # n vector containing cluster id of each subject 
clusterNo     = matrix(0, nrow = C , ncol = 1)  # gives number of subjects in each clusters 
clusterMember = matrix(0, nrow = C, ncol = n)   # gives the rows are clusters, columns subjects   
maxClusterID  = matrix(0, nrow = n , ncol = 1)  # n vector containing the best cluster configuration


####  standardizing the columns of YY 
# for (ii in 1:p)
# {
#        X[,ii] = (X[,ii] - mean(X[,ii]))/sd(X[,ii])
# }
####  Form of G  and Z ## 
zz  = X*X
GG  = X%*%t(X)/p

######################### Initialize Gbar and Gbar_off (empirical diagonal and off-diagonal means)########################
Gbar         =   sum(diag(GG))/n
Gbar_off     =   sum(GG - diag(diag(GG)))/(n*(n-1))



C_start       = 2     # number of clusters. 
clusterID     = c(rep(1,5), rep(2,5))                  #cluster id 
mu_diag       = array(0, C_start)
sigma_diag    = array(1, C_start)
mu            = matrix(0, nrow = C_start, ncol = C_start)
sigma_start   = matrix(1, nrow = C_start, ncol = C_start)



####################### True G mu and G sigma ##################################
truemu_diag       = rep(0,2)
truesigma_diag    = rep(0,2)
true_mu           = matrix(0, nrow = C_start, ncol = C_start)
true_sigma        = matrix(1, nrow = C_start, ncol = C_start)
truemu_diag[1]    = sum(rep(2.5^2,2), rep(1.5^2,2), rep(0^2,(50-4)))/p+sum(rep(sigma^2,50))/p                                                                                         
truemu_diag[2]    = sum(rep(0^2,2), rep(1.5^2,2), rep(0^2,(50-4)))/p+sum(rep(sigma^2,50))/p
truesigma_diag[1] = 2*sum(rep(1 + (2*2.5^2),2),rep(1 + (2*1.5^2),2), rep(1 + (2*0^2),50-4))/p^2 
truesigma_diag[2] = 2*sum(rep(1 + (2*0^2),2),rep(1 + (2*1.5^2),2), rep(1 + (2*0^2),50-4))/p^2 
true_mu[1,1]       = sum(rep(2.5^2,2),rep(1.5^2,2),rep(0^2,(50-4)))/p 
true_mu[2,2]       = sum(rep(0^2,2),rep(1.5^2,2),rep(0^2,(50-4)))/p 
true_sigma[1,1]    = truemu_diag[1]/2 
true_sigma[2,2]    = truemu_diag[2]/2
true_mu[1,2]       = sum(rep(2.5*0,2),rep(1.5*1.5,2),rep(0*0,(50-4)))/p

##################################################### Likelihood function ############################################################# 
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


#############################  True_Likelihood  ####################################
uu = 1 
for(ii in 1:n)
{
  uu = uu*dnorm(diag(GG)[ii],truemu_diag[clusterID[ii]],truesigma_diag[clusterID[ii]]) 
}
ww = 1
for(ii in 1:(n-1))
{
  for(jj in (ii+1):n)
  {
    if(clusterID[ii] == clusterID[jj]) 
      ww = ww*dnorm(GG[ii,jj],true_mu[clusterID[ii],clusterID[ii]],true_sigma[clusterID[ii],clusterID[ii]]) 
    else
      ww = ww*dnorm(GG[ii,jj],true_mu[clusterID[ii],clusterID[jj]],true_sigma[clusterID[ii], clusterID[jj]]) 
  }
}
Likelihood = uu*ww


###################### Prior on mu with penalty ######################### 
mu_number = C(C+3)/2




###################### Prior on sigma inverse-gamma ##########################




################### Prior on 








