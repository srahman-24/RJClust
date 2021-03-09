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
    X[(ii+5), (jj+2)] =   rnorm(1, 1.5, sigma)
  }
}
n1 = n2 = 5

##### Initialize Cluster ID ##########
clusterID     = matrix(-2, nrow = n , ncol = 1) # n vector containing cluster id of each subject 
clusterNo     = matrix(0, nrow = C , ncol = 1)  # gives number of subjects in each clusters 
clusterMember = matrix(0, nrow = C, ncol = n)   # gives the rows are clusters, columns subjects   
maxClusterID  = matrix(0, nrow = n , ncol = 1)  # n vector containing the best cluster configuration


####  standardizing the columns of YY 
for (ii in 1:p)
{
       X[,ii] = (X[,ii] - mean(X[,ii]))/sd(X[,ii])
}
####  Form of G  and Z ## 
zz  = X*X
GG  = X%*%t(X)/p

######################### Initialize Gbar and Gbar_off (empirical diagonal and off-diagonal means)########################
Gbar         =   sum(diag(GG))/n
Gbar_off     =   sum(GG - diag(diag(GG)))/(n*(n-1))



C_start  = 2 
mu_diag  = array(0,C_start)
mu       = matrix(0,nrow = C_start, ncol = C_start)


####################### Likelihood function ############################ 



###################### Prior on mu with penalty ######################### 




###################### Prior on sigma inverse-gamma ##########################




################### Prior on 








