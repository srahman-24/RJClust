#Guo Simulation Setting 1 & 3  (Results without prescreening) 
#True number of groups = 2 

rm(list=ls(all=TRUE)) 
library(mclust)
library(reshape) 
library(ggplot2)
######### Bayesian Clustering Algorithm ############ 
#n    = c(20,20,50,50)        # Equal cluster size settings


n1    = 10000
n2    = 10000 
n3    = 10000
n4    = 10000 
#n     = c(n1,n2)
#n     = c(n1,n2,n3)
n     = c(n1,n2,n3,n4)         # Unequal Cluster size settings
p     = 500           # first 20 being informative and remaining ones are non-informative 
C     = 4                     
sigma = 1               # noise level (sigma level should be different for different clusters)
N     = sum(n)

set.seed(1044)
X     = matrix(rnorm(N*p,0,sigma), nrow = N, ncol = p, byrow = TRUE)
p0    = 20
p1    = 20

#Cluster 1: N(2.5, sigma)(1-p0), N(1.5, sigma)(p0-(p0+p1)) 

X[1:n[1],1:p0]                                               =   rnorm(n[1]*10, 2.5, sigma)
X[1:n[1],(1+p0):(p0+p1)]                                     =   rnorm(n[1]*10, 1.5, sigma)

#Cluster 2: N(0, sigma)(1-10), N(1.5, sigma) (11-20)
X[(n[1]+1):(n[1]+n[2]),1:p0]                                   =   rnorm(n[2]*10, 0, sigma)
X[(n[1]+1):(n[1]+n[2]),(1+p0):(p0+p1)]                         =   rnorm(n[2]*10, 1.5, sigma)

#Cluster 3: N(0, sigma)(1-10), N(-1.5,sigma)(11-20)
X[(n[1]+n[2]+1):(n[1]+n[2]+n[3]),1:p0]                         =   rnorm(n[3]*10,  0, sigma)
X[(n[1]+n[2]+1):(n[1]+n[2]+n[3]),(1+p0):(p0+p1)]               =   rnorm(n[3]*10, -1.5, sigma)
# 
#Cluster 4: N(-2.5,sigma)(1-10), N(-1.5, sigma)(11-20)
X[(n[1]+n[2]+n[3]+1):(n[1]+n[2]+n[3]+n[4]),1:p0]               =   rnorm(n[4]*10, -2.5, sigma)
X[(n[1]+n[2]+n[3]+1):(n[1]+n[2]+n[3]+n[4]),(1+p0):(p0+p1)]     =   rnorm(n[4]*10, -1.5, sigma)
# 

####  standardizing the columns of X 
#system.time({
####  Form of G  ##
#table(kmeans(X, centers = 2)$cluster, c(rep(1,n1), rep(2,n2)))
system.time({ 
Group_km     = kmeans(X, centers = 9)$cluster  
}) 
table(Group_km, c(rep(1,n1), rep(2,n2), rep(3,n3), rep(4,n4))) 
f_rez(Group_km, c(rep(1,n1), rep(2,n2), rep(3,n3), rep(4,n4))) 

###SparceKmeans  : Daniella Witten and Tibshirani 

# Z            = scale(X, center = T, scale =  T)
# 
# system.time({
# 
# GG           = Z%*%t(Z)/p
# 
# ############ Applying the Hierachical and EM Algorithm Package #############
# 
# gg_wodiag     =  GG - diag(diag(GG))
# GG_new        =  cbind(gg_wodiag + diag(colSums(gg_wodiag)/(N-1)), diag(GG))
# MclustGG      =  Mclust(GG_new, modelNames = "VVI")
# summary(MclustGG) 
# 
# 
# #table(MclustGG$classification, c(rep(1,n1), rep(2,n2)))
# #table(MclustGG$classification, c(rep(1,n1), rep(2,n2), rep(3,n3)))
# RJtable = table(MclustGG$classification, c(rep(1,n1), rep(2,n2), rep(3,n3), rep(4,n4))) })
# 
# RJtable
 
### Subsamples :  n1 = 25; n2 = 25; n3 = 25; n4 = 25
### Automated 

system.time({
  
cut            = 100      ######    ###### (user specified)
num.samp       = N/cut
ns             = sample(1:N)
n.samp         = list() 

    
d   = 0 
CC  = list()
for(ii in 1:num.samp)
{
  n.samp[[ii]] = ns[((ii-1)*cut + 1):(ii*cut)]
  c1           = length(intersect(n.samp[[ii]], 1:n[1]))    #counting each representation of clusters
  c2           = length(intersect(n.samp[[ii]], (n[1]+1):(n[1]+n[2])))
  c3           = length(intersect(n.samp[[ii]], (n[1]+n[2]+1):(n[1]+n[2]+n[3])))
  c4           = length(intersect(n.samp[[ii]], (n[1]+n[2]+n[3]+1):(n[1]+n[2]+n[3]+n[4])))
  print(c1)
  print(c2)
  print(c3)
  print(c4)
  #Z1           = Z[n.samp[[ii]], ]
  Z1            = scale(X[n.samp[[ii]], ], center = T, scale =  T)

GG1            =  Z1%*%t(Z1)/p
gg_wodiag      =  GG1 - diag(diag(GG1))
GG_new         =  cbind(gg_wodiag + diag(colSums(gg_wodiag)/(cut-1)), diag(GG1))
MclustGG1      =  Mclust(GG_new, modelNames = "VVI")


for (jj in 1: MclustGG1$G)
{
  CC[[d+jj]]     = n.samp[[ii]][which(MclustGG1$classification == jj)]
  
}
d = d + MclustGG1$G

}


 

######## Check if NN = N
NN = 0; nn = rep(0,d);  J.diag = J.homo.off = NULL 

Jmat = NULL
Jmat = matrix(0, nrow = d, ncol = d)

for(ii in 1:d)
{
  NN              = NN + length(CC[[ii]])
  Xnew[[ii]]      = X[CC[[ii]],]
  nn[ii]          = length(CC[[ii]])
  RR              = X[CC[[ii]],]%*%t(X[CC[[ii]],])/p
  J.diag          = c(J.diag, mean(diag(RR)))
  J.homo.off      = c(J.homo.off, (sum(RR)- sum(diag(RR)))/(nn[ii]^2 - nn[ii]))  
  if(ii == d){ break }
  for(jj in (ii+1):d)
  {
       Jmat[ii,jj]   =  mean(X[CC[[ii]],]%*%t(X[CC[[jj]],])/p)
  }
}


Jmat = Jmat + t(Jmat) 
for(ii in 1:d)  
{
  Jmat[ii,ii]        =  J.homo.off[ii]
}  
Jmat        =  cbind(Jmat, J.diag)



Clust.J = Mclust(Jmat, modelNames = "VVI")
Clust.J$G
Clust.J$classification
Cluster1 = Cluster2 = list()
Group = Group.kmeans = rep(0,N)
for(ii in 1: Clust.J$G)
{   
  Cluster1[[ii]]  = which(Clust.J$classification == ii)
  for(jj in 1: length(Cluster1[[ii]]))
  {
    Group[c(CC[[Cluster1[[ii]][[jj]]]])] = ii
  }
}


#table(Group, c(rep(1,n1), rep(2,n2)))
#table(Group, c(rep(1,n1), rep(2,n2), rep(3,n3)))
RJscale_T = table(Group, c(rep(1,n1), rep(2,n2), rep(3,n3), rep(4,n4)))

})

  
RJscale_T
f_rez(Group, c(rep(1,n1), rep(2,n2), rep(3,n3), rep(4,n4)))  



# RRnew = Xnew%*%t(Xnew)/p 
# 
# ############################ Goal is to make a "d by d"  matrix #######################
# RR_diag         = list()
# RR_off          = list()
# RR_off[[1]]     = list()
# RR_diag[[1]]    = RRnew[1:nn[[1]],1:nn[[1]]]
# dd_J            = (sum(RR_diag[[1]]) - sum(diag(RR_diag[[1]])))/(nn[1]^2 - nn[1])  
# ee_J            = mean(diag(RR_diag[[1]]))
# for(ii in 2:d)
# {
#   RR_diag[[ii]]        = RRnew[(1 + sum(nn[1:(ii-1)])):(sum(nn[1:ii])), (1+ sum(nn[1:(ii-1)])):(sum(nn[1:ii]))]
#   RR_off[[1]][[ii]]    = RRnew[1:nn[1], (1+sum(nn[1:(jj-1)])):sum(nn[1:ii])] 
#   dd_J                 = c(dd_J, (sum(RR_diag[[ii]]) - sum(diag(RR_diag[[ii]])))/(nn[ii]^2 - nn[ii]))
#   ee_J                 = c(ee_J, mean(diag(RR_diag[[ii]])))
# }
# 
# 
# for(ii in 2:(d-1))
# {
#   
#   RR_off[[ii]]   =  list()
#   for(jj in (ii+1):d)
#   {
#     RR_off[[ii]][[jj]]  = RRnew[(sum(nn[1:(ii-1)])+1):sum(nn[1:ii]), (sum(nn[1:ii])+1):sum(nn[1:jj])]  
#   }
#   
# }
# 
# ######### Constructing the d by (d+1) J-matrix  ###########
# JJ = diag(dd_J)
# for(ii in 1:(d-1))
# {
#   for(jj in (ii+1):d)
#   {
#     JJ[ii,jj] = mean(RR_off[[ii]][[jj]]) 
#   }
# }
# 
# JJ = JJ + t(JJ) - diag(dd_J)
# JJ = cbind(JJ, ee_J)
# 
# 



################## kmeans implementation on RJ #####################
# kmean.clust = 2
# KM = kmeans(JJ, centers = kmean.clust)
# for(ii in 1:kmean.clust)
# {
#   Cluster2[[ii]]  = which(KM$cluster == ii)
#   for(jj in 1: length(Cluster1[[ii]]))
#   {
#     Group.kmeans[c(CC[[Cluster2[[ii]][[jj]]]])] = ii
#   }
# }
# table(Group.kmeans, c(rep(1,n1), rep(2,n2), rep(3,n3), rep(4,n4)))
# 




