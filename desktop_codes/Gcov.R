### function input: G, C, z, N
### parameters estimates: ss.diag, ss.off, sss, mean.diag, mean.off, ss.bound 
### function output : List of Cov Matrix 1,....,C.   
Gcov = function(GG,C,z,N,...)  
{
  
  ss.diag   = rep(1,C) 
  mean.diag = rep(1,C)
  ss.off    = ss.bound = nn = matrix(0, nrow = C, ncol = C)
  mean.off  = matrix(0, nrow = C, ncol = C)
  sss       = nnn = array(0,dim = c(C, C, C))
  gamma     = list()
  Cov       = list()
  
  

  labels    = rep(1,N) 
  for(ii in 1:C)
  {
    gamma[[ii]]         =  which(z[,ii] == 1) 
    labels[gamma[[ii]]] =  ii
  }
  
###Diagonal Covariance Estimates:
###Total number of Diagonal Parameters: (C+1)  
for(ii in 1:C)
{ 
  for(jj in ii:C)
    {
       if(ii == jj)
      { 
      #homo-off diagonals
      GG_1                = GG[gamma[[ii]], gamma[[jj]]]
      nn[ii,jj]           = (length(GG_1) - nrow(GG_1))
      mean.diag[ii]       = mean(diag(GG_1))   #diagonal mean
      ss.diag[ii]         = var(diag(GG_1))    #diagonal variance
      M = mean.off[ii,jj] = (sum(GG_1) - sum(diag(GG_1)))/nn[ii,jj]  #offdiagonals mean
      ss.off[ii,jj]       = (sum((GG_1 - M)^2) - sum(diag((GG_1 - M)^2)))/nn[ii,jj] #offdiagonals variance  
      }
    else
    {
      #hetero-off diagonals 
      GG_2                   = GG[gamma[[ii]], gamma[[jj]]]
      nn[ii,jj]              = length(GG_2)
      M = mean.off[ii,jj]    = mean(GG_2)
      ss.off[ii,jj]          = sum((GG_2 - M)^2)/nn[ii,jj]
    }
  }
}

ss.off =  ss.off + t(ss.off) - diag(diag(ss.off))

####### ss.diag and ss.off are completely cluster specific and symmetric 
###Boundary Covariance Estimates:  non-symmetric
###Total number of Boundary parameters: C  

for(ii in 1:C)
{
  GG_1                = GG[gamma[[ii]], gamma[[ii]]]
  for(jj in 1:C)
  {
    if(ii == jj)
    {
      nn[ii,jj]           = (length(GG_1) - nrow(GG_1))
      S.homo = 0
      for(kk in 1:nrow(GG_1))
      {
        S.homo            =  S.homo + sum((GG_1[kk,] - mean.off[ii,jj])*(GG_1[kk,kk] - mean.diag[ii])) 
        - (GG_1[kk,kk] - mean.off[ii,jj])*(GG_1[kk,kk] - mean.diag[ii])
      }
      ss.bound[ii,jj]     =  S.homo/nn[ii,jj]                                 
    }
    else
    {
      GG_2               = GG[gamma[[ii]], gamma[[jj]]]
      GG_2               = GG_2 - mean.off[ii,jj]   
      nn[ii,jj]          = length(GG_2)
      S.hetero = 0
      for(kk in 1:nrow(GG_2))
      {
        S.hetero   =   S.hetero + sum(GG_2[kk,]*(GG_1[kk,kk] - mean.diag[ii])) 
      }
      ss.bound[ii,jj]     =   S.hetero/nn[ii,jj]    
    }
  }
}
#Offdiagonals Estimates: symmetric 
## 1<=ii <=jj <= kk <= C 
## Total number of parameters = C choose 3 + 2*C choose 2 + C choose 1. 
## Hardest part. 

for(ii in 1:C)
{
  for(jj in ii:C)
  { 
    if(jj == ii)
    {
      GG_1      =    GG[gamma[[ii]], gamma[[jj]]]
      GG_1      =    (GG_1 - mean.off[ii,jj]) + diag(diag(GG_1)-mean.diag[ii]) - diag(diag(GG_1)-mean.off[ii,jj])  
    }
    else 
    {
      GG_1      =   GG[gamma[[ii]], gamma[[jj]]]
      GG_1      =   GG_1 - mean.off[ii,jj]
    }  
    for(kk in jj:C)
    { 
      GG_2      =   GG[gamma[[ii]], gamma[[kk]]]
      GG_2      =   GG_2 - mean.off[ii,kk]
      if(ii == jj && jj == kk )
      {
        GG_2      =    GG[gamma[[ii]], gamma[[kk]]]
        GG_2      =    (GG_2 - mean.off[ii,kk]) + diag(diag(GG_2) - mean.diag[ii]) - diag(diag(GG_2) - mean.off[ii,kk])  
        S=0
        nnn[ii,jj,kk]   =   nrow(GG_1)*ncol(GG_1)*ncol(GG_2) - length(GG_1) - length(GG_2) + nrow(GG_1)
        for(r1 in 1:nrow(GG_1))
        {
          S = S + sum(kronecker(GG_1[r1,],GG_2[r1,])) - sum(GG_1*GG_2) - sum(GG_1[r1,r1]*GG_2[r1,]) + 
            + sum((GG_1[r1,r1])^2) 
        }
        sss[ii,jj,kk] = S/nnn[ii,jj,kk]
      }
      if(ii == jj && jj < kk)
      {
        S = 0
        nnn[ii,jj,kk]  =  nrow(GG_1)*ncol(GG_1)*ncol(GG_2) - length(GG_2)
        for(r1 in 1:nrow(GG_1))
        {
          S = S + sum(kronecker(GG_1[r1,],GG_2[r1,])) - sum(GG_1[r1,r1]*GG_2[r1,])
        }
        sss[ii,jj,kk] = S/nnn[ii,jj,kk]
      }
      if(ii < jj && jj == kk)
      {
        S=0
        nnn[ii,jj,kk] = nrow(GG_1)*ncol(GG_1)*ncol(GG_2) - length(GG_1)
        for(r1 in 1:nrow(GG_1))
        {
          S = S + sum(kronecker(GG_1[r1,],GG_2[r1,])) - sum(GG_1[r1,]*GG_2[r1,])
        }
        sss[ii,jj,kk] = S/nnn[ii,jj,kk]
      }
      if(ii < jj && jj < kk)
      {
        S=0
        nnn[ii,jj,kk] = nrow(GG_1)*ncol(GG_1)*ncol(GG_2)
        for(r1 in 1:nrow(GG_1))
        {
          S = S +  sum(kronecker(GG_1[r1,],GG_2[r1,])) 
        }
        sss[ii,jj,kk] = S/nnn[ii,jj,kk]
      }
      sss[ii,kk,jj] = sss[jj,kk,ii] = sss[kk,jj,ii] = sss[kk,ii,jj] = sss[jj,ii,kk] = sss[ii,jj,kk]
    }
  }
}



### Estimating all the cluster Covariance Matrix: 
for(cc in 1: C)
{
  Cov[[cc]]               = matrix(0, nrow = (N+1), ncol = (N+1))
  Cov[[cc]][(N+1),(N+1)]  = ss.diag[cc]                 #boundary
  for(ii in 1:N)
  {
    if(ii != gamma[[cc]][1]) 
    {
      for(jj in ii:N)
      {
        if(jj != gamma[[cc]][1])
        {
          if(ii == jj)                              #diagonals 
          {
            Cov[[cc]][ii,jj] = ss.off[cc,labels[ii]]
          } 
          else
            Cov[[cc]][ii,jj] = sss[cc,labels[ii],labels[jj]]  #offdiagonals 
        }
        Cov[[cc]][jj,ii] = Cov[[cc]][ii,jj]
      }
      Cov[[cc]][(N+1), ii] = Cov[[cc]][ii, (N+1)]  =  ss.bound[cc, labels[ii]]                #boundary
    }
  }
  # Covariance on the Average Column 
  Cov[[cc]][gamma[[cc]][1], 1:N]              = colSums(Cov[[cc]][1:N, 1:N])/(N-1)
  Cov[[cc]][1:N, gamma[[cc]][1]]              = t(Cov[[cc]][gamma[[cc]][1], 1:N])
  Cov[[cc]][gamma[[cc]][1], gamma[[cc]][1]]   = sum(diag(Cov[[cc]])[1:N])/(N-1)
  Cov[[cc]][(N+1), gamma[[cc]][1]]            = Cov[[cc]][gamma[[cc]][1],(N+1)]  = sum(Cov[[cc]][1:N,(N+1)])/(N-1)
}


return(list(Cov=Cov, gamma = gamma))


}



