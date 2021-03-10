gamma     = list()
labels    = rep(1,N) 
K         = BIC_GG$G
mean.diag = rep(1,C)
nn        = matrix(0, nrow = C, ncol = C)
mean.off  = matrix(0, nrow = C, ncol = C)

for(ii in 1:K)
{
  gamma[[ii]]       =  which(BIC_GG$class == ii) 
  labels[gamma[[ii]]] =  ii
}


for(ii in 1:K)
{ 
  for(jj in ii:K)
  {
    if(ii == jj)
    { 
      #homo-off diagonals
      GG_1                = GG[gamma[[ii]], gamma[[jj]]]
      nn[ii,jj]           = (length(GG_1) - nrow(GG_1))
      mean.diag[ii]       = mean(diag(GG_1))   #diagonal mean
      M = mean.off[ii,jj] = (sum(GG_1) - sum(diag(GG_1)))/nn[ii,jj]  #offdiagonals mean
    }
    else
    {
      #hetero-off diagonals 
      GG_2                   = GG[gamma[[ii]], gamma[[jj]]]
      nn[ii,jj]              = length(GG_2)
      M = mean.off[ii,jj]    = mean(GG_2)
    }
  }
}


