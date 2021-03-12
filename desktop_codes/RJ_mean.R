RJ_mean = function(K, class, GG){

gamma     = list()
mean.diag = rep(1,K)
nn        = matrix(0, nrow = K, ncol = K)
mean.off  = matrix(0, nrow = K, ncol = K)
N         = length(class)

for(ii in 1:K)
{
  gamma[[ii]]         =  which(class == ii) 
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
      mean.off[ii,jj]     = (sum(GG_1) - sum(diag(GG_1)))/nn[ii,jj]  #offdiagonals mean
    }
    else
    {
      #hetero-off diagonals 
      GG_2                   = GG[gamma[[ii]], gamma[[jj]]]
      mean.off[ii,jj]        = mean(GG_2)
    }
  }
}

# mean.diag     # real diagonals 
# mean.off      # diagonals are homogeneous off diagonals
if(K > 1)
{
mean.off      =  mean.off + t(mean.off) - diag(diag(mean.off))
}


label = class
MU = matrix(0, nrow = K, ncol = N+1)
for(ii in 1: K)
{
  object  =  gamma[[ii]][1]
  for(jj in 1:N)
  {
    if(label[jj] == label[object])
    {
      MU[ii, jj] = diag(mean.off)[label[object]]
    }
    if(label[jj] != label[object])
    {
      MU[ii, jj] = mean.off[label[object], label[jj]]
    }
  }
      MU[ii, N+1] = mean.diag[label[object]]
}

return(MU)

}




