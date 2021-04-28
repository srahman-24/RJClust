library(Spectrum)
library(RJcluster)


Spec_func = function(X)
{
  res = Spectrum(t(X))
  return(Mutual_Information(res$assignments[1:40], group)$ami)
}

### 2016 TCGA Data. 

dim(brain[[1]])

dim(brain[[2]])

dim(brain[[3]])


boxplot(brain[[1]][1:10 ,])
