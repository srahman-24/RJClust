library(Spectrum)
library(RJcluster)

### write email to the maintainer or author. 

Spec_func = function(X)
{
  res = Spectrum(t(X))
  return(res)
}

Mutual_Information(res$assignments, group)$ami

### 2016 TCGA Data. 

dim(brain[[1]]) # RNA seq data

dim(brain[[2]]) #miRNA seq data

dim(brain[[3]]) # protein array data


boxplot(brain[[1]][1:10 ,])


brain[[1]]$TCGA.HT.7855
