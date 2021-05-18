library(Spectrum)
library(RJcluster)

### write email to the maintainer or author. 

source("/Users/srahman/Desktop/paper/RJClust/desktop_codes/sim_data.R")
X = sim_data(c(20,20,20,20), 220, 1,1,1)

Spec_func = function(X)
{
  res = Spectrum(t(X))
  return(res)
}

Mutual_Information(Spec_func(X$X)$assignments, X$group)$ami

'Error message: 
***Spectrum***
  detected views: 1
method: 1
kernel: density
calculating similarity matrix 1
Error in Rfast::Dist(t(mat)) : 
  Not compatible with requested type: [type=list; target=double].'

### 2016 TCGA Data. 

dim(brain[[1]]) # RNA seq data

dim(brain[[2]]) #miRNA seq data

dim(brain[[3]]) # protein array data


boxplot(brain[[3]][1:10 ,])


brain[[3]]
