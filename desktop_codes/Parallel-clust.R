library(mclust)
library(RJcluster)
library(clustvarsel)

### parallel
library(foreach)
library(parallel)
library(doParallel)
library(MASS)
numCores = detectCores()
#register a parallel backend using one of the packages that begin with do
registerDoParallel(numCores)

set.seed(44)
Seeds = sample(0:10000, 100, replace = FALSE)

system.time({
foreach(i = 1:10, .combine = rbind) %dopar% {
    set.seed(Seeds[i]) 
    rnorm(3)
}
})

foreach(i = 1:10, .combine = rbind) %dopar% {
  
sim   = sim_data(c(20,20,20,20), 220, 1,2, Seeds[i])
X     = sim$X 
group = sim$group

cl    = sparcl_fx(X, 4)
Mutual_Information(cl, group)$ami

cl    = Raftery_fx(X)
Mutual_Information(cl$class, group)$ami
cl$G

cl    = HDDC_fx(X)
Mutual_Information(cl$class, group)$ami
cl$K

}


