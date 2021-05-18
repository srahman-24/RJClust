library(parallel)
library(doParallel)
library(foreach)

### lapply()  and mclapply()  detectCores(), 

fx       = function(nstart) kmeans(Boston, 4, nstart = nstart)
starts   = rep(100, 40)
numCores = detectCores()
numCores

system.time({
  results = lapply(starts, fx)}
)

system.time({
  results <- mclapply(starts, fx, mc.cores = numCores)
})

x       =  iris[which(iris[,5] != "setosa"), c(1,5)]
res     =  NULL
trials  = seq(1, 10000)
boot_fx = function(trial) {
  ind     = sample(100, 100, replace=TRUE)
  result1 = glm(x[ind,2]~x[ind,1], family = binomial(logit))
  r  = coefficients(result1)
  res = cbind(res, r)
}

system.time({
  results <- lapply(trials, boot_fx)
})

system.time({
  results <- mclapply(trials, boot_fx, mc.cores = numCores)
})


##### foreach()
##### doParallel()
system.time({
  for (i in 1:3) {
    print(sqrt(i))
  }
})

system.time({
  lapply(1:3, sqrt)
})

system.time({
  mclapply(1:3, sqrt, mc.cores = numCores)
})

system.time({
  library(foreach)
  foreach(i = 1:3) %do% {
    sqrt(i)
  }})


foreach(i = 1:3, .combine = c) %do% {
  sqrt(i)
}


#### %dopar% is the doparallel vigenette
library(foreach)
library(doParallel)
registerDoParallel(numCores)

x <- iris[which(iris[,5] != "setosa"), c(1,5)]
trials <- 10000
system.time({
  r <- foreach(i = 1:trials, .combine = c ) %dopar% {
    ind <- sample(100, 100, replace = TRUE)
    result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
    coefficients(result1)
  }
})