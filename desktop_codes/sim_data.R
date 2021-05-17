sim_data = function(n,p,sigma1, sigma2, seed){

  set.seed(seed)
  group = c(rep(1, n[1]), rep(2, n[2]), rep(3, n[3]), rep(4, n[4]))
  # Set1:  sigma = 1 ( SNR : high level) , and n  = c(20,20,20,20) 
  # Set2:  sigma = 2 ( SNR : low level) ,  and n  = c(20,20,20,20) 
  # Set1:  sigma = 1 ( SNR : high level) , and n  = c(20,20,200,200) (unbalanced)
  
  X    = matrix(rnorm(sum(n)*p,0, sigma2), nrow = sum(n), ncol = p, byrow = TRUE)
  
  #Cluster 1: N(2.5, sigma)(1-10), N(1.5, sigma)(11-20) 
  
  X[1:n[1],1:10]                                       =   rnorm(n[1]*10, 2.5, sigma1)
  X[1:n[1],(1 + 10):(10 + 10)]                         =   rnorm(n[1]*10, 1.5, sigma1)
  
  #Cluster 2: N(0, sigma)(1-10), N(1.5, sigma) (11-20)
  X[(n[1] + 1):(n[1] + n[2]),1:10]                         =   rnorm(n[2]*10, 0, sigma1)
  X[(n[1] + 1):(n[1] + n[2]),(1 + 10):(10 + 10)]           =   rnorm(n[2]*10, 1.5, sigma1)
  
  #Cluster 3: N(0, sigma)(1-10), N(-1.5,sigma)(11-20)
  X[(n[1] + n[2] + 1):(n[1] + n[2] + n[3]),1:10]               =   rnorm(n[3]*10, 0, sigma1)
  X[(n[1] + n[2] + 1):(n[1] + n[2] + n[3]),(1 + 10):(10 + 10)] =   rnorm(n[3]*10, -1.5, sigma1)
  
  #Cluster 4: N(-2.5,sigma)(1-10), N(-1.5, sigma)(11-20)
  X[(n[1] + n[2] + n[3] + 1):(n[1] + n[2] + n[3] + n[4]),1:10]                  =   rnorm(n[4]*10, -2.5, sigma1)
  X[(n[1] + n[2] + n[3] + 1):(n[1] + n[2] + n[3] + n[4]),(1 + 10):(10 + 10)]    =   rnorm(n[4]*10, -1.5, sigma1)

  return(list(X = X, group = group))
    
}

