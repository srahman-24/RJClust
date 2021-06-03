library(mclust)
library(RJcluster)

file.choose()
data = read.table("/Users/srahman/Downloads/samples50.txt")
dim(data)
n = nrow(data)
p = ncol(data)
boxplot(log(data[, 1:10] + 1))
X = as.matrix(log(data + 1))
X = scale(X, center = T, scale = T)
boxplot(X[,1:10])

##### Ignore RJclust for n = 50, we will update the technical glitch ########
res = RJclust(X, penalty = "mclust", seed = 2)
res$class
table(res$class, true_labels = c(rep(0,25), rep(1,25)))

##### Implement the following if RJclust dont work ######
#### install.packages("mclust") 

library(mclust)

p           =  ncol(X)
N           =  nrow(X)
GG          =  tcrossprod(X, X)/p
gg_wodiag   =  GG - diag(diag(GG))
GG_new      =  cbind(gg_wodiag + diag(colSums(gg_wodiag)/(N-1)), diag(GG))

BIC_GG      =  Mclust(GG_new, modelNames = "VVI")
BIC_GG$G
plot(BIC_GG$BIC)
table(BIC_GG$classification, true_labels = c(rep(0,25), rep(1,25)))


blah blah 
