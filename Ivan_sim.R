
file.choose()
data = read.table("/Users/srahman/Downloads/samples300.txt")
dim(data)
n = nrow(data)
p = ncol(data)
boxplot(log(data[, 1:100]+1))
X = as.matrix(data) 

X = as.matrix(log(data +1))
X = scale(X, center = T, scale = T)
p           =  ncol(X)
N           =  nrow(X)
GG          =  tcrossprod(X, X)/p
gg_wodiag   =  GG - diag(diag(GG))
GG_new      =  cbind(gg_wodiag + diag(colSums(gg_wodiag)/(N-1)), diag(GG))

BIC_GG      =  Mclust(GG_new, modelNames = "VVI")
BIC_GG$G
plot(BIC_GG$BIC)
table(BIC_GG$classification, true_labels = c(rep(0,150), rep(1,150)))
#table(BIC_GG$classification, group)
#f_rez(BIC_GG$classification, group)
#plot(BIC_GG$BIC)