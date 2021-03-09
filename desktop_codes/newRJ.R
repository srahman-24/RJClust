
# toy example 
RR   = matrix(c(1:16), nrow = 4, ncol = 4, byrow = T)
RR
JJ   = NULL 
for(ii in 1:4)
{
  JJ = rbind(JJ, RR[ii,][-ii])
}
JJ   = cbind(JJ, diag(RR))
JJ

# Alizadeh 

Alizadeh= read.csv("/Users/srahman/Documents/Clustering/alizadeh-2000-v2_database1.csv", header = T, row.names = 1) 
D   = t(Alizadeh) 
Z1  = D[1:41,]
Z1  = D
p   = ncol(Z1)
n   = nrow(Z1)
RR  = Z1%*%t(Z1)/p
JJ   = NULL 
for(ii in 1:n)
{
  JJ = rbind(JJ, RR[ii,][-ii])
}
JJ   = cbind(JJ, diag(RR))
MclustJJ      = Mclust(JJ, modelNames = "VVI")
groupJ        = MclustJJ$class 
group = group_Ali1 = c(rep(2,4), rep(1,7), 2, 1, rep(2,3), rep(1,6), 2,1,2,1,rep(2,3),1,1,2,2,1,rep(2,4),1,2,2,1)
group = group_Ali = c(group[1:41], rep(3,9), rep(4,11),1)
table(group,groupJ)

GG  = Z1%*%t(Z1)/p
gg = GG 
gg_wodiag       = gg - diag(diag(gg))
gg_wdiag        = cbind(gg_wodiag, diag(gg))
GG_new          = cbind(gg_wodiag + diag(colSums(gg_wodiag))/(n-1), diag(gg))
Mclustgg        = Mclust(GG_new, modelNames = "VVI")
groupG          = Mclustgg$class
table(group,groupG)
