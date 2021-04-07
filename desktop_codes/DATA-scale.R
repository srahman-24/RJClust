##### Data #######

#### 1. Alizadeh 1 (no log, no center, no scale)
Alizadeh    = read.csv("/Users/srahman/Documents/Work/Clustering/alizadeh-2000-v1_database.csv", header = T, row.names = 1)
D  = t(Alizadeh)
boxplot(D[,1:100])
X = D
p = ncol(X)
N = nrow(X)
group = 1*grepl('DLBCL1.', colnames(Alizadeh)) 

#####2.  Alizadeh 2 
Alizadeh2          = read.csv("/Users/srahman/Documents/Work/Clustering/Microarray-data/alizadeh-2000-v2_database1.csv",header = T, row.names = 1) 
D                  = t(Alizadeh2)
boxplot(D[,1:100])
group = c(rep(1,41), rep(2,9), rep(3,11),1)

####3.  Alizadeh 3 
group1 = c(rep(2,4), rep(1,7), 2, 1, rep(2,3), rep(1,6), 2,1,2,1,rep(2,3),1,1,2,2,1,rep(2,4),1,2,2,1)
group  = group_Ali =  c(group1[1:41], rep(3,9), rep(4,11),1)

####4.  Armstrong 1
Armstrongv1   = read.csv("/Users/srahman/Documents/Work/Clustering/Microarray-data/armstrong-2002-v1_database.csv", header = T, row.names = 1) 
D             = t(Armstrongv1)
group = group_Arm1 = c(rep(1,24), rep(2,48))
boxplot(D[,1:100])
X = log(D)
boxplot(sqrt(colVars(X)))
X = scale(log(D), center = T, scale =  T)
boxplot(X[,1:100])


####5.  Armstrong 2 (log, center, no scale)
Armstrongv2   = read.csv("/Users/srahman/Documents/Work/Clustering/Microarray-data/armstrong-2002-v2_database.csv", header = T, row.names = 1)
D             = t(Armstrongv2)
boxplot(D[,1:100])
X             = scale(log(D), center = T, scale =  T)
boxplot(X[,1:100])
group         = group_Arm2 = c(rep(1,24), rep(2,20), rep(3,28))
X = scale(log(D), center = T, scale =  T)
boxplot(X[,1:100])

####6. Bhattacharya (log, center, scale)

Bhattacharjee = read.csv("/Users/srahman/Documents/Work/Clustering/Microarray-data/bhattacharjee-2001_database.csv", header = T, row.names = 1) #
D             = t(Bhattacharjee)
group         = group_Bhatt = c(rep(1,139), rep(2, 156-139), rep(3, 162-156), rep(4,183-162), rep(5, 203-183))
boxplot(D[,1:100])
X = log(D)
boxplot(sqrt(colvars(X)))
X = scale(log(D), center = T, scale =  T)
boxplot(X[,1:100])


####7. Bittner 
Bittner       = read.csv("/Users/srahman/Documents/Work/Clustering/Microarray-data/bittner-2000_database.csv", header = T, row.names = 1) #13
D             = t(Bittner)
boxplot(D[,1:100])
group         = c(1:nrow(D))
group         = group_Bittner = c(rep(2,12), rep(1,19), rep(2,7))



###8. Bredel
Bredel        = read.csv("/Users/srahman/Documents/Work/Clustering/Microarray-data/bredel-2005_database.csv", header = T, row.names = 1)  #11
D             = t(Bredel) 
boxplot(D[,1:100])
X = scale(log(D), center = T, scale =  T)
boxplot(X[,1:100])
group         = group_Bredel = c(rep(1,14), rep(2,14))


###9. Chowdary

Chowdary      = read.csv("/Users/srahman/Documents/Work/Clustering/Microarray-data/chowdary-2006_database.csv", header = T, row.names = 1) #15 
D             = t(Chowdary)
boxplot(D[,1:100])
X = scale(log(D), center = T, scale =  T)
boxplot(X[,1:100])
group         = c(rep(1,62), rep(2,42))


####10. 
Dyrskjot      = read.csv("/Users/srahman/Documents/Work/Clustering/Microarray-data/dyrskjot-2003_database.csv", header = T, row.names = 1)
D             = t(Dyrskjot)
boxplot(D[,1:100])
X = scale(log(D), center = T, scale =  T)
boxplot(X[,1:100])
group = NULL
group = c(1:nrow(D))
group1         = c(1:nrow(D));
group1[c(1,2,5,6,8,9,11,22,34)] = 1; #T2+
group1[c(3,4,7,10,12,13,14,16,17,18,19,23,26:33)] = 2; #Ta
group1[c(15,20,21,24,25,35:40)] = 3; #T1
group2 = group1;
group2[c(3,13,14,30)] = 4; #Ta2
group2[c(4,7,12,16,18,19,23,31,33)]= 5; #Ta3
group = group1; 
group = group2   #important 






