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


####10. Dyrskjot
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



####11. Garber 
Garber        = read.csv("/Users/srahman/Documents/Work/Clustering/Microarray-data/garber-2001_database.csv", header = T, row.names = 1)  
D             = t(Garber)
boxplot(D[,1:100])
X = D
group = NULL
group = group_Garber    = c(1:nrow(D))
group[c(1:3,8,9,13,15,16,22,23,25,27,36,41,48,63,66)] = 1
group[c(4,5,7,10,11,12,17,19,20,21,24,26,28,29,30,31,32,35,37:40,42:47,51,52,54,56:62,64,65)] = 2
group[c(6,14,34,50,55)] = 3
group[c(18,33,49,53)]   = 4


####12. Golub
Golub       = read.csv("/Users/srahman/Documents/Work/Clustering/Microarray-data/golub-1999-v2_database.csv", header = T, row.names = 1)  #12  
D             = t(Golub)
boxplot(D[,1:100])
X = scale(log(D), center = T, scale =  T)
boxplot(X[,1:100])
group        = NULL 
group = c(1:nrow(D))
group[c(1,4,5,7,8,12,13,15:22,24:27,39:58)] = 1
group[c(2,3,6,9,10,11,14,23,55)] = 2
group[c(28:38,59:72)] = 3
group = c(1:nrow(D))
group[c(1,4,5,7,8,12,13,15:22,24:27,39:58)] = 1
group[c(2,3,6,9,10,11,14,23,55)] = 2
group[c(28:38,59:72)]= 2



####13. Gordon 
Gordon        = read.csv("/Users/srahman/Documents/Work/Clustering/Microarray-data/gordon-2002_database.csv", header = T, row.names = 1)  #12  
D             = t(Gordon)
boxplot(D[,1:100])
X = scale(log(D), center = T, scale =  T)
boxplot(X[,1:100])
group = group_Gordon = c(rep(1,31), rep(2,150))


####14. Laiho
Laiho         = read.csv("/Users/srahman/Documents/Work/Clustering/Microarray-data/laiho-2007_database.csv", header = T, row.names = 1)   #4
D             = t(Laiho)
boxplot(D[,1:100])
X = scale(log(D), center = T, scale =  T)
boxplot(X[,1:100])
group = group_Laiho = c(rep(1,8), rep(2,29))



#####15. Lapointev1
Lapointev1    = read.csv("/Users/srahman/Documents/Work/Clustering/Microarray-data/lapointe-2004-v1_database.csv", header = T)
D             = t(Lapointev1)
boxplot(D[,1:100])
X = D
group  = NULL 
group  = c(1:nrow(D))
group[c(1,11,20,28,30,36,52,55,56,61,66)] = 1
group[c(2,3,8,9,10,12,16,17,18,19,21,24,26,29,31,32,35,37,39,40,41,42,43,44,45,46,49,50,51,53,54,58,59,60,62,64,65,68,69)] = 1
group[c(4,5,6,7,13,14,15,22,23,25,27,33,34,38,47,48,57,63,67)] = 2

#####16. Lapointev2
Lapointev2    = read.csv("/Users/srahman/Documents/Work/Clustering/Microarray-data/lapointe-2004-v2_database.csv", header = T)
D             = t(Lapointev2)
boxplot(D[,1:100])
X = D
group = NULL 
group = group_La2 = c(1:nrow(D))
group[c(1,17,32,44,46,55,81,85,86,96,103)] = 1
group[c(2,4,10,12,14,16,22,24,26,28,29,30,35,36,42,43,48,50,52,58,59,66,68,69,70,72,75,76,80,84,88,89,90,93,94,99,100,
        104,105,106,109)] = 2
group[c(3,5,11,13, 15,18,23,25,27,31,33,38,40,45,47,49,54,56,60:65,67,71,77,78,79,82,83,91,92,95,97,101,102,108,110)] = 3
group[c(6:9, 19:21, 34, 37, 39, 41, 51, 53, 57, 73, 74, 87, 98, 107)] = 4



Liang         = read.csv("C:\\Users\\srahman\\Documents\\Data Requested\\My Documents\\Clustering\\Real Data\\liang-2005_database.csv", header = T, row.names = 1)  
D             = t(Liang)  
group = NULL 
group = c(1:nrow(D))
group[c(1,4,5,6,10,11)] = 1
group[c(2,3,7,8,9,12:26,30:37)] = 2
group[c(27:29)] = 3