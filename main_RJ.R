###### main function ###### 
##### author: Shahina Rahman ######


rm(list = ls())

library(mclust)

source("RJMclust.R")
source("RJclust.R")

Dyrskjot      = read.csv("dyrskjot-2003_database.csv", header = T, row.names = 1)
D             = t(Dyrskjot)

# Uncomment the label information from Source1/Source2 of the dataset Dyrskjot from README.md file here. 

#group1         = c(1:nrow(D));
#group1[c(1,2,5,6,8,9,11,22,34)] = 1;  #T2+
#group1[c(3,4,7,10,12,13,14,16,17,18,19,23,26:33)] = 2; #Ta
#group1[c(15,20,21,24,25,35:40)] = 3;  #T1

#group2 = group1
#group2[c(3,13,14,30)] = 4; #Ta2
#group2[c(4,7,12,16,18,19,23,31,33)]= 5; #Ta3

group = group2;


Z  = scale(log(D), center = T, scale =T)

RJM = RJMclust(Z, group)
RJM$ami
RJM$n.clust
table(group, RJM$groupG)

RJ = RJclust(Z, C_max, group)
RJ$ami 
RJ$n.clust 
table(group, RJM$groupG)




