###### main function ###### 
##### author: Shahina Rahman ######


rm(list = ls())

library(mclust)
library(devtools)
#library(e1071)
source("RJMclust.R")
source("RJclust.R")

Dyrskjot      = read.csv("dyrskjot-2003_database.csv", header = T, row.names = 1)
D             = t(Dyrskjot)


group         = c(1:nrow(D))
group[c(1,2,5,6,8,9,11,22,34)] = 1
group[c(3,4,7,10,12,13,14,16,17,18,19,23,26:33)] = 2
group[c(15,20,21,24,25,35:40)] = 3
group[c(3,13,14,30)] = 4
group[c(4,7,12,16,18,19,23,31,33)]= 5

Z  = scale(log(D), center = T, scale =T)

RJM = RJMclust(Z)
RJM$ami
RJM$n.clust
table(group, RJM$groupG)

RJ = RJclust(Z, C_max)
RJ$ami 
RJ$n.clust 
table(group, RJM$groupG)




