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

# group information of the dataset Dyrskot can be obtained from README.md file. 

Z  = scale(log(D), center = T, scale =T)

RJM = RJMclust(Z)
RJM$ami
RJM$n.clust
table(group, RJM$groupG)

RJ = RJclust(Z, C_max)
RJ$ami 
RJ$n.clust 
table(group, RJM$groupG)




