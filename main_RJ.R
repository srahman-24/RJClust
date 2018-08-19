###### main function ###### 
##### author: Shahina Rahman ######


rm(list = ls())

library(mclust)

source("RJMclust.R")
source("RJclust.R")

Dyrskjot      = read.csv("Dyrskjot-2003.csv", header = T, row.names = 1)
D             = t(Dyrskjot)

# Paste the label information Source1/Source2 of the dataset Dyrskot from README.md file here. 

Z  = scale(log(D), center = T, scale =T)

RJM = RJMclust(Z)
RJM$ami
RJM$n.clust
table(group, RJM$groupG)

RJ = RJclust(Z, C_max)
RJ$ami 
RJ$n.clust 
table(group, RJM$groupG)




