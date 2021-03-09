comp_time = c(0.092, 7.72, 13.69, 526.23, 616.82, 689.36)
barplot(log(comp_time), ylab = "Average log(Secs)", col = c("blue", "gray", "orange", "red", "magenta", "purple"), 
        names.arg = c("RJClust", "EC", "GAP", "CC-HC", "CC-KM", "CLEST"), cex.lab = 1.3,
        cex.axis = 1.5, cex.names = 1.5)





Gclust = c(0,0,0,0,1)
GAP  =   c(1,0,3,3,1)
CC_ESM = c(1,0,1,2,1)
CC_HC =  c(1,1,3,3,1)
CC_KM =  c(1,1,1,2,0)
CLEST =  c(3,1,3,3,2)


counts = rbind(Gclust,GAP,CC_ESM,CC_HC,CC_KM,CLEST)
barplot(t(counts), beside = T,  col= terrain.colors(5), 
        cex.names = 1.1, cex.axis = 1.7, ylab = " |  C.est  -  C.true  |",cex.lab = 1.3)
legend("topleft", c("Data1","Data2","Data3","Data4","Data5"), cex=1.5, 
       fill=terrain.colors(5))

