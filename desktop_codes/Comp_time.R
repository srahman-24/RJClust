bic_Ali    = c(5544.825, bic_Ali)
#bic_Shipp  = c(6214.05, bic_eval)
bic_Arm    = c(5763.56, bic_Arm)
bic_Golub  = c(9202.624, bic_Golub)
bic_Chow   = c(4616.859, bic_Chow)
bic_Dyrs   = c(855.9858, bic_Dyrs)
bic_Pom    = c(2418.830, bic_Pom)  
bic_Su     = c(-8037.413, bic_Su)
bic_Nutt   = 

############ All in 1 graph  ##########
BIC.plot = data.frame(Data = c(rep("Alizadeh",6), rep("Armstorng", 7), rep("Dyrskot", 10),
                          rep("Chowdary",3), rep("Golub", 10), rep("Su", 19)), 
                 Clusters = c(1:6, 1:7, 1:10, 1:3, 1:10, 1:19), 
                 BIC.plot = c(bic_Ali, bic_Arm, bic_Dyrs+6000, log(bic_Chow)*1000, 
                         bic_Golub, bic_Shipp, bic_Su) )

ggplot(data = BIC.plot, aes(x = Clusters, y = BIC.plot, group = Data)) + 
  geom_line(aes(color = Data), size = 1.4) +
  geom_point(aes(color = Data), size = 3)

# Best performance  Datasets 

BIC = data.frame(Data = c(rep("Alizadeh",6), rep("Armstorng", 7), 
                          rep("Chowdary",3), rep("Golub", 10), rep("Dyrskot", 10), rep("Pomeroy",7)), 
                 Clusters = c(1:6, 1:7, 1:3, 1:10, 1:10, 1:7), 
                 BIC = c(bic_Ali, bic_Arm, log(bic_Chow)*700, 
                         bic_Golub, bic_Dyrs+6000, bic_Pom + 3300) )

ggplot(data = BIC, aes(x = Clusters, y = BIC, group = Data)) + 
  geom_line(aes(color = Data), size = 1.4) +
  geom_point(aes(color = Data), size = 3) + 
  xlim(1, 8) + ylim(5300, 11500) + labs(x = "C", title = "Optimal choice of Clusters") + 
  theme(plot.title   = element_text(hjust = 0.5, size=30, face="bold"),
        axis.title.x = element_text(size=30, face="bold"), 
        axis.title.y = element_text(size=30, face="bold"), 
        legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=12, face = "bold"))
  


#Computation time 
Ali_time  = c(0.146, 4.85, 19.21, 30.94, 28.33, 768.64, 699.91, 3136.17 )
Arm_time  = c(0.222, 5.76, 34.30, 23.91, 22.91, 723.9, 686.56, 1419.23)
Dyrs_time = c(0.047, 5.52, 5.110, 11.56, 10.32, 144.4, 129.13, 756.32)
Chow_time = c(0.130, 3.23, 5.14, 21.53, 20.42, 41.48, 35.96, 215.90)
Gol_time  = c(0.123, 6.57, 16.20, 22.44, 21.31, 916.06, 881.12, 5310.05) 
Pom_time  = c(0.110, 4.63, 15.83, 25.69, 22.02, 196.06, 182.45, 857.42) 
#Su_time   = c(2.898, 627.04, 306.64, 287.63, 280.24, 1164.06, 1482.45, 2857.42) 

Comp_time = data.frame(Data = c(rep("Chowdary",8), rep("Dyrskot", 8), rep("Pomeroy", 8),
                       rep("Armstrong", 8), rep("Alizadeh",8), rep("Golub", 8)),
                       Methods   = rep(c("RJ-Diagonal","RJ-Full", "GAP", "CC-km", "CC-hc", "CC-gmm", "CLEST", "CC-som"), 6), 
                       Comp.time = c(log(Chow_time), log(Dyrs_time), log(Pom_time), log(Arm_time), log(Ali_time), log(Gol_time)
                              ))

ggplot(data = Comp_time, aes(x = Data, y = Comp.time, group = Methods)) + 
  geom_line(aes(color = Methods), size = 1.4) +
  geom_point(aes(color = Methods), size = 3) + 
  labs(x = "Datasets", y = "log(secs)", title = "Computation Time") + 
  theme(plot.title   = element_text(hjust = 0.5, size=30, face="bold"),
        axis.title.x = element_text(size=30, face="bold"), 
        axis.title.y = element_text(size=30, face="bold"),
        axis.text.x  = element_text(size=10, face = "bold"),
        legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=12, face = "bold"))

Comp.time.mat = rbind(Ali_time, Arm_time, Dyrs_time, Chow_time, Gol_time, Pom_time)
colMeans(Comp.time.mat)

Comp_time1 = data.frame(Methods = c("RJ-Diagonal","RJ-Full", "GAP", "CC-km", "CC-hc", "CC-gmm", "CLEST", "CC-som"), 
                        Comp.mean = log(colMeans(Comp.time.mat)))

ggplot(Comp_time1, aes(x=Methods, y = Comp.mean, fill = Methods)) + geom_bar(stat= "identity")+
    coord_flip() + labs(y = "log(secs)", title = "Mean computation time") + 
  theme(plot.title   = element_text(hjust = 0.5, size=30, face="bold"),
        axis.title.x = element_text(size=30, face="bold"), 
        axis.title.y = element_text(size=30, face="bold"), 
        axis.text.x  = element_text(size=20, face = "bold"),
        axis.text.y  = element_text(size=12, face = "bold"),
        legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=12, face = "bold"))

