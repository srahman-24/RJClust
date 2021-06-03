
Runtime = read.csv("/Users/srahman/Documents/Work/Clustering/RJclust-RCode/Runtimes.csv")
Runtime

data_num = nrow(Runtime)
times = rbind(as.matrix(Runtime[,2]), as.matrix(Runtime[,3]), as.matrix(Runtime[,4]),
                  as.matrix(Runtime[,5]), as.matrix(Runtime[,6]), as.matrix(Runtime[,7]),
                  as.matrix(Runtime[,8]))

methods = c(rep("SPEC", 32), rep("AP", 32), rep("GAP", 32), rep("HDDC", 32),
            rep("RCC", 32), rep("RJ", 32), rep("Cvarsl", 32))


secs = times
Run = data.frame(log(secs), methods)

Run$times
Run$methods

library(ggplot2)
ggplot(Run, aes(x = methods, y = log(secs), fill = methods)) + geom_boxplot() + 
  theme(
    legend.position="none",
    plot.title   = element_text(color = "black", size = 18, face="bold"), 
    axis.title.x = element_text(color = "black", size = 18, face="bold"),
    axis.title.y = element_text(color = "black", size = 18, face="bold"), 
    axis.text.x  = element_text(face  = "bold",  size = 18),
    axis.text.y  = element_text(face  = "bold",  size = 18),
    axis.line    = element_line(colour = "darkblue", size = 2, linetype = "solid"))  +
     labs(x=" ", y = "Time (log secs)")  + coord_flip()

  #axis.line    = element_line(colour = "darkblue", size = 1, linetype = "solid"),
  #panel.grid.major.x = element_blank(),
  #panel.grid.minor.x = element_blank()

theme(
  plot.title   = element_text(color = "black", size = 24, face="bold"), 
  axis.title.x = element_text(color = "black", size = 24, face="bold"),
  axis.title.y = element_text(color = "black", size = 24, face="bold"), 
  axis.text.x  = element_text(face  = "bold",  size = 24),
  axis.text.y  = element_text(face  = "bold",  size = 24))  +

