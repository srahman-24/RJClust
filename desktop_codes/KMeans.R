library(datasets)
data(iris)    #Petal.length  #Petal.width

######## K-Means Clustering ######
install.packages("ggplot2")
library(ggplot2)
ggplot(iris, aes(Petal.Length, Petal.Width, color = Species)) + geom_point()
#Clustering
iris.cluster = kmeans(iris[,3:4],3, nstart = 20)
iris.cluster
table(iris.cluster$cluster, iris$Species)


######## Hierarchical Clustering ####### 
clusters <- hclust(dist(iris[, 3:4])) #diag = FALSE, upper = FALSE
plot(clusters)
clusterCut <- cutree(clusters, 3)
table(clusterCut, iris$Species)
######  using a different linkage method #####
clusters <- hclust(dist(iris[, 3:4]), method = 'average')
plot(clusters)
clusterCut <- cutree(clusters, 3)
table(clusterCut, iris$Species)



