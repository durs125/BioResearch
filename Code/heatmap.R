set.seed(123)                                                     # Set seed for reproducibility
data <- matrix(rnorm(100, 0, 10), nrow = 10, ncol = 10)           # Create example data
colnames(data) <- paste0("col", 1:10)                             # Column names
rownames(data) <- paste0("row", 1:10)     

#install.packages("ggplot2")                                       # Install ggplot2 package
#library("ggplot2")    
par(mfrow =2, mfcol=2)

 
setwd("/home/david/BioResearch/tmp/2Mech32PtRollAvg/")
#setwd("/home/david/BioResearch/tmp/2Mech62PtRollAvg/")
AmCV<-read.csv("heatmeanDelayTime.csv")

heatmap(as.matrix(mapply(as.numeric, AmCV),ncol=dim(AmCV)[2])[,-1],labRow = AmCV[,1] )
legend()
#legend(location, title, legend, ...)
AmCV<-(read.csv("heatAmplitude.CV.csv"))
x <- log(as.matrix(mapply(as.numeric, AmCV)))
rc <- rainbow(nrow(x), start = 0, end = .3)
cc <- rainbow(ncol(x[,-1]), start = 0, end = .3)


heatmap(log(as.matrix(mapply(as.numeric, AmCV),ncol=dim(AmCV)[2])[,-1]),Colv = NA, Rowv = NA, scale="column",labRow = AmCV[,1] ,col = cm.colors(256),     main = "Amplitude.CV")



AmM<-read.csv("heatMean.Amplitude.csv")

heatmap(as.matrix(mapply(as.numeric, AmM),ncol=dim(AmM)[2])[,-1])

PerCV<-read.csv("heatPeriod.CV.csv")

heatmap(as.matrix(mapply(as.numeric, PerCV),ncol=dim(PerCV)[2])[,-1])

PerM<-read.csv("heatPeriod.Mean.csv")

heatmap(as.matrix(mapply(as.numeric, PerM),ncol=dim(PerM)[2])[,-1])



library("ggplot2")
ggplot(nba.m, aes(variable, Name)) + geom_tile(aes(fill = rescale),
+     colour = "white") + scale_fill_gradient(low = "white",
+     high = "steelblue"))
