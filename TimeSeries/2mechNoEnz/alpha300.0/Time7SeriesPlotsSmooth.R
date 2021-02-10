

i <- 1
 library("data.table")
#setwd("C:\\Users\\David\ Infortunio.DavidInfortunio\\Documents\\OttResearch\\GAMA")

#setwd("/home/david/BioResearch/t2")
#setwd("/home/david/BioResearch/t3")
#setwd("/home/david/BioResearch/t1")
setwd("/home/david/BioResearch/python/PostProcessing/Simulations/alpha300.0/plots")
#setwd("/tmp/Infort/Gamma2")

j <-  1

#   var(datas[,3])


fiile <- list.files()
fiile <- fiile[grepl(".csv", fiile)]
fiille <- fiile[j]


numFile <- length(fiile)

datas <-
  read.table(
    fiille,
    sep = ",",
    header = TRUE,
    check.names = FALSE,
    col.names = c("Time", "Proteins", "peak")
  )
datan <- dim(datas)[1]
library(stringr)


fiillee <- fiille



plotBio <- function(fiillee) {
library(stringr)
  fiillee2  <-   str_replace(fiillee, pattern   =  "Time.*mma.", replacement = "Time")
    fiillee2  <-   str_replace(fiillee, pattern   =  "zz*t", replacement = "T")
    fiillee2  <-   str_replace(fiillee, pattern   =  "zz*.*t", replacement = "")
  fiillee2  <-   str_replace( pattern   =  "const",fiillee2,replacement = "")
  fiillee2  <-   str_replace( pattern   =  ".csv",fiillee2,replacement = "")
  pdf(paste("/home/david/BioResearch/Plots/", "GTS,",fiillee2, date(),".pdf",sep   =   "")) # To make PDFs
  datas <-
    read.table(
      fiillee,
      sep = ",",
      header = TRUE,
      check.names = FALSE,
      col.names = c("Time", "Protein1", "Protein2")
    )
  file2  <-   str_replace( pattern   =  "00t",fiillee,replacement = "")
  file2  <-   str_replace( pattern   =  "000",file2,replacement = "")
  file2  <-   str_replace( pattern   =  "RunNum",file2,replacement = "")
  file2  <-   str_replace( pattern   =  "Time*Gamma",file2,replacement = "Time")
  file2  <-   str_replace( pattern   =  ".csv",file2,replacement = "")
#datas <- datas[datas[,1]<=200,]#4320
  #ifelse(grepl("G",file2),file2 <- paste(file2,"",sep = ""), file2 <- paste(file2,"ernoulli SD",sep = ""))

  plot(
    datas[, 1],
    datas[, 2],
    type   =   "l",
    xlab   =   "Time ",
    ylab   = "Protein1", col = "Blue",

      xlim = c(100,1900), ylim = c(0, 1000),
    main   =   file2
  )
    lines(
    datas[, 1],
    datas[, 3],col = 2)
  
  
  arrays<-read.csv("/home/david/BioResearch/python/PostProcessing/Simulations/alpha300.0/peaks.csv",header=FALSE)
  
   points(x=(arrays[23,1:95]+100),y=arrays[24,1:95],col=2, pch =3) #100 is burn in time, only do first 1500 minutes of peaks max
  points(x=(arrays[55,1:95]+100),y=arrays[56,1:95], col = 3, pch =2 ) #100 is burn in time
  if (FALSE){
  avgs<-16
  avgProtein <- frollmean(datas$Protein,avgs)[-floor(avgs/2):-1]
  avgProtein[1:floor(avgs/2)] <- datas$Protein[0:floor(avgs/2)]
  avgProtein <- append(avgProtein, tail(datas$Protein, ceiling(avgs/2)))

  #plot.window(        xlim = c(5,195), ylim = c(0, 10+max( datas[, 2])))
  
  
  lines(
    datas[, 1],??
    avgProtein,
col = "blue"
  )
  
  x11 <- t(c(50,200))
  y11 <- t(c(9,9))
  lines(x11, y11,col =5)


    y8p5 <- t(c(12,12))
  lines(x11, y8p5,col =6)
  datas3<- datas
  datas3[,3] <-avgProtein
#  Set Ratio 1055 and 661
  #dev.off()
  write.csv(file = paste(fiillee2,"fileRollAvg62.csv"),datas3)

}
  dev.off()
}
wddds <- getwd()
setwd("/home/david/BioResearch")

#pdf(paste("GamaTimeSeries,",date(),".pdf",sep   =   ""))
#jpeg(paste("GAMA,",date(),".jpeg",sep   =   ""))

#lapply( fiile, plotBio)

setwd(wddds)
lapply( fiile, plotBio)
#plotBio(fiile[2])
#dev.off()

