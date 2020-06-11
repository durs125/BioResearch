#First I take a 9 point roll mean then I use the 2 line method, I do not use the rolling mean to change the peak height
# To make things right it should be f60,62, not f8,11
# Fourier transfourm via fft not posible since the spacing is wrong
  library("data.table")

i<-1
if (i){#initialize the system}
if ( length(dir("/scratch/Infor"))) {
  library(stringr, lib.loc  =  "/scratch/Infor" )
  #library(dotCall64, lib.loc  =  "/scratch/Infor" )
  library("Rcpp", lib.loc = "/scratch/Infor")
  #library(tidyverse, lib.loc  = "/home/infor/Infor" )

}else{
  library("stringr")
  library("Rcpp")
}
str_remove<-function(text,pattern) { str_replace(text,pattern,"") }

cleanName <- function(sList,pattern){
str_remove<-function(text,pattern) { str_replace(text,pattern,"") }
 if ( length(dir("/scratch/Infor"))) {
  library('stringr', lib.loc  =  "/scratch/Infor" )
 }  else{ #Make sure the libraries are available
  library("stringr")
 }#Make sure the libraries are available
 sList[1] <- str_remove(sList[1],pattern)
 return(sList)
}

#library(tidyverse, lib.loc  = "/home/infor/Infor" )

#install.packages("dotCall64", lib  =  "/scratch/Infor" )

library(tools)
library(parallel)
library(compiler)




#something I use when debugging

i <- 2

# if ( length(dir("/tmp/Infort/Gamma2"))) {

# dir.create("/tmp/")
# dir.create("/tmp/Infort")
# dir.create("/tmp/Infort/Gamma2")
# }
#setwd("C:\\Users\\David\ Infortunio.DavidInfortunio\\Documents\\OttResearch\\Gamma")
setwd("/tmp/Infort/Gamma2")


fiile <- list.files()
fiile <- fiile[grepl("csv",fiile)]
j <-  2
fiille <- fiile[j]
datas <-
  read.table(fiile[j],    sep = ",",    header = TRUE,    check.names = FALSE,    col.names = c("Time", "Proteins", "peak")  )
#   var(datas[,3])
datan <- dim(datas)[1]



outerLoop2 <- function(fiille) {



  datas <-
    read.table(
      fiille,
      sep = ",",
      header = TRUE,
      check.names = FALSE,
      col.names = c("Time", "Proteins", "peak")
    )

  names(datas)[3] <- "peak"

  datas[, 3] <- 0

  len <- dim(datas)[1]
  #  datas<-datas[datas$Time< = 200,]#make it remove values over the 90 minute time max
    # RollingAverage Rolling Average rolling mean
  avgs<-32
  avgProtein <- frollmean(datas$Protein,avgs)[-floor(avgs/2):-1]
  avgProtein[0:floor(avgs/2)] <- datas$Protein[0:floor(avgs/2)]
  avgProtein <- append(avgProtein, tail(datas$Protein, ceiling(avgs/2)))
  datas[, 3] <- findPeaks(len, avgProtein, datas$peak)
  
  # want to take the time at the highest value for a given peak value
  

Creast <- aggregate(datas[,c(2,1)],
              by  =  list(datas$peak),
              max,
              na.rm  =  TRUE)
  period <-
    aggregate(datas$Time,
              by  =  list(datas$peak),
              min,
              na.rm  =  TRUE)
  if (dim(Creast)[1] - 1 > 0) {
        period2 <- list(matrix(NA, nrow  =  dim(Creast)[1] - 1))
  }else{
        warning("There seems to be no oscilation")
  }
  difk <- as.numeric(Creast[c(-1, -2), 3] - Creast[c(-1, -(length(Creast[, 2]))), 3])
  period2[[1]] <- difk
  amplitudes <- aggregate(datas$Protein,
                          by  =  list(datas$peak),max, na.rm  =  TRUE)
  amplitudes <- amplitudes[-1, -1]  # get rid of first peak
   MeanAmps <- mean(amplitudes[-length(amplitudes)], na.rm  = TRUE)

   Everest <- sort(amplitudes,partial=length(amplitudes)-1)[length(amplitudes)-1]
   SDAmps <- sd(amplitudes[-length(amplitudes)], na.rm  = TRUE)
   MeanPeriod <- mean(period2[[1]], na.rm  = TRUE)
   SDPeriod <- sd(period2[[1]], na.rm  = TRUE)
  if ( length(dir("/scratch/Infor"))) {
    library(stringr, lib.loc  =  "/scratch/Infor" )


  }else{
    library("stringr")

  }
  
  str_remove<-function(text,pattern) { str_replace(text,pattern,"") }
   Var<-str_remove(fiille,"zz.*Time")
   #       Var<-aa[[1]][1];Var<-str_remove(Var,"zz.*Time")
   #            
   

     # Var <- str_remove(Var,"zz...0*0")
  Var <-  str_remove(Var,"=Var*.*")
CV <- str_remove(fiille,".*mma.")
CV <- str_remove(CV,"zz..0*0")


CV <- str_remove(fiille,".*mma.")

  yr <- "1?"
#ifelse(grepl(fiille,"yr"),yr = str_remove(CV,"=CV*.*")
       CV <- as.numeric(str_remove(CV,"=CV*.*"))
       t<- str_remove(fiille,"DelTime*.*")
       t<- str_remove(t,"t")

  return(c(fiille, MeanPeriod, SDPeriod, MeanAmps, SDAmps, Var,CV,yr,t,Everest))
}
Loops2 <- cmpfun(outerLoop2)

cppFunction(
       '
       Rcpp::NumericVector findPeaks(int len,  Rcpp::NumericVector Proteins, Rcpp::NumericVector peak ){
       int height = 0;

       for(int i = 2; i<len;  ++i){
// this is f 60,62
       if (height <= 62 ) {
       height = Proteins[i];
       peak[i] = peak[i - 1];
       }else{
       if (Proteins[i] < 61  ) {
       peak[i] = peak[i - 1] + 1;
       height = 0;
       } else{peak[i] = peak[i-1]; }
       }
       }
       return peak;
       }
       '
)


AnalyzeBio <- function(durs){

 library(tools)
 library(parallel)
 library(compiler)
 str_remove<-function(text,pattern) { str_replace(text,pattern,"") }
  if ( length(dir("/scratch/Infor"))) {
  library('stringr', lib.loc  =  "/scratch/Infor" )
 }  else{ #Make sure the libraries are available
  library("stringr")
 }#Make sure the libraries are available
 
 #setwd("/tmp/Infort/Gamma2")
  setwd("~/")
 setwd(durs)
 fiile <- list.files()
 fiile <- fiile[grepl("csv",fiile)]
 #numFile <- length(fiile)

 # get a list of csv files
 if ( length(dir("/scratch/Infor"))) {
  cl2 <- makePSOCKcluster(getOption("cl.cores", 19)) #cores is 3, sets the cluster
  clusterEvalQ(cl2,library("Rcpp",lib.loc = "/scratch/Infor"))

 }else{
  cl2 <- makePSOCKcluster(getOption("cl.cores", 3)) #cores is 3, sets the cluster
  clusterEvalQ(cl2,library("Rcpp"))

 }
  clusterEvalQ(cl2,library("data.table"))
 clusterEvalQ(cl2,
              cppFunction('
                                Rcpp::NumericVector findPeaks(int len, Rcpp::NumericVector Proteins, Rcpp::NumericVector peak ){
                                int height = 0;

                                for(int i = 2; i<len;  ++i){

                                if (height <=  62) {
                                height = Proteins[i];
                                peak[i] = peak[i - 1];
                                }else{
                                if (Proteins[i] < 61  ) {
                                peak[i] = peak[i - 1] + 1;
                                height = 61;
                                } else{peak[i] = peak[i-1]; }
                                }
                                }
                                return peak;
                                }
                                '))
 #error capital X vs lower case x
 breakOut <- parLapply(cl2, fiile, fun  =  Loops2)
 
 closeAllConnections()
 return(breakOut)      }





CAB <- cmpfun(AnalyzeBio)

#Function
extractFileName <- function(output2) {      ###Pull information from the file name

 output  <-  lapply(output2, cleanName, pattern  = "0000GRunNum...*")
 output  <-  lapply(output, cleanName, pattern  = "t*.*DelTime")
 output  <-  lapply(output, cleanName, pattern  = "0*.csv")
 output  <-  lapply(output, cleanName, pattern  = "e")
 output  <-  lapply(output, cleanName, pattern  = "Varienc...*")
 output  <-  lapply(output, cleanName, pattern  = "=Varienc...*")
   output  <-  lapply(output, cleanName, pattern  = "=Varianc...*")
 output  <-  lapply(output, cleanName, pattern  = "Varince...*")
 output  <-  lapply(output, cleanName, pattern  = "Varianc...*")
 output  <-  lapply(output, cleanName, pattern  = "00.*")
 return(output)
}#Pull information from the file name to make it useful


#Function
listToNumFrame <- function(output, output2){### Turn things into an Data Frame and make them numeric + add CV
 output3 <- t(mapply(rbind, output2))
 dfo3  <-  as.data.frame(output3)
 dfo  <-  t(mapply(rbind, output))
 dfo  <-  as.data.frame(dfo)
 CVAmp <-
  as.numeric(as.vector(dfo$V5)) / as.numeric(as.vector(dfo$V4))
 CVPeriod <-
  as.numeric(as.vector(dfo$V3)) / as.numeric(as.vector(dfo$V2))
 dfo2 <- cbind(as.vector(dfo$V1), CVAmp, CVPeriod,as.vector(dfo3$V1))
 dfoNA <- dfo[is.na(dfo$V3),]
 dfo$V1  <-  as.numeric(as.vector(dfo$V1))
 dfo  <-  dfo[order(dfo$V1),]
 dfo <- dfo[!is.na(dfo$V2), ]
 return(dfo)
}### Turn things into an Data Frame and make them numeric + add CV



#Function
dfoNAN  <- function(output, output2){### Turn things into an Data Frame and make them numeric + add CV
  output3 <- t(mapply(rbind, output2))
  dfo3  <-  as.data.frame(output3)
  dfo  <-  t(mapply(rbind, output))
  dfo  <-  as.data.frame(dfo)

  dfoNA <- dfo[is.na(dfo$V3),]

  return(dfoNA)
 }### Turn things into an Data Frame and make them numeric + add CV



#Function
dfo2FNC <- function(dfo, output, output2){
 output3 <- t(mapply(rbind, output2))
 dfo3  <-  as.data.frame(output3)
 CVAmp <-
  as.numeric(as.vector(dfo$V5)) / as.numeric(as.vector(dfo$V4))
 CVPeriod <-
  as.numeric(as.vector(dfo$V3)) / as.numeric(as.vector(dfo$V2))
 dfo2 <- cbind(as.vector(dfo$V1), CVAmp, CVPeriod,as.vector(dfo3$V1))

}




#Function
Plots <- function(output2){
 i <- 1
      ##Make sure the libraries are available
 if ( length(dir("/scratch/Infor"))) {
  library('stringr', lib.loc  =  "/scratch/Infor" )
  #library(dotCall64, lib.loc  =  "/scratch/Infor" )
  library("Rcpp", lib.loc = "/scratch/Infor")
  #library(tidyverse, lib.loc  = "/home/infor/Infor" )

 } else{#Make sure the libraries are available
  library("stringr")
  library("Rcpp")
 }#Make sure the libraries are available

 output <- extractFileName(output2)
 dfo <- listToNumFrame(output,output2) #last number is second largest peak


 dfo2 <- dfo2FNC(dfo, output,output2)
 dfoNA <- dfoNAN(output,output2)

  if(i) {#Find the mean and the variance within cells and accross cells
  mmPeriod  <-  aggregate(as.numeric(as.vector(dfo$V2)), by  =  list(dfo$V1), mean, na.rm = TRUE)
  msPeriod  <-  aggregate(as.numeric(as.vector(dfo$V3)), by  =  list(dfo$V1), mean, na.rm = TRUE)
  smPeriod <- aggregate(as.numeric(as.vector(dfo$V2)), by  =  list(dfo$V1), sd)
  ssPeriod  <-  aggregate(as.numeric(as.vector(dfo$V3)), by  =  list(dfo$V1), sd, na.rm = FALSE)
  mmAmp  <- aggregate(as.numeric(as.vector(dfo$V4)), by  =  list(dfo$V1), mean, na.rm = TRUE)
  msAmp  <- aggregate(as.numeric(as.vector(dfo$V5)), by  =  list(dfo$V1), mean, na.rm = TRUE)
  smAmp  <- aggregate(as.numeric(as.vector(dfo$V4)), by  =  list(dfo$V1), sd, na.rm = TRUE)
  ssAmp  <- aggregate(as.numeric(as.vector(dfo$V5)), by  =  list(dfo$V1), sd, na.rm = TRUE)
  
  #CVAmp <- as.vector(dfo$V5)/as.vector(dfo$V4)
  mcAmp  <- aggregate(as.numeric(as.vector(dfo2[,2])), by  =  list(dfo2[,1]), mean)
  mcPeriod  <-  aggregate(as.numeric(as.vector(dfo2[,3])), by  =  list(dfo2[,1]), mean)

  
    Everest<- aggregate(as.numeric(as.vector(dfo$V10)), by  =  list(dfo$V1), mean, na.rm = TRUE) #second largest peak

  tryCatch	(mmCountG  <-  aggregate(as.numeric(as.vector(dfo$V2)), by  =  list(dfo$V1), length))
  if(dim(dfoNA)[1]){tryCatch(mmCountB  <-  aggregate(dfoNA$V2, by  =  list(dfoNA$V1), length))}
 }  #Find the mean and the variance within cells and accross cells

       ###Set the directory and write the pdfs.
 ifelse( length(dir("/scratch/Infor")) , setwd("/scratch/Infor"),setwd("//home/david/BioResearch"))
 #runsN  <-  runsN
 runsN <- mmCountG[2,2]

 if(i){##Plot}
 pdf(paste(output2[[1]][1],"DelayedOscilationGAMMA,",runsN,"runs.pdf",sep  =  ""))

 plot(dfo2[,1],dfo2[,2],type  =  "p", xlab  =  "Gamma Function Variance ", ylab  =  paste("Coeficient of Variation(",runsN,")",sep = ""), main  =  " Coefficient of Variation Amp ")
 plot(dfo2[,1],dfo2[,3],type  =  "p", xlab  =  "Gamma Function Variance ", ylab  =  paste("CV of(",runsN,") runs",sep = ""), main  =  " Coefficient of Variation Period")
 #  CVPeriod <- dfo$V3/dfo$V2
 #       plot(ssAmp[,1],CVAmp,type  =  "p", xlab  =  "Gamma Function Variance ", ylab  =  paste("Std(",runsN,"-",sep = ""), main  =  "Average Coeficient of Variation Amp ")
 #        plot(smAmp[,1],CVPeriod,type  =  "p", xlab  =  "Gamma Function Variance ", ylab  =  paste("Std(",runsN," ",sep = ""), main  =  "Average Coeficient of Variation Period")

 #plot(ssAmp[,1],ssAmp[,2],type  =  "p", xlab  =  "Gamma Function Variance ", ylab  =  paste("StanDev(",runsN,"-run StanDev(amplitude in a run))",sep = ""), main  =  "StanDev in amplitude  StanDev ")
 #plot(smAmp[,1],smAmp[,2],type  =  "p", xlab  =  "Gamma Function Variance ", ylab  =  paste("StanDev(",runsN,"-run Mean(amplitude in a run))",sep = ""), main  =  "StanDev of amplitude mean within a run")
 
 plot(Everest[,1],Everest[,2],type  =  "p", xlab  =  "Gamma Function Variance ", ylab  =  paste("Mean(",runsN,"-run SecondLargest(amplitude in a run))",sep = ""), main  =  "Second Largest Amplitude ")
 
 plot(mmAmp[,1],mmAmp[,2],type  =  "p", xlab  =  "Gamma Function Variance ", ylab  =  paste("Mean(",runsN,"-run Mean(amplitude in a run))",sep = ""), main  =  "mean amplitude ")
 plot(msAmp[,1],msAmp[,2],type  =  "p", xlab  =  "Gamma Function Variance ", ylab  =  paste("Mean(",runsN,"-run StanDev(amplitude in a run))",sep = ""), main  =  "mean of amplitude StanDev within a run")
 
 plot(mmPeriod[,1],mmPeriod[,2],type  =  "p", xlab  =  "Gamma Function Variance ", ylab  =  paste("Mean(",runsN,"-run Mean(period in a run))",sep = ""), main  =  "mean period across runs")
 plot(msPeriod[,1],msPeriod[,2],type  =  "p", xlab  =  "Gamma Function Variance ", ylab  =  paste("Mean(",runsN,"-run StanDev (period in a run))",sep = ""), main  =  "mean of (period StanDev within a run)")
 
 # plot(meanMaxPeriod[,1],meanMaxPeriod[,2],type  =  "p", xlab  =  "Gamma Function Variance ", ylab  =  paste("Mean(",runsN,"-run Max distance between 2 peaks (period in a run))",sep = ""), main  =  "mean of (Max time between peaks within a run)")
 
  plot(msPeriod[,1],msPeriod[,2],type  =  "p", xlab  =  "Gamma Function Variance ", ylab  =  paste("Mean(",runsN,"-run StanDev (period in a run))",sep = ""), main  =  "mean of (period StanDev within a run)")
  
  
 #plot(smPeriod[,1],smPeriod[,2],type  =  "p", xlab  =  "Gamma Function Variance ", ylab  =  paste("StanDev(",runsN,"-run Mean (period in a run))",sep = ""), main  =  "StanDev of (period mean in a run)")
 #plot(ssPeriod[,1],ssPeriod[,2],type  =  "p", xlab  =  "Gamma Function Variance ", ylab  =  paste("StanDev(",runsN,"-run StanDev (period in a run))",sep = ""), main  =  "StanDev of (period StanDev in a run)")

 tryCatch(plot(mmCountG[,1],mmCountG[,2],type = "p", xlab = "Variance of delay ", ylab = "Total Runs in the data" , main = "Total Runs in the data") )

 if(dim(dfoNA)[1]){	tryCatch(plot(mmCountB[,1],mmCountB[,2],type = "p", xlab = "Max distance of  delay from 1 ", ylab = "NA's in data" , main = "unusable Data"))}

 #closeAllConnections()

 #closeAllConnections()
 }#first plot
 dev.off()

 if(TRUE){###LogPlots}

  pdf(paste("LogScalDelayedOscilationGamma32,", runsN, "runs.pdf", sep = ""))
  #plot(ssAmp[,1],log(ssAmp[,2]),type = "p", xlab = "Gamma Function Variance ", ylab = paste("log(Variance(",runsN,"-run Variance(amplitude in a run)))",sep=""), main = "variance in amplitude  variance ")
  plot(
   smAmp[, 1],
   log(smAmp[, 2]),
   type = "p",
   xlab = "Gamma Function Variance ",
   ylab = paste("log Variance(", runsN, "-run Mean(amplitude in a run))", sep =
                 ""),
   main = "variance of amplitude mean within a run"
  )
 # plot(    svAmp[, 1],   log(svAmp[, 2]),   type = "p",   xlab = "Gamma Function Variance ",   ylab = paste("log Mean(", runsN, "-run Variance(amplitude in a run))", #sep =                 ""),    main = "mean of amplitude variance  within a run"  )
  plot(
   mmAmp[, 1],
   log(mmAmp[, 2]),
   type = "p",
   xlab = "Gamma Function Variance ",
   ylab = paste("log Mean(", runsN, "-run Mean(amplitude in a run))", sep =
                 ""),
   main = "mean amplitude "
  )
  plot(
   mmPeriod[, 1],
   log(mmPeriod[, 2]),
   type = "p",
   xlab = "Gamma Function Variance ",
   ylab = paste("log Mean(", runsN, "-run Mean(period in a run))", sep = ""),
   main = "mean period across runs"
  )
  #plot(   svPeriod[, 1],   log(svPeriod[, 2]),   type = "p",   xlab = "Gamma Function Variance ",   ylab = paste("log Mean(", runsN, "-run Variance (period in a run))", sep =    ""),   main = "mean of (period variance  within a run)"  )
  plot(
   smPeriod[, 1],
   log(smPeriod[, 2]),
   type = "p",
   xlab = "Gamma Function Variance ",
   ylab = paste("log variance(", runsN, "-run Mean (period in a run))", sep =
                 ""),
   main = "variance of (period mean in a run)"
  )
  # plot(ssPeriod[,1],log(ssPeriod[,2]),type = "p", xlab = "Gamma Function Variance ", ylab = paste("log(variance(",runsN,"-run variance (period in a run))",sep=""), main = "variance of (period variance  in a run)")
  #closeAllConnections()
 }   ###Log Plots
 dev.off()
}
#pdfPlot <- cmpfun(Plots)




}#Initialize the system

ifelse((length(dir("/scratch/Infor"))) ,
       result <- tryCatch( output6 <- CAB(durs <- "/tmp/Infort/Gamma2") )
       ,
       result <-  tryCatch(   output6 <-    CAB(durs <-         "/tmp/Infort/Gamma2")  ))

   closeAllConnections()
   ifelse( length(dir("/scratch/Infor")),  setwd("/scratch/Infor"), setwd("/home/david/BioResearch/tmp"))
   save.image(paste(output6[[1]][1],"Gamma4day.RData"))
   output7<-lapply(output6,str_remove,pattern =",")
      output7<-lapply(output7,str_remove,pattern =",")
   head(output7)
   
   rnam<-c("File Name",	"Period Mean",	"Standard deviation Period",	"Mean Amplitude",	"Standard Deviation Amplitude",	"Variance in delay","CV of Delay",	"Filler for backwards compatibility","Filler for backwards compatibility",	"Second Highest Peak")
   fileN<-paste("stats1RA",output6[[1]][1],sep="")
write.csv(output7,file = fileN,row.names = rnam   )
rotate<-read.csv(fileN)
write.csv(file = paste("stats2","RA32",output6[[1]][1],sep=""),t(rotate))
#    graphPlease<-read.csv(paste("stats2",output6[[1]][1],sep=""))
#    graphPlease<- as.data.frame(graphPlease)
 #   graphPlease <- graphPlease[2:dim(graphPlease)[1],]
  #  Please<- graphPlease[,1:3]
   # Please[,1]<- graphPlease[,10]
  #  Please[,2]<-graphPlease[,4]/graphPlease[,3]
  #  Please[,2]<-graphPlease[,6]/graphPlease[,5]


#pdf(paste(output6[[1]][1],"MeanTimeDelayVsPeriodAndAmp.pdf",sep  =  ""))
#plot(Please[,1],graphPlease[,3],type  =  "p", xlab  =  "Gamma Function Mean Delay ", ylab  =  paste("Mean Period"," ",sep = ""), main  =  " Mean Period ")
#plot(Please[,1],graphPlease[,4],type  =  "p", xlab  =  "Gamma Function Mean Delay ", ylab  =  paste("StDev Period"," ",sep = ""), main  =  " StDev Period ")
#plot(Please[,1],graphPlease[,5],type  =  "p", xlab  =  "Gamma Function Mean Delay ", ylab  =  paste("Mean AMP"), main  =  " Mean Amp")
#plot(Please[,1],graphPlease[,6],type  =  "p", xlab  =  "Gamma Function Mean Delay ", ylab  =  paste("StDev AMP"), main  =  " StDev Amp")

#plot(Please[,1],Please[,2],type  =  "p", xlab  =  "Gamma Function Mean Delay ", ylab  =  paste("Coeficient of VariationPeriod"," ",sep = ""), main  =  " Coefficient of Varieation Period ")
#plot(Please[,1],Please[,3],type  =  "p", xlab  =  "Gamma Function Mean Delay ", ylab  =  paste("CV AMP"), main  =  " Coefficient of Varieation Amp")
#dev.off()
#pdfPlot <- cmpfun(Plots); tryCatch(pdfPlot(output2 = output6))

   #q("yes")

 closeAllConnections()
