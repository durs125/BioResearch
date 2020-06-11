
#Set the number of processors in the cluster, the initialization vector, and the number of runs vector
# install.packages("Rcpp", "~/R/x86_64-redhat-linux-gnu-library/3.6");  install.packages("Rcpp11", "~/R/x86_64-redhat-linux-gnu-library/3.6") ; install.packages("foreach", lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.6") install.packages("foreach", "/scratch/Infor"); install.packages("snow","/scratch/Infor")
   #install.packages("iterators","/scratch/Infor"); install.packages("doSNOW","/scratch/Infor");install.packages("snow","/scratch/Infor");install.packages("Rcpp","/scratch/Infor"); install.packages("stringr","/scratch/Infor");install.packages("Rcpp11","/scratch/Infor"); 

#   
#   foreach
# install.packages("snow","/scratch/Infor")install.packages("foreach","/scratch/Infor")
# between
 #initialization <- c(.1,60,80.5 , 1.45,.05,1.55)
 #initialization <- c(1,5,10.5 , 1.49,.01,1.51)
 # double startVar,  double  stepsize,   double maxVar,  double startzz,  double addzz, const double maxzz  )

# Do Fourier analysis, make sure thinggs go to zero, do more long runs, see if you can decouple the period and the delay time
 pLoopEXP <- function(sampleSize, initialization) {
   stepsize <- initialization[1]
   startVar <- initialization[2]
   maxVar <- initialization[3]
   startzz <- initialization[4]
   addzz <- initialization[5]
   maxzz <- initialization[6]
   setwd("..")
   setwd("..")
     setwd("..")
ifelse( length(dir("/scratch/Infor")), setwd("/scratch/Infor"),setwd("/home/david/BioResearch/"))
ifelse( length(dir("/scratch/Infor")), setwd("/scratch/Infor"),setwd("/home/david/BioResearch/"))
ifelse( length(dir("/scratch/Infor")), sourceCpp(paste( getwd(),'/NEnz7GammaD.cpp',sep = "")),sourceCpp(paste( getwd(),'/CPP/NEnz7GammaD.cpp',sep = "")))
ifelse( length(dir("/scratch/Infor")),   clusterEvalQ(cl2,Rcpp::sourceCpp(paste( getwd(),'/NEnz7GammaD.cpp',sep = ""))),  clusterEvalQ(cl2,Rcpp::sourceCpp(paste( getwd(),'/CPP/NEnz7GammaD.cpp',sep = ""))))

  # clusterEvalQ(cl2,Rcpp::sourceCpp(paste( getwd(),'/DF7ParallelBernouli.cpp',sep = "")))
  #foreach(i = 1:length(sampleSize) , .noexport = c("MakeGamaFiles")) %dopar%{
  foreach(i = 1:sampleSize[2] ) %dopar% {
  #foreach(i = 2:3 ) %dopar%{
   #MakeGamaFiles(number files needed,stepsize variance,startVariance, startzz, addzz, maxzz );
   #zz is multiple of yr
   # MakeGammaFiles2( NumberRuns, double stepsize,double startVar,  double startzz,  double addzz,  double maxzz )
   #zz is multiple of CV
   MakeGammaFiles2(sampleSize[1] + i,stepsize,startVar, maxVar ,startzz,addzz,maxzz)
   # MakeBernFiles(sampleSize[i])
  }
 }#end pLoop
 #end ploopExp
 
 qLoopEXP <- function(sampleSize, initialization) {
  stepsize <- initialization[1]
  startVar <- initialization[2]
  maxVar <- initialization[3]
  startzz <- initialization[4]
  addzz <- initialization[5]
  maxzz <- initialization[6]

  clusterEvalQ(cl2,Rcpp::sourceCpp(paste( getwd(),'/NEnz7GammaD.cpp',sep = "")))
  # clusterEvalQ(cl2,Rcpp::sourceCpp(paste( getwd(),'/DF7ParallelBernouli.cpp',sep = "")))
  #foreach(i = 1:length(sampleSize) , .noexport = c("MakeGamaFiles")) %dopar%{

 foreach(i = 0:sampleSize[2] ) %dopar% {
   #foreach(i = 2:3 ) %dopar%{
   #MakeGamaFiles(number files needed,stepsize variance,startVariance, startzz, addzz, maxzz );
   #zz is multiple of yr
  #  MakeGammaFiles2( NumberRuns, double stepsize,double startVar,  double startzz,  double addzz,  double maxzz )
   #zz is multiple of CV
   MakeGammaFiles2(sampleSize[1]+i,stepsize,startVar, maxVar ,startzz,addzz,maxzz)
   # MakeBernFiles(sampleSize[i])
  }
 }# End qLoop
 # End qLoop

 if ( length(dir("/scratch/Infor"))) {
  ifelse(length(dir("/tmp/Infort/Gamma2")), closeAllConnections(), dir.create("/tmp/Infort"))
  ifelse(length(dir("/tmp/Infort/Gamma2")), closeAllConnections(), dir.create("/tmp/Infort/Gamma2"))

  sampleSize <- c(1,40)
  closeAllConnections()
  library(compiler)
  library(parallel)
  library("foreach", lib.loc = "/scratch/Infor")
  library("snow", lib.loc = "/scratch/Infor")
  library("iterators", lib.loc = "/scratch/Infor")
  library("doSNOW", lib.loc = "/scratch/Infor")
  library("Rcpp", lib.loc = "/scratch/Infor")
  library("Rcpp11", lib.loc = "/home/infor/R/x86_64-redhat-linux-gnu-library/3.6")
 
 
 
  library(tools)
   b<-tryCatch(unlink("/tmp/Infort/Gamma2",recursive = TRUE)) # prevents the files from being reused in multiple runs

setwd("~/")
    setwd("/scratch/Infor")
   
   



 }else {
 setwd("~/")
  setwd("/home/david/BioResearch")
  sampleSize <-c(1,4)
  closeAllConnections()
  #dir.create("/tmp/Infort")
   #dir.create("/tmp/Infort/Bern")
   
   
   
   b<-tryCatch(unlink("/tmp/Infort/Gamma2",recursive = TRUE)) # prevents the files from being reused in multiple runs

   
   
   
  dir.create("/tmp/Infort/Gamma2",recursive = TRUE)
 # for (zz in initialization[5]:initialization[7]:initialization[6] {
setwd("~/")
    setwd("/home/david/BioResearch")
  # dir.create(cat("/tmp/Infort/Gamma2/t",zz,sep =""))
  #}
  library(compiler)
  library(parallel)
  library("snow")
  library("iterators")
  library("doSNOW")
  library("foreach" )
  library(tools)

  library("Rcpp")

 
 }



xx<-1
yy<-1


        for(xx in 1:1){
    for(yy in (0:15)){ 
     b<-tryCatch(unlink("/tmp/Infort/Gamma2",recursive = TRUE)) # prevents the files from being reused in multiple runs
     closeAllConnections()
 if ( length(dir("/scratch/Infor"))) {

  ifelse(length(dir("/tmp/Infort/Gamma2")), closeAllConnections(), dir.create("/tmp/Infort/Gamma2", recursive = TRUE))
  sampleSize <- c(1,40)
  closeAllConnections()
 
  cl2 <- makeCluster(18, type = "SOCK") #Change me for number processors on computer -1
  registerDoSNOW(cl2)
  clusterEvalQ(cl2,library("Rcpp", lib.loc = "/scratch/Infor"))
  clusterEvalQ(cl2,.packages(all.available = TRUE, lib.loc = NULL))

  clusterEvalQ(cl2,library("Rcpp11", lib.loc = "/home/infor/R/x86_64-redhat-linux-gnu-library/3.6")) #Do I need this???
  library(tools)




  clusterEvalQ(cl2,setwd("/scratch/Infor"))
  
 clusterEvalQ(cl2,Rcpp::sourceCpp(paste( getwd(),'/NEnz7GammaD.cpp',sep = "")))
 }else{
   dir.create("/tmp/Infort/Gamma2", recursive=TRUE)
 setwd("~/")
  setwd("/home/david/BioResearch")
  sampleSize <-c(1,1)
  closeAllConnections()
  #dir.create("/tmp/Infort")
   #dir.create("/tmp/Infort/Bern")
   
   
   
  
 # for (zz in initialization[5]:initialization[7]:initialization[6] {

  # dir.create(cat("/tmp/Infort/Gamma2/t",zz,sep =""))
  #}



  cl2 <- makeCluster(3, type = "SOCK") #ChangeMe for number processors on computer -1
  registerDoSNOW(cl2)
    clusterEvalQ(cl2,setwd("~/"))
  clusterEvalQ(cl2,setwd("/home/david/BioResearch/"))
  clusterEvalQ(cl2,library("Rcpp"))
  clusterEvalQ(cl2,library("Rcpp11"))
   clusterEvalQ(cl2,Rcpp::sourceCpp(paste( getwd(),'/CPP/NEnz7GammaD.cpp',sep = "")))
 #clusterEvalQ(cl2,Rcpp::sourceCpp(paste( getwd(),'/DF7ParallelBernouli.cpp',sep = "")))
 }

 i <- 1


 #cl2 <- makeCluster(2, type = "SOCK") #Change me for number processors on computer -1
 # registerDoSNOW(cl2)
 # clusterEvalQ(cl2,library("Rcpp"))
 #clusterEvalQ(cl2,Rcpp::sourceCpp('DF7ParallelGama.cpp'))
 # clusterEvalQ(cl2,.packages(all.available = TRUE, lib.loc = NULL))
  #       clusterEvalQ(cl2,setwd("~/"))
   #      clusterEvalQ(cl2,setwd("/home/david/BioResearch/"))



    
    #endfor yy
    
    

  #sampleSize <-c(1,2)
 #  initialization <- c(.1,1,15 , 2.1,.2,2.1)
  #ideal initialization <- c(.01,.2,10.1,.5,.5,8.1 )
  



       startVar<-.5+.5*yy
      # startzz <- .01+.25*xx
  initialization <- c(startVar,2.5,startVar , .05 ,.05,.8) #test run formerly to .9
  # double startVar,  double  stepsize,   double maxVar,  double startzz,  double addzz, const double maxzz  )

    # double startMean Delay,  double  stepsize,   double maxMean Delay,  double start CV Delay,  double CV Delay, const double CV delay  )
 timeTaken <- tryCatch(system.time(pLoopEXP(sampleSize, initialization)))
 timeTaken



 closeAllConnections()
 if (length(dir("/scratch/Infor"))) {
source("/scratch/Infor/R16PeakToPeak.R")
 }else{
 #setwd("/home/david/BioResearch/plotting")
 # source("/home/david/BioResearch/plotting/TimeSeriesPlots5Short.R")
 source("/home/david/BioResearch/plotting/R16PeakToPeak62.R")
  source("/home/david/BioResearch/plotting/R16PeakToPeak32.R")
 closeAllConnections()
}

 closeAllConnections()
 #install.packages("foreach", lib = "/scratch/Infor")
 

}#endfor yy
}#endfor xx
source("/home/david/BioResearch/plotting/summeryStats2.R")

