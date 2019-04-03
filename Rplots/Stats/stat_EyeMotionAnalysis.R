### Stat Analysis of Eye Movement Data - Distance and duration of hunt events vs Max Eye Vergence Etc ##


source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")


#### Load the Tracked Hunts Register ###
strRegisterDataFileName <- paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_SetC",".rds",sep="") #Processed Registry on which we add 
message(paste(" Importing Retracked HuntEvents from:",strRegisterDataFileName))
datTrackedEventsRegister <- readRDS(strRegisterDataFileName) ## THis is the Processed Register File On 
lMotionBoutDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_SetC.rds",sep="") ) #Processed Registry on which we add )
lEyeMotionDat <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData_SetC",".rds",sep="")) #Processed Registry on which we add )


dev.off()

idx <- 10
plot(lEyeMotionDat[[idx]][,"t"],lEyeMotionDat[[idx]][,"LEyeAngle"]-lEyeMotionDat[[idx]][,"REyeAngle"],type="l")
plot(lEyeMotionDat[[idx]][,"DistToPrey"])
plot(lEyeMotionDat[[idx]][,"t"],lEyeMotionDat[[idx]][,"DistToPrey"])
plot(lEyeMotionDat[[idx]][,"DistToPrey"],lEyeMotionDat[[idx]][,"LEyeAngle"]-lEyeMotionDat[[idx]][,"REyeAngle"] )
## Does max EyeV Correlate with Hunt Duration Or With Distance To Prey Travelled ?
### Create A DAta Struct with EyeV before Strike/Capture - Duration of hunt Event-Eye Vergence - bUt also Distance from Prey travelled

##Get the Eye Vergence before the strike
lMaxEyeV <- list()
for (idx in 1:NROW(lEyeMotionDat) )
{
  vEye <- lEyeMotionDat[[idx]][,"LEyeAngle"]-lEyeMotionDat[[idx]][,"REyeAngle"]
  idx_Start <- min(which ( vEye > G_THRESHUNTVERGENCEANGLE))
  ##Closest To Prey in far time 
  idx_End <- max(which(lEyeMotionDat[[idx]][,"DistToPrey"] < median(lEyeMotionDat[[idx]][,"DistToPrey"]))) ##max(which ( vEye > G_THRESHUNTVERGENCEANGLE))
  maxV <- max(vEye[idx_End-100:idx_End ])
  
  huntDuration_ms <- lEyeMotionDat[[idx]][,"t"][idx_End]-lEyeMotionDat[[idx]][,"t"][idx_Start]
  huntDistance_mm <- lEyeMotionDat[[idx]][,"DistToPreyInit"][idx_End]-lEyeMotionDat[[idx]][,"DistToPrey"][idx_End]
  
  lMaxEyeV[[idx]] <- list(maxEyeV = maxV,
                          duration_ms= huntDuration_ms,
                          distance_mm=huntDistance_mm,
                          RegIdx=lEyeMotionDat[[idx]][,"RegistarIdx"][idx_End],
                          groupID=lEyeMotionDat[[idx]][,"groupID"][idx_End],
                          CaptureStrike=lEyeMotionDat[[idx]][,"doesCaptureStrike"][idx_End] )

}


datMaxEyeV <- data.frame(do.call(rbind,lMaxEyeV))

duration_cor <- cor(unlist(datMaxEyeV$maxEyeV), unlist(datMaxEyeV$duration_ms), method = "pearson") 
duration_cov <- cov(unlist(datMaxEyeV$maxEyeV), unlist(datMaxEyeV$duration_ms))
distance_cor <- cor(unlist(datMaxEyeV$maxEyeV), unlist(datMaxEyeV$distance_mm) , method = "pearson")
distance_cov <- cov(unlist(datMaxEyeV$maxEyeV), unlist(datMaxEyeV$distance_mm)) 


##Does Duration Explain Eye Vergence ? ##
strPlotFileName <- paste(strPlotExportPath,"/stat/EyeVToDuration_corr.pdf",sep="")
pdf(strPlotFileName,width = 16,height = 18 ,paper = "a4",onefile = TRUE );
plot(datMaxEyeV$duration_ms,datMaxEyeV$maxEyeV,ylim=c(0,100),main=paste("P cor=",prettyNum(duration_cor) ),ylab="Max Eye Vergence",xlab="Hunt duration (ms)" )
abline(lm(unlist(datMaxEyeV$maxEyeV)~unlist(datMaxEyeV$duration_ms)), col="red") # regression line (y~x) 
## Robust locally weighted regression is  a method for smoothing by Ezekiel (1941, p. 51). The points are grouped accord- a scatterplot, (xi,yi), i = 1, . . . , n, in which the fitted value at xk ing to xi, and for each group the mean of the yi is plotted is the value of a  polynomial fit to the data using weighted least against the mean of the xi. More recently, Stone (1977) squares,
lines(lowess(unlist(datMaxEyeV$maxEyeV)~unlist(datMaxEyeV$duration_ms) ))
dev.off()

##Distance??
strPlotFileName <- paste(strPlotExportPath,"/stat/EyeVToDistance_corr.pdf",sep="")
pdf(strPlotFileName,width = 16,height = 18 ,paper = "a4",onefile = TRUE );
  plot(datMaxEyeV$distance_mm,datMaxEyeV$maxEyeV,ylim=c(0,100))
  abline(lm(unlist(datMaxEyeV$maxEyeV)~unlist(datMaxEyeV$distance_mm)), col="red") # regression line (y~x) 
  lines(lowess(unlist(datMaxEyeV$maxEyeV)~unlist(datMaxEyeV$distance_mm) ))
dev.off()

##How does Distance to target correlate to time
strPlotFileName <- paste(strPlotExportPath,"/stat/DistanceToDuration_corr.pdf",sep="")
pdf(strPlotFileName,width = 16,height = 18 ,paper = "a4",onefile = TRUE );

  plot(datMaxEyeV$distance_mm,datMaxEyeV$duration_ms,xlab="Distance (mm)",ylab="Hunt duration (ms)" )
  abline(lm(unlist(datMaxEyeV$duration_ms)~unlist(datMaxEyeV$distance_mm)), col="red") # regression line (y~x) 
  lines(lowess(unlist(datMaxEyeV$duration_ms)~unlist(datMaxEyeV$distance_mm) ))

dev.off()




### Other Aux Stuff ##
## Duration Of Hunt Events
strGroupID <- unique(datTrackedEventsRegister[unlist(datMaxEyeV$RegIdx),]$groupID)

boxplot(unlist(datMaxEyeV[datMaxEyeV$groupID==1,]$duration_ms ),
        unlist(datMaxEyeV[datMaxEyeV$groupID==2,]$duration_ms ) ,
        unlist(datMaxEyeV[datMaxEyeV$groupID==3,]$duration_ms ),main="Duration of Hunt Event",names=strGroupID  )

boxplot(unlist(datMaxEyeV[datMaxEyeV$groupID==1,]$distance_mm ),
        unlist(datMaxEyeV[datMaxEyeV$groupID==2,]$distance_mm ) ,
        unlist(datMaxEyeV[datMaxEyeV$groupID==3,]$distance_mm ),main="Distance of Hunt Event",names=strGroupID  )


## Check Efficiency Of Prey Distance Over Hunt Time  ##
vDistOverTime_A <- unlist(datMaxEyeV[datMaxEyeV$groupID==1,]$distance_mm) / unlist(datMaxEyeV[datMaxEyeV$groupID==1,]$duration_ms)
vDistOverTime_B <- unlist(datMaxEyeV[datMaxEyeV$groupID==2,]$distance_mm) / unlist(datMaxEyeV[datMaxEyeV$groupID==2,]$duration_ms)
vDistOverTime_C <- unlist(datMaxEyeV[datMaxEyeV$groupID==3,]$distance_mm) / unlist(datMaxEyeV[datMaxEyeV$groupID==3,]$duration_ms)

layout(matrix(c(1,2,3), 3, 1 ,byrow=TRUE))
ptbreaks <- seq(from=-0.001,to=0.005,by=1/5000)
hist(vDistOverTime_A,xlim=c(0,0.004),breaks=ptbreaks,main=strGroupID[1])
hist(vDistOverTime_B,xlim=c(0,0.004),breaks=ptbreaks,main=strGroupID[2])
hist(vDistOverTime_C,xlim=c(0,0.004),breaks=ptbreaks,main=strGroupID[3])


