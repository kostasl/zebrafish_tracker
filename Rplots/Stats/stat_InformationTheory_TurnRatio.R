## Calculate Mutual Information Between turn ratio and Distance to Prey

library(tools)
library(RColorBrewer);


source("config_lib.R")



#lEyeMotionDat <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData_SetC",".rds",sep="")) #
lFirstBoutPoints <-readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_Validated",".rds",sep="")) 



### Capture Speed vs Distance to prey ###
datCapture_NL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$NL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$NL[,"CaptureSpeed"],Undershoot=lFirstBoutPoints$NL[,"Turn"]/lFirstBoutPoints$NL[,"OnSetAngleToPrey"],RegistarIdx=lFirstBoutPoints$NL[,"RegistarIdx"],Validated= lFirstBoutPoints$NL[,"Validated"] ) )
datCapture_LL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$LL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$LL[,"CaptureSpeed"]),Undershoot=lFirstBoutPoints$LL[,"Turn"]/lFirstBoutPoints$LL[,"OnSetAngleToPrey"],RegistarIdx=lFirstBoutPoints$LL[,"RegistarIdx"],Validated= lFirstBoutPoints$LL[,"Validated"] )
datCapture_DL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$DL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$DL[,"CaptureSpeed"]),Undershoot=lFirstBoutPoints$DL[,"Turn"]/lFirstBoutPoints$DL[,"OnSetAngleToPrey"],RegistarIdx=lFirstBoutPoints$DL[,"RegistarIdx"],Validated= lFirstBoutPoints$DL[,"Validated"] )
##Select Validated Only
datCapture_NL <- datCapture_NL[datCapture_NL$Validated == 1, ]
datCapture_LL <- datCapture_LL[datCapture_LL$Validated == 1, ]
datCapture_DL <- datCapture_DL[datCapture_DL$Validated == 1, ]





##
## Calc Info In single Sample & Hunt Event
#InfoCalc <- function(DistRange,minDist,maxDist,Ulist)
  
  DistRange  <- seq(0,1,0.1)
  SpeedRange <- seq(0,100,1) ##We limit The information Obtained To Reasonable Ranges Of Phi (Vergence Angle)
  
  
  #Generates a list with all possible pairs of in the range of Phi And Distance X 
  Grid <- expand.grid(SpeedRange,DistRange)
  
  ### Tally 2D data points in the grid
  
  
  ## Normalize Grid To Sum of Points to obtain Joint Prob P(Capt Speed| Distance)
  
 
  PVec=rep(0,NROW(Grid)) 
  ##Calc Density for all input Space
  for (i in 1:NROW(Grid) )
  {
    PVec[i] <- phiDens(Grid[i,1],Grid[i,2],Ulist)
  }
  ##Normalize to Probability
  PVec=PVec/sum(PVec)
  
  # Convert Pvec to a matrix / 
  PMatrix=matrix(PVec,nrow=length(PhiRange),ncol=length(DistRange),byrow = FALSE)
  
  ##Image shows a transpose Of the PMatrix, - In reality Rows contain Phi, and Cols are X 
  #image(t(PMatrix),y=PhiRange,x=DistRange)
  MargVec=rowSums(PMatrix) ### Marginalize across X to obtain P(Response/Phi)
  
  Iloc=PMatrix/MargVec*length(DistRange) ##Information On Local x/For Each X - Assume X is unif. and so Prob[X]=1/Length(X)
  
  ###row sum
  sel=PMatrix>0
  #INFO=sums(PMatrix[sel]*log2(Iloc[sel]) )
  ### Return Marginals I_xPhi
  mInfo <- PMatrix*log2(Iloc)
  
  ##Calc Marginals only on relevant Region - up to min distance from Target
  ##Make Binary Array indicating OutOf Bound regions
  vIntRegion<-as.numeric( DistRange >= minDist)
  mIntRegion <- t(matrix(vIntRegion,dim(mInfo)[2],dim(mInfo)[1]))
  
  INFO=colSums(mInfo*mIntRegion,na.rm=TRUE )
  
#  return(INFO)
#}
