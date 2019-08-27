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



DistRange  <- seq(0,0.8,0.1)
SpeedRange <- seq(0,70,1) ##We limit The information Obtained To Reasonable Ranges Of Phi (Vergence Angle)

freqM_DF <- InfoCalc_get2DFreq(datCapture_DL$DistanceToPrey,datCapture_DL$CaptureSpeed,DistRange,SpeedRange)
freqM_LF <- InfoCalc_get2DFreq(datCapture_LL$DistanceToPrey,datCapture_LL$CaptureSpeed,DistRange,SpeedRange)
freqM_NF <- InfoCalc_get2DFreq(datCapture_NL$DistanceToPrey,datCapture_NL$CaptureSpeed,DistRange,SpeedRange)


H_Y <- H_entropy(rowSums(freqM_DF) )
H_X <- H_entropy(colSums(freqM_DF))
H_XY <- H_entropy(freqM_DF)
## Calc Mutual Info
H_X + H_Y - H_XY



H_Y <- H_entropy(rowSums(freqM_LF) )
H_X <- H_entropy(colSums(freqM_LF))
H_XY <- H_entropy(freqM_LF)
## Calc Mutual Info
H_X + H_Y - H_XY


##
## Calc Info In single Sample & Hunt Event
InfoCalc_get2DFreq <- function(datX,datY,XRange,yRange)
{
  
  ### Tally 2D data points in the grid
  nbins <- 10
  
  x.bin <- seq(floor(min(XRange)), ceiling(max(XRange)), length=nbins)
  y.bin <- seq(floor(min(yRange)), ceiling(max(yRange)), length=nbins)
  
  
  freq <-  as.data.frame(table(findInterval(datX, x.bin),findInterval(datY, y.bin)))
  freq[,1] <- as.numeric(freq[,1])
  freq[,2] <- as.numeric(freq[,2])
  
  ##Place frequencies on 2D matrix 
  freq2D <- diag(nbins)*0
  freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]

  ##Draw Matrix
  par(mfrow=c(1,2))
  image(x.bin, y.bin, freq2D, col=topo.colors(max(freq2D)))
  contour(x.bin, y.bin, freq2D, add=TRUE, col=rgb(1,1,1,.7))
  
  palette(rainbow(max(freq2D)))
  cols <- (freq2D[-1,-1] + freq2D[-1,-(nbins-1)] + freq2D[-(nbins-1),-(nbins-1)] + freq2D[-(nbins-1),-1])/4
  persp(freq2D, col=cols)
  
  
  
  return(freq2D)
}
 

H_entropy<- function(in_freq) 
{
  in_freq_norm <- in_freq/sum(in_freq)
  in_freq_norm = in_freq_norm[in_freq_norm >0 ]
  
  H = -sum(in_freq_norm* log2(in_freq_norm))  
  
  return(H)
}
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
