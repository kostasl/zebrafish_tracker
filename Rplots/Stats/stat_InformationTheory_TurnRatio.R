## Calculate Mutual Information Between turn ratio and Distance to Prey

library(tools)
library(RColorBrewer);


source("config_lib.R")



##
## Calc Info In single Sample & Hunt Event
InfoCalc_get2DFreq <- function(datX,datY,XRange,YRange)
{
  
  ### Tally 2D data points in the grid
  nXbins <- 16 ##ie 2^3 - x bit encoding of Distance
  nYbins <- 16 ##ie 2^3 - x bit encoding of Distance
  
  x.bin <- seq(floor(min(XRange)), (max(XRange)), length=nXbins)
  y.bin <- seq(floor(min(yRange)), (max(YRange)), length=nYbins)
  
  
  freq <-  as.data.frame(table(X=findInterval(datX, x.bin),Y=findInterval(datY, (y.bin) )))
  freq[,1] <- as.numeric(as.character(freq[,1]))
  freq[,2] <- as.numeric(as.character(freq[,2]))
  
  ##Place frequencies on 2D matrix 
  freq2D <- diag(nrow=nYbins,ncol=nXbins)*0
  ##Reverse the Y, so matrix pos matches speed labels
  freq2D[cbind(as.numeric( as.character(freq[,2]) ), as.numeric(as.character(freq[,1] )  ) ) ] <- freq[,3]
  
  ##Draw Matrix
  par(mfrow=c(1,3),pty="s")
  plot(datX,datY,xlim=XRange,ylim=YRange)
  image(x.bin, (y.bin), t(freq2D), col=topo.colors(max(freq2D)+1))
  contour(x.bin, y.bin, t(freq2D), add=TRUE,nlevels=5, col=rgb(1,1,1,.7))
  
  palette(rainbow(max(freq2D)))
  #cols <- (freq2D[-1,-1] + freq2D[-1,-(nYbins-1)] + freq2D[-(nXbins-1),-(nYbins-1)] + freq2D[-(nXbins-1),-1])/4
  persp(freq2D ,theta=45) #col=cols
  
  
  
  return(freq2D)
}##
##Debug ##
datX = datCapture_DL$Undershoot
datY = datCapture_DL$CaptureSpeed
XRange <- c(0,2)
YRange <- c(0,70)


##Shannon Entropy
H_entropy<- function(in_freq) 
{
  in_freq_norm <- in_freq/sum(in_freq)
  in_freq_norm = in_freq_norm[in_freq_norm >0 ]
  
  H = -sum(in_freq_norm* log2(in_freq_norm))  
  
  return(H)
}
## Normalize Grid To Sum of Points to obtain Joint Prob P(Capt Speed| Distance)

## Calc Empririca Marginal Entropies and Mutual information in 2D frequency/Count matrix - 
##Measures average reduction in uncertainty about X that results from learning Y, or vice-versa - information of X on Y
## I(X;Y) = I(Y;X)
calcMIEntropy <- function(freqM)
{
  H_Y <- H_entropy(rowSums(freqM) )
  H_X <- H_entropy(colSums(freqM))
  H_XY <- H_entropy(freqM)
  ## Calc Mutual Info
  H_X + H_Y - H_XY
  ##Conditional Entropy is related to Joint Entropy
  H_XCY <- H_XY - H_Y
  ## Calc Mutual Info between X and Y
  I_XY <- H_X -  H_XCY
  
  return (list(H_X=H_X, H_Y=H_Y, H_XY=H_XY, H_XGivenY=H_XCY, MutualInf_XY=I_XY) )
}


#lEyeMotionDat <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData_SetC",".rds",sep="")) #
lFirstBoutPoints <-readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_Validated",".rds",sep="")) 


#### Load  hunting stats- Generated in main_GenerateMSFigures.r - now including the cluster classification -
datCapture_NL <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_NL_clustered",".rds",sep="")) 
datCapture_LL <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_LL_clustered",".rds",sep="")) 
datCapture_DL <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_DL_clustered",".rds",sep="")) 

datCapture_NL <- datCapture_NL[datCapture_NL$Cluster == "fast",]
datCapture_LL <- datCapture_LL[datCapture_LL$Cluster == "fast",]
datCapture_DL <- datCapture_DL[datCapture_DL$Cluster == "fast",]

### Capture Speed vs Distance to prey ###
#datCapture_NL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$NL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$NL[,"CaptureSpeed"],Undershoot=lFirstBoutPoints$NL[,"Turn"]/lFirstBoutPoints$NL[,"OnSetAngleToPrey"],RegistarIdx=lFirstBoutPoints$NL[,"RegistarIdx"],Validated= lFirstBoutPoints$NL[,"Validated"] ) )
#datCapture_LL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$LL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$LL[,"CaptureSpeed"]),Undershoot=lFirstBoutPoints$LL[,"Turn"]/lFirstBoutPoints$LL[,"OnSetAngleToPrey"],RegistarIdx=lFirstBoutPoints$LL[,"RegistarIdx"],Validated= lFirstBoutPoints$LL[,"Validated"] )
#datCapture_DL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$DL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$DL[,"CaptureSpeed"]),Undershoot=lFirstBoutPoints$DL[,"Turn"]/lFirstBoutPoints$DL[,"OnSetAngleToPrey"],RegistarIdx=lFirstBoutPoints$DL[,"RegistarIdx"],Validated= lFirstBoutPoints$DL[,"Validated"] )
##Select Validated Only
#datCapture_NL <- datCapture_NL[datCapture_NL$Validated == 1, ]
#datCapture_LL <- datCapture_LL[datCapture_LL$Validated == 1, ]
#datCapture_DL <- datCapture_DL[datCapture_DL$Validated == 1, ]


XRange_NL  <- range(datCapture_NL$Undershoot) #seq(0,2,0.2)
YRange_NL <- range(datCapture_NL$CaptureSpeed) ##We limit The information Obtained To Reasonable Ranges Of Phi (Vergence Angle)

XRange_LL  <- range(datCapture_LL$Undershoot) #seq(0,2,0.2)
YRange_LL <- range(datCapture_LL$CaptureSpeed) ##We limit The information Obtained To Reasonable Ranges Of Phi (Vergence Angle)

XRange_DL  <- range(datCapture_DL$Undershoot) #seq(0,2,0.2)
YRange_DL <- range(datCapture_DL$CaptureSpeed) ##We limit The information Obtained To Reasonable Ranges Of Phi (Vergence Angle)


XRange  <- seq(0,2,0.2) #
YRange <- seq(0,60,5) ##We limit The information Obtained To Reasonable Ranges Of Phi (Vergence Angle)

freqM_DF <- InfoCalc_get2DFreq(datCapture_DL$Undershoot,datCapture_DL$CaptureSpeed,XRange,YRange)
freqM_LF <- InfoCalc_get2DFreq(datCapture_LL$Undershoot,datCapture_LL$CaptureSpeed,XRange,YRange)
freqM_NF <- InfoCalc_get2DFreq(datCapture_NL$Undershoot,datCapture_NL$CaptureSpeed,XRange,YRange)

calcMIEntropy(freqM_NF)
calcMIEntropy(freqM_LF)
calcMIEntropy(freqM_DF)

##Correlation Undershoot to Speed
cancor(datCapture_NL$Undershoot,datCapture_NL$CaptureSpeed)
cancor(datCapture_LL$Undershoot,datCapture_LL$CaptureSpeed)
cancor(datCapture_DL$Undershoot,datCapture_DL$CaptureSpeed)

sd(datCapture_NL$Undershoot)
sd(datCapture_LL$Undershoot)
sd(datCapture_DL$Undershoot)

sd(datCapture_NL$CaptureSpeed)
sd(datCapture_LL$CaptureSpeed)
sd(datCapture_DL$CaptureSpeed)

#### ### Speed Vs Distance

XRange  <- c(0,0.8) #
YRange <- x(0,60) ##We limit The information Obtained To Reasonable Ranges Of Phi (Vergence Angle)

freqM_DF <- InfoCalc_get2DFreq(datCapture_DL$DistanceToPrey,datCapture_DL$CaptureSpeed,XRange,YRange)
freqM_LF <- InfoCalc_get2DFreq(datCapture_LL$DistanceToPrey,datCapture_LL$CaptureSpeed,XRange,YRange)
freqM_NF <- InfoCalc_get2DFreq(datCapture_NL$DistanceToPrey,datCapture_NL$CaptureSpeed,XRange,YRange)

calcMIEntropy(freqM_NF)
calcMIEntropy(freqM_LF)
calcMIEntropy(freqM_DF)

##Correlation Speed to Distance
cancor(datCapture_NL$DistanceToPrey,datCapture_NL$CaptureSpeed)
cancor(datCapture_LL$DistanceToPrey,datCapture_LL$CaptureSpeed)
cancor(datCapture_DL$DistanceToPrey,datCapture_DL$CaptureSpeed)


##redundancy
#1-H_X/2^3
#1-H_Y/2^3

##Verify
library(entropy)

mi.empirical(freqM_NF,unit="log2" ) 
mi.empirical(freqM_LF,unit="log2" ) 
mi.empirical(freqM_DF,unit="log2" )  

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
