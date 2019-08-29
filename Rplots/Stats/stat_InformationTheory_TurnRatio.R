## Calculate Mutual Information Between turn ratio and Distance to Prey

library(tools)
library(RColorBrewer);


source("config_lib.R")



##
## Calc Info In single Sample & Hunt Event
InfoCalc_get2DFreq <- function(datX,datY,XRange,YRange,plot=FALSE)
{
  
  ### Tally 2D data points in the grid
  nXbins <- 16 ##ie 2^3 - x bit encoding of Distance
  nYbins <- 16 ##ie 2^3 - x bit encoding of Distance
  
  x.bin <- seq(floor(min(XRange)), (max(XRange)), length=nXbins)
  y.bin <- seq(floor(min(YRange)), (max(YRange)), length=nYbins)
  
  
  freq <-  as.data.frame(table(X=findInterval(datX, x.bin),Y=findInterval(datY, (y.bin) )))
  freq[,1] <- as.numeric(as.character(freq[,1]))
  freq[,2] <- as.numeric(as.character(freq[,2]))
  
  ##Place frequencies on 2D matrix 
  freq2D <- diag(nrow=nYbins,ncol=nXbins)*0
  ##Reverse the Y, so matrix pos matches speed labels
  freq2D[cbind(as.numeric( as.character(freq[,2]) ), as.numeric(as.character(freq[,1] )  ) ) ] <- freq[,3]
  
  ##Draw Matrix
  if (plot)
  {
    par(mfrow=c(1,3),pty="s")
    plot(datX,datY,xlim=XRange,ylim=YRange)
    image(x.bin, (y.bin), t(freq2D), col=topo.colors(max(freq2D)+1))
    contour(x.bin, y.bin, t(freq2D), add=TRUE,nlevels=5, col=rgb(1,1,1,.7))
    palette(rainbow(max(freq2D)))
    #cols <- (freq2D[-1,-1] + freq2D[-1,-(nYbins-1)] + freq2D[-(nXbins-1),-(nYbins-1)] + freq2D[-(nXbins-1),-1])/4
    persp(freq2D ,theta=45) #col=cols
  }
  
  
  return(freq2D)
}##
##Debug ##
#datX = datCapture_DL$Undershoot
#datY = datCapture_DL$CaptureSpeed


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



##Bootstrap Data analysis from 2 chosen columns of datCapture_NL (lFirstBout)  to get stats on correlations and Mutual Information
bootStrap_stat <- function(datCapture_X,datCapture_Y,N,XRange,YRange)
{
  datCapture <- data.frame(cbind(datCapture_X,datCapture_Y))
  l_sampleXYAnalysis <- list()
  for (i in 1:N)
  {
    
    #freqM_DF <- InfoCalc_get2DFreq(datCapture_DL$Undershoot,datCapture_DL$CaptureSpeed,XRange,YRange)
    #freqM_LF <- InfoCalc_get2DFreq(datCapture_LL$Undershoot,datCapture_LL$CaptureSpeed,XRange,YRange)
    idxSample <- sample(1:NROW(datCapture),size=floor(NROW(datCapture)*0.80))
    datSub <- datCapture[idxSample,]
    freqM_NF <- InfoCalc_get2DFreq(datSub[,1],datSub[,2],XRange,YRange)
    
    infC <- calcMIEntropy(freqM_NF)
    corrXY <- cor(datSub[,1],datSub[,2],method="pearson")
    l_sampleXYAnalysis[[i]] <- data.frame(MI = infC$MutualInf_XY,entropy_X = infC$H_X,entropy_Y = infC$H_Y,corr=corrXY)
    #inf_LF <-calcMIEntropy(freqM_LF)
    #inf_DF <-calcMIEntropy(freqM_DF)
  }
  datXYAnalysis <- do.call(rbind,l_sampleXYAnalysis)
  
  return(datXYAnalysis)
}


#lEyeMotionDat <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData_SetC",".rds",sep="")) #
lFirstBoutPoints <-readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_Validated",".rds",sep="")) 


#### Load  hunting stats- Generated in main_GenerateMSFigures.r - now including the cluster classification -
datCapture_NL <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_NL_clustered",".rds",sep="")) 
datCapture_LL <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_LL_clustered",".rds",sep="")) 
datCapture_DL <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_DL_clustered",".rds",sep="")) 

#datCapture_NL <- datCapture_NL[datCapture_NL$Cluster == "fast",]
#datCapture_LL <- datCapture_LL[datCapture_LL$Cluster == "fast",]
#datCapture_DL <- datCapture_DL[datCapture_DL$Cluster == "fast",]

### Capture Speed vs Distance to prey ###
#datCapture_NL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$NL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$NL[,"CaptureSpeed"],Undershoot=lFirstBoutPoints$NL[,"Turn"]/lFirstBoutPoints$NL[,"OnSetAngleToPrey"],RegistarIdx=lFirstBoutPoints$NL[,"RegistarIdx"],Validated= lFirstBoutPoints$NL[,"Validated"] ) )
#datCapture_LL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$LL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$LL[,"CaptureSpeed"]),Undershoot=lFirstBoutPoints$LL[,"Turn"]/lFirstBoutPoints$LL[,"OnSetAngleToPrey"],RegistarIdx=lFirstBoutPoints$LL[,"RegistarIdx"],Validated= lFirstBoutPoints$LL[,"Validated"] )
#datCapture_DL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$DL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$DL[,"CaptureSpeed"]),Undershoot=lFirstBoutPoints$DL[,"Turn"]/lFirstBoutPoints$DL[,"OnSetAngleToPrey"],RegistarIdx=lFirstBoutPoints$DL[,"RegistarIdx"],Validated= lFirstBoutPoints$DL[,"Validated"] )
##Select Validated Only
#datCapture_NL <- datCapture_NL[datCapture_NL$Validated == 1, ]
#datCapture_LL <- datCapture_LL[datCapture_LL$Validated == 1, ]
#datCapture_DL <- datCapture_DL[datCapture_DL$Validated == 1, ]

# XRange_NL  <- range(datCapture_NL$Undershoot) #seq(0,2,0.2)
# YRange_NL <- range(datCapture_NL$CaptureSpeed) ##We limit The information Obtained To Reasonable Ranges Of Phi (Vergence Angle)
# 
# XRange_LL  <- range(datCapture_LL$Undershoot) #seq(0,2,0.2)
# YRange_LL <- range(datCapture_LL$CaptureSpeed) ##We limit The information Obtained To Reasonable Ranges Of Phi (Vergence Angle)
# 
# XRange_DL  <- range(datCapture_DL$Undershoot) #seq(0,2,0.2)
# YRange_DL <- range(datCapture_DL$CaptureSpeed) ##We limit The information Obtained To Reasonable Ranges Of Phi (Vergence Angle)
# 

XRange  <- c(0,0.8) #
YRange <- c(0,60) ##We limit The information Obtained To Reasonable Ranges Of Phi (Vergence Angle)

stat_Cap_NF <- bootStrap_stat(datCapture_NL$DistanceToPrey,datCapture_NL$CaptureSpeed,1000,XRange,YRange)
stat_Cap_LF <- bootStrap_stat(datCapture_LL$DistanceToPrey,datCapture_LL$CaptureSpeed,1000,XRange,YRange)
stat_Cap_DF <- bootStrap_stat(datCapture_DL$DistanceToPrey,datCapture_DL$CaptureSpeed,1000,XRange,YRange)

## Distance Vs Capture Speed Mututal INformation 
bkSeq <- seq(0,1,0.02)
hist(stat_Cap_NF$MI,xlim=c(0,1),ylim=c(0,300),col=colourL[2],breaks = bkSeq ,
     xlab="Mutual information between capture speed and distance to prey", main="Bootstrapped Mutual information")
hist(stat_Cap_LF$MI,xlim=c(0,1),col=colourL[1],add=TRUE ,breaks = bkSeq)
hist(stat_Cap_DF$MI,xlim=c(0,1),col=colourL[3],add=TRUE,breaks = bkSeq )


##Correlation
bkSeq <- seq(-0.1,0.8,0.02)
hist(stat_Cap_NF$corr,xlim=c(-0.1,0.8),ylim=c(0,300),col=colourL[2],breaks = bkSeq,xlab="Pearson's correlation speed vs distance",main="Bootstraped 0.80" )
hist(stat_Cap_LF$corr,xlim=c(-0.1,0.8),col=colourL[1],add=TRUE ,breaks = bkSeq)
hist(stat_Cap_DF$corr,xlim=c(-0.1,0.8),col=colourL[3],add=TRUE,breaks = bkSeq )

##LF Explores More Speeds - Higher Speed Entropy
bkSeq <- seq(0,4,0.05)
hist(stat_Cap_NF$entropy_Y,xlim=c(1,4),col=colourL[2], breaks=bkSeq,xlab="Capture speed entropy  ",main=NA  )
hist(stat_Cap_LF$entropy_Y,xlim=c(1,4),col=colourL[1],add=TRUE, breaks=bkSeq)
hist(stat_Cap_DF$entropy_Y,xlim=c(1,4),col=colourL[3],add=TRUE, breaks=bkSeq)

##LF Explores More Distances - Higher Prey-Distance Entropy
hist(stat_Cap_NF$entropy_X,xlim=c(1,4),col=colourL[2], breaks=bkSeq,xlab="Distance to prey entropy  " ,main=NA )
hist(stat_Cap_LF$entropy_X,xlim=c(1,4),col=colourL[1], breaks=bkSeq,add=TRUE )
hist(stat_Cap_DF$entropy_X,xlim=c(1,4),col=colourL[3], breaks=bkSeq,add=TRUE )

### DENSITIES ####
pBw <- 0.02

# Plot Speed Vs Distance Correlation - bootstraped Stat ##
strPlotName = paste(strPlotExportPath,"/stat/fig5_statbootstrap_correlation_SpeedVsDistance.pdf",sep="")
pdf(strPlotName,width=7,height=7,title="Correlations In hunt variables",onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))
par(mar = c(3.9,4.7,1,1))

  plot(density(stat_Cap_NF$corr,kernel="gaussian",bw=pBw),
       col=colourLegL[1],xlim=c(0,1),lwd=3,lty=1,ylim=c(0,10),main=NA, xlab=NA,ylab=NA,cex=cex,cex.axis=cex) #expression(paste("slope ",gamma) ) )
  lines(density(stat_Cap_LF$corr,kernel="gaussian",bw=pBw),col=colourLegL[2],lwd=3,lty=2)
  lines(density(stat_Cap_DF$corr,kernel="gaussian",bw=pBw),col=colourLegL[3],lwd=3,lty=3)
  
# legend("topright",         legend=c(  expression (),
#                    bquote(NF~ ''  ),
#                    bquote(LF ~ '' ),
#                    bquote(DF ~ '' )  ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
#         col=colourLegL,lty=c(1,2,3),lwd=3,cex=cex)
  mtext(side = 1,cex=cex,cex.main=cex, line = lineXAxis, expression(paste("Correlation of capture speed to prey distance  ") ))
  mtext(side = 2,cex=cex,cex.main=cex, line = lineAxis, expression("Density function"))

dev.off()





############# BOOTSTRAP 


meanMI <- list(MI_NF=mean(stat_Cap_NF$MI),MI_LF=mean(stat_Cap_LF$MI),MI_DF=mean(stat_Cap_DF$MI))
barplot(c(meanMI$MI_NF,meanMI$MI_LF,meanMI$MI_DF),ylim=c(0,1) )

# library(ggplot)
# dodge <- position_dodge(width = 0.9)
# limits <- aes(ymax = mean(datXYAnalysis$MI) + sd(datXYAnalysis$MI)/sqrt(NROW((datXYAnalysis$MI))),
#               ymin = mean(datXYAnalysis$MI) - 2*sd(datXYAnalysis$MI)/sqrt(NROW((datXYAnalysis$MI))))
# 
# p <- ggplot(data = datXYAnalysis, aes(y = MI ))
XRange  <- c(0,2) #
YRange <- c(0,60) ##We limit The information Obtained To Reasonable Ranges Of Phi (Vergence Angle)

###### UNDERSHOOT TO SPEED
  stat_CapTurnVsSpeed_NF <- bootStrap_stat(datCapture_NL$Undershoot,datCapture_NL$CaptureSpeed,10000,XRange,YRange)
  stat_CapTurnVsSpeed_LF <- bootStrap_stat(datCapture_LL$Undershoot,datCapture_LL$CaptureSpeed,10000,XRange,YRange)
  stat_CapTurnVsSpeed_DF <- bootStrap_stat(datCapture_DL$Undershoot,datCapture_DL$CaptureSpeed,10000,XRange,YRange)
  
  
  # TURN Vs Capture Speed Mututal INformation 
  bkSeq <- seq(0,4,0.02)
  hist(stat_CapTurnVsSpeed_NF$MI,xlim=c(0,4),ylim=c(0,300),col=colourL[2],breaks = bkSeq ,
       xlab="Mutual information between capture speed and distance to prey", main="Bootstrapped Mutual information")
  hist(stat_CapTurnVsSpeed_LF$MI,xlim=c(0,4),col=colourL[1],add=TRUE ,breaks = bkSeq)
  hist(stat_CapTurnVsSpeed_DF$MI,xlim=c(0,4),col=colourL[3],add=TRUE,breaks = bkSeq )
  
  
  
  # TURN Correlation to Speed / For LF undershooting is combined with faster captures (and more distal) - Not for NF, or DF
  bkSeq <- seq(-0.8,0.8,0.02)
  hist(stat_CapTurnVsSpeed_NF$corr,xlim=c(-0.8,0.8),ylim=c(0,300),col=colourL[2],breaks = bkSeq,xlab="Pearson's correlation turn-ratio vs speed",main="Bootstraped 0.80" )
  hist(stat_CapTurnVsSpeed_LF$corr,xlim=c(-0.8,0.8),col=colourL[1],add=TRUE ,breaks = bkSeq)
  hist(stat_CapTurnVsSpeed_DF$corr,xlim=c(-0.8,0.8),col=colourL[3],add=TRUE,breaks = bkSeq )

  #  PLot Density Turn Vs Speed
  strPlotName = paste(strPlotExportPath,"/stat/fig6_statbootstrap_correlation_TurnVsSpeed.pdf",sep="")
  pdf(strPlotName,width=7,height=7,title="Correlations In hunt variables - turn-ratio vs capture Speed",onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))
  par(mar = c(3.9,4.7,1,1))
  
    pBw <- 0.02
    plot(density(stat_CapTurnVsSpeed_NF$corr,kernel="gaussian",bw=pBw),
         col=colourLegL[1],xlim=c(-0.5,0.5),lwd=3,lty=1,ylim=c(0,10),main=NA, xlab=NA,ylab=NA,cex=cex,cex.axis=cex) #expression(paste("slope ",gamma) ) )
    lines(density(stat_CapTurnVsSpeed_LF$corr,kernel="gaussian",bw=pBw),col=colourLegL[2],lwd=3,lty=2)
    lines(density(stat_CapTurnVsSpeed_DF$corr,kernel="gaussian",bw=pBw),col=colourLegL[3],lwd=3,lty=3)
    mtext(side = 1,cex=cex,cex.main=cex, line = lineXAxis, expression(paste("Correlation of turn-ratio to capture speed  ") ))
    mtext(side = 2,cex=cex,cex.main=cex, line = lineAxis, expression("Density function"))
  
  dev.off()  
  
  
    
  
  
  
#### UNDERSHOOT TO DISTANCE FROM PREY 
  XRange  <- c(0,2) #
  YRange <- c(0,0.8) ##We limit The information Obtained To Reasonable Ranges Of Phi (Vergence Angle)
  
  stat_CapTurnVsDist_NF <- bootStrap_stat(datCapture_NL$Undershoot,datCapture_NL$DistanceToPrey,10000,XRange,YRange)
  stat_CapTurnVsDist_LF <- bootStrap_stat(datCapture_LL$Undershoot,datCapture_LL$DistanceToPrey,10000,XRange,YRange)
  stat_CapTurnVsDist_DF <- bootStrap_stat(datCapture_DL$Undershoot,datCapture_DL$DistanceToPrey,10000,XRange,YRange)
  
    # TURN Vs Capture Speed Mutual INformation  
  bkSeq <- seq(0,4,0.02)
  hist(stat_CapTurnVsDist_NF$MI,xlim=c(0,4),ylim=c(0,300),col=colourL[2],breaks = bkSeq ,
       xlab="Mutual information between turn-ratio and distance to prey", main="Bootstrapped Mutual information")
  hist(stat_CapTurnVsDist_LF$MI,xlim=c(0,4),col=colourL[1],add=TRUE ,breaks = bkSeq)
  hist(stat_CapTurnVsDist_DF$MI,xlim=c(0,4),col=colourL[3],add=TRUE,breaks = bkSeq )
  
  # TURN Correlation to Speed - DF/LF UNdershoot correlates with distance increase -  NF, the opposite correlation arises
  #require( tikzDevice )
  #strPlotName = paste(strPlotExportPath,"/Correlations_HuntTurnVsDist.tex",sep="")
  #tikz( strPlotName )
  
    bkSeq <- seq(-0.8,0.8,0.02)
    hist(stat_CapTurnVsDist_NF$corr,xlim=c(-0.8,0.8),ylim=c(0,300),col=colourL[2],breaks = bkSeq,xlab="Pearson's correlation Turn-ratio vs distance",main="Bootstraped 0.80" )
    hist(stat_CapTurnVsDist_LF$corr,xlim=c(-0.8,0.8),col=colourL[1],add=TRUE ,breaks = bkSeq)
    hist(stat_CapTurnVsDist_DF$corr,xlim=c(-0.8,0.8),col=colourL[3],add=TRUE,breaks = bkSeq )
    
#  dev.off()


  ### DENSITY PLOT  
  strPlotName = paste(strPlotExportPath,"/stat/fig6_statbootstrap_correlation_TurnVsDistance.pdf",sep="")
  pdf(strPlotName,width=7,height=7,title="Correlations In hunt variables - turn-ratio vs capture Speed",onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))
  par(mar = c(3.9,4.7,1,1))
  
  pBw <- 0.02
  plot(density(stat_CapTurnVsDist_NF$corr,kernel="gaussian",bw=pBw),
       col=colourLegL[1],xlim=c(-0.5,0.5),lwd=3,lty=1,ylim=c(0,10),main=NA, xlab=NA,ylab=NA,cex=cex,cex.axis=cex) #expression(paste("slope ",gamma) ) )
  lines(density(stat_CapTurnVsDist_LF$corr,kernel="gaussian",bw=pBw),col=colourLegL[2],lwd=3,lty=2)
  lines(density(stat_CapTurnVsDist_DF$corr,kernel="gaussian",bw=pBw),col=colourLegL[3],lwd=3,lty=3)
  mtext(side = 1,cex=cex,cex.main=cex, line = lineXAxis, expression(paste("Correlation of turn-ratio to distance to prey  ") ))
  mtext(side = 2,cex=cex,cex.main=cex, line = lineAxis, expression("Density function"))
  
  dev.off()  
  
  
  
  
  
  
  ### CORRELOLAGRAM ###

library(corrgram)
layout(matrix(c(1,2,3),1,3, byrow = FALSE))
#Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.7,12,1))

strPlotName = paste(strPlotExportPath,"/Correlations_Huntvariables_NF.pdf",sep="")
pdf(strPlotName,width=8,height=8,title="Correlations In hunt variables",onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))

corrgram(cbind(Speed=datCapture_NL$CaptureSpeed,Dist=datCapture_NL$DistanceToPrey,Turnratio=datCapture_NL$Undershoot,FastCluster=datCapture_NL$Cluster)
               , order=FALSE, lower.panel=panel.pie ,
         upper.panel=NULL, text.panel=panel.txt,
         main="NF Hunt variable correlations")
dev.off()

strPlotName = paste(strPlotExportPath,"/Correlations_Huntvariables_LF.pdf",sep="")
pdf(strPlotName,width=8,height=8,title="Correlations In hunt variables",onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))
  corrgram(cbind(Speed=datCapture_LL$CaptureSpeed,Dist=datCapture_LL$DistanceToPrey,Turnratio=datCapture_LL$Undershoot,FastCluster=datCapture_LL$Cluster)
         , order=FALSE, lower.panel=panel.pie ,
         upper.panel=NULL, text.panel=panel.txt,
         main="LF Hunt variable correlations")
dev.off()

strPlotName = paste(strPlotExportPath,"/Correlations_Huntvariables_DF.pdf",sep="")
pdf(strPlotName,width=8,height=8,title="Correlations In hunt variables",onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))

  corrgram(cbind(Speed=datCapture_DL$CaptureSpeed,Dist=datCapture_DL$DistanceToPrey,Turnratio=datCapture_DL$Undershoot,FastCluster=datCapture_DL$Cluster)
           , order=FALSE, lower.panel=panel.pie ,
           upper.panel=NULL, text.panel=panel.txt,
           main="DF Hunt variable correlations")
dev.off()

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
