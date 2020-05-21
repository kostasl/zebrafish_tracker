### Kostas Lagogiannis 2019-06-24 

## 3D Gaussian Model for each group, to discover covariance structure in Undershoot to Distance/Speed 
## ******** No clustering Between slow and Fast Swims ****** 
## I made this to complement the Clustering Method, so as to characterize the overall covariance structure
## Update 17/10/19 : To Model Group Behaviour we need to model individual larvae behaviour and estimate differences between group behaviour
##                  ie do not use events to infer changes in group behaviour


library(rjags)
library(runjags)

source("config_lib.R")
source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")



## Returns estaimes of each larvae behaviour from the Model
## Added Estimates Of Covariance 
getEstimatesPerLarva <- function(drawG,stail)
{
  ##Find perdio of  pattern of LarvaId / So we can infer hunt events per larva
  maxIdx <- head(which(drawG$Lid == max(drawG$Lid) ),1 )
  tblEventPerLarva <- table(drawG$Lid[1:maxIdx])
  
  ldist <- list()
  lSpeed <- list()
  lTurnRatio <- list()
  lCovar_SpeedDist <- list()
  lCovar_SpeedTurn <- list()
  lCovar_DistTurn <- list()
  nsam <- NROW(drawG$mu[1,1,,])
  ##Iterate Through each Larva and get a mean estimate of behaviour according to model
  for ( i in (1:head(as.numeric(drawG$NLarv),1)) )
  {
    lTurnRatio[[i]]  <- sapply(tail(drawG$mu[i,1,,],stail),mean)
    ldist[[i]]       <- sapply(tail(drawG$mu[i,3,,],stail),mean)
    lSpeed[[i]]      <- sapply(tail(drawG$mu[i,2,,],stail),mean)
    Ci <- 2;Cj <- 3; ##Speed Distance to Prey Covar Per Larva
    lCovar_SpeedDist[[i]] <- sapply(  drawG$cov[i,Ci,Cj,(nsam-stail):nsam,]/(sqrt(drawG$cov[i,Ci,Ci,(nsam-stail):nsam,]*drawG$cov[i,Cj,Cj,(nsam-stail):nsam,]) )   ,mean )
    Ci <- 2;Cj <- 1; ##Speed Turn-Ratio  Covar Per Larva
    lCovar_SpeedTurn[[i]] <- sapply(  drawG$cov[i,Ci,Cj,(nsam-stail):nsam,]/(sqrt(drawG$cov[i,Ci,Ci,(nsam-stail):nsam,]*drawG$cov[i,Cj,Cj,(nsam-stail):nsam,]) )   ,mean )
    Ci <- 3;Cj <- 1; ##Speed Turn-Ratio  Covar Per Larva
    lCovar_DistTurn[[i]] <- sapply(  drawG$cov[i,Ci,Cj,(nsam-stail):nsam,]/(sqrt(drawG$cov[i,Ci,Ci,(nsam-stail):nsam,]*drawG$cov[i,Cj,Cj,(nsam-stail):nsam,]) )   ,mean )
    
    
    #Ci <- 2;Cj <- 3;
    #lCovar_SpeedDist[[i]] <- (drawG$cov[,Ci,Cj,(nsam-stail):nsam,nchains]/( sqrt(drawG$cov[,Ci,Ci,(nsam-ntail):nsam,nchains])*sqrt(drawG$cov[,Cj,Cj,(nsam-ntail):nsam,nchains]) )   )
  }
  ##Overlay The Density From The Estimated Mean Overshoot Of Each Larva
  mEstDistToPrey <- unlist(lapply(ldist,mean) )
  mEstSpeed      <- unlist(lapply(lSpeed,mean) ) 
  mEstTurnRatio  <- unlist(lapply(lTurnRatio,mean) ) 
  Covar_SpeedDist <- unlist(lapply(lCovar_SpeedDist,mean) ) 
  Covar_SpeedTurn <- unlist(lapply(lCovar_SpeedTurn,mean) ) 
  Covar_DistTurn <- unlist(lapply(lCovar_DistTurn,mean) ) 
  return(cbind(HuntEvents=tblEventPerLarva,Undershoot=mEstTurnRatio, DistanceToPrey=mEstDistToPrey,CaptureSpeed=mEstSpeed,
               Covar_SpeedDist,Covar_SpeedTurn,Covar_DistTurn))  
  
}

## Processes The Draw Samples, extracts covariance coeficient, from larva that have more than minEventCount Hunt Events
getCov_Coeff <- function(Ci,Cj,draw_LF,draw_NF,draw_DF,ntail,minEventCount=0)
{
  
  ##Covariance 
  nsam <- NROW(draw_LF$muG[,3,,1])
  ##Select Those with More Than N Hunt Events
  
  tblEventPerLarva <- list()
  maxIdx <- head(which(draw_LF$Lid == max(draw_LF$Lid) ),1 )
  tblEventPerLarva$LF <- table(draw_LF$Lid[1:maxIdx])
  maxIdx <- head(which(draw_NF$Lid == max(draw_NF$Lid) ),1 )
  tblEventPerLarva$NF <- table(draw_NF$Lid[1:maxIdx])
  maxIdx <- head(which(draw_DF$Lid == max(draw_DF$Lid) ),1 )
  tblEventPerLarva$DF <- table(draw_DF$Lid[1:maxIdx])

  ##If Called with Requirement to filter for Larvae that have at least  n Data points   
  if (minEventCount > 0)
  {
    Lid_SubSet <- list(LF=as.numeric(names(tblEventPerLarva$LF[tblEventPerLarva$LF > minEventCount])),
                     NF=as.numeric(names(tblEventPerLarva$NF[tblEventPerLarva$NF > minEventCount])),
                     DF=as.numeric(names(tblEventPerLarva$DF[tblEventPerLarva$DF > minEventCount]))
                    )
  }else ##Alll Estimates so Unscertainty in Estimates can be fully evaluated
  {
    
    Lid_SubSet <- list(LF=1:NROW(draw_LF$mu),
                       NF=1:NROW(draw_NF$mu),
                       DF=1:NROW(draw_DF$mu)
                      )
  }
  
  
  
  nchains <- 1:7
  cov_coeff <- list()
  cov_coeff$LF <- (draw_LF$cov[Lid_SubSet$LF,Ci,Cj,(nsam-ntail):nsam,nchains]/( sqrt(draw_LF$cov[Lid_SubSet$LF,Ci,Ci,(nsam-ntail):nsam,nchains])*sqrt(draw_LF$cov[Lid_SubSet$LF,Cj,Cj,(nsam-ntail):nsam,nchains]) )   )
  cov_coeff$NF <- (draw_NF$cov[Lid_SubSet$NF,Ci,Cj,(nsam-ntail):nsam,nchains]/sqrt(draw_NF$cov[Lid_SubSet$NF,Ci,Ci,(nsam-ntail):nsam,nchains]*draw_NF$cov[Lid_SubSet$NF,Cj,Cj,(nsam-ntail):nsam,nchains]) )
  cov_coeff$DF <- (draw_DF$cov[Lid_SubSet$DF,Ci,Cj,(nsam-ntail):nsam,nchains]/sqrt(draw_DF$cov[Lid_SubSet$DF,Ci,Ci,(nsam-ntail):nsam,nchains]*draw_DF$cov[Lid_SubSet$DF,Cj,Cj,(nsam-ntail):nsam,nchains]) )
  
  ## Get Posterior Sample Differences Between Groups
  cov_coeff_LFVsNF <- tail(cov_coeff$LF,ntail*NROW(nchains)) - tail(cov_coeff$NF,ntail*NROW(nchains))
  cov_coeff_LFVsDF <- tail(cov_coeff$LF,ntail*NROW(nchains)) - tail(cov_coeff$DF,ntail*NROW(nchains))
  cov_coeff_NFVsDF <- tail(cov_coeff$NF,ntail*NROW(nchains)) - tail(cov_coeff$DF,ntail*NROW(nchains))
  #cov_LF$ij <- t(data.frame(draw_LF$cov[,Ci,Cj,,] ))   cov_LF$ii <- t(data.frame(draw_LF$cov[,Ci,Ci,,]))   cov_LF$jj <- t(data.frame(draw_LF$cov[,Cj,Cj,,]))
  ##Calculate Posterior Density Of Covar Coefficient rho
  #cov_LF$coeff <- cov_LF$ij/sqrt(cov_LF$ii*cov_LF$jj)
  
  ##Average Over Columns/Samples - Produce Distribution of Mean Cov Of Group -  Across Samples (Vector With n mean points)
  ##What is the Members Cov On Average?Distibution of E[rho_g] = Sum(rho_i,NLarvae)/NLarvae
  message("E[rho_g] Covar LF:",mean(apply(cov_coeff$LF,1,"mean") ) )
  message("E[rho_g] Covar NF:",mean(apply(cov_coeff$NF,1,"mean") ) )
  message("E[rho_g] Covar DF:",mean(apply(cov_coeff$DF,1,"mean") ) )
  message("Prob that we observe a positive covariance in a group :")
  P_LFGtThanZero <- length(cov_coeff$LF[cov_coeff$NF > 0])/length(cov_coeff$LF)
  P_NFGtThanZero <- length(cov_coeff$DF[cov_coeff$NF > 0])/length(cov_coeff$NF)
  P_DFGtThanZero <- length(cov_coeff$DF[cov_coeff$DF > 0])/length(cov_coeff$DF)
  message("Prob Observe +ve covar in LF > 0 =",P_LFGtThanZero)
  message("Prob Observe +ve covar in NF > 0 =",P_NFGtThanZero)
  message("Prob Observe +ve covar in DF > 0 =",P_DFGtThanZero)
  
  message("Calculate Probabilities that *mean* covariance of group is positive : ")
  ##Using length(cov_coeff$NF[cov_coeff$NF > 0])/length(cov_coeff$NF) can give estimate of samling a covariance
  ##that is above zero in the population'
  
  muCovarPerSample <- list()
  muCovarPerSample$LF <- apply(cov_coeff$LF,2,"mean")
  muCovarPerSample$NF <- apply(cov_coeff$NF,2,"mean")
  muCovarPerSample$DF <- apply(cov_coeff$DF,2,"mean")
  muCovarPerSample$LFvsNF <- muCovarPerSample$LF - muCovarPerSample$NF
  muCovarPerSample$LFvsDF <- muCovarPerSample$LF - muCovarPerSample$DF
  P_LFMeanGreaterThanZero <- length(muCovarPerSample$LF[muCovarPerSample$LF >0] )/length(muCovarPerSample$LF)
  P_NFMeanGreaterThanZero <- length(muCovarPerSample$NF[muCovarPerSample$NF >0] )/length(muCovarPerSample$NF)
  P_DFMeanGreaterThanZero <- length(muCovarPerSample$DF[muCovarPerSample$DF >0] )/length(muCovarPerSample$DF)
  P_LFMeanGreaterThanNF <- length(  muCovarPerSample$LFvsNF[  muCovarPerSample$LFvsNF > 0])/length(  muCovarPerSample$LFvsNF)
  P_LFMeanGreaterThanDF <- length(  muCovarPerSample$LFvsDF[  muCovarPerSample$LFvsDF > 0])/length(  muCovarPerSample$LFvsDF)
  message("Prob Mean group covar LF > 0= ",P_LFMeanGreaterThanZero)
  message("Prob Mean group covar NF > 0= ",P_NFMeanGreaterThanZero)
  message("Prob Mean group covar DF > 0= ",P_DFMeanGreaterThanZero)
  message("Prob Mean group covar LF > NF= ",P_LFMeanGreaterThanNF)
  message("Prob Mean group covar DF > LF= ",1-P_LFMeanGreaterThanDF)
  
  return(cov_coeff)
}

## Plots the Data Density and the 2 Gaussians fititng high and low speed capture swims
plotCaptureSpeedFit <- function(datSpeed,drawMCMC,colourIdx,nchain = 1)
{
  xquant <- seq(0,70,1)
  XLIM <- c(0,60)
  YLIM <- c(0,0.15)
  pdistBW <- 5 ## mm/sec
  strKern <- "gaussian"
  #ntail <- NROW(drawMCMC$mu[1,2,,nchain])*0.10
  ntail <- min(50,NROW(drawMCMC$mu[1,1,,1])*0.10)
  
  plot(density(datSpeed$CaptureSpeed,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM,ylim=YLIM ,main=NA)
  for (i in 1:(ntail-1) )
  {
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[1,2,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[1,2,ntail-i,nchain],1)),type='l',col=colourHLine[colourIdx],lty=1 )
    #lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[2,2,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[2,2,ntail-i,nchain],1)),type='l',col=colourH[colourIdx],lty=2 )
  }
  ##Data
  lines(density(datSpeed$CaptureSpeed,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM )
  legend("topright",title="",
         legend=c( paste("Data Density "), #(Bw:",prettyNum(digits=2, pdistBW ),")" ) ,
                   paste("Model low speed " ),
                   paste("Model high speed " )),
         col=c("black",colourR[4],colourLegL[colourIdx]),lwd=c(3,1,1),lty=c(1,1,2) ) 
  
}


## Plots the Data Density and the 2 Gaussians fititng high and low speed capture swims
plotUndeshootClusterFit <- function(datTurn,drawMCMC,colourIdx,nchain = 1)
{
  xquant <- seq(0,2,0.02)
  XLIM <- c(0,2)
  YLIM <- c(0,3)
  pdistBW <- 0.1 ## mm/sec
  strKern <- "gaussian"
  ntail <- min(50,NROW(drawMCMC$mu[1,1,,1])*0.10)
  
  plot(density(datTurn$Undershoot,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM,ylim=YLIM ,main=NA)
  for (i in 1:(ntail-1) )
  {
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[1,1,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[1,1,ntail-i,nchain],1)),type='l',col=colourHLine[colourIdx],lty=1 )
    #lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[2,1,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[2,1,ntail-i,nchain],1)),type='l',col=colourLegL[colourIdx],lty=2 )
  }
  
  lines(density(datTurn$Undershoot,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM )
  legend("topright",title="",
         legend=c( paste("Data Density "), #(Bw:",prettyNum(digits=2, pdistBW ),")" ) ,
                   paste("Model low speed " ),
                   paste("Model high speed " )),
         col=c("black",colourR[4],colourLegL[colourIdx]),lwd=c(3,1,1),lty=c(1,1,2) ) 
  
}


plotDistanceClustFit <- function(datDist,drawMCMC,colourIdx,nchain = 1)
{
  xquant <- seq(-0.1,0.8,0.05)
  pdistBW <- DIM_MMPERPX ## Manuall annotation  error is at least 1 px error , so smoothing with this bw is relevant
  strKern <- "gaussian"
  ntail <- min(50,NROW(drawMCMC$mu[1,1,,1])*0.10)
  plot(density(datDist$DistanceToPrey,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=c(0,0.8),ylim=c(0,5) ,
       main=NA)
  for (i in 1:(ntail-1) )
  {
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[1,3,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[1,3,ntail-i,nchain],1)),type='l',col=colourHLine[colourIdx],lty=1 )
    #lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[2,3,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[2,3,ntail-i,nchain],1)),type='l',col=colourLegL[colourIdx],lty=2 )
  }
  lines(density(datDist$DistanceToPrey,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=c(0,0.8) )
  legend("topright",title=NA,
         legend=c( paste("Data Density "), #(Bw:",prettyNum(digits=2, pdistBW ),")" ) ,
                   paste("Model low speed " ),
                   paste("Model high speed " )),
         col=c("black",colourR[4],colourLegL[colourIdx]),lwd=c(3,1,1),lty=c(1,1,2) ) 

}


plotCovar_SpeedVsTurnRatio <- function(draw_LF,draw_NF,draw_DF,ntail)
{
  ##Speed TO Distance Covariance Coeff      
  ### Show Speed Fit ###
  outer = FALSE
  line = 1 ## SubFig Label Params
  lineAxis = 4.7
  lineXAxis = 3.5
  lineTitle = 0.5
  
  cex = 1.4
  adj  = 3.5
  padj <- 0.9
  las <- 1
  ####       BOTTOM,LEFT,TOP,RIGHT
  par(mar = c(6.9,5.5,4.5,1))
  layout(matrix(c(1,1,2,3,4),1,5, byrow = TRUE))
  Ci <- 1
  Cj <- 2
  ##Calculate mean cov. coeff based on larvae that have >1 hunt events 
  #(For Larvae with n<2 Cov Is meaningless/ and just adds to 0 group mean)
  cov_coeff <- plotModelCovCoeff(Ci,Cj,draw_LF,draw_NF,draw_DF,ntail,0)
  mtext(side = 2,cex=cex,padj=padj, line = lineAxis, expression("Density function") )
  mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Group covariance of turn-ratio to capture speed" ) )  )
  mtext("A",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=-15,adj=adj,cex.main=cex,cex=cex)
  
  nNF <- dim(cov_coeff$NF)[1]
  nLF <- dim(cov_coeff$LF)[1]
  nDF <- dim(cov_coeff$DF)[1]
  
  legend("topleft",
         legend=c(  expression (),
                    bquote(NF[""] ~'/'~.(nNF)),
                    bquote(LF[""] ~'/'~.(nLF)),
                    bquote(DF[""] ~'/'~.(nDF))
                    #, bquote(All ~ '#' ~ .(ldata_ALL$N)  )
         ),
         lty=c(2,1,3), col=colourLegL,cex=cex,lwd=3)
  #Get Cov coeff for All larvae regardless of min number of events
  cov_coeff_SpeedTurnRatio <- getCov_Coeff(Ci,Cj,draw_LF,draw_NF,draw_DF,ntail,0)
  
  ## plot ecdf / with Confidence Interval About Mean (SEM)
  par(mar = c(6.9,5.0,4.5,1))
  plotECDF_withCI(cov_coeff_SpeedTurnRatio$LF,lModelEst_LF[,"HuntEvents"],colourLegL[2],pchL[4],colourHLine[2],NewPlot=TRUE)
  mtext(side = 2,cex=cex, line = lineAxis,padj=padj, expression("Cumulative function") )
  mtext("B",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=-15,adj=adj,cex.main=cex,cex=cex)
  ####       BOTTOM,LEFT,TOP,RIGHT
  #par(mar = c(6.9,2.5,4.5,1))
  plotECDF_withCI(cov_coeff_SpeedTurnRatio$NF,lModelEst_NF[,"HuntEvents"],colourLegL[1],pchL[6],colourHLine[1],NewPlot=TRUE)
  mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Covariance of turn-ratio to capture speed per larva" ) )  )
  plotECDF_withCI(cov_coeff_SpeedTurnRatio$DF,lModelEst_DF[,"HuntEvents"],colourLegL[3],pchL[5],colourHLine[3],NewPlot=TRUE)
  
  
  
  ##Filtered With Larvae having N>X hunt Events
  #  Covar_SpeedTurnNH <-  HuntModelStat_LF[HuntModelStat_LF[,"HuntEvents"] > 2 ,"Covar_SpeedTurn"]
  #plot(ecdf(Covar_SpeedTurnNH ),col=colourLegL[2],main=NA,pch=pchL[4],xlab=NA,ylab=NA,cex=cex,cex.axis=cex,xlim=c(-1,1))
  
  
  nNF <- dim(cov_coeff_SpeedTurnRatio$NF)[1]
  nLF <- dim(cov_coeff_SpeedTurnRatio$LF)[1]
  nDF <- dim(cov_coeff_SpeedTurnRatio$DF)[1]
  
  legend(x=-1.1,y=1.0,
         legend=c(  expression (),
                    bquote(NF[""] ~'/'~.(nNF) ),
                    bquote(LF[""] ~'/'~.(nLF) ),
                    bquote(DF[""] ~'/'~.(nDF) )
                    #, bquote(All ~ '#' ~ .(ldata_ALL$N)  )
         ),title="Group/#Larvae",
         pch=c(pchL[6],pchL[4],pchL[5]), col=colourLegL,cex=cex)
  
}

plotCovar_SpeedVsDistance <- function(draw_LF,draw_NF,draw_DF,ntail)
{
  
  par(mar = c(6.9,5.7,4.5,1))
  layout(matrix(c(1,1,2,3,4),1,5, byrow = TRUE))
  outer = FALSE
  line = 1 ## SubFig Label Params
  lineAxis = 4.7
  lineXAxis = 3.5
  lineTitle = 0.5
  cex = 1.4
  adj  = 3.5
  las <- 1
  padj <- 0.5
  
  Ci <- 2
  Cj <- 3
  cov_coeff <- plotModelCovCoeff(Ci,Cj,draw_LF,draw_NF,draw_DF,ntail,0,XRange=c(-0.4,0.4))
  mtext(side = 2,cex=cex, line = lineAxis,padj=padj, expression("Density function") )
  mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Group covariance of capture speed to distance" ) )  )
  mtext("A",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=-15,adj=adj,cex.main=cex,cex=cex)
  
  nNF <- dim(cov_coeff$NF)[1]
  nLF <- dim(cov_coeff$LF)[1]
  nDF <- dim(cov_coeff$DF)[1]
  
  legend("topleft",
         legend=c(  expression (),
                    bquote(NF[""] ~'/'~.(nNF)),
                    bquote(LF[""] ~'/'~.(nLF)),
                    bquote(DF[""] ~'/'~.(nDF))
                    #, bquote(All ~ '#' ~ .(ldata_ALL$N)  )
         ),
         lty=c(2,1,3), col=colourLegL,cex=cex,lwd=3)
  
  ## plot ecdf
  #plot(ecdf(apply(cov_coeff$LF,1,"mean")),col=colourLegL[2],main=NA,pch=pchL[4],xlab=NA,ylab=NA,cex=cex,cex.axis=cex)
  #lines(ecdf(apply(cov_coeff$NF,1,"mean")),col=colourLegL[1],pch=pchL[6],cex=cex)
  #lines(ecdf(apply(cov_coeff$DF,1,"mean")),col=colourLegL[3],pch=pchL[5],cex=cex)
  #mtext(side = 2,cex=cex, line = lineAxis,padj=padj, expression("Cumulative function") )
  #mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Est. capture speed-distance covariance per larva" ) )  )
  #mtext("B",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=-15,adj=adj,cex.main=cex,cex=cex)
  
  #Get Cov coeff for All larvae regardless of min number of events
  cov_coeff <- getCov_Coeff(Ci,Cj,draw_LF,draw_NF,draw_DF,ntail,0)
  
  par(mar = c(6.9,5.0,4.5,1))
  ## plot ecdf / with Confidence Interval About Mean (SEM)
  plotECDF_withCI(cov_coeff$LF,lModelEst_LF[,"HuntEvents"],colourLegL[2],pchL[4],colourHLine[2],NewPlot=TRUE)
  mtext(side = 2,cex=cex, line = lineAxis,padj=padj, expression("Cumulative function") )
  mtext("B",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=-15,adj=adj,cex.main=cex,cex=cex)
  
  plotECDF_withCI(cov_coeff$NF,lModelEst_NF[,"HuntEvents"],colourLegL[1],pchL[6],colourHLine[1],NewPlot=TRUE)
  mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Larva covariance of capture speed to distance " ) )  )
  plotECDF_withCI(cov_coeff$DF,lModelEst_DF[,"HuntEvents"],colourLegL[3],pchL[5],colourHLine[3],NewPlot=TRUE)
  
  
  
  nNF <- dim(cov_coeff$NF)[1]
  nLF <- dim(cov_coeff$LF)[1]
  nDF <- dim(cov_coeff$DF)[1]
  
  legend(x=-1.1,y=1.0,
         legend=c(  expression (),
                    bquote(NF[""] ~'/'~.(nNF) ),
                    bquote(LF[""] ~'/'~.(nLF) ),
                    bquote(DF[""] ~'/'~.(nDF) )
                    #, bquote(All ~ '#' ~ .(ldata_ALL$N)  )
         ),title="Group/#Larvae",
         pch=c(pchL[6],pchL[4],pchL[5]), col=colourLegL,cex=cex)
  
}

plotCovar_DistanceVsTurnRatio <- function(draw_LF,draw_NF,draw_DF,ntail)
{
  
  par(mar = c(6.9,5.7,4.5,1))
  layout(matrix(c(1,1,2,3,4),1,5, byrow = TRUE))
  outer = FALSE
  line = 1 ## SubFig Label Params
  lineAxis = 4.7
  lineXAxis = 3.5
  lineTitle = 0.5
  
  cex = 1.4
  adj  = 3.5
  padj <- 0.9
  las <- 1
  
  Ci <- 1
  Cj <- 3
  cov_coeff_TurnRatio <-  plotModelCovCoeff(Ci,Cj,draw_LF,draw_NF,draw_DF,ntail,0,XRange=c(-0.4,0.4))
  mtext(side = 2,cex=cex, line = lineAxis,padj=padj, expression("Density function") )
  mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Group covariance of turn-ratio to capture distance" ) )  )
  mtext("A",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=-15,adj=adj,cex.main=cex,cex=cex)
  
  nNF <- dim(cov_coeff_TurnRatio$NF)[1]
  nLF <- dim(cov_coeff_TurnRatio$LF)[1]
  nDF <- dim(cov_coeff_TurnRatio$DF)[1]
  
  legend("topleft",
         legend=c(  expression (),
                    bquote(NF[""] ~'/'~.(nNF)),
                    bquote(LF[""] ~'/'~.(nLF)),
                    bquote(DF[""] ~'/'~.(nDF))
                    #, bquote(All ~ '#' ~ .(ldata_ALL$N)  )
         ),
         lty=c(2,1,3), col=colourLegL,cex=cex,lwd=3)
  
  cov_coeff_TurnRatio <- getCov_Coeff(Ci,Cj,draw_LF,draw_NF,draw_DF,ntail,0)
  ## plot ecdf
  par(mar = c(6.9,5.0,4.5,1))
  ## plot ecdf / with Confidence Interval About Mean (SEM)
  plotECDF_withCI(cov_coeff_TurnRatio$LF,lModelEst_LF[,"HuntEvents"],colourLegL[2],pchL[4],colourHLine[2],NewPlot=TRUE)
  mtext(side = 2,cex=cex, line = lineAxis,padj=padj, expression("Cumulative function") )
  mtext("B",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=-15,adj=adj,cex.main=cex,cex=cex)
  
  plotECDF_withCI(cov_coeff_TurnRatio$NF,lModelEst_NF[,"HuntEvents"],colourLegL[1],pchL[6],colourHLine[1],NewPlot=TRUE)
  mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Larva covariance of capture distance to turn-ratio " ) )  )
  plotECDF_withCI(cov_coeff_TurnRatio$DF,lModelEst_DF[,"HuntEvents"],colourLegL[3],pchL[5],colourHLine[3],NewPlot=TRUE)
  
  nNF <- dim(cov_coeff_TurnRatio$NF)[1]
  nLF <- dim(cov_coeff_TurnRatio$LF)[1]
  nDF <- dim(cov_coeff_TurnRatio$DF)[1]
  
  legend(x=-1.1,y=1.0,
         legend=c(  expression (),
                    bquote(NF[""] ~'/'~.(nNF) ),
                    bquote(LF[""] ~'/'~.(nLF) ),
                    bquote(DF[""] ~'/'~.(nDF) )
                    #, bquote(All ~ '#' ~ .(ldata_ALL$N)  )
         ),title="Group/#Larvae",
         pch=c(pchL[6],pchL[4],pchL[5]), col=colourLegL,cex=cex)
  
}

initfunct <- function(nchains,N)
{
  initlist <- replicate(nchains,list(#mID=c(rbinom(N,1,0.5)), 
#                                     sigma = matrix(c (  c(runif(1,min=0,max=0.1),runif(1,min=0,max=2)),
#s                                                         c(runif(1,min=0,max=0.1),runif(1,min=0,max=15))  ),nrow=2,byrow=T  ),
#                                     mu  = matrix(c (  c( rnorm(1,mean=1,sd=sqrt(1/10) ), rnorm(1,mean=8,sd=sqrt(1/2) ) ),
#                                                        c( rnorm(1,mean=1, sd=sqrt(1/10) ) , rnorm(1,mean=30, sd=sqrt(1/0.1) )    ) )
#                                                     ,nrow=2,byrow = T  ),
                                     ".RNG.name"="base::Super-Duper",
                                     ".RNG.seed"=round(runif(1,0,60000)) ),
                                     simplify=FALSE)
  return(initlist)
}


plotChk_Undershootfit <- function (draw_F)
{
  ##Undershoot Posterior Group Vs INdividual Density
  with(draw_F,{
    plot(density(tail(muG[,1,,], stail) ),ylim=c(0,16),xlim=c(0,2),lwd=2,col="red")
    for ( i in (1:NLarv[1] ) )
      lines( density( tail( mu[i,1,,],stail)),lty=2)
    ###Show Inferred Distribution
    lines(seq(0,2,0.1),dnorm(seq(0,2,0.1),mean=mean( tail(muG[,1,,],stail)),sd=sqrt(mean( tail( 1/tG[,1,,],stail))) ),col="purple",lwd=4)
    
  })
  x <- seq(0,2,by=0.05)
  points(x, 10*dnorm(x,mean=1,sd=1/sqrt(8)),col="red",cex=2 )
}

##Plot Covariance Density of Mean Group Behaviour, given by averaging the covariance of individual larvae
plotModelCovCoeff <- function(Ci,Cj,draw_LF,draw_NF,draw_DF,ntail,minEventCount=0,XRange=c(-0.6,0.6))
{
  cov_coeff <-getCov_Coeff(Ci,Cj,draw_LF,draw_NF,draw_DF,ntail,minEventCount)
  ###Distribution Of Estimated mean Covariances Per Group
  ## Apply mean per sample / across larvae
  pBW <- 0.005
  
  ylimR <- c(0,15)
  plot( density(apply(cov_coeff$LF,2,"mean"),
                from=-1,to=1,n=300,bw=pBW),xlim=XRange ,col=colourLegL[2],lwd=3,main=NA,xlab=NA,ylab=NA,ylim=ylimR,lty=1,cex=cex,cex.axis=cex)
  lines(density( apply(cov_coeff$NF,2,"mean"),
                 from=-1,to=1,n=300,bw=pBW),xlim=XRange,col=colourLegL[1],lwd=3,lty=2 )
  lines(density( apply(cov_coeff$DF,2,"mean"),
                 from=-1,to=1,n=300,bw=pBW),xlim=XRange,col=colourLegL[3],lwd=3,lty=3 )
  
  
  ###Distribution Of Estimated Covariances In Group
  ## Apply mean per sample / across larvae
  #pBW <- 0.02
  #XRange <- c(-0.5,0.5)
  #plot( density(tail(cov_coeff$LF,ntail*NROW(nchains)),
  #              from=-1,to=1,n=200,bw=pBW),col=colourLegL[2],lwd=3,main=NA,xlab=NA,ylab=NA,ylim=ylimR,lty=1,cex=cex,cex.axis=cex)
  #lines(density(tail(cov_coeff$NF,ntail*NROW(nchains) ) ,
  #               from=-1,to=1,n=200,bw=pBW),xlim=XRange,col=colourLegL[1],lwd=3,lty=2 )
  #lines(density( tail(cov_coeff$DF,ntail*NROW(nchains) ),
   #              from=-1,to=1,n=200,bw=pBW),xlim=XRange,col=colourLegL[3],lwd=3,lty=3 )

  return(cov_coeff)
}

## Plots ECDF Of Estimated Covariance Values Per Larva, Along with a shaded area of Confidence - Using SD/NEvents Per Larva
# The e.c.d.f. (empirical cumulative distribution function) Fn is a step function with jumps i/n at observation 
# values, where i is the number of tied observations at that value.
plotECDF_withCI <- function(cov_coeff,NEventsPerEstimate,colourPoint,pCh=1,colourShade,NewPlot=FALSE,minHEvents=0)
{

  
  ## plot ecdf / with Confidence Interval About Mean (SEM)
  if (NewPlot)
    plot(ecdf(apply(cov_coeff,1,"mean")),col=colourPoint,main=NA,pch=pCh,xlab=NA,ylab=NA,cex=cex,cex.axis=cex,xlim=c(-1,1) )
  edcdf_covar_SemH <- ecdf( apply(cov_coeff,1,"mean") + apply(cov_coeff,1,"sd")/sqrt(NEventsPerEstimate) )
  edcdf_covar_SemL <- ecdf(apply(cov_coeff,1,"mean") - apply(cov_coeff,1,"sd")/sqrt(NEventsPerEstimate) )
  # the unique data values (if there were no ties)
  xxH <- knots(edcdf_covar_SemH);  xxL <- knots(edcdf_covar_SemL) 
  ##Join Shaded Areas to Show Std Error About Mean
  X0 <- c(rev(xxH),xxL)
  Y0 <- c( edcdf_covar_SemH(rev(xxH) ), edcdf_covar_SemL(xxL) )
  polygon(X0,Y0,col=colourShade,border=colourPoint)
  lines(ecdf(apply(cov_coeff,1,"mean")),col=colourPoint,main=NA,pch=pCh,xlab=NA,ylab=NA,cex=cex,cex.axis=cex)
  #lines(edcdf_covar_SemL,pch=pCh,cex=0.6,col=colourPoint)
  #lines(edcdf_covar_SemH,pch=pCh,cex=0.6,col=colourPoint)
  
  
}

##  3D Gaussian Hierarchical  Model of Larvae Hunt Behaviour for (N=60) larva, given k Data points of capture swims
## Some Larvgae had produced no data 
## Estimating Hunt Behaviour per Larvae before inferring mean group behaviour
strmodel3Variables_LarvaHuntBehaviour <- "
var x_rand[NLarv,3];

model {

##Draw  speed,turn-ratio,distance from 3d gaussian
for (i in 1:N)
{
  ##Draw from gaussian model of Each Larva
  c[i,1:3] ~ dmnorm(mu[Lid[i],],prec[Lid[i], , ]) ## data in column 1 and 2
  
}

  

## for each Gaussian in the mixture - Single Gaussian  Here -
for  (l in 1:NLarv)
{
  ##Covariance matrix and its inverse -> the precision matrix
  prec[l,1:3,1:3] ~ dwish(R,3)
  cov[l,1:3,1:3]  <- inverse(prec[l,1:3,1:3])  
  
  ## Larva priors Are linked to the Group's Priors
  mu[l,1] ~ dnorm(muG[1,1], tG[1,1]) ##turn ratio
  mu[l,2] ~ dnorm(muG[1,2], tG[1,2]) ##cap speed
  mu[l,3] ~ dnorm(muG[1,3], tG[1,3]) ##Distance prey
  
  ## Synthesize data from the distribution for This Larva
  x_rand[l,] ~ dmnorm(mu[l,],prec[l,,])
  
}

### Make Group Priors 
for  (g in 1:1)
{
  muG[g,1] ~ dnorm(1, 8)T(0.0,2) ##turn ratio
  muG[g,2] ~ dnorm(25,0.001)T(0,) ##cap speed
  muG[g,3] ~ dnorm(0.1,0.1)T(0,) ##Distance prey
  
  tG[g,1] ~ dgamma(10,3)
  tG[g,2] ~ dgamma(4,25) ##sd < 20
  tG[g,3] ~ dgamma(10,3)
}

for(i in 1:3){
  for(j in 1:3){
      R[i,j] <- equals(i,j)*1e-4
    }
}
  ## Possible to establish Wishart Prior? 
  #CG[1:3,1:3] ~ dwish(R,3)
  #covG[1:3,1:3] <- inverse(CG)
  
} "



strModelPDFFileName <- "/stat/fig7S1-stat_3Dmodel_TurnRatioVsSpeedAndDistance.pdf"
strDataPDFFileName <- "/stat/fig7-UndershootCaptureSpeedCV_scatter_Valid.pdf"
strCaptSpeedDensityPDFFileName <- "/stat/fig7-stat_modelCaptureSpeed_Valid.pdf"
strUndershootDensityPDFFileName <- "/stat/fig7-stat_modelUndershoot_Valid.pdf"
strDistanceDensityPDFFileName <- "/stat/stat_modelDistance_Valid.pdf"
strModelCovarPDFFileName <- "/stat/fig7-stat_modelCaptureSpeedVsUndershootAndDistance_COVar.pdf"



### Load PreCalculated Model Results ###
load(paste0(strDataExportDir,"stat_Larval3DGaussianBehaviouModel_RJags.RData"))

datTrackedEventsRegister <- readRDS( paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate.rds","",sep="") ) ## THis is the Processed Register File On 
#lMotionBoutDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_SetC.rds",sep="") ) #Processed Registry on which we add )
#lEyeMotionDat <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData_SetC",".rds",sep="")) #
lFirstBoutPoints <-readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_wCapFrame_Validated.rds",sep="")) 

#### LOAD Capture First-Last Bout hunting that include the cluster classification - (made in stat_CaptureSpeedVsDistanceToPrey)
##22/10/19- Updated with Time To get To prey INfo 
datCapture_NL <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_wCapFrame_NL_clustered.rds",sep="")) 
datCapture_LL <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_wCapFrame_LL_clustered.rds",sep="")) 
datCapture_DL <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_wCapFrame_DL_clustered.rds",sep="")) 

#datCapturePointsClustered <- list()
#datCapturePointsClustered$NF <- merge(datCapture_NL, lFirstBoutPoints$NL,by="RegistarIdx",no.dups=FALSE,suffixes=c(NA,".bt"))
#datCapturePointsClustered$LF <- merge(datCapture_LL, lFirstBoutPoints$LL,by="RegistarIdx",no.dups=FALSE,suffixes=c(NA,".bt"))
#datCapturePointsClustered$DF <- merge(datCapture_DL, lFirstBoutPoints$DL,by="RegistarIdx",no.dups=FALSE,suffixes=c(NA,".bt"))

datTurnVsStrikeSpeed_NL <- data.frame( cbind(TurnRatio=lFirstBoutPoints$NL[,"Turn"]/lFirstBoutPoints$NL[,"OnSetAngleToPrey"],CaptureSpeed=lFirstBoutPoints$NL[,"CaptureSpeed"]),DistanceToPrey=lFirstBoutPoints$NL[,"DistanceToPrey"],Validated= lFirstBoutPoints$NL[,"Validated"] )
datTurnVsStrikeSpeed_LL <- data.frame( cbind(TurnRatio=lFirstBoutPoints$LL[,"Turn"]/lFirstBoutPoints$LL[,"OnSetAngleToPrey"],CaptureSpeed=lFirstBoutPoints$LL[,"CaptureSpeed"]),DistanceToPrey=lFirstBoutPoints$LL[,"DistanceToPrey"],Validated= lFirstBoutPoints$LL[,"Validated"] )
datTurnVsStrikeSpeed_DL <- data.frame( cbind(TurnRatio=lFirstBoutPoints$DL[,"Turn"]/lFirstBoutPoints$DL[,"OnSetAngleToPrey"],CaptureSpeed=lFirstBoutPoints$DL[,"CaptureSpeed"]),DistanceToPrey=lFirstBoutPoints$DL[,"DistanceToPrey"],Validated= lFirstBoutPoints$DL[,"Validated"] )

###Filter For Hunting Where Prey Is approached into strike distance, rather than Initial Prey Distance being within strike Distance
#datTurnVsStrikeSpeed_NL <- datTurnVsStrikeSpeed_NL[datTurnVsStrikeSpeed_NL$OnSetDistance > 0.6,]
#datTurnVsStrikeSpeed_LL <- datTurnVsStrikeSpeed_LL[datTurnVsStrikeSpeed_LL$OnSetDistance > 0.6,]
#datTurnVsStrikeSpeed_DL <- datTurnVsStrikeSpeed_DL[datTurnVsStrikeSpeed_DL$OnSetDistance > 0.6,]
###Validated Only
replace(datTurnVsStrikeSpeed_NL$Validated, is.na(datTurnVsStrikeSpeed_NL$Validated), 0)
replace(datTurnVsStrikeSpeed_LL$Validated, is.na(datTurnVsStrikeSpeed_LL$Validated), 0)
replace(datTurnVsStrikeSpeed_DL$Validated, is.na(datTurnVsStrikeSpeed_DL$Validated), 0) 

###Validated Only
datTurnVsStrikeSpeed_NL <- datTurnVsStrikeSpeed_NL[datTurnVsStrikeSpeed_NL$Validated == 1, ]
datTurnVsStrikeSpeed_LL <- datTurnVsStrikeSpeed_LL[datTurnVsStrikeSpeed_LL$Validated == 1, ]
datTurnVsStrikeSpeed_DL <- datTurnVsStrikeSpeed_DL[datTurnVsStrikeSpeed_DL$Validated == 1, ]

datTurnVsStrikeSpeed_ALL <- rbind(datTurnVsStrikeSpeed_NL,datTurnVsStrikeSpeed_LL,datTurnVsStrikeSpeed_DL)



##Get Hunt Success
datHuntLabelledEventsSB <- getLabelledHuntEventsSet()
datFishSuccessRate <- getHuntSuccessPerFish(datHuntLabelledEventsSB)


##Merge Exp IDs - to identify events of individuals
##Take all expID from the successful hunt Events we have extracted hunt variables from 
vexpID <- list(LF = datTrackedEventsRegister[datCapture_LL$RegistarIdx,]$expID,
               NF=datTrackedEventsRegister[datCapture_NL$RegistarIdx,]$expID,
               DF=datTrackedEventsRegister[datCapture_DL$RegistarIdx,]$expID)


## Merge EXP ID
## Add Exp ID Column - Signifying Which Larvae Executed the Capture Success Hunt- 
datCapture_LF_wExpID <- cbind(datCapture_LL,expID=vexpID$LF,groupID=2)
datCapture_NF_wExpID <- cbind(datCapture_NL,expID=vexpID$NF,groupID=3)
datCapture_DF_wExpID <- cbind(datCapture_DL,expID=vexpID$DF,groupID=1)
datCapture_ALL_wExpID <- rbind(datCapture_LF_wExpID,datCapture_NF_wExpID,datCapture_DF_wExpID)


##Merge Hunt Power To Hunt-Capture Variables 
datMergedCapAndSuccess_LF <- merge(x=datCapture_LF_wExpID,y=datFishSuccessRate,by="expID",all.x=TRUE)
datMergedCapAndSuccess_NF <- merge(x=datCapture_NF_wExpID,y=datFishSuccessRate,by="expID",all.x=TRUE)
datMergedCapAndSuccess_DF <- merge(x=datCapture_DF_wExpID,y=datFishSuccessRate,by="expID",all.x=TRUE)

###Empirical Distribution
datHuntLarvaStat <- aggregate(datCapture_ALL_wExpID,by=list(datCapture_ALL_wExpID$expID),mean)

table(datCapture_ALL_wExpID$expID)

##
##
steps <-50000
nchains <- 7
nthin <- 10
#str_vars <- c("mu","rho","sigma","x_rand") #Basic model 
str_vars <- c("mu","cov","x_rand","muG","tG","NLarv","Lid","R") #Mixture Model
##Make Serial Larvae ID, that links each hunt event to an individual larva 
## Maintain RegIDx so we trace Back
groupSizeLF <-groupSizeNF <- groupSizeDF  <- 60#NROW(unique(expID)
ldata_LF <- with(datMergedCapAndSuccess_LF, {list(c=cbind(Undershoot,CaptureSpeed,DistanceToPrey,Cluster),Efficiency=Efficiency,RegIdx=RegistarIdx,Lid=as.numeric(as.factor(as.numeric(expID)) ) ,N=NROW(expID), NLarv=groupSizeLF)  }) ##Live fed
ldata_NF <- with(datMergedCapAndSuccess_NF, {list(c=cbind(Undershoot,CaptureSpeed,DistanceToPrey,Cluster),Efficiency=Efficiency,RegIdx=RegistarIdx,Lid=as.numeric(as.factor(as.numeric(expID)) ) ,N=NROW(expID), NLarv=groupSizeNF)   }) ##Live fed
ldata_DF <- with(datMergedCapAndSuccess_DF, {list(c=cbind(Undershoot,CaptureSpeed,DistanceToPrey,Cluster),Efficiency=Efficiency,RegIdx=RegistarIdx,Lid=as.numeric(as.factor(as.numeric(expID)) ) ,N=NROW(expID), NLarv=groupSizeDF)   }) ##Live fed

#ldata_ALL <-with(datCapture_ALL_wExpID, {list(c=cbind(Undershoot,CaptureSpeed,DistanceToPrey),Efficiency=Efficiency,RegIdx=RegistarIdx,Lid=as.numeric(as.factor(as.numeric(expID)) ),Gid=groupID ,N=NROW(expID),NLarv=NROW(unique(expID))  ) }) ##Live fed list(c=datTurnVsStrikeSpeed_ALL,N=NROW(datTurnVsStrikeSpeed_ALL)) ##Dry fed
ldata_LF_fast <- with(datMergedCapAndSuccess_LF[datMergedCapAndSuccess_LF$Cluster == "fast",], {list(c=cbind(Undershoot,CaptureSpeed,DistanceToPrey,Cluster),Efficiency=Efficiency,RegIdx=RegistarIdx,Lid=as.numeric(as.factor(as.numeric(expID)) ) ,N=NROW(expID), NLarv=groupSizeLF ) }) ##Live fed
ldata_NF_fast <- with(datMergedCapAndSuccess_NF[datMergedCapAndSuccess_NF$Cluster == "fast",], {list(c=cbind(Undershoot,CaptureSpeed,DistanceToPrey,Cluster),Efficiency=Efficiency,RegIdx=RegistarIdx,Lid=as.numeric(as.factor(as.numeric(expID)) ) ,N=NROW(expID), NLarv=groupSizeNF ) }) ##Live fed
ldata_DF_fast <- with(datMergedCapAndSuccess_DF[datMergedCapAndSuccess_DF$Cluster == "fast",], {list(c=cbind(Undershoot,CaptureSpeed,DistanceToPrey,Cluster),Efficiency=Efficiency,RegIdx=RegistarIdx,Lid=as.numeric(as.factor(as.numeric(expID)) ) ,N=NROW(expID), NLarv=groupSizeDF  ) }) ##Live fed

#ldata_ALL <-with(datCapture_ALL_wExpID, {list(c=cbind(Undershoot,CaptureSpeed,DistanceToPrey),Efficiency=Efficiency,RegIdx=RegistarIdx,Lid=as.numeric(as.factor(as.numeric(expID)) ),Gid=groupID ,N=NROW(expID),NLarv=NROW(unique(expID))  ) }) ##Live fed list(c=datTurnVsStrikeSpeed_ALL,N=NROW(datTurnVsStrikeSpeed_ALL)) ##Dry fed
ldata_LF_slow <- with(datMergedCapAndSuccess_LF[datMergedCapAndSuccess_LF$Cluster == "slow",], {list(c=cbind(Undershoot,CaptureSpeed,DistanceToPrey,Cluster),Efficiency=Efficiency,RegIdx=RegistarIdx,Lid=as.numeric(as.factor(as.numeric(expID)) ) ,N=NROW(expID), NLarv=groupSizeLF ) }) ##Live fed
ldata_NF_slow <- with(datMergedCapAndSuccess_NF[datMergedCapAndSuccess_NF$Cluster == "slow",], {list(c=cbind(Undershoot,CaptureSpeed,DistanceToPrey,Cluster),Efficiency=Efficiency,RegIdx=RegistarIdx,Lid=as.numeric(as.factor(as.numeric(expID)) ) ,N=NROW(expID), NLarv=groupSizeNF) }) ##Live fed
ldata_DF_slow <- with(datMergedCapAndSuccess_DF[datMergedCapAndSuccess_DF$Cluster == "slow",], {list(c=cbind(Undershoot,CaptureSpeed,DistanceToPrey,Cluster),Efficiency=Efficiency,RegIdx=RegistarIdx,Lid=as.numeric(as.factor(as.numeric(expID)) ) ,N=NROW(expID), NLarv=groupSizeDF ) }) ##Live fed


##Save for Ready to Eat Public Comsumption
saveRDS(ldata_LF,file=paste0(strDataExportDir,"pubDat/huntEpisodeDataMergedWithLarvalSuccess_LF.rds") )
saveRDS(ldata_NF,file=paste0(strDataExportDir,"pubDat/huntEpisodeDataMergedWithLarvalSuccess_NF.rds") )
saveRDS(ldata_DF,file=paste0(strDataExportDir,"pubDat/huntEpisodeDataMergedWithLarvalSuccess_DF.rds") )
saveRDS(datHuntLarvaStat,file=paste0(strDataExportDir,"pubDat/LarvaEmpiricalMeanHuntBehaviour.rds"))

### RUN JAGS MODEL ###
    jags_model_LF <- jags.model(textConnection(strmodel3Variables_LarvaHuntBehaviour), data = ldata_LF, 
                             n.adapt = 100, n.chains = nchains, quiet = F,inits=initfunct(nchains,ldata_LF$N))
    update(jags_model_LF, 300)
    draw_LF=jags.samples(jags_model_LF,steps,thin=nthin,variable.names=str_vars)
    
    ## Not Fed
    jags_model_NF <- jags.model(textConnection(strmodel3Variables_LarvaHuntBehaviour), data = ldata_NF, 
                             n.adapt = 500, n.chains = nchains, quiet = F,inits=initfunct(nchains,ldata_NF$N)) 
    update(jags_model_NF,300)
    draw_NF=jags.samples(jags_model_NF,steps,thin=nthin,variable.names=str_vars)
    
    ## Dry  Fed
    jags_model_DF <- jags.model(textConnection(strmodel3Variables_LarvaHuntBehaviour), data = ldata_DF, 
                             n.adapt = 500, n.chains = nchains, quiet = F,inits=initfunct(nchains,ldata_DF$N))
    update(jags_model_DF, 300)
    draw_DF=jags.samples(jags_model_DF,steps,thin=nthin,variable.names=str_vars)
####### END OF RUN MODELS ##

message("Mean LF Und:", prettyNum( mean(draw_LF$muG[,1,,]) , digits=3),
        " Speed : ",prettyNum( mean(draw_LF$muG[,2,,1]), digits=3),
        " Distance : ",prettyNum(mean(draw_LF$muG[,3,,1]), digits=3)
)


message("Mean NF Und:", prettyNum( mean(draw_NF$muG[,1,,]) , digits=3),
        " Speed : ",prettyNum( mean(draw_NF$muG[,2,,1]), digits=3),
        " Distance : ",prettyNum(mean(draw_NF$muG[,3,,1]), digits=3)
)

message("Mean DF Und:", prettyNum( mean(draw_DF$muG[,1,,]) , digits=3),
        " Speed : ",prettyNum( mean(draw_DF$muG[,2,,1]), digits=3),
        " Distance : ",prettyNum(mean(draw_DF$muG[,3,,1]), digits=3)
)

save(draw_NF,draw_LF,draw_DF,file = paste0(strDataExportDir,"stat_Larval3DGaussianBehaviouMode_All60_CapturesOnly_RJags.RData"))
## ALL  groups
#jags_model_ALL <- jags.model(textConnection(strmodel_capspeedVsUndershoot_Mixture), data = ldata_ALL, 
                            #n.adapt = 500, n.chains = 3, quiet = F)
#update(jags_model_ALL, 300)
#draw_ALL=jags.samples(jags_model_ALL,steps,thin=2,variable.names=str_vars)
schain <- 1:7
stail <- 300


##Capt Speed
plot(density(tail(draw_LF$muG[,2,,schain], 150) ),ylim=c(0,0.5),xlim=c(0,50),main="Capt. Speed",col=colourLegL[2],lty=1)
lines(density(tail(draw_NF$muG[,2,,schain], 150)),col=colourLegL[1],lty=2)
lines(density(tail(draw_DF$muG[,2,,schain], 150)),col=colourLegL[3],lty=3)

plot(density(tail(draw_LF$muG[,1,,schain], 150) ),ylim=c(0,5),xlim=c(0,2),main="Turn-Ratio",col=colourLegL[2])
lines(density(tail(draw_NF$muG[,1,,schain], 150) ),ylim=c(0,5),xlim=c(0,2),col=colourLegL[1])
lines(density(tail(draw_DF$muG[,1,,schain], 150) ),ylim=c(0,5),xlim=c(0,2),col=colourLegL[3])


## Speed Posterior Group Vs INdividual Density
with(draw_LF,{
  plot(density(tail(muG[,2,,], 150) ),ylim=c(0,1),xlim=c(0,60),lwd=2,col="red",main="Speed LF")
  for ( i in (1:NLarv[1] ) )
    lines( density( tail( mu[i,2,,],stail)),lty=2)
  ###Show Inferred Distribution
  lines(1:100,dnorm(1:100,mean=mean( tail(muG[,2,,3],stail)),sd=sqrt(mean( tail( 1/tG[,2,,1],stail))) ),col="purple",lwd=4)
})
##Compare To Empirical - Change group DF,LF,NF-- V Good Match!
lines(density(datHuntLarvaStat[datHuntLarvaStat$groupID==2,]$CaptureSpeed,bw=2),col="blue",lwd=2)

##Speed Posterior Group Vs INdividual Density
with(draw_DF,{
  plot(density(tail(muG[,2,,], 500) ),ylim=c(0,1),xlim=c(0,60),lwd=2,col="red",main="Speed DF")
  for ( i in (1:NLarv[1] ) )
    lines( density( tail( mu[i,2,,],stail)),lty=2)
  ###Show Inferred Distribution
   lines(1:100,dnorm(1:100,mean=mean( tail(muG[,2,,3],stail)),sd=sqrt(mean( tail( 1/tG[,2,,1],stail)) )) ,col="purple",lwd=4)
})
##Compare To Empirical - Change group DF,LF,NF-- V Good Match!
lines(density(datHuntLarvaStat[datHuntLarvaStat$groupID==1,]$CaptureSpeed,bw=2),col="blue",lwd=2)

##Speed Posterior Group Vs INdividual Density
with(draw_NF,{
  plot(density(tail(muG[,2,,], 150) ),ylim=c(0,1),xlim=c(0,60),lwd=2,col="red",main="Speed NF")
  for ( i in (1:NLarv[1] ) )
    lines( density( tail( mu[i,2,,],stail)),lty=2)
  ###Show Inferred Distribution
  lines(1:100,dnorm(1:100,mean=mean( tail(muG[,2,,3],stail)),sd=sqrt(mean( tail( 1/tG[,2,,1],stail))) ),col="purple",lwd=4)
})
##Compare To Empirical - Change group DF,LF,NF-- V Good Match!
lines(density(datHuntLarvaStat[datHuntLarvaStat$groupID==3,]$CaptureSpeed,bw=2),col="blue",lwd=2)



##Distance Posterior Group Vs INdividual Density
with(draw_NF,{
  plot(density(tail(muG[,3,,1], 100) ),ylim=c(0,16),xlim=c(0,1),lwd=2,col="red")
  for ( i in (1:NLarv[1] ) )
    lines( density( tail( mu[i,3,,1],stail)),lty=2)
})

##CONVERGENCE - CHECK
#Idx: Lid,Variable (1Under),Sample Row,Chain - 
plot(density(tail(draw_LF$muG[,2,,1],1000) ),type='l')
for (c in schain)
  lines(density(tail(draw_LF$muG[,2,,c],1000) ),col="red")


##Make Model Estimate Per Larvae Including RegIDx For Back Reference 
expID_LF <- unique(datTrackedEventsRegister[ldata_LF$RegIdx[which(ldata_LF$Lid %in% 1:NROW(lModelEst_LF) )],]$expID )
expID_NF <- unique(datTrackedEventsRegister[ldata_NF$RegIdx[which(ldata_NF$Lid %in% 1:NROW(lModelEst_NF) )],]$expID )
expID_DF <- unique(datTrackedEventsRegister[ldata_DF$RegIdx[which(ldata_DF$Lid %in% 1:NROW(lModelEst_DF) )],]$expID )
datHunterStat_Model <- rbind(data.frame(lModelEst_LF,groupIDF="LL",expID=expID_LF),
                             data.frame(lModelEst_NF,groupIDF="NL",expID=expID_NF),
                             data.frame(lModelEst_DF,groupIDF="DL",expID=expID_DF))

save(datHunterStat_Model,file = paste0(strDataExportDir,"stat_Larval3DGaussianBehaviourModelPerLarva.RData"))


#### Compare Group Model To Density Obtain through Mean Estimated Behaviour For Each Larva (Dashed) ####
### Obtain Estimated Mean Values For Each Larva & Plot Group Population
## Plot Distance Density
plot(density(sapply(tail(draw_LF$mu[,3,,],stail),mean)),col=colourLegL[2],lty=1 ,lwd=3,main="Distance to Prey Group Post. Vs Mean Larv. post.",ylim=c(0,6)) ##Mean Group Undershoot From Mean Of Each Larva
lModelEst_LF <- getEstimatesPerLarva(draw_LF,stail)
lines(density( unlist(lapply(lModelEst_LF[,"DistanceToPrey"],mean) ) ),lty=2,col=colourLegL[2],lwd=2 )

lines(density(sapply(tail(draw_NF$mu[,3,,],stail),mean)),col=colourLegL[1] ,lty=1 ,lwd=3) ##Mean Group Undershoot From Mean Of Each Larva
lModelEst_NF <- getEstimatesPerLarva(draw_NF,stail)
lines(density( unlist(lapply(lModelEst_NF[,"DistanceToPrey"],mean) ) ),col=colourLegL[1],lty=2,lwd=2 )

lines(density(sapply(tail(draw_DF$mu[,3,,],stail) ,mean)),col=colourLegL[3] ,lty=1,lwd=3) ##Mean Group Undershoot From Mean Of Each Larva
lModelEst_DF <- getEstimatesPerLarva(draw_DF,stail)
lines(density( unlist(lapply(lModelEst_DF[,"DistanceToPrey"],mean) ) ),col=colourLegL[3],lty=2,lwd=2 )

## Plot Model Speed Density
plot(density(sapply(tail(draw_LF$mu[,2,,],stail),mean)),col=colourLegL[2] ,lwd=2,main="Capture Speed",ylim=c(0,0.1)) ##Mean Group Undershoot From Mean Of Each Larva
lines(density(sapply(tail(draw_NF$mu[,2,,],stail),mean)),col=colourLegL[3] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva
lines(density(sapply(tail(draw_DF$mu[,2,,],stail) ,mean)),col=colourLegL[1] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva

## Plot Model Undershoot Density / Mean Sample point Across larva 
plot(density(sapply(tail(draw_LF$mu[,1,,],stail),mean)),col=colourLegL[2] ,lwd=2,main="Turn ratio") ##Mean Group Undershoot From Mean Of Each Larva
lines(density(sapply(tail(draw_NF$mu[,1,,],stail),mean)),col=colourLegL[3] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva
lines(density(sapply(tail(draw_DF$mu[,1,,],stail) ,mean)),col=colourLegL[1] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva

#### Compare Models - Group Behaviour In Probs ####
message("Compare Group Behaviours - Model Based")#draw_DF$muG[,3,,1]

drawCapSpeedLFVsNF <- tail(draw_LF$muG[,2,,],stail)-tail(draw_NF$muG[,2,,],stail)
drawCapSpeedLFVsDF <- tail(draw_LF$muG[,2,,],stail)-tail(draw_DF$muG[,2,,],stail)
drawCapSpeedNFVsDF <- tail(draw_NF$muG[,2,,],stail)-tail(draw_DF$muG[,2,,],stail)
message("Mean Cap Speed LF:",prettyNum( mean(tail(draw_LF$muG[,2,,],stail) ), digits=3) )
message("Mean Cap Speed DF:",prettyNum( mean(tail(draw_DF$muG[,2,,],stail) ), digits=3) )
message("Mean Cap Speed NF:",prettyNum( mean(tail(draw_NF$muG[,2,,],stail) ), digits=3) )

message("Mean diff Cap Speed LF-NF:",prettyNum( mean(drawCapSpeedLFVsNF),digits=3) )
message("Mean diff Cap Speed LF-DF:",prettyNum( mean(drawCapSpeedLFVsDF),digits=3) )
message("Mean diff Cap Speed NF-DF:",prettyNum( mean(drawCapSpeedNFVsDF),digits=3) )

## Probability Of Differences ##
PCapSpeed_LFgtNF <- length(drawCapSpeedLFVsNF[drawCapSpeedLFVsNF > 0])/length(drawCapSpeedLFVsNF) ## ProbValLessThan(dLLbVsNF,0)
PCapSpeed_LFgtDF <- length(drawCapSpeedLFVsDF[drawCapSpeedLFVsDF > 0])/length(drawCapSpeedLFVsDF) ## ProbValLessThan(dLLbVsNF,0)
PCapSpeed_DFgtNF <- 1-length(drawCapSpeedNFVsDF[drawCapSpeedNFVsDF > 0])/length(drawCapSpeedNFVsDF) ## ProbValLessThan(dLLbVsNF,0)
message("Prob that LF Capt Speed > NF: ",prettyNum(PCapSpeed_LFgtNF,digits=3))
message("Prob that LF Capt Speed > DF: ",prettyNum(PCapSpeed_LFgtDF,digits=3))
message("Prob that DF Capt Speed > NF: ",prettyNum(PCapSpeed_DFgtNF,digits=3))

## Turn Ratio ##
drawTurnRatioLFVsNF <- tail(draw_LF$muG[,1,,],stail)-tail(draw_NF$muG[,1,,],stail)
drawTurnRatioLFVsDF <- tail(draw_LF$muG[,1,,],stail)-tail(draw_DF$muG[,1,,],stail)
drawTurnRatioNFVsDF <- tail(draw_NF$muG[,1,,],stail)-tail(draw_DF$muG[,1,,],stail)
message("Mean TurRatio LF:",prettyNum( mean(tail(draw_LF$muG[,1,,],stail) ), digits=3) )
message("Mean TurRatio DF:",prettyNum( mean(tail(draw_DF$muG[,1,,],stail) ), digits=3) )
message("Mean TurRatio NF:",prettyNum( mean(tail(draw_NF$muG[,1,,],stail) ), digits=3) )

message("Mean diff TurRatio LF-NF:",prettyNum( mean(drawTurnRatioLFVsNF),digits=3) )
message("Mean diff TurRatio LF-DF:",prettyNum( mean(drawTurnRatioLFVsDF),digits=3) )
message("Mean diff TurRatio NF-DF:",prettyNum( mean(drawTurnRatioNFVsDF),digits=3) )

## Prob - Diff in Turnration -ve means Higher Undershoot For LF
PUndershoot_LFgtNF <- length(drawTurnRatioLFVsNF[drawTurnRatioLFVsNF < 0])/length(drawTurnRatioLFVsNF) ## ProbValLessThan(dLLbVsNF,0)
PUndershoot_LFgtDF <- length(drawTurnRatioLFVsDF[drawTurnRatioLFVsDF < 0])/length(drawTurnRatioLFVsDF) ## ProbValLessThan(dLLbVsNF,0)
PUndershoot_DFgtNF <- 1-length(drawTurnRatioNFVsDF[drawTurnRatioNFVsDF < 0])/length(drawTurnRatioNFVsDF) ## ProbValLessThan(dLLbVsNF,0)

message("Prob that LF Undershoot > NF: ",prettyNum(PUndershoot_LFgtNF,digits=3))
message("Prob that LF Undershoot > DF: ",prettyNum(PUndershoot_LFgtDF,digits=3))
message("Prob that DF Undershoot > NF: ",prettyNum(PUndershoot_DFgtNF,digits=3))

## Compare Distance 
drawCapDistLFVsNF <- tail(draw_LF$muG[,3,,],stail)-tail(draw_NF$muG[,3,,],stail)
drawCapDistLFVsDF <- tail(draw_LF$muG[,3,,],stail)-tail(draw_DF$muG[,3,,],stail)
drawCapDistNFVsDF <- tail(draw_NF$muG[,3,,],stail)-tail(draw_DF$muG[,3,,],stail)
message("Mean Cap Dist LF:",prettyNum( mean(tail(draw_LF$muG[,3,,],stail) ), digits=3) )
message("Mean Cap Dist DF:",prettyNum( mean(tail(draw_DF$muG[,3,,],stail) ), digits=3) )
message("Mean Cap Dist NF:",prettyNum( mean(tail(draw_NF$muG[,3,,],stail) ), digits=3) )
message("Mean diff Cap Dist LF-NF:",prettyNum( mean(drawCapDistLFVsNF),digits=3) )
message("Mean diff Cap Dist LF-DF:",prettyNum( mean(drawCapDistLFVsDF),digits=3) )
message("Mean diff Cap Dist NF-DF:",prettyNum( mean(drawCapDistNFVsDF),digits=3) )

## Prob - Diff in Turnration -ve means Higher Undershoot For LF
PCapDist_LFgtNF <- length(drawCapDistLFVsNF[drawCapDistLFVsNF > 0])/length(drawCapDistLFVsNF) ## ProbValLessThan(dLLbVsNF,0)
PCapDist_LFgtDF <- length(drawCapDistLFVsDF[drawCapDistLFVsDF > 0])/length(drawCapDistLFVsDF) ## ProbValLessThan(dLLbVsNF,0)
PCapDist_NFgtDF <- length(drawCapDistNFVsDF[drawCapDistNFVsDF > 0])/length(drawCapDistNFVsDF) ## ProbValLessThan(dLLbVsNF,0)
message("Prob that LF Cap Dist > NF: ",prettyNum(PCapDist_LFgtNF,digits=3))
message("Prob that LF Cap Dist > DF: ",prettyNum(PCapDist_LFgtDF,digits=3))
message("Prob that NF Cap Dist > DF: ",prettyNum(PCapDist_NFgtDF,digits=3))


#Idx: Lid,Variable (2=SpeedDist Covar),Sample Row,Chain
## 


### Estimate  densities  ###
nContours <- 5
ntail <- 1200 #NROW(draw_NF$mu[1,1,,1])*0.20
#load(paste0(strDataExportDir,"stat_CaptSpeedVsUndershootAndDistance_RJags.RData"))
lModelEst_LF <- getEstimatesPerLarva(draw_LF,ntail)
lModelEst_NF <- getEstimatesPerLarva(draw_NF,ntail)
lModelEst_DF <- getEstimatesPerLarva(draw_DF,ntail)


## Check out the covar coeffient , compare estimated densities
pBw   <- 0.02
#dALLb_rho<-density(tail(draw_ALL$rho[,,1],ntail),kernel="gaussian",bw=pBw)


#load(paste0(strDataExportDir,"stat_CaptSpeedVsDistance_Covariance_RJags.RData"))
###Check COnv
draw <- draw_DF
varIdx <- 1
plot(draw$muG[,varIdx,,1],type='l',ylim=c(0,2),col=rfc(nchains)[1] )
for (i in 2:nchains)
  lines(draw$muG[,varIdx,,i],type='l',ylim=c(0,2),col=rfc(nchains)[i] )

##Get the synthesized data:
#plot(tail(draw_NF$x_rand[1,,1],ntail ),tail(draw_NF$x_rand[2,,1],ntail ),col=colourH[1])
#points(tail(draw_LF$x_rand[1,,1],ntail ),tail(draw_LF$x_rand[2,,1],ntail ),col=colourH[2])
#points(tail(draw_DF$x_rand[1,,1],ntail ),tail(draw_DF$x_rand[2,,1],ntail ),col=colourH[3])


#### MAIN COVARIANCE PLOT  ##
#### COVARIANCE PLOTS OVER ALL Capture Bouts Fast/Slow ####
load(file = paste0(strDataExportDir,"stat_Larval3DGaussianBehaviouMode_All60_CapturesOnly_RJags.RData"))

## Turn-Ratio(1)xSpeed(2) Covariance Coeff: Calc as rho=Cij/(sigmai*sigmaj)
pdf(file= paste0(strPlotExportPath,"/stat/stat_3dmodel_SpeedVsTurn_Covar.pdf"),width=14,height=7,
    title="Covariance in 3D statistical model for Capture Strike speed / Undershoot Ratio / Distance to Prey")
  plotCovar_SpeedVsTurnRatio(draw_LF,draw_NF,draw_DF,ntail)
dev.off()

## Distance(3)xSpeed(2) Covariance Coeff: Calc as rho=Cij/(sigmai*sigmaj)


pdf(file= paste0(strPlotExportPath,"/stat/stat_3dmodel_SpeedVsDistance_Covar.pdf"),width=14,height=7,
    title="Covariance in 3D statistical model for Capture Strike speed / Undershoot Ratio / Distance to Prey")
  plotCovar_SpeedVsDistance(draw_LF,draw_NF,draw_DF,ntail)
dev.off()


## TurnRatio(1)xDistance(3) Covariance Coeff: Calc as rho=Cij/(sigmai*sigmaj)

pdf(file= paste0(strPlotExportPath,"/stat/stat_3dmodel_TurnVsDistance_Covar.pdf"),width=14,height=7,
    title="Covariance in 3D statistical model for Capture Strike speed / Undershoot Ratio / Distance to Prey")
  plotCovar_DistanceVsTurnRatio(draw_LF,draw_NF,draw_DF,ntail)
dev.off()
#### END OF ALL Bouts Capture Covariance ##



### Show covariance In the High Speed Capture Cluster ##
#### COVARIANCE PLOTS OVER FAST Capture Bouts  ####
ntail <- 5000
### Lid,Matrix i,Matrix j,Sample,chain

load(file = paste0(strDataExportDir,"stat_Larval3DGaussianBehaviouModel_FastCapturesOnly_RJags.RData"))

## Turn-Ratio(1)xSpeed(2) Covariance Coeff: Calc as rho=Cij/(sigmai*sigmaj)
pdf(file= paste0(strPlotExportPath,"/stat/stat_3dmodel_Fast_SpeedVsTurn_Covar.pdf"),width=14,height=7,
    title="Covariance in 3D statistical model for FAST Capture Strike speed / Undershoot Ratio / Distance to Prey")
  plotCovar_SpeedVsTurnRatio(draw_LF,draw_NF,draw_DF,ntail)
dev.off()


pdf(file= paste0(strPlotExportPath,"/stat/stat_3dmodel_Fast_SpeedVsDistance_Covar.pdf"),width=14,height=7,
    title="Covariance in 3D statistical model for Capture Strike speed / Undershoot Ratio / Distance to Prey")
  plotCovar_SpeedVsDistance(draw_LF,draw_NF,draw_DF,ntail)
dev.off()


## TurnRatio(1)xDistance(3) Covariance Coeff: Calc as rho=Cij/(sigmai*sigmaj)
pdf(file= paste0(strPlotExportPath,"/stat/stat_3dmodel_Fast_TurnVsDistance_Covar.pdf"),width=14,height=7,
    title="Covariance in 3D statistical model for FAST Capture Strike speed / Undershoot Ratio / Distance to Prey")
  plotCovar_DistanceVsTurnRatio(draw_LF,draw_NF,draw_DF,ntail)
dev.off()
#### END OF COVAR Over Fast Capture Only ##


#### COVARIANCE PLOTS OVER SLOW Capture Bouts  ####
load(file = paste0(strDataExportDir,"stat_Larval3DGaussianBehaviouMode_SLOWCapturesOnly_RJags.RData"))

## Turn-Ratio(1)xSpeed(2) Covariance Coeff: Calc as rho=Cij/(sigmai*sigmaj)
pdf(file= paste0(strPlotExportPath,"/stat/stat_3dmodel_Slow_SpeedVsTurn_Covar.pdf"),width=14,height=7,
    title="Covariance in 3D statistical model for SLOW Capture Strike speed / Undershoot Ratio / Distance to Prey")
  plotCovar_SpeedVsTurnRatio(draw_LF,draw_NF,draw_DF,ntail)
dev.off()

## Distance(3)xSpeed(2) Covariance Coeff: Calc as rho=Cij/(sigmai*sigmaj)


pdf(file= paste0(strPlotExportPath,"/stat/stat_3dmodel_Slow_SpeedVsDistance_Covar.pdf"),width=14,height=7,
    title="Covariance in 3D statistical model for SLOW Capture Strike speed / Undershoot Ratio / Distance to Prey")
  plotCovar_SpeedVsDistance(draw_LF,draw_NF,draw_DF,ntail)
dev.off()


## TurnRatio(1)xDistance(3) Covariance Coeff: Calc as rho=Cij/(sigmai*sigmaj)

pdf(file= paste0(strPlotExportPath,"/stat/stat_3dmodel_Slow_TurnVsDistance_Covar.pdf"),width=14,height=7,
    title="Covariance in 3D statistical model for SLOW Capture Strike speed / Undershoot Ratio / Distance to Prey")
  plotCovar_DistanceVsTurnRatio(draw_LF,draw_NF,draw_DF,ntail)
dev.off()
#### END OF SLOW Capture Covariance ##



##
##### 3D OPENGL plot - Balls Model of Group Mean behaviour #####
##
library( rgl )
ntail <- 700

##Prepare Data
datMu3D <-  data.frame( cbind.data.frame(
                        TurnR=tail(c(draw_LF$muG[,1,,]), ntail) ,
                        CSpeed=tail(c(draw_LF$muG[,2,,]), ntail) ,
                        Dist=tail(c(draw_LF$muG[,3,,]), ntail), col=colourLegL[2],pch=pchL[4],group="LF" )  )

datMu3D <- rbind(datMu3D,
                  data.frame( cbind.data.frame(
                   TurnR=tail(c(draw_NF$muG[,1,,]), ntail) ,
                   CSpeed=tail(c(draw_NF$muG[,2,,]), ntail),
                   Dist=tail(c(draw_NF$muG[,3,,]), ntail), col=colourLegL[1],pch=pchL[6]),group="NF"  )
              )

datMu3D <- rbind(datMu3D,
                 data.frame( cbind.data.frame(
                   TurnR=tail(c(draw_DF$muG[,1,,1]), ntail) ,
                   CSpeed=tail(c(draw_DF$muG[,2,,1]), ntail),
                   Dist=tail(c(draw_DF$muG[,3,,1]), ntail), col=colourLegL[3],pch=pchL[5],group="DF")  )
                 )
datMu3D$col  <- as.character( datMu3D$col)
  
  ##Open Window And Plot
  open3d()
  bbox <- par3d('bbox') 
  rgl::plot3d( x=datMu3D$TurnR, y=datMu3D$CSpeed, z=datMu3D$Dist, col = datMu3D$col, type = "s", radius = 0.5,
               #xlab="Turn Ratio", ylab="Capture Speed (mm/sec)",zlab="Distance to prey (mm)",
               xlab="", ylab="",zlab="",
               xlim=c(0.5,1.5), ylim=c(10,50), zlim=c(0,0.8),
               box = TRUE ,aspect = TRUE,axes=FALSE
               #,expand = 1.5
               )
  box3d()
  title3d(main=NULL)
  rgl::axis3d('x+-',at=c(0.5,0.8,1,1.2,1.5))
  rgl::axis3d('z-+',at=seq(0.0,0.7,len=8))
  rgl::axis3d('y+-',at=seq(50,10,len=5),labels=rev(seq(10,50,len=5))) 
  
  
  #mtext3d("Turn Ratio", "x+-", line = 2, at = NULL, pos = NA) 
  #mtext3d("Capture Speed (mm/sec)", "y+-", line = 2, at = NULL, pos = NA) 
  #mtext3d("Distance (mm)", "z-+", line = 4, at = NULL, pos = NA,angle=90) 
  
  
  rgl::rgl.viewpoint(0,-60,fov=35,type = c("userviewpoint") )
  rgl::rgl.viewpoint(0,0)

#decorate3d(
#           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
#           top = TRUE, aspect = FALSE, expand = 1.03,cex=cex)

###Warning External Editing of PDF may fail. a Ghostscript conversion can fix this: 
#Use : gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/screen -dNOPAUSE -dBATCH  -dQUIET -sOutputFile=output.pdf fig7_Modelballs3D_SpeedVsTurn_view.pdf
#rgl::rgl.postscript( paste0(strPlotExportPath,"/fig7_Modelballs3D_Perspective_TurnVsSpeed_view.pdf"),fmt="pdf",drawText = FALSE )

rgl::rgl.postscript( paste0(strPlotExportPath,"/fig7_Modelballs3D_Perspective_TurnVsSpeed_view5.pdf"),fmt="pdf",drawText = FALSE )
rgl::rgl.postscript( paste0(strPlotExportPath,"/fig7_Modelballs3D_Perspective_TurnVsSpeed_doc5.tex"),fmt="tex",drawText = TRUE )
rgl::rgl.snapshot( paste0(strPlotExportPath,"/fig7_Modelballs3D__Perspective_TurnVsSpeed_view5.png"),fmt="png" )

## I use the tex Axis doc to combine 3d Fig with axis text. I then compile a pdf, which I import into inkScape (using Cairo to rasterize image), so I can adjust size, add axis labels etc etc.
## END OF 3D plot Messing with exporting
## 


zLL <- kde2d(c(tail(draw_LF$muG[,1,,1],ntail)), c(tail(draw_LF$muG[,2,,1],ntail)),n=180)
zNL <- kde2d(c(tail(draw_NF$muG[,1,,1],ntail)), c(tail(draw_NF$muG[,2,,1],ntail)),n=180)
zDL <- kde2d(c(tail(draw_DF$muG[,1,,1],ntail)), c(tail(draw_DF$muG[,2,,1],ntail)),n=180)
#zALL <- kde2d(c(tail(draw_ALL$mu[,1,,1],ntail)), c(tail(draw_ALL$mu[,2,,1],ntail)),n=80)
# 
# zLLD <- kde2d(c(tail(draw_LF$mu[,1,,],ntail)), c(tail(draw_LF$mu[,3,,],ntail)),n=180)
# zNLD <- kde2d(c(tail(draw_NF$mu[,1,,],ntail)), c(tail(draw_NF$mu[,3,,],ntail)),n=180)
# zDLD <- kde2d(c(tail(draw_DF$mu[,1,,],ntail)), c(tail(draw_DF$mu[,3,,],ntail)),n=180)
# 
# 
zLLS <- kde2d(c(tail(draw_LF$muG[,3,,1],ntail)), c(tail(draw_LF$muG[,2,,1],ntail)),n=180)
zNLS <- kde2d(c(tail(draw_NF$muG[,3,,1],ntail)), c(tail(draw_NF$muG[,2,,1],ntail)),n=180)
zDLS <- kde2d(c(tail(draw_DF$muG[,3,,1],ntail)), c(tail(draw_DF$muG[,2,,1],ntail)),n=180)

####  Detailed 2D Section FIGURE Of Each Pair   #####
## PLot Model 2D Section of Estimated Means ##
## Open Output PDF 
pdf(file= paste(strPlotExportPath,strModelPDFFileName,sep=""),width=14,height=7,
    title="A 3D statistical model for Capture Strike speed / Turn Ratio / Distance to Prey")
    
    ### Show Speed Fit ###
    outer = FALSE
    line = 1 ## SubFig Label Params
    lineAxis = 2.7
    lineXAxis = 3.0
    cex = 1.4
    adj  = 3.5
    padj <- -8.0
    las <- 1
    
    
    layout(matrix(c(1,2),1,2, byrow = TRUE))
    ##Margin: (Bottom,Left,Top,Right )
    par(mar = c(3.9,4.7,2,1))
    
    ## Plot the mean of the 2D Models ##
    ##Collect Draws from all chains
    plot(datMu3D$TurnR ,datMu3D$CSpeed,col=datMu3D$col,pch=datMu3D$pch, xlim=c(0.5,1.5),ylim=c(10,50),ylab=NA,xlab=NA,cex=cex,cex.axis=cex  )
    
    mtext(side = 1,cex=cex, line = lineAxis, expression("Turn ratio" )) #["~gamma~"]"
    mtext(side = 2,cex=cex, line = lineAxis, expression("Capture Speed (mm/sec)  " ))
    #mtext("A",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)
    
    contour(zDL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1,lty=2)
    contour(zLL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1,lty=2)
    contour(zNL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1,lty=2)
    
    legend("topright",
           legend=c(  expression (),
                      bquote(NF[""]  ),
                      bquote(LF[""]  ),
                      bquote(DF[""]  )
                      #, bquote(All ~ '#' ~ .(ldata_ALL$N)  )
                      ),
           pch=c(pchL[6],pchL[4],pchL[5]), col=colourLegL,cex=cex)

    ## Distance To Prey Vs Speed ##
    plot(datMu3D$Dist,datMu3D$CSpeed,col=datMu3D$col,pch=datMu3D$pch,  xlim=c(0,0.5),ylim=c(10,50),ylab=NA,xlab=NA ,cex=cex,cex.axis=cex )
  
    mtext(side = 1,cex=cex, line = lineAxis, expression("Distance to prey (mm)  " ))
    mtext(side = 2,cex=cex, line = lineAxis, expression(" Capture Speed (mm/sec)" ))
    
    contour(zDLS, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1,lty=2)
    contour(zLLS, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1,lty=2)
    contour(zNLS, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1,lty=2)

dev.off()
       
 






#### Show Covar Of Undershoot to Speed  Membership
plot(dNLb_rhoUS,col=colourLegL[1],xlim=c(-1,1),ylim=c(0.4,10),lwd=4,lty=1,main=NA,xlab=NA,ylab=NA,cex=cex,cex.axis=cex )
lines(dLLb_rhoUS,col=colourLegL[2],lwd=3,lty=2)
lines(dDLb_rhoUS,col=colourLegL[3],lwd=3,lty=3)


##### Individual Rand Vars Fit ###
## Capture Speeds ##
pdf(file= paste(strPlotExportPath,strCaptSpeedDensityPDFFileName ,sep=""))
par(mar = c(3.9,4.3,1,1))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
npchain<-3
plotCaptureSpeedFit(datTurnVsStrikeSpeed_NL,draw_NF,1,npchain)
title(main="Model capture Speed")
plotCaptureSpeedFit(datTurnVsStrikeSpeed_LL,draw_LF,2,npchain)
plotCaptureSpeedFit(datTurnVsStrikeSpeed_DL,draw_DF,3,npchain)

dev.off()

## Undershoot  ##
pdf(file= paste(strPlotExportPath,strUndershootDensityPDFFileName ,sep=""))
par(mar = c(3.9,4.3,1,1))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
plotUndeshootClusterFit(datTurnVsStrikeSpeed_NL,draw_NF,1)
title(main="Model undershoot on 1st turn to prey")
plotUndeshootClusterFit(datTurnVsStrikeSpeed_LL,draw_LF,2)
plotUndeshootClusterFit(datTurnVsStrikeSpeed_DL,draw_DF,3)
dev.off()

## Distance ##
pdf(file= paste(strPlotExportPath,strDistanceDensityPDFFileName ,sep=""))
par(mar = c(3.9,4.3,1,1))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
plotDistanceClustFit(datTurnVsStrikeSpeed_NL,draw_NF,1)
title(main="Model distance from prey prior to capture")
plotDistanceClustFit(datTurnVsStrikeSpeed_LL,draw_LF,2)
plotDistanceClustFit(datTurnVsStrikeSpeed_DL,draw_DF,3)

dev.off()
## plot 
##plot(xquant,dnorm(xquant,mean=tail(draw_NF$mu[2,2,,1],1),sd=tail(draw_NF$sigma[2,2,,1],1)),type='l',col=colourH[1],lty=1 )



## SLow Clust
clust <- 2
###Undershoot-Speed Covar
plot(density(draw_LF$sigma[clust,1,,1]*draw_LF$sigma[clust,2,,1]*draw_LF$rho[clust,1,,1]),
     col=colourLegL[2],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4))
lines(density(draw_NF$sigma[clust,1,,1]*draw_NF$sigma[clust,2,,1]*draw_NF$rho[clust,1,,1]),
     col=colourLegL[1],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4))
lines(density(draw_DF$sigma[clust,1,,1]*draw_DF$sigma[clust,2,,1]*draw_DF$rho[clust,1,,1]),
      col=colourLegL[3],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4))

###Speed Distance
plot(density(draw_LF$sigma[clust,3,,1]*draw_LF$sigma[clust,2,,1]*draw_LF$rho[clust,2,,1]),
     col=colourLegL[2],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4))
lines(density(draw_NF$sigma[clust,3,,1]*draw_NF$sigma[clust,2,,1]*draw_NF$rho[clust,2,,1]),
      col=colourLegL[1],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4))
lines(density(draw_DF$sigma[clust,3,,1]*draw_DF$sigma[clust,2,,1]*draw_DF$rho[clust,2,,1]),
      col=colourLegL[3],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4))


#mcmc_samples <- coda.samples(jags_model, c("mu", "rho", "sigma", "x_rand"),                             n.iter = 5000)


####         PLOT EMPIRICAL                              ##############
###        UNdershoot Vs Capture speed               ###
densNL <-  kde2d(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed,n=80)
densLL <-  kde2d(datTurnVsStrikeSpeed_LL$Undershoot, datTurnVsStrikeSpeed_LL$CaptureSpeed,n=80)
densDL <-  kde2d(datTurnVsStrikeSpeed_DL$Undershoot, datTurnVsStrikeSpeed_DL$CaptureSpeed,n=80)

covLL <- cov( 1/datTurnVsStrikeSpeed_LL$Undershoot,datTurnVsStrikeSpeed_LL$CaptureSpeed)
covDL <- cov( 1/datTurnVsStrikeSpeed_DL$Undershoot,datTurnVsStrikeSpeed_DL$CaptureSpeed)
covNL  <- cov( 1/datTurnVsStrikeSpeed_NL$Undershoot,datTurnVsStrikeSpeed_NL$CaptureSpeed)


pdf(file= paste(strPlotExportPath,"distal",strDataPDFFileName,sep=""))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.5,1,1))

plot(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed,col=colourLegL[1],
     xlab=NA,ylab=NA,ylim=c(0,60),xlim=c(0,2),main=NA,cex=cex)
lFit <- lm(datTurnVsStrikeSpeed_NL$CaptureSpeed ~ datTurnVsStrikeSpeed_NL$Undershoot)
abline(lFit,col=colourLegL[1],lwd=3.0) ##Fit Line / Regression
contour(densNL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[4],lty=2,lwd=1)
legend("topright",
       legend=paste("NF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ) ,cex=cex)  #prettyNum(digits=3, cov(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed)

plot(datTurnVsStrikeSpeed_LL$Undershoot, datTurnVsStrikeSpeed_LL$CaptureSpeed,col=colourLegL[2],
     ylim=c(0,60),xlim=c(0,2),xlab=NA,ylab=NA,cex=cex)
lFit <- lm(datTurnVsStrikeSpeed_LL$CaptureSpeed ~ datTurnVsStrikeSpeed_LL$Undershoot)
abline(lFit,col=colourLegL[2],lwd=3.0) ##Fit Line / Regression
contour(densLL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[4],lty=2,lwd=1)
mtext(side = 2,cex=cex, line = lineAxis-0.7, expression("Capture Speed (mm/sec) " ))
legend("topright",
       legend=paste("LF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex ) 


plot(datTurnVsStrikeSpeed_DL$Undershoot, datTurnVsStrikeSpeed_DL$CaptureSpeed,col=colourLegL[3],ylim=c(0,60),xlim=c(0,2),
     xlab=NA,ylab=NA,main=NA,cex=cex)
lFit <- lm(datTurnVsStrikeSpeed_DL$CaptureSpeed ~ datTurnVsStrikeSpeed_DL$Undershoot)
abline(lFit,col=colourLegL[3],lwd=3.0) ##Fit Line / Regression
contour(densDL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[4],lty=2,lwd=1)
mtext(side = 1,cex=cex, line = lineXAxis, expression("Turn ratio ["~gamma~"]" ))
legend("topright",
       legend=paste("DF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex ) 


dev.off()



## EMPIRICAL - UNdeshoot vs Prey Distance 
pdf(file= paste(strPlotExportPath,"/stat/UndershootAnalysis/fig7-UndershootDistanceCV_Distal_scatter.pdf",sep=""))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(4.5,4.3,0.5,1))

plot(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$DistanceToPrey,col=colourLegL[1],
     xlab=NA,ylab=NA,ylim=c(0,1.0),xlim=c(0,2),main=NA,cex=cex)
lFit <- lm(datTurnVsStrikeSpeed_NL$DistanceToPrey ~ datTurnVsStrikeSpeed_NL$Undershoot)
abline(lFit,col=colourLegL[1],lwd=3.0) ##Fit Line / Regression
legend("topright",
       legend=paste("NF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex )  #prettyNum(digits=3, cov(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed)

plot(datTurnVsStrikeSpeed_LL$Undershoot, datTurnVsStrikeSpeed_LL$DistanceToPrey,col=colourLegL[2],
     ylim=c(0,1),xlim=c(0,2.0),xlab=NA,ylab=NA,cex=cex)
lFit <- lm(datTurnVsStrikeSpeed_LL$DistanceToPrey ~ datTurnVsStrikeSpeed_LL$Undershoot)
abline(lFit,col=colourLegL[2],lwd=3.0) ##Fit Line / Regression
mtext(side = 2,cex=cex, line = 2.2, expression("Distance to prey  (mm) " ))
legend("topright",
       legend=paste("LF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex ) 


plot(datTurnVsStrikeSpeed_DL$Undershoot, datTurnVsStrikeSpeed_DL$DistanceToPrey,col=colourLegL[3],
     ylim=c(0,1.0),xlim=c(0,2),   xlab=NA,ylab=NA,main=NA,cex=cex)
lFit <- lm(datTurnVsStrikeSpeed_DL$DistanceToPrey ~ datTurnVsStrikeSpeed_DL$Undershoot)
abline(lFit,col=colourLegL[3],lwd=3.0) ##Fit Line / Regression
mtext(side = 1,cex=cex, line = lineXAxis, expression("Turn ratio ["~gamma~"]" ))
legend("topright",
       legend=paste("DF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ) ,cex=cex) 

dev.off()


### Onset/ DETECTION Angle Density supplementary angle figure
pdf(file= paste(strPlotExportPath,"/stat/fig7S1-stat_3Dmodel_TurnRatioVsSpeedAndDistance.pdf",sep=""))
  ### Show Speed Fit ###
  outer = FALSE
  line = 1 ## SubFig Label Params
  lineAxis = 2.7
  lineXAxis = 3.0
  cex = 1.4
  adj  = 3.5
  padj <- -8.0
  las <- 1
  
  ##Margin: (Bottom,Left,Top,Right )
  par(mar = c(4.5,4.3,0.5,1))
  
  plot(density(lFirstBoutPoints$NL[,"OnSetAngleToPrey"],bw=10),col=colourLegL[1],xlim=c(-120.0,120),lwd=3,lty=1,main=NA,xlab=NA,ylab=NA)
  lines(density(lFirstBoutPoints$LL[,"OnSetAngleToPrey"],bw=10),col=colourLegL[2],xlim=c(-120.0,120),lwd=3,lty=2)
  lines(density(lFirstBoutPoints$DL[,"OnSetAngleToPrey"],bw=10),col=colourLegL[3],xlim=c(-120.0,120),lwd=3,lty=3,main=NA)
  
  legend("topleft",
         legend=c(  expression (),
                    bquote(NF[""] ~ '#' ~ .(NROW(lFirstBoutPoints$NL))  ),
                    bquote(LF[""] ~ '#' ~ .(NROW(lFirstBoutPoints$LL))  ),
                    bquote(DF[""] ~ '#' ~ .(NROW(lFirstBoutPoints$DL))  )
         ), 
         col=colourLegL,lty=c(1,2,3,4),lwd=3,cex=cex)
  
  mtext(side = 2,cex=cex, line = lineAxis, expression("Density function") )
  mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Prey azimuth upon detection  " ) )  )
  #mtext("B",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=c
dev.off()


### Plot Initial Distance to Prey Vs Final Distance Scatter
plot(lFirstBoutPoints$LL[,"OnSetDistanceToPrey"],lFirstBoutPoints$LL[,"DistanceToPrey"],col=colourLegL[2],pch=pchL[2])
points(lFirstBoutPoints$NL[,"OnSetDistanceToPrey"],lFirstBoutPoints$NL[,"DistanceToPrey"],col=colourLegL[1],pch=pchL[1])
points(lFirstBoutPoints$DL[,"OnSetDistanceToPrey"],lFirstBoutPoints$DL[,"DistanceToPrey"],col=colourLegL[3],pch=pchL[3])



### Onset/ DETECTION Angle Density supplementary angle figure
pdf(file= paste(strPlotExportPath,"/stat/UndershootAnalysis/fig5S2-DetectionAngleVsDistance_scatter.pdf",sep=""))
  plot(lFirstBoutPoints$LL[,"OnSetAngleToPrey"],lFirstBoutPoints$LL[,"OnSetDistanceToPrey"],col=colourLegL[2],pch=pchL[2],xlim=c(-120.0,120),lwd=2,lty=1,main=NA,xlab=NA,ylab=NA)
  points(lFirstBoutPoints$NL[,"OnSetAngleToPrey"],lFirstBoutPoints$NL[,"OnSetDistanceToPrey"],col=colourLegL[1],pch=pchL[1],lwd=2)
  points(lFirstBoutPoints$DL[,"OnSetAngleToPrey"],lFirstBoutPoints$DL[,"OnSetDistanceToPrey"],col=colourLegL[3],pch=pchL[3],lwd=2)
  mtext(side = 2,cex=cex, line = lineAxis, expression("Prey distance upon detection (mm)") )
  mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Prey azimuth upon detection (deg)  " ) )  )
dev.off()
##




############# Plot Position Of Prey Prior Capture Bout 

pdf(file= paste(strPlotExportPath,"/PreyPositionPriorCapture_Validated.pdf",sep=""))

plotCaptureBoutPreyPositions()
dev.off()