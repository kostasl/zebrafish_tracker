### Kostas Lagogiannis 2019-06-24 
## Learning steers the ontogeny of an efficient hunting sequence in zebrafish larvae.

## 3D Multivariate Gaussian Model capturing group behaviour
## ******** No clustering Between slow and Fast Swims ****** 
## I made this to complement the Clustering Method, so as to characterize the overall covariance structure


library(rjags)
library(here) ###Used ti Set correct paths

#library(runjags)
setwd(here())
#set_here(path="./pub") 
source("common_lib.R")
#source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
#source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
#source("HuntingEventAnalysis_lib.r")




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
  ## Fig Label Params
  outer = FALSE
  line = 1 
  lineAxis = 4.7
  lineXAxis = 3.5
  lineTitle = 0.5
  
  cex = 1.4
  adj  = 3.5
  padj <- 0.5
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
  ## Fig Label Params
  outer = FALSE
  line = 1 
  lineAxis = 4.7
  lineXAxis = 3.5
  lineTitle = 0.5
  
  cex = 1.4
  adj  = 3.5
  padj <- 0.5
  las <- 1
  
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
  ## Fig Label Params
  outer = FALSE
  line = 1 
  lineAxis = 4.7
  lineXAxis = 3.5
  lineTitle = 0.5
  cex = 1.4
  adj  = 3.5
  padj <- 0.5
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


##
##
steps <-15000
nchains <- 7
nthin <- 5
#str_vars <- c("mu","rho","sigma","x_rand") #Basic model 
str_vars <- c("mu","cov","x_rand","muG","tG","NLarv","Lid") #Mixture Model

getwd()
setwd("/home/kostasl/workspace/zebrafishtrack/Rplots")
##Contains Serial Larvae ID, that links each hunt event to particular larva 
ldata_LF  <-  readRDS(file=paste0("dat/huntEpisodeDataMergedWithLarvalSuccess_LF.rds") )
ldata_NF  <-  readRDS(file=paste0("dat/huntEpisodeDataMergedWithLarvalSuccess_NF.rds") )
ldata_DF  <-  readRDS(file=paste0("dat/huntEpisodeDataMergedWithLarvalSuccess_DF.rds") )
datHuntLarvaStat <- readRDS(file=paste0("dat/LarvaEmpiricalMeanHuntBehaviour.rds"))


lModelEst_LF <- getEstimatesPerLarva(draw_LF,stail)
lModelEst_NF <- getEstimatesPerLarva(draw_NF,stail)
lModelEst_DF <- getEstimatesPerLarva(draw_DF,stail)


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

save(draw_NF,draw_LF,draw_DF,file = paste0(strDataExportDir,"stat_Larval3DGaussianBehaviouModel_RJags.RData"))



schain <- 1:7 ## Chains used for data visualiazation
stail <- 500 ## Number Of  Chain Samples to Use for Plots - from the end of the chain


##CONVERGENCE - CHECK
##Each of draw_XX arrays is structured  as
# Idx: Lid,Variable (1 turnratio,2-speed,3-distance),Sample Row,Chain
plot(density(tail(draw_LF$muG[,3,,1],1000) ),type='l',main="Distance - Across Chains")
for (c in schain)
  lines(density(tail(draw_LF$muG[,3,,c],1000) ),col="red")

plot(density(tail(draw_LF$muG[,2,,1],1000) ),type='l',main="Speed- Across Chains")
for (c in schain)
  lines(density(tail(draw_LF$muG[,2,,c],1000) ),col="red")

plot(density(tail(draw_LF$muG[,1,,1],1000) ),type='l',main="Turn Ratio- Across Chains")
for (c in schain)
  lines(density(tail(draw_LF$muG[,1,,c],1000) ),col="red")



## ############################################### ##
#### 3D Display of Group Behaviour MODEL          ####
## ############################################### ##
  library( rgl )
  ntail <- 300
  
  ##Prepare Data
  datMu3D <-  data.frame( cbind.data.frame(
    TurnR=tail(draw_LF$muG[,1,,1], ntail) ,
    CSpeed=tail(draw_LF$muG[,2,,1], ntail),
    Dist=tail(draw_LF$muG[,3,,1], ntail), col=colourLegL[2],pch=pchL[4],group="LF" )  )
  
  datMu3D <- rbind(datMu3D,
                   data.frame( cbind.data.frame(
                     TurnR=tail(draw_NF$muG[,1,,1], ntail) ,
                     CSpeed=tail(draw_NF$muG[,2,,1], ntail),
                     Dist=tail(draw_NF$muG[,3,,1], ntail), col=colourLegL[1],pch=pchL[6]),group="NF"  )
  )
  
  datMu3D <- rbind(datMu3D,
                   data.frame( cbind.data.frame(
                     TurnR=tail(draw_DF$muG[,1,,1], ntail) ,
                     CSpeed=tail(draw_DF$muG[,2,,1], ntail),
                     Dist=tail(draw_DF$muG[,3,,1], ntail), col=colourLegL[3],pch=pchL[5],group="DF")  )
  )
  datMu3D$col  <- as.character( datMu3D$col)
  
  ##Open Window And Plot
  open3d()
  bbox <- par3d('bbox') 
  rgl::plot3d( x=datMu3D$TurnR, y=datMu3D$CSpeed, z=datMu3D$Dist, col = datMu3D$col, type = "s", radius = 0.5,
               #xlab="Turn Ratio", ylab="Capture Speed (mm/sec)",zlab="Distance to prey (mm)",
               xlab="Turn Ratio", ylab="Speed",zlab="Distance",
               xlim=c(0.5,1.5), ylim=c(10,50), zlim=c(0,0.8),
               box = TRUE ,aspect = TRUE,axes=FALSE
               #,expand = 1.5
  )
  box3d()
  title3d(main=NULL)
  rgl::axis3d('x+-',at=c(0.5,0.8,1,1.2,1.5))
  rgl::axis3d('z-+',at=seq(0.0,0.7,len=8))
  rgl::axis3d('y+-',at=seq(50,10,len=5),labels=rev(seq(10,50,len=5))) 
  
  
  rgl::rgl.viewpoint(0,-60,fov=35,type = c("userviewpoint") )
  rgl::rgl.viewpoint(0,0)
  


## ############################################# ##
## Speed Posterior Group Vs Individual Density
## ########################################## ##
## Live Fed
with(draw_LF,{
  plot(density(tail(muG[,2,,], 150) ),ylim=c(0,1),xlim=c(0,60),lwd=2,col="red",main="Compare Data to Model - Speed LF")
  for ( i in (1:NLarv[1] ) )
    lines( density( tail( mu[i,2,,],stail)),lty=2)
  ###Show Inferred Distribution
  lines(1:100,dnorm(1:100,mean=mean( tail(muG[,2,,3],stail)),sd=sqrt(mean( tail( 1/tG[,2,,1],stail))) ),col="purple",lwd=4)
})
##Compare To Empirical - Change group DF,LF,NF-- V Good Match!
lines(density(datHuntLarvaStat[datHuntLarvaStat$groupID==2,]$CaptureSpeed,bw=2),col="blue",lwd=3,lty=2)

##Dry Fed DF
##Speed Posterior Group Vs INdividual Density
with(draw_DF,{
  plot(density(tail(muG[,2,,], 150) ),ylim=c(0,1),xlim=c(0,60),lwd=2,col="red",main="Compare Data to Model  DF")
  for ( i in (1:NLarv[1] ) )
    lines( density( tail( mu[i,2,,],stail)),lty=2)
  ###Show Inferred Distribution
   lines(1:100,dnorm(1:100,mean=mean( tail(muG[,2,,3],stail)),sd=sqrt(mean( tail( 1/tG[,2,,1],stail)) )) ,col="purple",lwd=4)
})
##Compare To Empirical - Change group DF,LF,NF-- V Good Match!
lines(density(datHuntLarvaStat[datHuntLarvaStat$groupID==1,]$CaptureSpeed,bw=2),col="blue",lwd=2,lty=2)

##Not Fed
##Speed Posterior Group Vs INdividual Density
with(draw_NF,{
  plot(density(tail(muG[,2,,], 150) ),ylim=c(0,1),xlim=c(0,60),lwd=2,col="red",main="Speed NF")
  for ( i in (1:NLarv[1] ) )
    lines( density( tail( mu[i,2,,],stail)),lty=2)
  ###Show Inferred Distribution
  lines(1:100,dnorm(1:100,mean=mean( tail(muG[,2,,3],stail)),sd=sqrt(mean( tail( 1/tG[,2,,1],stail))) ),col="purple",lwd=4)
})
##Compare To Empirical - Change group DF,LF,NF-- V Good Match!
lines(density(datHuntLarvaStat[datHuntLarvaStat$groupID==3,]$CaptureSpeed,bw=2),col="blue",lwd=2,lty=2)



#### Distance Posterior Group Vs INdividual Density
with(draw_NF,{
  plot(density(tail(muG[,3,,1], 100) ),ylim=c(0,16),xlim=c(0,1),lwd=2,col="red")
  for ( i in (1:NLarv[1] ) )
    lines( density( tail( mu[i,3,,1],stail)),lty=2)
})

### Compare Group Model To Density Obtain through Mean Estimated Behaviour For Each Larva
### Obtain Estimated Mean Values For Each Larva & Plot Group Population
## Plot Distance Density
plot(density(sapply(tail(draw_LF$mu[,3,,],stail),mean)),col=colourLegL[2] ,lwd=2,main="Distance to Prey",ylim=c(0,6)) ##Mean Group Undershoot From Mean Of Each Larva
lines(density( unlist(lapply(lModelEst_LF[,"DistanceToPrey"],mean) ) ) )
lines(density(sapply(tail(draw_NF$mu[,3,,],stail),mean)),col=colourLegL[1] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva
lines(density( unlist(lapply(lModelEst_NF[,"DistanceToPrey"],mean) ) ) )
lines(density(sapply(tail(draw_DF$mu[1,3,,],stail) ,mean)),col=colourLegL[3] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva
lines(density( unlist(lapply(lModelEst_DF[,"DistanceToPrey"],mean) ) ) )


## Plot Speed Density
plot(density(sapply(tail(draw_LF$mu[,2,,],stail),mean)),col=colourLegL[2] ,lwd=2,main="Capture Speed",ylim=c(0,0.1)) ##Mean Group Undershoot From Mean Of Each Larva
lines(density(sapply(tail(draw_NF$mu[,2,,],stail),mean)),col=colourLegL[3] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva
lines(density(sapply(tail(draw_DF$mu[,2,,],stail) ,mean)),col=colourLegL[1] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva

## Plot Undershoot Density / Mean Sample point Across larva 
plot(density(sapply(tail(draw_LF$mu[,1,,],stail),mean)),col=colourLegL[2] ,lwd=2,main="Turn ratio") ##Mean Group Undershoot From Mean Of Each Larva
lines(density(sapply(tail(draw_NF$mu[,1,,],stail),mean)),col=colourLegL[3] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva
lines(density(sapply(tail(draw_DF$mu[,1,,],stail) ,mean)),col=colourLegL[1] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva


#### Plot Covariance Analysis ####
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

