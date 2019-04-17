### Kostas Lagogiannis 2019-04-15 
## Discovered relationship between the last bout speed - ie Capture speed and the undershoot ratio -
## higher undershoot predicts higher capture speeds, while undershoot also seems to predict higher distance from prey 
##  suggesting that LF stays further away from prey, so it does stronger capture bouts and it is undershoot that allows it to do it.
## Their ability to judge distance is also revealed in the the eye vergence prior to capture, where there is a relationship between EyeV and distance to prey   is shown 

### Stat Model on Capture speed vs undershoot


source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")


strRegisterDataFileName <- paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_SetC",".rds",sep="") #Processed Registry on which we add 


##



### PLOT ####
########################################################
###        UNdershoot Vs Capture speed               ###
datTurnVsStrikeSpeed_NL <- data.frame( cbind(Undershoot=lFirstBoutPoints$NL[,"Turn"]/lFirstBoutPoints$NL[,"OnSetAngleToPrey"],CaptureSpeed=lFirstBoutPoints$NL[,"CaptureSpeed"]) )
datTurnVsStrikeSpeed_LL <- data.frame( cbind(Undershoot=lFirstBoutPoints$LL[,"Turn"]/lFirstBoutPoints$LL[,"OnSetAngleToPrey"],CaptureSpeed=lFirstBoutPoints$LL[,"CaptureSpeed"]) )
datTurnVsStrikeSpeed_DL <- data.frame( cbind(Undershoot=lFirstBoutPoints$DL[,"Turn"]/lFirstBoutPoints$DL[,"OnSetAngleToPrey"],CaptureSpeed=lFirstBoutPoints$DL[,"CaptureSpeed"]) )
densNL <-  kde2d(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed,n=80)
densLL <-  kde2d(datTurnVsStrikeSpeed_LL$Undershoot, datTurnVsStrikeSpeed_LL$CaptureSpeed,n=80)
densDL <-  kde2d(datTurnVsStrikeSpeed_DL$Undershoot, datTurnVsStrikeSpeed_DL$CaptureSpeed,n=80)

covLL <- cov( 1/datTurnVsStrikeSpeed_LL$Undershoot,datTurnVsStrikeSpeed_LL$CaptureSpeed)
covDL <- cov( 1/datTurnVsStrikeSpeed_DL$Undershoot,datTurnVsStrikeSpeed_DL$CaptureSpeed)
covNL  <- cov( 1/datTurnVsStrikeSpeed_NL$Undershoot,datTurnVsStrikeSpeed_NL$CaptureSpeed)


pdf(file= paste(strPlotExportPath,"/stat/UndershootAnalysis/UndershootCaptureSpeed_scatter.pdf",sep=""))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.3,1,1))

plot(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed,col=colourP[1],
     xlab=NA,ylab=NA,ylim=c(0,60),xlim=c(0,2),main=NA)
lFit <- lm(datTurnVsStrikeSpeed_NL$CaptureSpeed ~ datTurnVsStrikeSpeed_NL$Undershoot)
abline(lFit,col=colourH[1],lwd=3.0) ##Fit Line / Regression
contour(densNL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[1],lty=2,lwd=3)
legend("topright",
       legend=paste("NF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ) )  #prettyNum(digits=3, cov(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed)

plot(datTurnVsStrikeSpeed_LL$Undershoot, datTurnVsStrikeSpeed_LL$CaptureSpeed,col=colourP[2],
     ylim=c(0,60),xlim=c(0,2),xlab=NA,ylab=NA)
lFit <- lm(datTurnVsStrikeSpeed_LL$CaptureSpeed ~ datTurnVsStrikeSpeed_LL$Undershoot)
abline(lFit,col=colourH[2],lwd=3.0) ##Fit Line / Regression
contour(densLL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[2],lty=2,lwd=3)
mtext(side = 2,cex=0.8, line = 2.2, expression("Capture Speed (mm/sec) " ))
legend("topright",
       legend=paste("LF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ) ) 


plot(datTurnVsStrikeSpeed_DL$Undershoot, datTurnVsStrikeSpeed_DL$CaptureSpeed,col=colourP[3],ylim=c(0,60),xlim=c(0,2),
     xlab=NA,ylab=NA,main=NA)
lFit <- lm(datTurnVsStrikeSpeed_DL$CaptureSpeed ~ datTurnVsStrikeSpeed_DL$Undershoot)
abline(lFit,col=colourH[3],lwd=3.0) ##Fit Line / Regression
contour(densDL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[3],lty=2,lwd=3)
mtext(side = 1,cex=0.8, line = 2.2, expression("Undershoot "~(gamma) ))
legend("topright",
       legend=paste("DF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ) ) 


dev.off()
