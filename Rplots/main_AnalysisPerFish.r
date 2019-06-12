### Success Analysis Per Fish - 
## Using Labelled Data from SB - with sparse corrections by KL where wrong labels have been spotted
### KOstasl 13-08-18
source("DataLabelling/labelHuntEvents_lib.r")
source("config_lib.R")
#strProcDataFileName <-paste("setn-12","-HuntEvents-SB-ALL",sep="") ##To Which To Save After Loading
#strProcDataFileName <- paste("setn14-HuntEventsFixExpID-SB-Updated-Merged",sep="") ##To Which To Save After Loading
strProcDataFileName <- paste("setn15-HuntEvents-SB-Updated-Merged") ##To Which To Save After Loading
strProcDataFileName <- "setn15-HuntEvents-SB-Updated-Merged2"
message(paste(" Loading Hunt Event List to Process... "))
#load(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
datHuntEventAllGroupToLabel <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))
##<- datHuntEvent
groupsList <- c("DL","NL","LL") ##unique(datHuntEventAllGroupToLabel$groupID)
str_FilterLabel <- "UnLabelled"

datFishSuccessRate <- getHuntSuccessPerFish(datHuntEventAllGroupToLabel)
tblEventsTracked <- table(datHuntEventAllGroupToLabel$expID, datHuntEventAllGroupToLabel$markTracked,useNA="always" )
remove(datFishSuccessRateMerged)

##Histogram of Success Doesn t Show Any obvious differences
#layout(matrix(c(1,2,3), 3, 1 ,byrow=TRUE))
#hist(datFishSuccessRate[datFishSuccessRate$groupID == "LL",]$Success,breaks = seq(0,45,3),ylim=c(0,50),main="LL",col="#0000FFAA")
#hist(datFishSuccessRate[datFishSuccessRate$groupID == "LL",]$Fails,breaks = seq(0,45,3),ylim=c(0,50),main="LL",col="#FF0000AA",add=T)
#hist(datFishSuccessRate[datFishSuccessRate$groupID == "NL",]$Success,breaks = seq(0,45,3),ylim=c(0,50),main="NL",col="#0000FFAA")
#hist(datFishSuccessRate[datFishSuccessRate$groupID == "NL",]$Fails,breaks = seq(0,45,3),ylim=c(0,50),main="NL",col="#FF0000AA",add=T)
#hist(datFishSuccessRate[datFishSuccessRate$groupID == "DL",]$Success,breaks = seq(0,45,3),ylim=c(0,50),main="DL",col="#0000FFAA")
#hist(datFishSuccessRate[datFishSuccessRate$groupID == "DL",]$Fails,breaks = seq(0,45,3),ylim=c(0,50),main="DL",col="#FF0000AA",add=T)
#datFishSuccessRate[,"Success"] <- as.numeric(datFishSuccessRate[,"Success"])
#


datFishSuccessRateMerged <- cbind(datFishSuccessRate,markUnTrackable=data.frame(tblEventsTracked[row.names(datFishSuccessRate),1]),
      markTracked=data.frame(tblEventsTracked[row.names(datFishSuccessRate),2]),
      notTracked=data.frame(tblEventsTracked[row.names(datFishSuccessRate),3]))

names(datFishSuccessRateMerged)[(NCOL(datFishSuccessRate)+1):NCOL(datFishSuccessRateMerged)] <- c("markUnTrackable","markTracked","notTracked") ##Set Field Names - notTracked : Have not been detailed retracked yet

vScoreIdx        <- ((datFishSuccessRate[,"Success"]*datFishSuccessRate[,"Success"])/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails"]))
vEfficiencyRatio <- (datFishSuccessRate[,"Success"]/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails"]))
vEfficiencyRatio_Strike <- (datFishSuccessRate[,"Success"]/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails_WS"]))
vEfficiencyRatio_NStrike <- (datFishSuccessRate[,"Success"]/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails_NS"]))
#vScoreIdx[is.nan(vScoreIdx) ] <- 0
datFishSuccessRateMerged <- cbind(datFishSuccessRateMerged,vScoreIdx,vEfficiencyRatio,vEfficiencyRatio_Strike,vEfficiencyRatio_NStrike)

## Subset only the active/Larvae - ones that have hunted 
datFishSuccessRateActive <- datFishSuccessRateMerged[!is.nan(datFishSuccessRateMerged$vScoreIdx),]
datFishSuccessRateActive_WS <- datFishSuccessRateMerged[!is.nan(datFishSuccessRateMerged$vEfficiencyRatio_Strike),]
datFishSuccessRateActive_NS <- datFishSuccessRateMerged[!is.nan(datFishSuccessRateMerged$vEfficiencyRatio_NStrike),]






## Plot Density of Hunting POWER S^2/(S+F)
densDLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vScoreIdx)
densNLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vScoreIdx)
densLLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vScoreIdx)

cdfDLScore <- ecdf(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vScoreIdx)
cdfNLScore <- ecdf(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vScoreIdx)
cdfLLScore <- ecdf(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vScoreIdx)

dev.off() ##Clear Old plot

## Plot CDF HUNT POWER ##
pdf(file= paste(strPlotExportPath,"/stat/efficiency/fig3-ecdf_huntpower.pdf",sep=""))
par(mar = c(3.9,4.3,1,1))

plotHuntPowerData(datHuntEventAllGroupToLabel)

# plot(cdfNLScore,lty=2,lwd=3,col=colourLegL [1],cex=1.2,cex.axis=1.3,xlim=c(0,12),pch=pchL[1],ylim=c(0.03,1.01),
#      main=NA,ylab=NA,  xlab=NA)
# plot(cdfLLScore,add=T,lty=1,lwd=3,col=colourLegL[2],pch=pchL[2],ylim=c(0,1.01),cex=1.2)
# plot(cdfDLScore,add=T,lty=1,lwd=3,col=colourLegL[3],pch=pchL[3],ylim=c(0,1.01),cex=1.2)
# 
# mtext(side = 1,cex=1.2, line = 2.8, expression( "Hunt power " ~ N[S]^2/(N[S]+N[F]) ,paste("") )   )
# mtext(side = 2,cex=1.2, line = 2.8, expression("Cumulative function " ))
# 
# legend("bottomright",legend=paste(c("NL #","LL #","DL #"),c(densNLScore$n,densLLScore$n,densDLScore$n) ),
#        col = colourLegL,pch=pchL,cex=1.2)
dev.off()




strPlotFileName <- paste(strPlotExportPath,"/stat/HuntEfficiency_hist.pdf",sep="")
pdf(strPlotFileName,width = 16,height = 18 ,paper = "a4",onefile = TRUE );

## Plot Histogram of efficiency ##
layout(matrix(c(1,2,3), 3, 1 ,byrow=TRUE))
ptbreaks <- seq(from=0,to=1,by=1/10)
layout(matrix(c(1,2,3), 3, 1 ,byrow=TRUE))
hist(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vEfficiencyRatio,col=colourR[3],main=paste("DL #",NROW(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]))
     ,xlab="",breaks=ptbreaks,ylim=c(0,20))
hist(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vEfficiencyRatio,col=colourR[2],main=paste("LL #",NROW(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",])),
     xlab="",breaks=ptbreaks,ylim=c(0,20))
hist(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vEfficiencyRatio,col=colourR[1],main=paste("NL #",NROW(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]))
     ,xlab="Hunt  Efficiency (Success/(Fail+Succ.) )score",breaks=ptbreaks,ylim=c(0,20))
mtext("Efficiency","top",side=3)
dev.off()

###Hunt Power Histogram

strPlotFileName <- paste(strPlotExportPath,"/stat/HuntPower_hist.pdf",sep="")
pdf(strPlotFileName,width = 16,height = 18 ,paper = "a4",onefile = TRUE );
layout(matrix(c(1,2,3), 3, 1 ,byrow=TRUE))
hist(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vEfficiencyRatio*datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$HuntEvents,
     col=colourR[3],main=paste("DL #",NROW(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",])),
     xlab="",breaks=10,ylim=c(0,40))

hist(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vEfficiencyRatio*datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$HuntEvents,
     col=colourR[2],main=paste("LL #",NROW(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",])),
     xlab="",breaks=10,ylim=c(0,40))

hist(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vEfficiencyRatio*datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$HuntEvents,
     col=colourR[1],main=paste("NL #",NROW(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",])),
     xlab="",breaks=10,ylim=c(0,40))
####
dev.off()



###Scatter of hunt efficiency vs hunt count
strPlotFileName <- paste(strPlotExportPath,"/stat/HuntEfficiencyVsHuntCount_scatter.pdf",sep="")
pdf(strPlotFileName,width = 16,height = 18 ,paper = "a4",onefile = TRUE );


  plot(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vEfficiencyRatio,
       datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$HuntEvents,col=colourR[2],pch=3,
       ylab="Hunt Event Count",xlab="Efficiency")
  points(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vEfficiencyRatio,
              datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$HuntEvents,col=colourR[1],pch=16 )
  points(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vEfficiencyRatio,
         datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$HuntEvents,col=colourR[3],pch=21 )

vDL_expID <- rownames(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",])
vLL_expID <- rownames(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",])  

text(unlist(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vEfficiencyRatio)*1.01,
     unlist( datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$HuntEvents)*1.01,vDL_expID
    ,cex=0.7)
text(unlist(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vEfficiencyRatio)*1.01,
     unlist( datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$HuntEvents)*1.01,vLL_expID
     ,cex=0.7)

dev.off()



## Plot Efficiency  Density 
densDLEffScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vEfficiencyRatio,bw=0.05)
densNLEffScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vEfficiencyRatio,bw=0.05)
densLLEffScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vEfficiencyRatio,bw=0.05)


dev.off() ##Clear Old plot
strPlotFileName <- paste(strPlotExportPath,"/stat/HuntEfficiency_Density.pdf",sep="")
pdf(strPlotFileName,width = 16,height = 16 ,paper = "a4",onefile = TRUE );

plot(densNLEffScore,col=colourH[1],main="Hunt efficiency density",type="l",lwd=3,lty=1,ylim=c(0,3.0),xlim=c(0,1) )
lines(densLLEffScore,col=colourH[2],lwd=3,lty=2)
lines(densDLEffScore,col=colourH[3],lwd=3,lty=3,xlab="Hunting efficiency score" )
legend("topright",legend=paste(c("NL #","LL #","DL #"),
                               c(densNLEffScore$n,densLLEffScore$n,densDLEffScore$n) ),col = colourH,lty=c(1,2,3,4),lwd=3)

dev.off()

## Plot Efficiency based on Fails_ With A strike Density 
densDLEffScore_WS <- density(datFishSuccessRateActive_WS[datFishSuccessRateActive_WS$groupID == "DL",]$vEfficiencyRatio_Strike)
densNLEffScore_WS <- density(datFishSuccessRateActive_WS[datFishSuccessRateActive_WS$groupID == "NL",]$vEfficiencyRatio_Strike)
densLLEffScore_WS <- density(datFishSuccessRateActive_WS[datFishSuccessRateActive_WS$groupID == "LL",]$vEfficiencyRatio_Strike)

dev.off() ##
plot(densLLEffScore_WS,col=colourH[2],main="Hunt efficiency With Strike density",type="l",lwd=2,ylim=c(0,2.0),xlim=c(0,1))
lines(densNLEffScore_WS,col=colourH[3],lwd=2)
lines(densDLEffScore_WS,col=colourH[1],lwd=2,xlab="Hunting efficiency score" )
legend("topright",legend=paste(c("DL #","LL #","NL #"),c(densDLEffScore$n,densLLEffScore$n,densNLEffScore$n) ),fill = colourH)


## Plot Efficiency based on Fails_ With NO Strike Density 
densDLEffScore_NS <- density(datFishSuccessRateActive_NS[datFishSuccessRateActive_NS$groupID == "DL",]$vEfficiencyRatio_NStrike)
densNLEffScore_NS <- density(datFishSuccessRateActive_NS[datFishSuccessRateActive_NS$groupID == "NL",]$vEfficiencyRatio_NStrike)
densLLEffScore_NS <- density(datFishSuccessRateActive_NS[datFishSuccessRateActive_NS$groupID == "LL",]$vEfficiencyRatio_NStrike)

dev.off() ##
plot(densLLEffScore_NS,col=colourH[2],main="Hunt efficiency WithOUT a Strike ",type="l",lwd=2,ylim=c(0,2.0),xlim=c(0,1))
lines(densNLEffScore_NS,col=colourH[3],lwd=2)
lines(densDLEffScore_NS,col=colourH[1],lwd=2,xlab="Hunting efficiency score" )
legend("topright",legend=paste(c("DL #","LL #","NL #"),c(densDLEffScore_NS$n,densLLEffScore_NS$n,densNLEffScore_NS$n) ),fill = colourH)

### Plot No strike Failure
layout(matrix(c(1,2,3), 3, 1 ,byrow=TRUE))
hist(datFishSuccessRateActive_NS[datFishSuccessRateActive_NS$groupID == "DL",]$Fails_NS  )
hist(datFishSuccessRateActive_NS[datFishSuccessRateActive_NS$groupID == "LL",]$Fails_NS  )
hist(datFishSuccessRateActive_NS[datFishSuccessRateActive_NS$groupID == "NL",]$Fails_NS  )

## Plot Histograp of scor/Gaine ##
ptbreaks <- seq(from=0,to=max(vScoreIdx,na.rm=TRUE)+1,by=1)
layout(matrix(c(1,2,3), 3, 1 ,byrow=TRUE))
hist(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vScoreIdx,col=colourR[1],main=paste("DL #",NROW(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]))
     ,xlab="",breaks=ptbreaks,ylim=c(0,45))
hist(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vScoreIdx,col=colourR[2],main=paste("LL #",NROW(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",])),
     xlab="",breaks=ptbreaks,ylim=c(0,45))
hist(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vScoreIdx,col=colourR[3],main=paste("NL #",NROW(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]))
     ,xlab="Hunt efficiency score",breaks=ptbreaks,ylim=c(0,45))



plot(densLLScore,col=colourH[2],main="Hunt Power",
     xlab=expression(S^2/(S+F),paste("")) ,type="l",lwd=2,pch=pchL[2],ylim=c(0,0.7))
lines(densNLScore,col=colourH[3],lwd=2)
lines(densDLScore,col=colourH[1],lwd=2,pch=pchL[3])
legend("topright",legend=paste(c("DL #","LL #","NL #"),c(densDLScore$n,densLLScore$n,densNLScore$n) ),fill = colourH)


###### Plot Power Under  Strike Fails 
## Plot Density of Hunting POWER S^2/(S+F)
densDLScore_WS <- density(datFishSuccessRateActive_WS[datFishSuccessRateActive_WS$groupID == "DL",]$vEfficiencyRatio_Strike*datFishSuccessRateActive_WS[datFishSuccessRateActive_WS$groupID == "DL",]$Success)
densNLScore_WS <- density(datFishSuccessRateActive_WS[datFishSuccessRateActive_WS$groupID == "NL",]$vEfficiencyRatio_Strike*datFishSuccessRateActive_WS[datFishSuccessRateActive_WS$groupID == "NL",]$Success)
densLLScore_WS <- density(datFishSuccessRateActive_WS[datFishSuccessRateActive_WS$groupID == "LL",]$vEfficiencyRatio_Strike*datFishSuccessRateActive_WS[datFishSuccessRateActive_WS$groupID == "LL",]$Success)

cdfDLScore_WS <- ecdf(datFishSuccessRateActive_WS[datFishSuccessRateActive_WS$groupID == "DL",]$vEfficiencyRatio_Strike*datFishSuccessRateActive_WS[datFishSuccessRateActive_WS$groupID == "DL",]$Success)
cdfNLScore_WS <- ecdf(datFishSuccessRateActive_WS[datFishSuccessRateActive_WS$groupID == "NL",]$vEfficiencyRatio_Strike*datFishSuccessRateActive_WS[datFishSuccessRateActive_WS$groupID == "NL",]$Success)
cdfLLScore_WS <- ecdf(datFishSuccessRateActive_WS[datFishSuccessRateActive_WS$groupID == "LL",]$vEfficiencyRatio_Strike*datFishSuccessRateActive_WS[datFishSuccessRateActive_WS$groupID == "LL",]$Success)


dev.off() ##Plot How Hunt Power is modified when comparing to Failures With A strike Only
plot(densLLScore_WS,col=colourH[2],main="Hunt Power On Fail With Strike",
     xlab=expression(S^2/(S+F),paste("")) ,type="l",lwd=2,ylim=c(0,0.7))
lines(densNLScore_WS,col=colourH[3],lwd=2)
lines(densDLScore_WS,col=colourH[1],lwd=2)
legend("topright",legend=paste(c("DL #","LL #","NL #"),c(densDLScore$n,densLLScore$n,densNLScore$n) ),fill = colourH)

## Plot CDF ##
strPlotName = paste(strPlotExportPath,"/stat/efficiency/ecdf_huntpower_WS.pdf",sep="")
pdf(strPlotName,width=8,height=8,title="Compare Hunt Power Between Groups With Strike",onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))
plot(cdfDLScore_WS,lty=2,lwd=1,col=colourH[1],xlim=c(0,12),xlab=expression(S^2/(S+F[WS]),paste("")),cex.axis=1.5,main="CDF Hunt power with Strike" )
plot(cdfLLScore_WS,add=T,lty=1,lwd=2,col=colourH[2])
plot(cdfNLScore_WS,add=T,lty=1,lwd=2,col=colourH[3])
legend("bottomright",legend=paste(c("DL #","LL #","NL #"),c(densDLScore_WS$n,densLLScore_WS$n,densNLScore_WS$n) ),fill = colourH)
dev.off()

###### Plot Power Under NO Strike Fails 
## Plot Density of Hunting POWER S^2/(S+F)
densDLScore_NS <- density(datFishSuccessRateActive_NS[datFishSuccessRateActive_NS$groupID == "DL",]$vEfficiencyRatio_NStrike*datFishSuccessRateActive_NS[datFishSuccessRateActive_NS$groupID == "DL",]$Success)
densNLScore_NS <- density(datFishSuccessRateActive_NS[datFishSuccessRateActive_NS$groupID == "NL",]$vEfficiencyRatio_NStrike*datFishSuccessRateActive_NS[datFishSuccessRateActive_NS$groupID == "NL",]$Success)
densLLScore_NS <- density(datFishSuccessRateActive_NS[datFishSuccessRateActive_NS$groupID == "LL",]$vEfficiencyRatio_NStrike*datFishSuccessRateActive_NS[datFishSuccessRateActive_NS$groupID == "LL",]$Success)

cdfDLScore_NS <- ecdf(datFishSuccessRateActive_NS[datFishSuccessRateActive_NS$groupID == "DL",]$vEfficiencyRatio_NStrike*datFishSuccessRateActive_NS[datFishSuccessRateActive_NS$groupID == "DL",]$Success)
cdfNLScore_NS <- ecdf(datFishSuccessRateActive_NS[datFishSuccessRateActive_NS$groupID == "NL",]$vEfficiencyRatio_NStrike*datFishSuccessRateActive_NS[datFishSuccessRateActive_NS$groupID == "NL",]$Success)
cdfLLScore_NS <- ecdf(datFishSuccessRateActive_NS[datFishSuccessRateActive_NS$groupID == "LL",]$vEfficiencyRatio_NStrike*datFishSuccessRateActive_NS[datFishSuccessRateActive_NS$groupID == "LL",]$Success)

dev.off() ##Plot How Hunt Power is modified when comparing to Failures With A strike Only
plot(densLLScore_NS,col=colourH[2],main="Hunt Power On Fail With NO Strike",
     xlab=expression(S^2/(S+F),paste("")) ,type="l",lwd=2,ylim=c(0,0.7))
lines(densNLScore_NS,col=colourH[3],lwd=2)
lines(densDLScore_NS,col=colourH[1],lwd=2)
legend("topright",legend=paste(c("DL #","LL #","NL #"),c(densDLScore_NS$n,densLLScore_NS$n,densNLScore_NS$n) ),fill = colourH)

## Plot CDF POWER not Strike##
strPlotName = paste(strPlotExportPath,"/stat/efficiency/ecdf_huntpower_NS.pdf",sep="")
pdf(strPlotName,width=8,height=8,title="Compare Hunt Power Between Groups With NO Strike",onefile = TRUE) #c
plot(cdfDLScore_NS,lty=2,lwd=1,col=colourH[1],xlim=c(0,12),xlab=expression(S^2/(S+F[NS]),paste("")),cex.axis=1.5,main="CDF Hunt power with NO Strike" )
plot(cdfLLScore_NS,add=T,lty=1,lwd=2,col=colourH[2])
plot(cdfNLScore_NS,add=T,lty=1,lwd=2,col=colourH[3])
legend("bottomright",legend=paste(c("DL #","LL #","NL #"),c(cdfDLScore_NS$n,cdfLLScore_NS$n,cdfNLScore_NS$n) ),fill = colourH)
dev.off()
###\note Difference in Hunt Power is greater when comparing against no strike Fails
## Motivation
densDLMotivation <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$Success+datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$Fails)
densNLMotivation <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$Success+datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$Fails)
densLLMotivation <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$Success+datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$Fails)
plot(densLLMotivation,col=colourH[2],main="Hunt score density",
     xlab="Hunting score  " ,type="l",lwd=2,ylim=c(0,0.1))
lines(densNLMotivation,col=colourH[3],lwd=2)
lines(densDLMotivation,col=colourH[1],lwd=2)
legend("topright",legend=paste(c("DL #","LL #","NL #"),c(densDLMotivation$n,densLLMotivation$n,densNLMotivation$n) ),fill = colourH)

##How Many Fish From Each Group Have A Score Higher Than :
tblSuccessDist <- table(datFishSuccessRate[datFishSuccessRate$Success > 0, ]$groupID,datFishSuccessRate[datFishSuccessRate$Success > 0, ]$Success )

## Examine Table of Hunt Power
table(datFishSuccessRateActive[datFishSuccessRateActive$vEfficiencyRatio > 0.0,]$groupID,
      round(datFishSuccessRateActive[datFishSuccessRateActive$vEfficiencyRatio > 0.0,]$Success*datFishSuccessRateActive[datFishSuccessRateActive$vEfficiencyRatio > 0.0,]$vEfficiencyRatio*10)/10 )


hist(round(datFishSuccessRateActive[datFishSuccessRateActive$vEfficiencyRatio > 0.0,]$Success*datFishSuccessRateActive[datFishSuccessRateActive$vEfficiencyRatio > 0.0,]$vEfficiencyRatio*10)/10)