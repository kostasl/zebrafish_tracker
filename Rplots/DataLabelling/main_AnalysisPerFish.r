### Success Analysis Per Fish - 
## Using Labelled Data from SB - with sparse corrections by KL where wrong labels have been spotted
### KOstasl 13-08-18
source("DataLabelling/labelHuntEvents_lib.r")

#strProcDataFileName <-paste("setn-12","-HuntEvents-SB-ALL",sep="") ##To Which To Save After Loading
#strProcDataFileName <- paste("setn14-HuntEventsFixExpID-SB-Updated-Merged",sep="") ##To Which To Save After Loading
strProcDataFileName <- paste("setn15-HuntEvents-SB-Updated-Merged") ##To Which To Save After Loading
message(paste(" Loading Hunt Event List to Process... "))
#load(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
datHuntEventAllGroupToLabel <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))
##<- datHuntEvent
groupsList <- c("DL","NL","LL") ##unique(datHuntEventAllGroupToLabel$groupID)
str_FilterLabel <- "UnLabelled"

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

datFishSuccessRate <- getHuntSuccessPerFish(datHuntEventAllGroupToLabel)
tblEventsTracked <- table(datHuntEventAllGroupToLabel$expID, datHuntEventAllGroupToLabel$markTracked,useNA="always" )
remove(datFishSuccessRateMerged)

datFishSuccessRateMerged <- cbind(datFishSuccessRate,markUnTrackable=data.frame(tblEventsTracked[row.names(datFishSuccessRate),1]),
      markTracked=data.frame(tblEventsTracked[row.names(datFishSuccessRate),2]),
      notTracked=data.frame(tblEventsTracked[row.names(datFishSuccessRate),3]))

names(datFishSuccessRateMerged)[5:7] <- c("markUnTrackable","markTracked","notTracked") ##Set Field Names - notTracked : Have not been detailed retracked yet

vScoreIdx        <- ((datFishSuccessRate[,"Success"]-datFishSuccessRate[,"Fails"])/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails"]))
vEfficiencyRatio <- (datFishSuccessRate[,"Success"]/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails"]))
#vScoreIdx[is.nan(vScoreIdx) ] <- 0
datFishSuccessRateMerged <- cbind(datFishSuccessRateMerged,vScoreIdx,vEfficiencyRatio)

## Subset only the active/Larvae - ones that have hunted 
datFishSuccessRateActive <- datFishSuccessRateMerged[!is.nan(datFishSuccessRateMerged$vScoreIdx),]


## Plot Histogram of efficiency ##
ptbreaks <- seq(from=0,to=1,by=1/10)
layout(matrix(c(1,2,3), 3, 1 ,byrow=TRUE))
hist(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vEfficiencyRatio,col=colourR[1],main=paste("DL #",NROW(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]))
     ,xlab="",breaks=ptbreaks,ylim=c(0,20))
hist(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vEfficiencyRatio,col=colourR[2],main=paste("LL #",NROW(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",])),
     xlab="",breaks=ptbreaks,ylim=c(0,20))
hist(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vEfficiencyRatio,col=colourR[3],main=paste("NL #",NROW(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]))
     ,xlab="Hunt  Efficiency (Success/(Fail+Succ.) )score",breaks=ptbreaks,ylim=c(0,20))


## Plot Efficiency  Density 
densDLEffScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vEfficiencyRatio)
densNLEffScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vEfficiencyRatio)
densLLEffScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vEfficiencyRatio)

dev.off() ##Clear Old plot
plot(densLLEffScore,col=colourH[2],main="Hunt efficiency density",type="l",lwd=2,ylim=c(0,2.0))
lines(densNLEffScore,col=colourH[3],lwd=2)
lines(densDLEffScore,col=colourH[1],lwd=2,xlab="Hunting efficiency score" )


## Plot Histograp of scor/Gaine ##

ptbreaks <- seq(from=-1,to=1,by=2/10)
layout(matrix(c(1,2,3), 3, 1 ,byrow=TRUE))
hist(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vScoreIdx,col=colourR[1],main=paste("DL #",NROW(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]))
     ,xlab="",breaks=ptbreaks,ylim=c(0,20))
hist(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vScoreIdx,col=colourR[2],main=paste("LL #",NROW(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",])),
     xlab="",breaks=ptbreaks,ylim=c(0,20))
hist(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vScoreIdx,col=colourR[3],main=paste("NL #",NROW(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]))
     ,xlab="Hunt efficiency score",breaks=ptbreaks,ylim=c(0,20))

## Plot Density 
densDLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vScoreIdx)
densNLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vScoreIdx)
densLLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vScoreIdx)

dev.off() ##Clear Old plot
plot(densLLScore,col=colourH[2],main="Hunt score density",
     xlab="Hunting score  " ,type="l",lwd=2,ylim=c(0,1.0))
lines(densNLScore,col=colourH[3],lwd=2)
lines(densDLScore,col=colourH[1],lwd=2)


##How Many Fish From Each Group Have A Score Higher Than :
tblSuccessDist <- table(datFishSuccessRate[datFishSuccessRate$Success > 0, ]$groupID,datFishSuccessRate[datFishSuccessRate$Success > 0, ]$Success )

tblSuccessDist