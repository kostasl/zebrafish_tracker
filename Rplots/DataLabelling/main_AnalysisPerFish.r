### Success Analysis Per Fish - 
## Using Labelled Data from SB - with sparse corrections by KL where wrong labels have been spotted
### KOstasl 13-08-18


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

vScoreIdx <- ((datFishSuccessRate[,"Success"]-datFishSuccessRate[,"Fails"])/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails"]))
vScoreIdx[is.nan(vScoreIdx) ] <- 0

##How Many Fish From Each Group Have A Score Higher Than :
tblSuccessDist <- table(datFishSuccessRate[datFishSuccessRate$Success > 0, ]$groupID,datFishSuccessRate[datFishSuccessRate$Success > 0, ]$Success )
