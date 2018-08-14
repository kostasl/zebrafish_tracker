### Success Analysis Per Fish - 
## Using Labelled Data from SB - with sparse corrections by KL where wrong labels have been spotted
### KOstasl 13-08-18



strProcDataFileName <-paste("setn12-HuntEvents-SB-Updated",sep="") ## Latest Updated HuntEvent Labelled data that integrates new COming Labels
strProcDataFileName <- "setn14-HuntEventsFixExpID-SB-Updated"
#strProcDataFileName <-paste("setn-12-HuntEvents-SB-ALL_19-07-18",sep="") ## Latest Updated HuntEvent Labelled data
message(paste(" Loading Hunt Event List to Analyse... "))
#load(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
datHuntLabelledEventsSB <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))

datFishSuccessRate <- getHuntSuccessPerFish(getHuntSuccessPerFish)
##Return Data Frame Of Counts For Success /Failures According to Labelled Data - Indicating Fish Group membership
getHuntSuccessPerFish <- function(datHuntLabelledEvents)
{
  tblResSB <- table(convertToScoreLabel(datHuntLabelledEventsSB$huntScore),datHuntLabelledEventsSB$groupID)
  tblFishScores <- table(datHuntLabelledEventsSB$expID, convertToScoreLabel(datHuntLabelledEventsSB$huntScore) )
  tblFishScoresLabelled<- tblFishScores[tblFishScores[,1] < 1, ] ##Pick Only THose ExpId (Fish) Whose Labelling Has (almost!) Finished
  datFishSuccessRate <- data.frame( cbind("Success" = tblFishScoresLabelled[,3]+tblFishScoresLabelled[,12],
                                        "Fails"= tblFishScoresLabelled[,4]+tblFishScoresLabelled[,10]+tblFishScoresLabelled[,11],
                                        "groupID"=NA) ) #
  for (e in row.names(tblFishScoresLabelled) )
    datFishSuccessRate[e,"groupID"] <- unique( datHuntLabelledEventsSB[datHuntLabelledEventsSB$expID == e,"groupID"] )
  
return (datFishSuccessRate)
  
}

###PLOT
layout()
plot(datFishSuccessRate[datFishSuccessRate$groupID == "LL",]$Fails,datFishSuccessRate[datFishSuccessRate$groupID == "LL",]$Success,
     pch=16,
     xlab="#Capture Failures",ylab="#Capture Success")
points(datFishSuccessRate[datFishSuccessRate$groupID == "NL",]$Fails,datFishSuccessRate[datFishSuccessRate$groupID == "NL",]$Success,pch=2)
points(datFishSuccessRate[datFishSuccessRate$groupID == "DL",]$Fails,datFishSuccessRate[datFishSuccessRate$groupID == "DL",]$Success,pch=4)
legend(2,15,legend=c("LL","NL","DL"), pch=c(16,2,4))

##Histogram of Succe
layout(matrix(c(1,2,3), 3, 1 ,byrow=TRUE))
hist(datFishSuccessRate[datFishSuccessRate$groupID == "LL",]$Success,breaks = seq(0,45,3),ylim=c(0,50),main="LL",col="#0000FFAA")
hist(datFishSuccessRate[datFishSuccessRate$groupID == "LL",]$Fails,breaks = seq(0,45,3),ylim=c(0,50),main="LL",col="#FF0000AA",add=T)

hist(datFishSuccessRate[datFishSuccessRate$groupID == "NL",]$Success,breaks = seq(0,45,3),ylim=c(0,50),main="NL",col="#0000FFAA")
hist(datFishSuccessRate[datFishSuccessRate$groupID == "NL",]$Fails,breaks = seq(0,45,3),ylim=c(0,50),main="NL",col="#FF0000AA",add=T)

hist(datFishSuccessRate[datFishSuccessRate$groupID == "DL",]$Success,breaks = seq(0,45,3),ylim=c(0,50),main="DL",col="#0000FFAA")
hist(datFishSuccessRate[datFishSuccessRate$groupID == "DL",]$Fails,breaks = seq(0,45,3),ylim=c(0,50),main="DL",col="#FF0000AA",add=T)
#datFishSuccessRate[,"Success"] <- as.numeric(datFishSuccessRate[,"Success"])
#
vScoreIdx <- ((datFishSuccessRate[,"Success"]-datFishSuccessRate[,"Fails"])/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails"]))
vScoreIdx[is.nan(vScoreIdx) ] <- 0

##How Many Fish From Each Group Have A Score Higher Than :
tblSuccessDist <- table(datFishSuccessRate[datFishSuccessRate$Success > 0, ]$groupID,datFishSuccessRate[datFishSuccessRate$Success > 0, ]$Success )
