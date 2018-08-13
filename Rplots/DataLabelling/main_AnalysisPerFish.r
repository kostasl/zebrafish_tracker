### Success Analysis Per Fish - 
## Using Labelled Data from SB - with sparse corrections by KL where wrong labels have been spotted
### KOstasl 13-08-18



strProcDataFileName <-paste("setn12-HuntEvents-SB-Updated",sep="") ## Latest Updated HuntEvent Labelled data that integrates new COming Labels
strProcDataFileName <- "setn12-HuntEventsFixExpID-SB-Updated"
#strProcDataFileName <-paste("setn-12-HuntEvents-SB-ALL_19-07-18",sep="") ## Latest Updated HuntEvent Labelled data
message(paste(" Loading Hunt Event List to Analyse... "))
#load(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
datHuntLabelledEventsSB <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))
tblResSB <- table(convertToScoreLabel(datHuntLabelledEventsSB$huntScore),datHuntLabelledEventsSB$groupID)


tblFishScores <- table(datHuntLabelledEventsSB$expID, convertToScoreLabel(datHuntLabelledEventsSB$huntScore) )

datFishSuccessRate <- data.frame( cbind("Success" = tblFishScores[,3]+tblFishScores[,12],"Fails"= tblFishScores[,4]+tblFishScores[,10]+tblFishScores[,11],"groupID"=NA) ) #
for (e in row.names(datFishSuccessRate) )
  datFishSuccessRate[e,"groupID"] <- unique( datHuntLabelledEventsSB[datHuntLabelledEventsSB$expID == e,"groupID"] )

#datFishSuccessRate[,"Success"] <- as.numeric(datFishSuccessRate[,"Success"])
#
vScoreIdx <- ((datFishSuccessRate[,"Success"]-datFishSuccessRate[,"Fails"])/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails"]))
vScoreIdx[is.nan(vScoreIdx) ] <- 0

##How Many Fish From Each Group Have A Score Higher Than :
tblSuccessDist <- table(datFishSuccessRate[datFishSuccessRate$Success > 0, ]$groupID,datFishSuccessRate[datFishSuccessRate$Success > 0, ]$Success )
