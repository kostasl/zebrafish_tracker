### Auxilliry functions involving Labelling data on hunt events ###  #
###
##   Can COMPARE Two Labelling Sets Using : ########
huntComp <- compareLabelledEvents(datHuntLabelledEventsSB,datHuntLabelledEventsSB2)
huntComp$huntScore <- convertToScoreLabel(huntComp$huntScore) ##Convert to Labels
huntComp$huntScoreB <- convertToScoreLabel(huntComp$huntScoreB)
datLabelClash <- huntComp[huntComp$huntScore != huntComp$huntScoreB & !is.na(huntComp$huntScoreB) & !is.na(huntComp$huntScore) 
                          & huntComp$huntScore!="UnLabelled",] ##Bring Out The labelled Mismatches##Compare:
tblLabelCompare <- table(huntComp$huntScore, huntComp$huntScoreB) ##COlumns is HuntEventB scores
write.csv(tblLabelCompare,file=paste(strDatDir,"/LabelledSet/","tblCompareLabellingSummary.csv",sep="") )
write.csv(datLabelClash,file=paste(strDatDir,"/LabelledSet/","tblLabelClash.csv",sep="") )
##########################

##################################### IMPORT / UPDATE LABELS ####
#### IMport New Labels From Set ###
strDataFileName <-paste("setn-12","-HuntEvents-SB-ALL",sep="") ##To Which To Save After Loading
datHuntLabelledEventsTarget <-readRDS(file=paste(strDatDir,"/LabelledSet/",strDataFileName,".rds",sep="" )) ##Save With Dataset Idx Identifier
strDataFileName <-paste("setn-12-HuntEvents-SB-ALL-09-07-18",sep="") ##To Which To Save After Loading
datHuntLabelledEventsSource <-readRDS(file=paste(strDatDir,"/LabelledSet/",strDataFileName,".rds",sep="" )) ##Save With Dataset Idx Identifier

## Check For Diferencess ##
#comp2 <- compareLabelledEvents(datHuntLabelledEventsTarget,datHuntLabelledEventsSource)
#datLabelClash <- comp2[comp2$huntScore != comp2$huntScoreB & !is.na(comp2$huntScoreB) & !is.na(comp2$huntScore) 
#                       & comp2$huntScore!=0,] ##Bring Out The labelled Mismatches##Compare:
datHuntLabelledEvents <- digestHuntLabels(datHuntLabelledEventsTarget,datHuntLabelledEventsSource)
##Save Onto New File 
strOutDataFileName <- "setn12-HuntEvents-SB-Updated"
saveRDS(datHuntLabelledEvents,file=paste(strDatDir,"/LabelledSet/",strOutDataFileName,".rds",sep="" ))
####### #############################
#######



