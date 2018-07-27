
### Auxilliary functions involving Labelling data on hunt events / And Merging Of Ongoing Labelling records ###  #
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

##################################### IMPORT / UPDATE LABELS ###################
####      IMport New Labels From Set                ###
strDataFileName <-paste("setn12-HuntEvents-SB-Updated",sep="") ##On Which To Add To
datHuntLabelledEventsTarget <-readRDS(file=paste(strDatDir,"/LabelledSet/",strDataFileName,".rds",sep="" )) ##Save With Dataset Idx Identifier
tblResT <- table(convertToScoreLabel(datHuntLabelledEventsTarget$huntScore),datHuntLabelledEventsTarget$groupID)

strDataFileName <-paste("setn-12-HuntEvents-SB-ALL_19-07-18",sep="") ##From Which To Obtain New Labelled HuntEvent records
datHuntLabelledEventsSource <-readRDS(file=paste(strDatDir,"/LabelledSet/",strDataFileName,".rds",sep="" )) ##Save With Dataset Idx Identifier
tblResS2 <- table(convertToScoreLabel(datHuntLabelledEventsSource$huntScore),datHuntLabelledEventsSource$groupID)

## Do Merging and return New merged data Frame
datHuntLabelledEvents <- digestHuntLabels(datHuntLabelledEventsTarget,datHuntLabelledEventsSource)
tblResM <- table(convertToScoreLabel(datHuntLabelledEvents$huntScore),datHuntLabelledEvents$groupID)

## Check For Diferences ##
## Compare Changes Betwee Previous And New Labelled Data Set##
## Check Were New File and Merged File Differ
huntCompM <- compareLabelledEvents(datHuntLabelledEvents,datHuntLabelledEventsSource)
huntCompM$huntScore <- convertToScoreLabel(huntCompM$huntScore) ##Convert to Labels
huntCompM$huntScoreB <- convertToScoreLabel(huntCompM$huntScoreB)
##List the Labelled Records where the Labels do not Agree
datDiff <- huntCompM[huntCompM$huntScore != huntCompM$huntScoreB 
          & huntCompM$huntScore != "UnLabelled"
          & !is.na(huntCompM$huntScore)
          & !is.na(huntCompM$huntScoreB),]

# Save to Output Table
tblLabelCompare <- table(huntCompM$huntScore, huntCompM$huntScoreB) ##COlumns is HuntEventB scores
tblLabelClash <-  table(datDiff$huntScore, datDiff$huntScoreB) ##Show Which Labels Have been Converted To Which
write.csv(tblLabelCompare,file=paste(strDatDir,"/LabelledSet/","tblCompareLabellingSummary.csv",sep="") )
write.csv(datDiff,file=paste(strDatDir,"/LabelledSet/","tblLabelClash.csv",sep="") )

##Save Onto New File 
strOutDataFileName <- "setn12-HuntEvents-SB-Updated"
saveRDS(datHuntLabelledEvents,file=paste(strDatDir,"/LabelledSet/",strOutDataFileName,".rds",sep="" ))
####### #############################
#######

###Use it To Locate One Of the Detail Retracked HuntEvents In the Labelled Group
findLabelledEvent <- function (EventRegisterRec)
{
  strDataFileName <-paste("setn12-HuntEvents-SB-Updated",sep="") ##On Which To Add To
  datLabelledHuntEventAllGroups <-readRDS(file=paste(strDatDir,"/LabelledSet/",strDataFileName,".rds",sep="" )) ##Save With Dataset Idx Identifier
  
  
  
  recs<- datLabelledHuntEventAllGroups[as.character(datLabelledHuntEventAllGroups$groupID) == as.character(EventRegisterRec$groupID) &
                                  as.character(datLabelledHuntEventAllGroups$eventID) == as.character(EventRegisterRec$eventID) &
                                  as.character(datLabelledHuntEventAllGroups$expID) == as.character(EventRegisterRec$expID)
                                ,]
  ##If Start Frame Is there - Check For Closest Match 
  if (any(names(EventRegisterRec) == "startFrame"))
  {
    d<-(recs$startFrame - EventRegisterRec$startFrame)
    recs <- recs[ which( abs(d) == min(abs(d) )), ] 
    
    }
  
  return(recs)
}


