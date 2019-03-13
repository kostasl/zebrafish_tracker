source("DataLabelling/labelHuntEvents_lib.r")
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

strDataFileName <-paste("setn-12-HuntEvents-SB-ALL_13-08-18",sep="") ##From Which To Obtain New Labelled HuntEvent records

datHuntLabelledEventsSource <-readRDS(file=paste(strDatDir,"/LabelledSet/",strDataFileName,".rds",sep="" )) ##Save With Dataset Idx Identifier
tblResS2 <- table(convertToScoreLabel(datHuntLabelledEventsSource$huntScore),datHuntLabelledEventsSource$groupID)

## Do Merging and return New merged data Frame
datHuntLabelledEvents <- digestHuntLabels(datHuntLabelledEventsTarget,datHuntLabelledEventsSource)


### Once some events get labelled their FrameRange May Change and So Cannot Be Matched by digest -
#Assume Row Names Have been Retained - Merge Unlabelled Using Row Names 
unlabelledRows <- row.names(datHuntLabelledEventsTarget[convertToScoreLabel(datHuntLabelledEventsTarget$huntScore) == "UnLabelled", ])
datHuntLabelledEvents[unlabelledRows,"huntScore"] <- datHuntLabelledEventsSource[unlabelledRows,]$huntScore

NARows <- row.names(datHuntLabelledEventsTarget[convertToScoreLabel(datHuntLabelledEventsTarget$huntScore) == "NA", ])
datHuntLabelledEvents[NARows,"huntScore"] <- datHuntLabelledEventsSource[NARows,]$huntScore
##Show Summary Of Labelling
tblResM <- table(convertToScoreLabel(datHuntLabelledEvents$huntScore),datHuntLabelledEvents$groupID)

## Check For Diferences ##
## Compare Changes Betwee Previous And New Labelled Data Set##
## Check Were New File and Merged File Differ
huntCompM <- compareLabelledEvents(datHuntLabelledEvents,datHuntLabelledEventsSource)
huntCompM$huntScore <- convertToScoreLabel(huntCompM$huntScore) ##Convert to Labels
huntCompM$huntScoreB <- convertToScoreLabel(huntCompM$huntScoreB)
##List the Labelled Records where the Labels do not Agree
datDiff <- huntCompM[huntCompM$huntScore != huntCompM$huntScoreB 
          #& huntCompM$huntScore != "UnLabelled"
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

## PATCH DATASET / Found Error with Dublicate ExpID ###
##13-08-18 Found Issue With Dublicate ExpID - Patch Code to Fix Issue:
strProcDataFileName <-paste("setn12-HuntEvents-SB-Updated",sep="") ## Latest Updated HuntEvent Labelled data that integrates new COming Labels
message(paste(" Loading Hunt Event List to Patcg... "))

datHuntLabelledEventsSB <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))

datHuntLabelledEventsSB$expID <- as.numeric(as.character( datHuntLabelledEventsSB$expID) )
datHuntLabelledEventsSB[datHuntLabelledEventsSB$expID == 3871 & datHuntLabelledEventsSB$groupID == "DL",]$expID <- 3881

strOutDataFileName <- "setn12-HuntEventsFixExpID-SB-Updated"
saveRDS(datHuntLabelledEventsSB,file=paste(strDatDir,"/LabelledSet/",strOutDataFileName,".rds",sep="" ))

########## END OF PATCH ##########

### APPEND  Latest Detected Events To Labelled Dataset and Save On New File set14
lastDEvents <- datHuntLabelledEventsKL[datHuntLabelledEventsKL$dataSetID %in% c(17,18),]

datHuntLabelledEventsM <- rbind(datHuntLabelledEvents,lastDEvents[ lastDEvents$groupID %in% c("NL","LL","DL") ,])
datHuntLabelledEventsM$expID <- as.numeric(as.character( datHuntLabelledEventsM$expID) )
datHuntLabelledEventsM[datHuntLabelledEventsM$expID == 3871 & datHuntLabelledEventsM$groupID == "DL",]$expID <- 3881

##Check Table Summary of Labels
tblResM2 <- table(convertToScoreLabel(datHuntLabelledEventsM$huntScore),datHuntLabelledEventsM$groupID)

# Save On NEW File ##
strOutDataFileName <- "setn14-HuntEventsFixExpID-SB-Updated"
saveRDS(datHuntLabelledEventsM,file=paste(strDatDir,"/LabelledSet/",strOutDataFileName,".rds",sep="" ))
#############

## 24/08/18 IMPORT NEW DATASET HUNT EVENTS ONTO LABELLED SET #####
##### Here I also Import / REPLACE - EMPTY Test Condition Records after I made HuntEvent Detection Stricter  ###

#strProcDataFileName <-paste("setn-12","-HuntEvents-SB-ALL",sep="") ##To Which To Save After Loading
strProcDataFileName <- paste("setn14-HuntEventsFixExpID-SB-Updated-Merged",sep="") ##To Which To Save After Loading
strMergedDataFileName <- paste("setn15-HuntEvents-SB-Updated-Merged",sep="") ##To Which To Save After Loading
#strProcDataFileName <- paste("setn14-D5-18-HuntEvents-Merged") ##To Which To Save After Loading
message(paste(" Loading Original Labelled Hunt Event List to Process... "))
#load(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
datHuntEventAllGroupToLabel <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))

## Subset The Ones We Will Keep
datHuntEventLiveOnly <- datHuntEventAllGroupToLabel[datHuntEventAllGroupToLabel$groupID %in% c("LL","DL","NL"),]

message(paste(" Loading Hunt Event List From Empty Test Condition (LE,DE,NE) to Process... "))
strProcDataFileName <- paste("setn15-D6-5-HuntEvents-Merged",sep="") ##To Which To Save After Loading
datHuntEventAllGroupToImport <- readRDS(file=paste(strDatDir,"/",strProcDataFileName,".rds",sep="" ))
##There Will Be Fields Missing - Correct
message("Check  For Missing Fields")
names(datHuntEventAllGroupToImport) 
names(datHuntEventLiveOnly)
datHuntEventAllMerged <- rbind(datHuntEventLiveOnly,datHuntEventAllGroupToImport)
message("Write Merged Data Having Both Labelled And New Dataset Events")
strfile <- paste(strDatDir,"/",strMergedDataFileName,".rds",sep="" )
saveRDS(datHuntEventAllMerged,file=strfile )
message(strfile)
############ END OF IMPORT NEW DATASET HUNT EVENTS ONTO LABELEED SET ### 

### FIND a HUNT Event Record in the Labelled Events from the HuntAnalusis Register -
## Use it To Locate One Of the Detail Retracked HuntEvents In the Labelled Group
## You can the Use mainLabellingBlind, and give the rowID so as to replay the Video in the tracker
findLabelledEvent <- function (EventRegisterRec)
{
  
 
  #strDataFileName <- paste("setn14-D5-18-HuntEvents-Merged") ##To Which To Save After Loading
  #strDataFileName <-paste("setn14-HuntEventsFixExpID-SB-Updated-Merged",sep="") ##To Which To Save After Loading
  strDataFileName <-paste("setn15-HuntEvents-SB-Updated-Merged2",sep="") ##To Which To Save After Loading
  
  message(paste(" Loading Hunt Event List to Validate : ","/LabelledSet/",strDataFileName,".rds" ))
  
  datLabelledHuntEventAllGroups <-readRDS(file=paste(strDatDir,"/LabelledSet/",strDataFileName,".rds",sep="" )) ##Save With Dataset Idx Identifier
  
  if (is.na(EventRegisterRec$startFrame))
      warning("Missing startFrame from Event Register")
  
  
  recs<- datLabelledHuntEventAllGroups[as.character(datLabelledHuntEventAllGroups$groupID) == as.character(EventRegisterRec$groupID) &
                                  as.character(datLabelledHuntEventAllGroups$eventID) == as.character(EventRegisterRec$eventID) &
                                  as.character(datLabelledHuntEventAllGroups$expID) == as.character(EventRegisterRec$expID)
                                ,]
  
  ##If Start Frame Is there - Check For Closest Match 
  if (any(names(EventRegisterRec) == "startFrame"))
  {
    d<-(recs$startFrame - EventRegisterRec$startFrame)
    if (!is.na(d)) ##Get The BEst Match FOr Start Frame- If startFrame is not NA
    {
      recs <- recs[ which( abs(d) == min(abs(d) )), ] 
    }
    
    }
  
  return(recs)
}

### 
##################################### FIX Missing EXPID / Larvae From Hunt Event List #####################
## Notes: I discovered that when Ploting hunt event counts, the total number of larvae in each group was lower than expected - ie, 51 instead of at least
## 60 larvae per group. In the eventList used for labelling I had added for this reason eventID = 0, as a record so all larvae appear even if they produced no Hunt events.
## However, there were still some missing. The code below runs HuntEvent Redection, looks for the missing records and Merges with the existing File and Adds the Hunt Events As required.
source("HuntingEventAnalysis_lib.r")
load(paste(strDatDir,"datAllFramesFix1_Ds-5-19.RData",sep="/")) ##Raw Data Tracker Frames 
load(paste(strDatDir,"groupsrcdatListPerDataSet_Ds-5-19.RData",sep="/"))
strProcDataFileName <- "setn15-HuntEvents-SB-Updated-Merged" ##Warning Set Includes Repeated Test For some LF fish - One In Different Food Density
message(paste(" Loading Hunt Event List to Analyse... "))
#load(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
datHuntLabelledEventsSBMerged <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))

##Make A copy 
datHuntLabelledEventsSBMerged_fixed <- datHuntLabelledEventsSBMerged[!is.na(datHuntLabelledEventsSBMerged$groupID),]



### LOAD Imported Data Sets - Starting From firstDataSet
##Alternatevelly Load The Complete Set From datAllFrames_Ds-5-16-.RData ##Avoids data.frame bug rbind
firstDataSet = NROW(strDataSetDirectories)-14
lastDataSet = NROW(strDataSetDirectories)
dataSetsToProcess = seq(from=firstDataSet,to=lastDataSet)

##For Some Reason Eventlist Does not have an event0 for all missing (non BHunting) Experiments 

##For Some Reason Eventlist Does not have an event0 for all missing (non BHunting) Experiments 
## I Needed to Fix With Adding 0 Events
vExpIDS <- list() ## The Source
vExpIDT <- list() ##The Target                                                                   
vExpIDS[["NE"]] <- levels(factor(datAllFrames[datAllFrames$groupID == "NE",]$expID))
vExpIDS[["NL"]] <- levels(factor(datAllFrames[datAllFrames$groupID == "NL",]$expID))
vExpIDS[["LE"]] <- levels(factor(datAllFrames[datAllFrames$groupID == "LE",]$expID))
vExpIDS[["LL"]] <- levels(factor(datAllFrames[datAllFrames$groupID == "LL",]$expID))
vExpIDS[["DE"]] <- levels(factor(datAllFrames[datAllFrames$groupID == "DE",]$expID))
vExpIDS[["DL"]] <- levels(factor(datAllFrames[datAllFrames$groupID == "DL",]$expID))

##Obtain ExpID Found In HuntEvent List For Each Group 
for (g in names(vExpIDS))
{
  
  datHuntEvent = detectHuntEvents(datAllFrames[datAllFrames$groupID == g,],vExpIDS[[g]],dataSetsToProcess)
  if (NROW(datHuntEvent[is.na(datHuntEvent$larvaID),]$larvaID ) > 0)
    datHuntEvent[is.na(datHuntEvent$larvaID),]$larvaID <- 5 ##Fix Unaccounted for LarvaID
  datHuntEvent$filenames = "." ##Add The Filename Field 
  vExpIDT[[g]] <- levels(factor(datHuntLabelledEventsSBMerged_fixed[datHuntLabelledEventsSBMerged_fixed$groupID == g,]$expID))  ##Get LIst Of ExpID already In EventList
  datMissingRecs <- datHuntEvent[!(datHuntEvent$expID %in%  vExpIDT[[g]]),] ##Get List Of Missing Records compared to the newly detected/processed HuntEvents
  
  message(paste("*Before Merge:Hunt Events List Has  N: ",NROW(table(datHuntLabelledEventsSBMerged_fixed[datHuntLabelledEventsSBMerged_fixed$groupID == g,]$expID)),
                " Larvae in ",g) )
  datHuntLabelledEventsSBMerged_fixed <- rbind(datHuntLabelledEventsSBMerged_fixed,datMissingRecs)
  
  message(paste("*After Merge:Hunt Events List Now Accounts for N: ",NROW(table(datHuntLabelledEventsSBMerged_fixed[datHuntLabelledEventsSBMerged_fixed$groupID == g,]$expID)),
                " Larvae in ",g) )
}

##Save the Update List
strProcDataFileName <- "setn15-HuntEvents-SB-Updated-Merged2"
saveRDS(datHuntLabelledEventsSBMerged_fixed,file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))

################ The

