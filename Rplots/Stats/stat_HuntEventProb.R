###Estimate/ Compare Hunt Event Prob. vs prey density for each group

vHuntEventLabels <- c("UnLabelled","NA","Success","Fail","No_Target","Not_HuntMode/Delete","Escape","Out_Of_Range","Duplicate/Overlapping","Fail-No Strike","Fail-With Strike",
                      "Success-SpitBackOut",
                      "Debri-Triggered")
huntLabels <- factor(x=5,levels=c(0,1,2,3,4,5,6,7,8,9,10,11,12),labels=vHuntEventLabels )##Set To NoTHuntMode


strCondTags = "DL"
strCondTagsE = "DE"
message(paste("#### ProcessGroup ",strCondTags," ###############"))
strDataFileName <- paste("./Stats/data/setn-12-HuntEvents",strCondTags,sep="-") ##To Which To Save After Loading
message(paste(" Loading Hunt Events: ",strDataFileName))
load(file=paste(strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier

datHuntEventFiltL <- datHuntEvent[datHuntEvent$huntScore != which(levels(huntLabels) == "NA") &
                                   datHuntEvent$huntScore != which(levels(huntLabels) == "Not_HuntMode/Delete") &
                                   datHuntEvent$huntScore != which(levels(huntLabels) == "Out_Of_Range") & 
                                   datHuntEvent$huntScore != which(levels(huntLabels) == "Duplicate/Overlapping") |
                                   datHuntEvent$eventID   == 0 , ] ##Keep THose EventID 0 so as to identify All experiments - even those with no events

strDataFileName <- paste("./Stats/data/setn-12-HuntEvents",strCondTagsE,sep="-") ##To Which To Save After Loading
message(paste(" Loading Hunt Events: ",strDataFileName))
load(file=paste(strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier

datHuntEventFiltE <- datHuntEvent[datHuntEvent$huntScore != which(levels(huntLabels) == "NA") &
                                    datHuntEvent$huntScore != which(levels(huntLabels) == "Not_HuntMode/Delete") &
                                    datHuntEvent$huntScore != which(levels(huntLabels) == "Out_Of_Range") & 
                                    datHuntEvent$huntScore != which(levels(huntLabels) == "Duplicate/Overlapping") |
                                    datHuntEvent$eventID   == 0 , ] ##Keep THose EventID 0 so as to identify All experiments - even those with no events
datHuntEventFiltComb <- rbind(datHuntEventFiltL,datHuntEventFiltE)



X11()
hist(datHuntEventFiltComb$PreyCount,breaks=seq(-1,80,length=60),main=paste(strCondTags,strCondTagsE,sep="+") )