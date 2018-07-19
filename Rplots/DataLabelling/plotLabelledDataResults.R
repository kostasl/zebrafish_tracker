###Plot Labelled Results ##
source("plotHuntStat_lib.r")

##source() Needs Plot Hunt Stat for PieChart
## We can update To The Latest DataLabelling by running auxFunctions 
## digestHuntLabels(datHuntLabelledEventsTarget,datHuntLabelledEventsSource)
############# PLOT LABELLED RESULTS ##########
strProcDataFileName <-paste("setn12-HuntEvents-SB-Updated",sep="") ## Latest Updated HuntEvent Labelled data that integrates new COming Labels
#strProcDataFileName <-paste("setn-12-HuntEvents-SB-ALL_19-07-18",sep="") ## Latest Updated HuntEvent Labelled data
message(paste(" Loading Hunt Event List to Process... "))
#load(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
datHuntLabelledEventsSB <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))
tblResSB <- table(convertToScoreLabel(datHuntLabelledEventsSB$huntScore),datHuntLabelledEventsSB$groupID)
##Load My KL Labelled File
strProcDataFileName <-paste("setn14-D5-18-HuntEvents-Merged",sep="") ##This is the KL labelled data set
datHuntLabelledEventsKL <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))
tblResKL <- table(convertToScoreLabel(datHuntLabelledEventsKL$huntScore),datHuntLabelledEventsKL$groupID)


###Label Summary  /No Target/Escapes /Success / Fail

#display.brewer.all() to see avaulable options
rfc <- colorRampPalette(rev(brewer.pal(11,'Set3')));
r <- c(rfc(30),"#FF0000");
strPlotName = paste(strPlotExportPath,"/HuntEventsLabellingSummary-Both.pdf",sep="")
pdf(strPlotName,width=16,height=8,title="Summary Of Manual Hunt Event Labelling for both SB and KL sets",onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))
#X11()
#plot.new()
layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
#layout(matrix(c(1,2,3), 2, 3, byrow = TRUE))
ScoreLabels <- c("Success","Fail","No Target","Escape")
colourH <-   c(rfc(NROW(ScoreLabels)),"#FF0000");
pieChartLabelledEvents(tblResSB,"DL") #pieChartLabelledEvents(tblResSB,"DL")
pieChartLabelledEvents(tblResSB,"NL")
pieChartLabelledEvents(tblResSB,"LL")
text(x=1.4,y=-0.8,labels = "SB",cex=1.5)  
pieChartLabelledEvents(tblResKL,"DL") #pieChartLabelledEvents(tblResSB,"DL")
pieChartLabelledEvents(tblResKL,"NL")
pieChartLabelledEvents(tblResKL,"LL")
text(x=1.4,y=0.8,labels = "KL",cex=1.5)

legend(-1.55,0.2,legend=ScoreLabels,
       fill=colourH,
       col = colourH,
       bg = "white",cex=1.0,
       merge=FALSE,horiz=FALSE)

dev.off()


###Label ##Sucess Fail Summary Only 
strPlotName = paste(strPlotExportPath,"/HuntEventsLabellingSuccessSummary-SB.pdf",sep="")
pdf(strPlotName,width=29,height=8,title="Summary Of Manual Hunt Event Labelling Conditioning On Prey Tracking Being On",onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))
#layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
layout(matrix(c(1,2,3), 2, 3, byrow = TRUE))
colourH <-  c("#B3B3B3","#66C2A5") #c(rfc(NROW(ScoreLabels)),"#FF0000");
pieChartLabelledSuccessVsFails(tblResSB,"DL") #pieChartLabelledEvents(tblResSB,"DL")
pieChartLabelledSuccessVsFails(tblResSB,"NL")
pieChartLabelledSuccessVsFails(tblResSB,"LL")
text(x=1.4,y=-0.8,labels = "SB",cex=1.5)  
# pieChartLabelledEvents(tblResKL,"DL")
# pieChartLabelledEvents(tblResKL,"NL")
# pieChartLabelledEvents(tblResKL,"LL")
# text(x=1.4,y=-0.8,labels = "KL",cex=1.5)
#plot.new()
legend(-1.40,1.0,legend=c("Fail","Success"),
       fill=colourH,
       col = colourH,
       bg = "white",cex=3.0,
       merge=FALSE,horiz=FALSE)

dev.off()

