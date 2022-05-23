##plot Hunting Event Statistics

##Requires Violin Plot vioplot * run : install.packages("vioplot")
#library(vioplot)
library(RColorBrewer);

source("plotHuntStat_lib.r")

### Load Merged Hunt Event Files datAllHuntEvent ###

message(paste("Load datAllHuntEvent :",strDataFileName ) )
load(file=paste(strDatDir,"/LabelledSet/",strDataFileName,sep="" )) ##Save With Dataset Idx Identifier
datAllHuntEvent <- datAllHuntEvents#datHuntEventAllGroupToLabel

## Load Hunt Stats from Merged Data
strDataFileName <- paste("setn",NROW(unique(datAllHuntEvent$dataSetID)),"-D",min(as.numeric(datAllHuntEvent$dataSetID) ),
                         "-",max(as.numeric(datAllHuntEvent$dataSetID) ),"-datHuntStat",sep="" )
message(paste("Load datHuntStat :",strDataFileName ) )
load(file=paste(strDataExportDir,"/",strDataFileName,".RData",sep="" ))
#save(datHuntStat,file=paste(strDataExportDir,"/",strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier



##par(bg="black")

colourH <- c(rgb(0.01,0.7,0.01,0.2),rgb(0.9,0.01,0.01,0.2),rgb(0.01,0.01,0.9,0.2))

################-HUNTING STAT PLOTS -###################
## Bar Plot Mean Hunting Events Per Animal #
## Common Subtitle info            ########
sampleSize = sum(unlist(datHuntStat[,"nLarva"],use.names = FALSE))
totalFrames = sum(unlist(datHuntStat[,"totalFrames"],use.names = FALSE))

FPS = 420;
strsub = paste("#n=", sampleSize, " #F:",totalFrames,
               "(",format(totalFrames/FPS/60,digits =3),"min)",
               " F_H/F:",format(sum(unlist(datHuntStat[,"totalHuntFrames"] ,use.names = FALSE) )/totalFrames,digits =3) ,
               " #Hunts:",sum(unlist(datHuntStat[,"groupHuntEvents"] ,use.names = FALSE) ),
               collapse=NULL)



strPlotName = paste("plots/HuntStat_N",sampleSize,".pdf",sep="")

## PDF OUTPUT ##
bonefile = FALSE

if (bonefile)
{
  pdf(strPlotName,paper = "a4",title="zfish Hunt Stat",width=16.27,height=22,onefile=TRUE) #col=(as.integer(filtereddatAllFrames$expID)) width=8,height=8
  par(mfrow=c(2,2)) ##MultiPlot Page
}



#huntFrames = sum(unlist(datMotionStat[,"huntframes"],use.names = FALSE))
##### Done Subtitle ##


##PLot Distribution Of Hunt Events Among Their Sampled Prey Counts - Not the initial Prey Counts##
strPlotName = "plots/HuntEventsVsPreyCount_Hist.pdf"

pdf(strPlotName,width=8,height=10,title="Hunt Vs the Prey Count they Occured under, for each Condition") #col=(as.integer(filtereddatAllFrames$expID))
#X11()
plotHuntEventPreyCountHist(strCondTags, dataSetsToProcess)
dev.off()


strPlotName = "plots/HuntEventsPerLarvaVsPreyCount_Hist.pdf"
pdf(strPlotName,width=8,height=10,title="Hunt/Larva Vs the Prey Count they Occured under, for each Condition") #col=(as.integer(filtereddatAllFrames$expID))
plotMeanHuntEventPerLarvaVsPreyCountHist(datAllHuntEvent)
dev.off()

### Interval Per Prey Count - Examine if there is a lag between hunt episodes that goes up
## in LL, so as to explain the diminished hunting rates with Increasing prey numbers 
strPlotName = "plots/HuntEventIntervalsPerLarvaVsPreyCount_Hist.pdf"
pdf(strPlotName,width=8,height=10,title="Time between Hunt Episodes of the Same event Vs the Prey Count they Occured under, for each Condition") #col=(as.integer(filtereddatAllFrames$expID))
plotMeanHuntIntervalPerLarvaVsPreyCountHist(datAllHuntEvent)
dev.off()

## Episode Duration Summary Box Plots
strPlotName = paste(strPlotExportPath,"/EpisodeDurationOnLabelledSet.pdf",sep="")
pdf(strPlotName,width=9,height=10,title="Episode Duration (T) compared across labels and Conditions (Labelled set)") #col=(as.integer(filtereddatAllFrames$expID))
boxPlotHuntEpisodeDuration(datAllHuntEvent)
dev.off()

########### MEAN and Distribution of Prey Count At Start of Hunt EVENTS ##### 
strPlotName = "plots/meanInitialPreyCountPerCond.pdf"
vDat <- datHuntStat[,"vHInitialPreyCount"]
vDatSetID <- datHuntStat[,"vDataSetID"]
vIDTable <- datHuntStat[,"vIDLookupTable"] ##vIDTable$DL <- vIDTable$DL[vIDTable$DL$expID!=3830,]


datmean <- unlist(datHuntStat[,"initPreyCount"],use.names = FALSE)
datse <- unlist(datHuntStat[,"initsePreyCount"],use.names = FALSE)
strtitle <- "Mean Initial Prey Count per Experiment"


yl <- c(0,sampleSize/5) ##Remove The 3830 from DL Hack
xl <- c(0,max(vDat$LL,vDat$NL,vDat$DL[!is.na(vDat$DL)])+10)


if (!bonefile)
  pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))

ylim <- max(unlist(vDat,use.names=FALSE),na.rm=TRUE)
xbarcenters <- barplotPerCondition(vDat,datmean,datse,strtitle,strsub,strPlotName,ylim)
plotConnectedPointsPairs(vIDTable,vDat,strCondTags,xbarcenters)

if (!bonefile)
dev.off()

##Do Initial Prey  Histogram for All Groups
strPlotName = "plots/meanInitialPreyCountPerCond_Hist.pdf"
vDat <- datHuntStat[,"vHInitialPreyCount"]
pdf(strPlotName,width=8,height=8,title="Initial Prey Count Per Experiment Histograms") #col=(as.integer(filtereddatAllFrames$expID))

par(mfrow=c(3,2)) ##MultiPlot Page
layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE))

yl <- c(0,ceiling(max(datHuntStat[,"vHLarvaEventCount"]$LL,datHuntStat[,"vHLarvaEventCount"]$DL,datHuntStat[,"vHLarvaEventCount"]$NL)/2)+5 )
xl <- c(0,70)
binSize <- 3
br <- seq(0,75,binSize)

hist(vDat$LE,col=colourH[1],ylim=yl,xlim=xl,breaks=br,
     main="LE",
     xlab = "Number of Prey at start of experiment",
     ylab = "# Experiments")
vDat$LL[vDat$LL > xl[2]] <- xl[2] ##Fix To Max
hist(vDat$LL,col=colourH[1],ylim=yl,xlim=xl,breaks=br,
     main="LL",
     xlab = "Number of Prey at start of experiment",
     ylab = "# Experiments")

hist(vDat$NE,col=colourH[2],ylim=yl,xlim=xl,breaks=br,
     main="NE",
     xlab = "Number of Prey at start of experiment",
     ylab = "# Experiments")

vDat$NL[vDat$NL > xl[2]] <- xl[2] ##Fix To Max
hist(vDat$NL,col=colourH[2],ylim=yl,xlim=xl,breaks=br,
     main="NL",
     xlab = "Number of Prey at start of experiment",
     ylab = "# Experiments")

hist(vDat$DE,col=colourH[3],ylim=yl,xlim=xl,breaks=br,
     main="DE",
     xlab = "Number of Prey at start of experiment",
     ylab = "# Experiments")
hist(vDat$DL,col=colourH[3],ylim=yl,xlim=xl,breaks=br,
     main="DL",
     xlab = "Number of Prey at start of experiment",
     ylab = "# Experiments")

dev.off()
###### # # # # # ## ## # ##  # # ## #  ## # # # 


#####  Scattter OF Hunt Events Vs Initial  Prey Counts ### ##
strPlotName = "plots/HuntEventsVsInitPreyCount_scatter.pdf"
pdf(strPlotName,width=8,height=8,title="Number of Prey In Live Test Conditions") 

yl <- c(0,ceiling(max(datHuntStat[,"vHLarvaEventCount"]$LL,datHuntStat[,"vHLarvaEventCount"]$DL,datHuntStat[,"vHLarvaEventCount"]$NL)/10)*10 )
xl <- c(0,60)
vB <- seq(0,xl[2],xl[2]/10)

hist(vDat$LL,col=colourH[1],ylim=yl,xlim=xl,breaks=vB,
     main="Number of Prey Vs Hunting Events",
     xlab = "Number of Prey at start of experiment",
     ylab = "Number of Hunt Events")
hist(vDat$LE,col=colourH[1],ylim=yl,xlim=xl,breaks=vB,add=T)

vDat$NL[vDat$NL > xl[2]] <- xl[2] ##Fix To Max
#vDat$NL[vDat$NL <= max(vDat$NE)] <- xl[2] ##Fix To Max
hist(vDat$NL,col=colourH[2],ylim=yl,xlim=xl,breaks=vB,add=T)
hist(vDat$NE,col=colourH[2],ylim=yl,xlim=xl,breaks=vB,add=T)
hist(vDat$DL,col=colourH[3],ylim=yl,xlim=xl,breaks=vB,add=T)
hist(vDat$DE,col=colourH[3],ylim=yl,xlim=xl,breaks=vB,add=T)
box()
par(new=TRUE)
##Add Scatter Of Hunt Events
plot(datHuntStat[,"vHInitialPreyCount"]$LL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL),
     pch=4,col=1,
     xlab="",ylab="",ylim=yl,xlim=xl,main="")
points(datHuntStat[,"vHInitialPreyCount"]$LE,as.numeric(datHuntStat[,"vHLarvaEventCount"]$LE),
     pch=4,col=1,
     xlab="",ylab="",ylim=yl,xlim=xl,main="")


points(datHuntStat[,"vHInitialPreyCount"]$NL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$NL),
       pch=2 ,col=2,type ="p",
       xlab="",ylab="",ylim=yl,xlim=xl)
points(datHuntStat[,"vHInitialPreyCount"]$NE,as.numeric(datHuntStat[,"vHLarvaEventCount"]$NE),
       pch=2 ,col=2,type ="p",
       xlab="",ylab="",ylim=yl,xlim=xl)


points(datHuntStat[,"vHInitialPreyCount"]$DL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$DL),
       pch=19 ,col=6,type ="p",
       xlab="",ylab="",ylim=yl,xlim=xl)
points(datHuntStat[,"vHInitialPreyCount"]$DE,as.numeric(datHuntStat[,"vHLarvaEventCount"]$DE),
       pch=19 ,col=6,type ="p",
       xlab="",ylab="",ylim=yl,xlim=xl)


legend(50,55,legend=c("Live Fed","Not Fed","Dry Fed"),
       fill=colourH,
       col = c(1, 2,6),pch = c(4,2,19),
       bg = "gray90",lty = c(2, -1, 1),
       merge=TRUE)

dev.off()
##########  # ## # # # # 

##### Number Of Events Recorded/Seen Vs Prey ###


#####  Scattter OF Hunt Events Vs Initial  Prey Counts ### ##
strPlotName = "plots/HuntEventsCountVsInitPreyCount_scatter.pdf"
pdf(strPlotName,width=8,height=8,title="Number of Prey In Live Test Conditions") 


yl <- c(0,ceiling(max(datHuntStat[,"vHLarvaEventCount"]$LL,datHuntStat[,"vHLarvaEventCount"]$DL,datHuntStat[,"vHLarvaEventCount"]$NL)/10)*10 )
xl <- c(0,60)


hist(vDat$LL,col=colourH[1],ylim=yl,xlim=xl,breaks=10,
     main="Number of Prey Vs Hunting Events",
     xlab = "Number of Prey at start of experiment",
     ylab = "Number of  Recorded Events")
hist(vDat$NL,col=colourH[2],ylim=yl,xlim=xl,breaks=10,add=T)
hist(vDat$DL,col=colourH[3],ylim=yl,xlim=xl,breaks=10,add=T)
box()
par(new=TRUE)
##Add Scatter Of Hunt Events
plot(datHuntStat[,"vHInitialPreyCount"]$LL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL),
     pch=4,col=1,
     xlab="",ylab="",ylim=yl,xlim=xl,main="")

points(datHuntStat[,"vHInitialPreyCount"]$LE,as.numeric(datHuntStat[,"vHLarvaEventCount"]$LE),
     pch=4,col=1,
     xlab="",ylab="",ylim=yl,xlim=xl,main="")


points(datHuntStat[,"vHInitialPreyCount"]$NL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$NL),
       pch=2 ,col=2,type ="p",
       xlab="",ylab="",ylim=yl,xlim=xl)

points(datHuntStat[,"vHInitialPreyCount"]$NE,as.numeric(datHuntStat[,"vHLarvaEventCount"]$NE),
       pch=2 ,col=2,type ="p",
       xlab="",ylab="",ylim=yl,xlim=xl)


points(datHuntStat[,"vHInitialPreyCount"]$DL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$DL),
       pch=19 ,col=6,type ="p",
       xlab="",ylab="",ylim=yl,xlim=xl)


points(datHuntStat[,"vHInitialPreyCount"]$DE,as.numeric(datHuntStat[,"vHLarvaEventCount"]$DE),
       pch=19 ,col=6,type ="p",
       xlab="",ylab="",ylim=yl,xlim=xl)

legend(50,55,legend=c("LL","NL","DL"),
       fill=colourH,
       col = c(1, 2,6),pch = c(4,2,19),
       bg = "gray90",lty = c(2, -1, 1),
       merge=TRUE)

dev.off()
##########  # ## # # # # 





########### MEAN and Distribution of FINAL PREY COUNT - On Last Hunt EVENTS ##### 
strPlotName = "plots/meanFinalPreyCountPerCond_Hist.pdf"
vDat <- datHuntStat[,"vHPreyReductionPerLarva"]
vDatSetID <- datHuntStat[,"vDataSetID"]
vIDTable <- datHuntStat[,"vIDLookupTable"] ##vIDTable$DL <- vIDTable$DL[vIDTable$DL$expID!=3830,]


datmean <- unlist(datHuntStat[,"meanPreyReductionPerLarva"],use.names = FALSE)
datse <- unlist(datHuntStat[,"sePreyReductionPerLarva"],use.names = FALSE)
strtitle <- "Mean Final Prey Count per Experiment"


pdf(strPlotName,width=8,height=8,title="Number of Prey at End of Experiment Histogram") 
yl <- c(0,ceiling(max(datHuntStat[,"vHLarvaEventCount"]$LL,datHuntStat[,"vHLarvaEventCount"]$DL,datHuntStat[,"vHLarvaEventCount"]$NL)/2) )
xl <- c(0,70)
binSize <- 3
br <- seq(0,75,binSize)

par(mfrow=c(3,2)) ##MultiPlot Page
layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE))

hist(vDat$LE,col=colourH[1],ylim=yl,xlim=xl,breaks=br,
     main="LE",
     xlab = "Number of Prey at End of experiment",
     ylab = "# Experiments")
hist(vDat$LL,col=colourH[1],ylim=yl,xlim=xl,breaks=br,
     main="LL",
     xlab = "Number of Prey at End of experiment",
     ylab = "# Experiments")

hist(vDat$NE,col=colourH[2],ylim=yl,xlim=xl,breaks=br,
     main="NE",
     xlab = "Number of Prey at End of experiment",
     ylab = "# Experiments")
hist(vDat$NL,col=colourH[2],ylim=yl,xlim=xl,breaks=br,
     main="NL",
     xlab = "Number of Prey at End of experiment",
     ylab = "# Experiments")

hist(vDat$DE,col=colourH[3],ylim=yl,xlim=xl,breaks=br,
     main="DE",
     xlab = "Number of Prey at End of experiment",
     ylab = "# Experiments")
hist(vDat$DL,col=colourH[3],ylim=yl,xlim=xl,breaks=br,
     main="DL",
     xlab = "Number of Prey at End of experiment",
     ylab = "# Experiments")

dev.off()

### Change In Prey Count Plot ###
strPlotName = "plots/meanPreyCountChangePerCond_Hist.pdf"
vDatF <- datHuntStat[,"vHPreyReductionPerLarva"]
vDatI <- datHuntStat[,"vHInitialPreyCount"]
vDatSetID <- datHuntStat[,"vDataSetID"]
vIDTable <- datHuntStat[,"vIDLookupTable"] ##vIDTable$DL <- vIDTable$DL[vIDTable$DL$expID!=3830,]


pdf(strPlotName,width=8,height=8,title="Number of Prey at End of Experiment Histogram") 
yl <- c(0,ceiling(max(datHuntStat[,"vHLarvaEventCount"]$LL,datHuntStat[,"vHLarvaEventCount"]$DL,datHuntStat[,"vHLarvaEventCount"]$NL)/1.8) )
xl <- c(-30,30)
binSize <- 3
br <- seq(-40,40,binSize)
#br <- br[br != 0]

par(mfrow=c(3,2)) ##MultiPlot Page
layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE))

blnc <- vDatI$LE-vDatF$LE
hist(blnc,col=colourH[1],ylim=yl,xlim=xl,breaks=br,
     main=paste("LE ",NROW(blnc[blnc<= -binSize])," vs ",NROW(blnc[blnc>= +binSize]),sep=""),
     xlab = "Reduction in Prey by End of experiment",
     ylab = "# Experiments")

blnc <- round(vDatI$LL-vDatF$LL)
hist(blnc,col=colourH[1],ylim=yl,xlim=xl,breaks=br,
     main=paste("LL ",NROW(blnc[blnc<= -binSize])," vs ",NROW(blnc[blnc>= +binSize]),sep=""),
     xlab = "Reduction in Prey by End of experiment",
     ylab = "# Experiments")

blnc <- round(vDatI$NE-vDatF$NE)
hist(blnc,col=colourH[2],ylim=yl,xlim=xl,breaks=br,
     main=paste("NE ",NROW(blnc[blnc<= -binSize])," vs ",NROW(blnc[blnc>= +binSize]),sep=""),
     xlab = "Reduction in Prey by End of experiment",
     ylab = "# Experiments")

blnc <- round(vDatI$NL-vDatF$NL)
hist(blnc,col=colourH[2],ylim=yl,xlim=xl,breaks=br,
     main=paste("NL ",NROW(blnc[blnc<= -binSize])," vs ",NROW(blnc[blnc>= +binSize]),sep=""),
     xlab = "Reduction in Prey by End of experiment",
     ylab = "# Experiments")

blnc <- round(vDatI$DE-vDatF$DE)
xx <- hist(blnc,col=colourH[3],ylim=yl,xlim=xl,breaks=br,
     main=paste("DE ",NROW(blnc[blnc<= -binSize])," vs ",NROW(blnc[blnc>= +binSize]),sep=""),
     xlab = "Reduction in Prey by End of experiment",
     ylab = "# Experiments")

blnc <- round(vDatI$DL-vDatF$DL)
hist(blnc,col=colourH[3],ylim=yl,xlim=xl,breaks=br,
     main=paste("DL ",NROW(blnc[blnc<= -binSize])," vs ",NROW(blnc[blnc>= +binSize]),sep=""),
     xlab = "Reduction in Prey by End of experiment",
     ylab = "# Experiments")

dev.off()


#############
###
### Compare FED to Non Fed Group ##
# strPlotName = "plots/NumberOfHuntEventsVsPreyCount-FedGroups.pdf"
# pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
# plot(datHuntStat[,"vHInitialPreyCount"]$NL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$NL),
#      pch=4,col=1,
#      xlab="Initial Prey Count",ylab="# hunt events",ylim=yl,xlim=xl,main=" Modulation of hunt rate by stimulus")
# 
# points(datHuntStat[,"vHInitialPreyCount"]$LL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL),
#        pch=2 ,col=2,type ="p",
#        xlab="Initial Prey Count",ylab="# hunt events",ylim=yl,xlim=xl)
# 
# points(datHuntStat[,"vHInitialPreyCount"]$DL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$DL),
#        pch=19 ,col=6,type ="p",
#        xlab="Initial Prey Count",ylab="# hunt events",ylim=yl,xlim=xl)
# 
# 
# legend(0.7, 55, c("NL","LL","DL"), col = c(1, 2,6),
#        text.col = "green4", lty = c(2, -1, 1), pch = c(4,2,19),
#        merge = TRUE, bg = "gray90")  
# dev.off()  
# 
########### END oF Prey Vs Hunt Events Plot ###########











#### SLICES PLot Mean Hunt Events Within Prey Cound Ranges / Slices ########
datHuntVsPreyLL <- cbind(datHuntStat[,"vHInitialPreyCount"]$LL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL) )
datHuntVsPreyNL <- cbind(datHuntStat[,"vHInitialPreyCount"]$NL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NL) )
datHuntVsPreyDL <- cbind(datHuntStat[,"vHInitialPreyCount"]$DL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$DL) )

checkRanges <- list(range(0,20),range(10,30),range(20,40),range(30,100))


strPlotName = paste("plots/meanHuntEventsPerGroupVsInitPreyCount_RangeSliced.pdf",sep="")
pdf(strPlotName,width=8,height=8,title="Event Counts Within Initial Prey Range ") 
par(mfrow=c(2,2)) ##MultiPlot Page
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
for (rng in checkRanges)
{
  print(rng)
  preyCntRange <- rng


  datSliceLL <- datHuntVsPreyLL[datHuntVsPreyLL[,1] > preyCntRange[1] & datHuntVsPreyLL[,1] < preyCntRange[2],2 ]
  datSliceNL <- datHuntVsPreyNL[datHuntVsPreyNL[,1] > preyCntRange[1] & datHuntVsPreyNL[,1] < preyCntRange[2],2 ]
  datSliceDL <- datHuntVsPreyDL[datHuntVsPreyDL[,1] > preyCntRange[1] & datHuntVsPreyDL[,1] < preyCntRange[2],2 ]

  datmean <- c(mean(datSliceLL,na.rm = TRUE),mean(datSliceNL,na.rm = TRUE),mean(datSliceDL,na.rm = TRUE))
  datse   <-  c(sd(datSliceLL,na.rm = TRUE)/sqrt(length(datSliceLL)),sd(datSliceNL,na.rm = TRUE)/sqrt(length(datSliceNL)),sd(datSliceDL,na.rm = TRUE)/sqrt(length(datSliceDL)) )

  strLocalsub = paste("")
  
  barCenters <- barplot(height = datmean,
        names.arg = c("LL","NL","DL"),
        beside = true, las = 2,
        ylim = c(0, 35),
        cex.names = 0.75, xaxt = "n",
        main = paste("#Hunts with InitPreyCount ",preyCntRange[1],preyCntRange[2]),
        sub = strLocalsub,
        ylab = "#",
        border = "black", axes = TRUE)
  
  datlbls <- c( paste("LL \nn=",length(datSliceLL) ),
                paste("NL \nn=",length(datSliceNL) ),
                paste("DL \nn=",length(datSliceDL) ) )
  text(x = barCenters+0.2, y = 0, srt = 45,
       adj = 1.8, labels = datlbls, xpd = TRUE)
  
  
  segments(barCenters, datmean - datse * 2, barCenters,
         datmean + datse * 2, lwd = 1.5)
}
dev.off()
#################################### # # # # # # # ######################




########### MEAN Prey Count REDUCTION from Initial During Hunt EVENTS ##### 

## Reduction Per Experiment/Larva # 
# ** \Note THat Currently it reports mean final prey count  ###
strPlotName = "plots/meanPreyCountReductionPerLarva.pdf"
vDat <- datHuntStat[,"vHPreyReductionPerLarva"] #INverse Sign And Denote Reduction
#vDat$DL <- vDat$DL[!is.na(vDat$DL)]
vDatSetID <- datHuntStat[,"vDataSetID"]
vIDTable <- datHuntStat[,"vIDLookupTable"]

###This Prey Count Reduction. Shows a little *increase* for NL, DL but not for LL 
datmean <- unlist(datHuntStat[,"meanPreyReductionPerLarva"],use.names = FALSE)
datse <- unlist(datHuntStat[,"sePreyReductionPerLarva"],use.names = FALSE)
strtitle <- "Mean Prey Count Reduction from Initial For each Larva"

if (!bonefile)
  pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))

ylim <- max(unlist(vDat,use.names=FALSE),na.rm=TRUE)
xbarcenters <- barplotPerCondition(vDat,datmean,datse,strtitle,strsub,strPlotName,ylim)
plotConnectedPointsPairs(vIDTable,vDat,strCondTags,xbarcenters)

if (!bonefile)
  dev.off()


##
strPlotName = "plots/meanPreyCountReductionPerLarvaHuntEpisode.pdf"
vDat <- datHuntStat[,"vHNablaPreyCount"] #INverse Sign And Denote Reduction
#vDat$DL <- vDat$DL[!is.na(vDat$DL)]
vDatSetID <- datHuntStat[,"vDataSetID"]
vIDTable <- datHuntStat[,"vIDLookupTable"]
#vIDTable$DL <- vIDTable$DL[vIDTable$DL$expID!=3830,]

###This Prey Count Reduction. Shows a little *increase* for NL, DL but not for LL 
datmean <- unlist(datHuntStat[,"meanNablaPreyCount"],use.names = FALSE)
datse <- unlist(datHuntStat[,"seNablaPreyCount"],use.names = FALSE)
strtitle <- "Mean Prey Count Reduction from Initial at each Hunt Event Per Larva"

if (!bonefile)
  pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))


ylim <- max(unlist(vDat,use.names=FALSE),na.rm=TRUE)
xbarcenters <- barplotPerCondition(vDat,datmean,datse,strtitle,strsub,strPlotName,ylim)

plotConnectedPointsPairs(vIDTable,vDat,strCondTags,xbarcenters)

if (!bonefile)
  dev.off()


###### Histogram ##########
strPlotName = "plots/meanReductionOfPreyCountPerLarvaHuntEpisode_Hist.pdf"
pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
lDeltas <- plotPairedChangeHistogram(vIDTable,vDat,strCondTags,uLim = 50,lLim = -10)
title("Prey Counts Events", sub = strsub, cex.main = 1.2,   font.main= 1.5, col.main= "black", cex.sub = 1.0, font.sub = 2, col.sub = "black")
dev.off()
## Plot Mean Change ######

strPlotName = "plots/meanPreyCountPerLarvaHuntEpisode_PairedChange.pdf"
meanDelta <- sapply(lDeltas,mean,na.rm =TRUE)
seDelta   <-  sapply(lDeltas,sd,na.rm =TRUE)/sqrt(sapply(lDeltas,NROW))
ylim <- 2*max(unlist(meanDelta,use.names=FALSE))
pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
plot.window(xlim=c(0,100), ylim=c(-10, ylim) )
barplotPerCondition(lDeltas,meanDelta,seDelta,"Mean Change in Hunt Prey Count ",strsub,strPlotName,ylim)
dev.off()

####


############################################
########### MEAN HUNTING EVENTS ##### 
#X11()
strPlotName = "plots/meanHuntEventsPerCond.pdf"
vDat <- datHuntStat[,"vHLarvaEventCount"]
vDatSetID <- datHuntStat[,"vDataSetID"]
vIDTable <- datHuntStat[,"vIDLookupTable"]
strCondTags <- names(datHuntStat[,1])
#vDat$DL <- vDat$DL[!is.na(vDat$DL)] 
#vIDTable$DL <- vIDTable$DL[vIDTable$DL$expID!=3830,]

datmean <- unlist(datHuntStat[,"meanHuntingEventsPerLarva"],use.names = FALSE)
datse <- unlist(datHuntStat[,"seHuntingEventsPerLarva"],use.names = FALSE)
strtitle <- "Mean Hunting Events Per Condition"

if (!bonefile)
  pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))


ylim <- max(unlist(vDat,use.names=FALSE))
xbarcenters <- barplotPerCondition(datHuntStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
plotConnectedPointsPairs(vIDTable,vDat,strCondTags,xbarcenters)

if (!bonefile)
  dev.off()

###### Histogram ##########
  strPlotName = "plots/meanHuntingEventsPerLarva_ChangeHist.pdf"
  pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
  lDeltas <- plotPairedChangeHistogram(vIDTable,vDat,strCondTags,uLim = max(vDat$LL,vDat$DL,vDat$NL),lLim = -10)
  title("Hunt Events", sub = strsub, cex.main = 1.2,   font.main= 1.5, col.main= "black", cex.sub = 1.0, font.sub = 2, col.sub = "black")
  dev.off()
  
  ## Plot Mean Change ######
  
  strPlotName = "plots/meanHuntingEventsPairedChange.pdf"
  meanDelta <- sapply(lDeltas,mean)
  seDelta   <-  sapply(lDeltas,sd)/sqrt(sapply(lDeltas,NROW))
  ylim <- 2*max(unlist(meanDelta,use.names=FALSE))
  pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
  #plot.window(xlim=c(0,100), ylim=c(-10, ylim) )
  barplotPerCondition(lDeltas,meanDelta,seDelta,"Mean Change in Hunt Event Count ",strsub,strPlotName,ylim)
  dev.off()
  
#### HUNT EVENTS HISTOGRAM ###
  uLim = max(vDat$LL,vDat$DL,vDat$NL)+10
  lLim = min(vDat$LL,vDat$DL,vDat$NL)
  res = 20
  strPlotName = "plots/HuntingEventsCount_Hist.pdf"
  pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
  layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE))
  
  hist(vDat$LE,ylim = c(0,NROW(vDat$LL)/2),breaks = seq(lLim,uLim,uLim/res), main ="LE ",xlab="# Hunt Events")
  hist(vDat$LL,ylim = c(0,NROW(vDat$LL)/2),breaks = seq(lLim,uLim,uLim/res), main ="LL ",xlab="# Hunt Events")
  hist(vDat$NE,ylim = c(0,NROW(vDat$LL)/2),breaks = seq(lLim,uLim,uLim/res), main ="NE ",xlab="# Hunt Events")
  hist(vDat$NL,ylim = c(0,NROW(vDat$LL)/2),breaks = seq(lLim,uLim,uLim/res), main ="NL ",xlab="# Hunt Events")
  hist(vDat$DE,ylim = c(0,NROW(vDat$LL)/2),breaks = seq(lLim,uLim,uLim/res), main ="DE ",xlab="# Hunt Events")
  hist(vDat$DL,ylim = c(0,NROW(vDat$LL)/2),breaks = seq(lLim,uLim,uLim/res), main ="DL ",xlab="# Hunt Events")
  
  dev.off()
  
  #######################################
  
#########HUnT EVENTS VS Initial Prey Count Scatter Plot 

  yl <- c(0, max(vDat$LL,vDat$DL,vDat$NL)+10)
  xl <- c(0,50)
strPlotName = "plots/NumberOfHuntEventsVsInitPreyCount.pdf"
pdf(strPlotName,width=8,height=8,title=strtitle,onefile=TRUE) #col=(as.integer(filtereddatAllFrames$expID))
par(bg="white",fg="black")
plot(datHuntStat[,"vHInitialPreyCount"]$LE,as.numeric(datHuntStat[,"vHLarvaEventCount"]$LE),
     pch=1,col=rDataset[as.numeric( datHuntStat[,"vDataSetID"]$LE ) ],
     xlab="Initial Prey Count",ylab="# hunt events",ylim=yl,xlim=xl,
     main="LE - Modulation of hunt rate by stimulus")

plot(datHuntStat[,"vHInitialPreyCount"]$LL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL),
     pch=2 ,col=rDataset[as.numeric( datHuntStat[,"vDataSetID"]$LL ) ],type ="p",
     xlab="Initial Prey Count",ylab="# hunt events",ylim=yl,xlim=xl,
     main="LL - Modulation of hunt rate by stimulus")

plot(datHuntStat[,"vHInitialPreyCount"]$NE,as.numeric(datHuntStat[,"vHLarvaEventCount"]$NE),
     pch=3,col=rDataset[as.numeric( datHuntStat[,"vDataSetID"]$NE ) ],
     xlab="Initial Prey Count",ylab="# hunt events",ylim=yl,xlim=xl,
     main="NE - Modulation of hunt rate by stimulus")

plot(datHuntStat[,"vHInitialPreyCount"]$NL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$NL),
     pch=4,col=rDataset[as.numeric( datHuntStat[,"vDataSetID"]$NL ) ],
     xlab="Initial Prey Count",ylab="# hunt events",ylim=yl,xlim=xl,
     main="NF Live - Modulation of hunt rate by stimulus")

plot(datHuntStat[,"vHInitialPreyCount"]$DE,as.numeric(datHuntStat[,"vHLarvaEventCount"]$DE),
     pch=5,col=rDataset[as.numeric( datHuntStat[,"vDataSetID"]$DE ) ],
     xlab="Initial Prey Count",ylab="# hunt events",ylim=yl,xlim=xl,
     main="DE  - Modulation of hunt rate by stimulus")  

plot(datHuntStat[,"vHInitialPreyCount"]$DL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$DL),
     pch=6,col=rDataset[as.numeric( datHuntStat[,"vDataSetID"]$DL ) ],
     xlab="Initial Prey Count",ylab="# hunt events",ylim=yl,xlim=xl,
     main="DL - Modulation of hunt rate by stimulus")  
dev.off()


######## EPISODE DURATION ############


strPlotName = "plots/meanEpisodeDurationOfGroup.pdf"
datmean <- unlist(datHuntStat[,"meanEpisodeDuration"],use.names = FALSE) #Of the Group
#vDat    <- datHuntStat[,"vmeanHLarvaDuration"] # Mean Episode Duration of Each LArva
vDat    <- datHuntStat[,"vmeanHEpisodeDurationPerLarva"]

vDatSetID <- datHuntStat[,"vDataSetID"]
datse   <- unlist(datHuntStat[,"seEpisodeDuration"],use.names = FALSE)
strtitle <- "Mean Duration of each Hunting Episode"

if (!bonefile)
  pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$expID))

lDat <- unlist(vDat,use.names=FALSE)
ylim <- max(replace(lDat,is.na(lDat),0) )
xbarcenters <- barplotPerCondition(datHuntStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
for (g in strCondTags)
{
  idx <- match(g,strCondTags)
  vpt <- unlist(vDat[g],use.names=TRUE)
  ##OPTIONAL - PLot NA points As Zero - 
  ##vpt <-replace(vpt,is.na(vpt),0) ##NA Appear Where for exp Where no Hunting Episodes Occured
  
  #points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=idx,col=r[idx] )
  #points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=as.numeric(vDatSetID[[g]]),col=rDataset[as.numeric(vDatSetID[[g]])] )
}

plotConnectedPointsPairs(vIDTable,vDat,strCondTags,xbarcenters)

if (!bonefile)
  dev.off()


###### Histogram ##########
  strPlotName = "plots/meanEpisodeDurationOfGroup_Hist.pdf"
  pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
  lDeltas <- plotPairedChangeHistogram(vIDTable,vDat,strCondTags,uLim = 500,lLim = -610)
  title("Episode Duration", sub = strsub, cex.main = 1.2,   font.main= 1.5, col.main= "black", cex.sub = 1.0, font.sub = 2, col.sub = "black")
  dev.off()
  
  ## Plot Mean Change ######
  strPlotName = "plots/meanEpisodeDurationChange.pdf"
  meanDelta <- sapply(lDeltas,mean, na.rm = TRUE)
  seDelta   <-  sapply(lDeltas,sd,na.rm = TRUE)/sqrt(sapply(lDeltas[is.na(lDeltas)==FALSE],NROW))
  ylim <- 4*max(unlist(meanDelta,use.names=FALSE))
  pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
  #plot.window(xlim=c(0,100), ylim=c(-10, ylim) )
  barplotPerCondition(lDeltas,meanDelta,seDelta,"Mean Change in Episode Duration ",strsub,strPlotName,ylim)
  dev.off()
  

#######################################




######## HUNTING DURATION PER LARVA ############
#X11()
strPlotName = "plots/meanHuntDurationOfGroup.pdf"
datmean <- unlist(datHuntStat[,"meanDuration"],use.names = FALSE)
datse   <- unlist(datHuntStat[,"seDuration"],use.names = FALSE)
vDat    <- datHuntStat[,"vHDurationPerLarva"] #Total H Duration Per Larva
vDatSetID <- datHuntStat[,"vDataSetID"]
strtitle <- "Duration of Hunting per Larva"

if (!bonefile)
  pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$expID))

ylim <- max(unlist(vDat,use.names=FALSE))
xbarcenters <- barplotPerCondition(datHuntStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
# for (g in strCondTags)
# {
#   idx <- match(g,strCondTags)
#   vpt = unlist(vDat[g],use.names=FALSE)
#   #points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=idx,col=r[idx] )
#   points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=as.numeric(vDatSetID[[g]]),col=rDataset[as.numeric(vDatSetID[[g]])] )
# }

plotConnectedPointsPairs(vIDTable,vDat,strCondTags,xbarcenters)

if (!bonefile)
  dev.off()


###### Histogram ##########
strPlotName = "plots/HuntDuration_Hist.pdf"
pdf(strPlotName,width=8,height=8,title=strtitle,onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))
lDeltas <- plotPairedChangeHistogram(vIDTable,vDat,strCondTags,uLim = 25000,lLim = -10000)
title("Hunt Duration", sub = strsub, cex.main = 1.2,   font.main= 1.5, col.main= "black", cex.sub = 1.0, font.sub = 2, col.sub = "black")
dev.off()

strPlotName = "plots/HuntDurationChange.pdf"
meanDelta <- sapply(lDeltas,mean, na.rm = TRUE)
seDelta   <-  sapply(lDeltas,sd,na.rm = TRUE)/sqrt(sapply(lDeltas[is.na(lDeltas)==FALSE],NROW))
ylim <- 4*max(unlist(meanDelta,use.names=FALSE))
pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
#plot.window(xlim=c(0,100), ylim=c(-10, ylim) )
barplotPerCondition(lDeltas,meanDelta,seDelta,"Mean Change in Episode Duration ",strsub,strPlotName,ylim)
dev.off()


#######################################


#######################################



############### HUNT RATIOn ###############
#X11()

strPlotName <- "plots/meanHuntRatioOfGroup.pdf"
datmean     <- unlist(datHuntStat[,"meanHuntRatioOfGroup"],use.names = FALSE)
datse       <- unlist(datHuntStat[,"seHuntRatioOfGroup"],use.names = FALSE)
vDat        <- datHuntStat[,"vLarvaHRatio"]
vDatSetID   <- datHuntStat[,"vDataSetID"]
strtitle    <- "Ratio of Time spent Hunting Over all Frames"

if (!bonefile)
  pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$expID))

ylim <- max(unlist(vDat,use.names=FALSE))
xbarcenters <- barplotPerCondition(datHuntStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
# for (g in strCondTags)
# {
#   idx <- match(g,strCondTags)
#   vpt = unlist(vDat[g],use.names=FALSE)
#   #points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=idx,col=r[idx] )
#   points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=as.numeric(vDatSetID[[g]]),col=rDataset[as.numeric(vDatSetID[[g]])] )
# }

plotConnectedPointsPairs(vIDTable,vDat,strCondTags,xbarcenters)


dev.off()


###### Histogram ##########
strPlotName = "plots/HuntRatio_Hist.pdf"
pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
lDeltas <- plotPairedChangeHistogram(vIDTable,vDat,strCondTags,uLim = 0.3,lLim = -0.3)
title("Hunt Ratio", sub = strsub, cex.main = 1.2,   font.main= 1.5, col.main= "black", cex.sub = 1.0, font.sub = 2, col.sub = "black")
dev.off()

## Plot Mean Change ######
strPlotName = "plots/HuntRatioChange.pdf"
meanDelta <- sapply(lDeltas,mean, na.rm = TRUE)
seDelta   <-  sapply(lDeltas,sd,na.rm = TRUE)/sqrt(sapply(lDeltas[is.na(lDeltas)==FALSE],NROW))
ylim <- 4*max(unlist(meanDelta,use.names=FALSE))
pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
#plot.window(xlim=c(0,100), ylim=c(-10, ylim) )
barplotPerCondition(lDeltas,meanDelta,seDelta,"Mean Change in Episode Duration ",strsub,strPlotName,ylim)
dev.off()

#######################################

###Find Hunt Rate Per  Prey Count



###############################################
