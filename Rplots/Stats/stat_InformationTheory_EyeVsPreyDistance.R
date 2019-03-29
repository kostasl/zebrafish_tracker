## 23-10-2018
### Consider Information content in eye vergence on the distance to prey
## We use Model and the sampled parameter values to obtain an estimate of the mean 
##  information content in each hunt episode 
## We calculate the mutual information using P(Phi | X), assuming X is uniform

source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")
source("Stats/stat_InformationTheory_EyeVsPreyDistance_lib.R")

fitseqNo <- 11
strTag <- "VarD" ##Str Appended to output to indicate calc Conditions : VarD : variable integral period depending on hunt event strike distance, FixD : All inf integrated until 0.6mm from prey

#### CalcInformation ##
load(file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RJags_LL",fitseqNo,".RData",sep="") )
load(file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RJags_NL",fitseqNo,".RData",sep="") )
load(file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RJags_DL",fitseqNo,".RData",sep="") )

## Can Load Last Inf Matrix 
#load(file=paste(strDataExportDir,"/stat_infoMat_EyeVergenceVsDistance_sigmoidFit5mm-5bit.RData",sep=""))
       

#### Load the Tracked Hunts Register ###
strRegisterDataFileName <- paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register",".rds",sep="") #Processed Registry on which we add 
message(paste(" Importing Retracked HuntEvents from:",strRegisterDataFileName))
datTrackedEventsRegister <- readRDS(strRegisterDataFileName) ## THis is the Processed Register File On 

strGroupID <- levels(datTrackedEventsRegister$groupID)

NSamples <-50 ## Number of Fit Lines (rows) to add in info matrrix


##For the 3 Groups 
colourH <- c(rgb(0.01,0.01,0.9,0.8),rgb(0.01,0.7,0.01,0.8),rgb(0.9,0.01,0.01,0.8),rgb(0.00,0.00,0.0,1.0)) ##Legend
colourP <- c(rgb(0.01,0.01,0.8,0.5),rgb(0.01,0.6,0.01,0.5),rgb(0.8,0.01,0.01,0.5),rgb(0.20,0.40,0.5,1.0)) ##points DL,LL,NL
colourR <- c(rgb(0.01,0.01,0.9,0.4),rgb(0.01,0.7,0.01,0.4),rgb(0.9,0.01,0.01,0.4),rgb(0.00,0.00,0.0,0.3)) ##Region (Transparency)
pchL <- c(16,2,4)
ltL <-  c(1,2,3) ##line Types

### Begin Script ###

## Sample Matrices Of Information / Retursn Struct Containing Mat and Id vectors ##
lInfStructLL <- calcInfoOfHuntEvent(drawLL,dataLL,groupID=2)
lInfStructNL <- calcInfoOfHuntEvent(drawNL,dataNL,groupID=3)
lInfStructDL <- calcInfoOfHuntEvent(drawDL,dataDL,groupID=1)

mInfMatrixLL <- lInfStructLL$infoMatrix
mInfMatrixNL <- lInfStructNL$infoMatrix
mInfMatrixDL <- lInfStructDL$infoMatrix


X11()
hist(mInfMatrixLL,col=colourH[2],xlim=c(0,3),breaks = seq(-1,3,1/20) )
X11()
hist(mInfMatrixDL,col=colourH[1],xlim=c(0,3),breaks = seq(-1,3,1/20))
X11()
hist(mInfMatrixNL,col=colourH[3],xlim=c(0,3),breaks = seq(-1,3,1/20))


X11()
hist(colMeans(mInfMatrixLL),col=colourH[2],xlim=c(0,3),breaks = seq(0,3,1/20),main="mean Inf Per Event ")
X11()
hist(colMeans(mInfMatrixDL),col=colourH[1],xlim=c(0,3),breaks = seq(0,3,1/20))
X11()
hist(colMeans(mInfMatrixNL),col=colourH[3],xlim=c(0,3),breaks = seq(0,3,1/20))


### Plot CDF ###
## Match the N 
pdf(file= paste(strPlotExportPath,"/stat/stat_InfSigmoidExp_EyeVsDistance_CDF_",strTag,".pdf",sep=""))
subset_mInfMatrixLL <- mInfMatrixLL[,sample(1:58,58)]
plot(ecdf(mInfMatrixDL),col=colourH[1],main="Information In Eye Vergence CDF",
     xlab="Information (bits) ",lty=1,lwd=2,xlim=c(0,2.5))
plot(ecdf(mInfMatrixLL),col=colourH[2],add=T,lty=2,lwd=2)
plot(ecdf(mInfMatrixNL),col=colourH[3],add=T,lty=3,lwd=2)
legend("topleft",legend=paste(c("DL n=","LL n=","NL n="),c(NCOL(mInfMatrixDL),NCOL(mInfMatrixLL) ,NCOL(mInfMatrixNL) ) ) 
       ,col=colourH,lty=c(1,2,3),lwd=2)
dev.off()


##Plot Density 
dLLphi<-density(mInfMatrixLL)
dDLphi<-density(mInfMatrixDL)
dNLphi<-density(mInfMatrixNL)

#X11()
pdf(file= paste(strPlotExportPath,"/stat/stat_InfSigmoidExp_EyeVsDistance_Density_",strTag,".pdf",sep=""))
plot(dLLphi,col=colourH[2],type="l",lwd=3,lty=2,
     ylim=c(0,2),main="Mutual Information Distance to Eye Vergence ",
     xlab=expression(paste(" (bits)" ) ) )
lines(dNLphi,col=colourH[3],lwd=3,lty=3)
lines(dDLphi,col=colourH[1],lwd=3,lty=1)
legend("topleft",legend=paste(c("DL n=","LL n=","NL n="),c(NCOL(mInfMatrixDL),NCOL(subset_mInfMatrixLL) ,NCOL(mInfMatrixNL) ) ) 
       ,col=colourH,lty=c(1,2,3),lwd=2)
dev.off()

save(lInfStructLL,lInfStructDL,lInfStructNL,drawLL,drawDL,drawNL,file=paste(strDataExportDir,"/stat_infoMat_EyeVergenceVsDistance_sigmoidFit5mm-5bit_",strTag,".RData",sep=""))      
####



########## INF VS UNDERSHOOT ANALYSIS #######################
####  Examine Correlation Of Information Vs Undershooting ###

### LOAD DATA ###
##Load List the information measured ##
load(file=paste(strDataExportDir,"/stat_infoMat_EyeVergenceVsDistance_sigmoidFit5mm-5bit_",strTag,".RData",sep=""))
### Load Regression Data ###
#load(paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RJAgsOUt_",fitseqNo,".RData",sep=""))
load(file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RJags_LL",fitseqNo,".RData",sep="") )
load(file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RJags_NL",fitseqNo,".RData",sep="") )
load(file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RJags_DL",fitseqNo,".RData",sep="") )

## Can Load Last Inf Matrix 
## Load the First Bout turn data --
lFirstBoutPoints <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData",".rds",sep="") ) #Processed Registry on which we add )
#####

##Convert to data frame 
datFirstBouts <- data.frame(do.call(rbind,lFirstBoutPoints )) ##data.frame(lFirstBoutPoints[["LL"]] )
##Select The First Bout Relevant Records - 
#Note we are subsetting the events at the intersection of having a 1st bout turn to prey AND a strike to prey

##Make Data Frame For Debug of Events #
datTrackedEventsRegisterWithFirstBout <- datTrackedEventsRegister[datFirstBouts$RegistarIdx,]

datFirstBoutVsInfLL <- mergeFirstTurnWithInformation(datFirstBouts,lInfStructLL  )
## Error Check That subset of data using RegIdx Indeed belongs to the intended group
stopifnot(unique(datTrackedEventsRegister[unlist(datFirstBoutVsInfLL$RegistarIdx),"groupID"]) == "LL" )  ##Check for Errors in Reg idx - Group should match registry

datFirstBoutVsInfNL <- mergeFirstTurnWithInformation(datFirstBouts,lInfStructNL  )
## Error Check That subset of data using RegIdx Indeed belongs to the intended group
stopifnot(unique(datTrackedEventsRegister[unlist(datFirstBoutVsInfNL$RegistarIdx),"groupID"]) == "NL" )  ##Check for Errors in Reg idx - Group should match registry

datFirstBoutVsInfDL <- mergeFirstTurnWithInformation(datFirstBouts,lInfStructDL  )
## Error Check That subset of data using RegIdx Indeed belongs to the intended group
stopifnot(unique(datTrackedEventsRegister[unlist(datFirstBoutVsInfDL$RegistarIdx),"groupID"]) == "DL" )  ##Check for Errors in Reg idx - Group should match registry

datVerifyNL <- cbind(lInfStructNL$vsampleRegisterIdx, unique(dataNL$RegistrarIdx)[ lInfStructNL$vsamplePSeqIdx])

### Plot Information Vs Undershoot Ratio ####
pdf(file= paste(strPlotExportPath,"/stat/stat_InfVsTurnRatio_",strTag,".pdf",sep=""))
plot(unlist(datFirstBoutVsInfLL$UnderShootRatio),unlist(datFirstBoutVsInfLL$MInf),
     ylim=c(0,2),xlim=c(0,2),
     xlab=( expression(paste(Phi,"/",theta," Turn Ratio  ") )  ),ylab="mutual Inf in Eye V",
     main="Information Vs Undershoot ",col=colourH[2],pch=pchL[2])
segments(1,-10,1,20);
points(unlist(datFirstBoutVsInfNL$UnderShootRatio),unlist(datFirstBoutVsInfNL$MInf),ylim=c(0,2),xlim=c(0,2),
       col=colourH[3],pch=pchL[3])
points(unlist(datFirstBoutVsInfDL$UnderShootRatio),unlist(datFirstBoutVsInfDL$MInf),ylim=c(0,2),xlim=c(0,2),
       col=colourH[1],pch=pchL[1])
legend("topleft",legend=paste(c("DL n=","LL n=","NL n="),c(NROW(datFirstBoutVsInfDL),NROW(datFirstBoutVsInfLL) ,NROW(datFirstBoutVsInfNL) ) ) 
       ,col=colourH,pch=pchL,lty=c(1,2,3),lwd=2)
dev.off()

boxplot(unlist(datFirstBoutVsInfLL$MInf),unlist(datFirstBoutVsInfNL$MInf),unlist(datFirstBoutVsInfDL$MInf),names=c("LL","NL","DL") )
boxplot(unlist(datFirstBoutVsInfLL$UnderShootRatio),unlist(datFirstBoutVsInfNL$UnderShootRatio),unlist(datFirstBoutVsInfDL$UnderShootRatio),names=c("LL","NL","DL") )

## plot Undershot Ratio  - Showing RegisterIDs ###
##dataLL$distMax[ lInfStructDL$vsamplePSeqIdx]
## Verification that distMax from Data belongs to the same event in lInfStructNL example : 
## 
pdf(file= paste(strPlotExportPath,"/stat/stat_InfVsTurnRatio_",strTag,"_WithIDs.pdf",sep=""))
plot(unlist(datFirstBoutVsInfLL$UnderShootRatio),unlist(datFirstBoutVsInfLL$MInf),
     ylim=c(0,2),xlim=c(0,2),
     xlab=( expression(paste(Phi,"/",theta," Turn Ratio  ") )  ),ylab="mutual Inf in Eye V",
     main="Information Vs Undershoot ",col=colourH[2],pch=pchL[2])
segments(1,-10,1,20);
points(unlist(datFirstBoutVsInfNL$UnderShootRatio),unlist(datFirstBoutVsInfNL$MInf),ylim=c(0,2),xlim=c(0,2),
       col=colourH[3],pch=pchL[3])
points(unlist(datFirstBoutVsInfDL$UnderShootRatio),unlist(datFirstBoutVsInfDL$MInf),ylim=c(0,2),xlim=c(0,2),
       col=colourH[1],pch=pchL[1])
text(unlist(datFirstBoutVsInfNL$UnderShootRatio)*1.01,unlist(datFirstBoutVsInfNL$MInf)*1.01,datTrackedEventsRegister[ unlist(datFirstBoutVsInfNL$RegistarIdx),"expID"],cex=0.7,col=colourP[3])
text(unlist(datFirstBoutVsInfLL$UnderShootRatio)*1.01,unlist(datFirstBoutVsInfLL$MInf)*1.01,datTrackedEventsRegister[ unlist(datFirstBoutVsInfLL$RegistarIdx),"expID"],cex=0.7)
text(unlist(datFirstBoutVsInfDL$UnderShootRatio)*1.01,unlist(datFirstBoutVsInfDL$MInf)*1.01,datTrackedEventsRegister[ unlist(datFirstBoutVsInfDL$RegistarIdx),"expID"],cex=0.7,col=colourP[1])

legend("topright",legend=paste(c("DL n=","LL n=","NL n="),c(NROW(datFirstBoutVsInfDL),NROW(datFirstBoutVsInfLL) ,NROW(datFirstBoutVsInfNL) ) ) 
       ,col=colourH,pch=pchL,lty=c(1,2,3),lwd=2)
dev.off()

## Plot UNDERSHOOT With ONSET Distance Information ####
##Make Colour Map Normalized from 0mm to 5mm 
rfcDark <- colorRampPalette(rev(brewer.pal(8,'YlOrRd')));
colMap5mm <- rfcDark(61)
distOnSet_NL <-  ( as.numeric(prettyNum( unlist( dataNL$distMax[ lInfStructNL$vsamplePSeqIdx]), digits=3) ) )
colMap5mm_NL <- colMap5mm[ (distOnSet_NL *10 ) ]
distOnSet_DL <- as.numeric(prettyNum( unlist( dataDL$distMax[ lInfStructDL$vsamplePSeqIdx]), digits=3) ) 
colMap5mm_DL <- colMap5mm[ (distOnSet_DL  *10 ) ]
distOnSet_LL <- as.numeric(prettyNum( unlist( dataLL$distMax[ lInfStructLL$vsamplePSeqIdx]), digits=3) ) 
colMap5mm_LL <- colMap5mm[ distOnset_LL*10  ]

#par(bg = "white")
#par(fg = "black")
pdf(file= paste(strPlotExportPath,"/stat/stat_InfVsTurnRatioAndDistance_",strTag,"_LF.pdf",sep=""))
plot(unlist(datFirstBoutVsInfLL$UnderShootRatio),unlist(datFirstBoutVsInfLL$MInf),
     ylim=c(0,2),xlim=c(0,2),
     xlab=( expression(paste(Phi,"/",theta," Turn Ratio  ") )  ),ylab="mutual Inf in Eye V",
     main="Information Vs Undershoot LF",bg="black",col=colMap5mm_LL,pch=17)
text(1.5,1.8,"0mm"); text(2.0,1.8,"6mm")
points(seq(1.5,2.0,(2.0-1.5)/( NROW(colMap5mm)-1)  ),rep(1.7,NROW(colMap5mm)),       pch=15,col=colMap5mm )
points(unlist(datFirstBoutVsInfLL$UnderShootRatio), unlist(datFirstBoutVsInfLL$MInf), pch= 2)
segments(1,-10,1,20);
plot(distOnset_LL,unlist(datFirstBoutVsInfLL$MInf),
     xlab=( expression(paste(tau," Onset distance ") )  ),ylab="mutual Inf in Eye V",
     xlim=c(0,6),ylim=c(0,2),
     main="Information Vs Onset Distance ",bg="black",col=colMap5mm_LL,pch=17) #col=colMap5mm_LL
dev.off()


pdf(file= paste(strPlotExportPath,"/stat/stat_InfVsTurnRatioAndDistance_",strTag,"_NF.pdf",sep=""))
plot(unlist(datFirstBoutVsInfNL$UnderShootRatio),unlist(datFirstBoutVsInfNL$MInf),
     ylim=c(0,2),xlim=c(0,2),
     xlab=( expression(paste(Phi,"/",theta," Turn Ratio  ") )  ),ylab="mutual Inf in Eye V",
     main="Information Vs Undershoot NF",bg="black",col=colMap5mm_LL,pch=17)
text(1.5,1.8,"0mm"); text(2.0,1.8,"6mm")
points(seq(1.5,2.0,(2.0-1.5)/( NROW(colMap5mm)-1)  ),rep(1.7,NROW(colMap5mm)),pch=15,col=colMap5mm )
points(unlist(datFirstBoutVsInfNL$UnderShootRatio), unlist(datFirstBoutVsInfNL$MInf), pch= 2)
segments(1,-10,1,20);
##2nd Page plot
plot(distOnSet_NL,unlist(datFirstBoutVsInfNL$MInf),
     xlim=c(0,6),ylim=c(0,2),
     xlab=( expression(paste(tau," Onset distance ") )  ),ylab="mutual Inf in Eye V",
     main="Information Vs Onset Distance ",bg="black",col=colMap5mm_NL,pch=15) # colMap5mm_NL
dev.off()


pdf(file= paste(strPlotExportPath,"/stat/stat_InfVsTurnRatioAndDistance_",strTag,"_DF.pdf",sep=""))
plot(unlist(datFirstBoutVsInfDL$UnderShootRatio),unlist(datFirstBoutVsInfDL$MInf),
     ylim=c(0,2),xlim=c(0,2),
     xlab=( expression(paste(Phi,"/",theta," Turn Ratio  ") )  ),ylab="mutual Inf in Eye V",
     main="Information Vs Undershoot DF",bg="black",col=colMap5mm_LL,pch=17)
text(1.5,1.8,"0mm"); text(2.0,1.8,"6mm")
points(seq(1.5,2.0,(2.0-1.5)/( NROW(colMap5mm)-1)  ),rep(1.7,NROW(colMap5mm)),pch=15,col=colMap5mm )
points(unlist(datFirstBoutVsInfDL$UnderShootRatio), unlist(datFirstBoutVsInfDL$MInf), pch= 2)
segments(1,-10,1,20);
##2nd Page plot
plot(distOnSet_DL,unlist(datFirstBoutVsInfDL$MInf),
     xlim=c(0,6),ylim=c(0,2),
     xlab=( expression(paste(tau," Onset distance ") )  ),ylab="mutual Inf in Eye V",
     main="Information Vs Onset Distance ",bg="black",col=colMap5mm_DL ,pch=19)#colMap5mm_DL
dev.off()


pdf(file= paste(strPlotExportPath,"/stat/stat_InfVsTurnRatioAndDistanceNorm_",strTag,"_LF.pdf",sep=""))
plot(unlist(datFirstBoutVsInfLL$UnderShootRatio),unlist(datFirstBoutVsInfLL$MInf)/(distOnset_LL),
     ylim=c(0,1),xlim=c(0,2),
     xlab=( expression(paste(Phi,"/",theta," Turn Ratio  ") )  ),ylab="mutual Inf in Eye V",
     main="Information/mm Vs Undershoot LF",bg="black",col=colMap5mm_LL,pch=17)
text(1.5,0.87,"0mm"); text(2.0,0.87,"6mm")
points(seq(1.5,2.0,(2.0-1.5)/( NROW(colMap5mm)-1)  ),rep(0.8,NROW(colMap5mm)),pch=15,col=colMap5mm )
points(unlist(datFirstBoutVsInfLL$UnderShootRatio), unlist(datFirstBoutVsInfLL$MInf)/(distOnset_LL), pch= 2)
segments(1,-10,1,20);
plot(distOnset_LL,unlist(datFirstBoutVsInfLL$MInf),
     xlab=( expression(paste(tau," Onset distance ") )  ),ylab="mutual Inf in Eye V",
     xlim=c(0,6),ylim=c(0,2),
     main="Information Vs Onset Distance ",bg="black",col=colMap5mm_LL,pch=17) #col=colMap5mm_LL
dev.off()


pdf(file= paste(strPlotExportPath,"/stat/stat_InfVsTurnRatioAndDistanceNorm_",strTag,"_DF.pdf",sep=""))
plot(unlist(datFirstBoutVsInfDL$UnderShootRatio),unlist(datFirstBoutVsInfDL$MInf)/(distOnSet_DL),
     ylim=c(0,1),xlim=c(0,2),
     xlab=( expression(paste(Phi,"/",theta," Turn Ratio  ") )  ),ylab="mutual Inf in Eye V",
     main="Information/mm Vs Undershoot DF",bg="black",col=colMap5mm_LL,pch=17)
text(1.5,0.87,"0mm"); text(2.0,0.87,"6mm")
points(seq(1.5,2.0,(2.0-1.5)/( NROW(colMap5mm)-1)  ),rep(0.8,NROW(colMap5mm)),pch=15,col=colMap5mm )
points(unlist(datFirstBoutVsInfDL$UnderShootRatio), unlist(datFirstBoutVsInfDL$MInf)/(distOnSet_DL), pch= 2)
segments(1,-10,1,20);
##2nd Page plot
plot(distOnSet_DL,unlist(datFirstBoutVsInfDL$MInf),
       xlim=c(0,6),ylim=c(0,2),
       xlab=( expression(paste(tau," Onset distance ") )  ),ylab="mutual Inf in Eye V",
       main="Information Vs Onset Distance ",bg="black",col=colMap5mm_DL ,pch=19)#colMap5mm_DL
dev.off()

pdf(file= paste(strPlotExportPath,"/stat/stat_InfVsTurnRatioAndDistanceNorm_",strTag,"_NF.pdf",sep=""))
plot(unlist(datFirstBoutVsInfNL$UnderShootRatio),unlist(datFirstBoutVsInfNL$MInf)/(distOnSet_NL),
     ylim=c(0,1),xlim=c(0,2),
     xlab=( expression(paste(Phi,"/",theta," Turn Ratio  ") )  ),ylab="mutual Inf in Eye V",
     main="Information/mm Vs Undershoot NF ",bg="black",col=colMap5mm_LL,pch=17)
text(1.5,0.87,"0mm"); text(2.0,0.87,"6mm")
points(seq(1.5,2.0,(2.0-1.5)/( NROW(colMap5mm)-1)  ),rep(0.8,NROW(colMap5mm)),pch=15,col=colMap5mm )
points(unlist(datFirstBoutVsInfNL$UnderShootRatio), unlist(datFirstBoutVsInfNL$MInf)/(distOnSet_NL), pch= 2)
segments(1,-10,1,20);
##2nd Page plot
plot(distOnSet_NL,unlist(datFirstBoutVsInfNL$MInf),
       xlim=c(0,6),ylim=c(0,2),
       xlab=( expression(paste(tau," Onset distance ") )  ),ylab="mutual Inf in Eye V",
       main="Information Vs Onset Distance ",bg="black",col=colMap5mm_NL,pch=15) # colMap5mm_NL
dev.off()

### Inf Vs DISTANCE Plot All Groups 

pdf(file= paste(strPlotExportPath,"/stat/stat_InfVsDistance_",strTag,"_All.pdf",sep=""))
plot(distOnset_LL,unlist(datFirstBoutVsInfLL$MInf),
     xlab=( expression(paste(tau," Onset distance ") )  ),ylab="mutual Inf in Eye V",
     xlim=c(0,6),ylim=c(0,2),
     main="Information Vs Onset Distance ",bg="black",col=colourP[2],pch=17) #col=colMap5mm_LL

points(distOnSet_DL,unlist(datFirstBoutVsInfDL$MInf),
       xlim=c(0,6),ylim=c(0,2),
       xlab=( expression(paste(tau," Onset distance ") )  ),ylab="mutual Inf in Eye V",
       main="Information Vs Onset Distance ",bg="black",col=colourP[1] ,pch=19)#colMap5mm_DL

points(distOnSet_NL,unlist(datFirstBoutVsInfNL$MInf),
       xlim=c(0,6),ylim=c(0,2),
       xlab=( expression(paste(tau," Onset distance ") )  ),ylab="mutual Inf in Eye V",
       main="Information Vs Onset Distance ",bg="black",col=colourP[3],pch=15) # colMap5mm_NL
legend("topright",legend=paste(c("DL n=","LL n=","NL n="),c(NROW(datFirstBoutVsInfDL),NROW(datFirstBoutVsInfLL) ,NROW(datFirstBoutVsInfNL) ) ) 
       ,col=colourP,pch=c(19,17,15),lty=c(1,2,3),lwd=2)
dev.off()





## Make Bar plot Of Undershoot Vs Inf - Bin Results
datFirstBoutVsInf <- rbind(datFirstBoutVsInfLL,datFirstBoutVsInfNL,datFirstBoutVsInfDL)
bins <-  20
vlevels <- seq(from=1,to=10,by=10/bins)
vUndershootBin <- cut(unlist(datFirstBoutVsInf$UnderShootRatio),breaks=seq(0.0,2,by=2/bins))
vInfVsUndershoot <- cbind(datFirstBoutVsInf$MInf,vUndershootBin)
##Do Mean For Each Level/Bin
lInfPerUndershoot = list()
vmeanInf = vector()
vsdInf = vector()
idx <- 0
for (p in vlevels)
{
  idx <- idx +1
  lInfPerUndershoot[[idx]] <- unlist(vInfVsUndershoot[vInfVsUndershoot[,2] == p,1])
  vmeanInf[idx] <- mean(unlist(vInfVsUndershoot[vInfVsUndershoot[,2] == p,1]))
  vsdInf[idx]   <- sd(unlist(vInfVsUndershoot[vInfVsUndershoot[,2] == p,1]))
}

pdf(file= paste(strPlotExportPath,"/stat/stat_boxplot_InfVsUndershoot_AllGroups_",strTag,".pdf",sep=""))
## Box Plot Results Undershoot Vs Information
boxplot(unlist(vInfVsUndershoot[,1])~unlist(vInfVsUndershoot[,2]),names=sort(unique(vUndershootBin) ),
        xlab="Overshoot ratio",
        ylab="Distance Information (bits)",
        main="Relating Undershoot to Larva's Ability to Infer Distance  ")

dev.off()


##Plot Vs Undershoot Angle ###
plot(datFirstBoutVsInfLL$UnderShootAngle,datFirstBoutVsInfLL$MInf,
     ylim=c(0,2),xlim=c(-60,60),xlab=("OnSetAngleToPrey - Turn Angle"),ylab="mutual Inf in Eye V",
     main="Information Vs Undershoot ",pch=pchL[2])
points(datFirstBoutVsInfNL$UnderShootAngle,datFirstBoutVsInfNL$MInf,ylim=c(0,2),
       col="red",pch=pchL[3])
points(datFirstBoutVsInfDL$UnderShootAngle,datFirstBoutVsInfDL$MInf,ylim=c(0,2),
       col="blue",pch=pchL[1])



########








## Entropy ##
#if you use this package please cite: Jean Hausser and Korbinian Strimmer. 2009. Entropy inference
#and the James-Stein estimator, with application to nonlinear gene association networks.  J. Mach.
#Learn.  Res.10:  1469-1484.  Available online fromhttp://jmlr.csail.mit.edu/papers/v10/  hausser09a.html.
library("entropy")

KL.empirical(binNLphi,binLLphi,unit="log2")


mp <- barplot(c(entropy(binDLphi,unit="log2"),entropy(binLLphi,unit="log2"),entropy(binNLphi,unit="log2")),col=colourR,xlab="Group ",
              main=expression(paste("Entropy of eye vergence angle during hunt" ) ),sub=expression(paste("Resolution: ", Phi,": 1 deg" ) ) ,
              ylab="(bits)", ylim=c(0,8) ) #legend.text = strGroupID
text(mp,y=c(-0.3,-0.3,-0.3),labels=strGroupID , xpd = TRUE, col = "black")


## Mutual Information ##
binLLphiVsDist <- discretize2d(dataLL$phi, dataLL$distP, numBins1=80, numBins2=10,r1=c(0,80),r2=c(0.5,2))
binNLphiVsDist <- discretize2d(dataNL$phi, dataNL$distP, numBins1=80, numBins2=10,r1=c(0,80),r2=c(0.5,2))
binDLphiVsDist <- discretize2d(dataDL$phi, dataDL$distP, numBins1=80, numBins2=10,r1=c(0,80),r2=c(0.5,2))

## Joint Entropy
H12_LL = entropy(binLLphiVsDist,unit="log2" )
H12_NL = entropy(binNLphiVsDist ,unit="log2")
H12_DL = entropy(binDLphiVsDist,unit="log2" )
##Plot Mutual Information
mp <- barplot(c(H12_DL,H12_LL,H12_NL),col=colourR,xlab="Group ",
              main=expression(paste("Joint entropy prey distance and eye vergence" ) ),sub=expression(paste("Resolution: ", Phi,": 1 deg, Distance : 0.5mm" ) ) ,
              ylab="(bits)", ylim=c(0,10) ) #legend.text = strGroupID
text(mp,y=c(-0.3,-0.3,-0.3),labels=strGroupID , xpd = TRUE, col = "black")

# mutual information
mi.empirical(binLLphiVsDist,unit="log2" ) # //Mutual Informations
mi.empirical(binNLphiVsDist,unit="log2" ) # approximately zero
mi.empirical(binDLphiVsDist,unit="log2" ) # approximately zero

# another way to compute mutual information
# compute marginal entropies
H12_LL = entropy(binLLphiVsDist,unit="log2" )
H1_LL = entropy(rowSums(binLLphiVsDist),unit="log2")
H2_LL = entropy(colSums(binLLphiVsDist),unit="log2")
H1_LL+H2_LL-H12_LL # mutual entropy

H12_NL = entropy(binNLphiVsDist,unit="log2" )
H1_NL = entropy(rowSums(binNLphiVsDist),unit="log2")
H2_NL = entropy(colSums(binNLphiVsDist),unit="log2")
H1_NL+H2_NL-H12_NL # mutual entropy


H12_DL = entropy(binDLphiVsDist,unit="log2" )
H1_DL = entropy(rowSums(binDLphiVsDist),unit="log2")
H2_DL = entropy(colSums(binDLphiVsDist),unit="log2")
H1_DL+H2_DL-H12_DL # mutual entropy



# sample from continuous uniform distribution
x1 = runif(10000)
X11()
hist(x1, xlim=c(0,1), freq=FALSE)
# discretize into 10 categories
y1 = discretize(x1, numBins=10, r=c(0,1))
# compute entropy from counts
entropy(y1) # empirical estimate near theoretical maximu
log(10)

#discretize2d(x1, x2, numBins1=10, numBins2=10)

# sample from a non-uniform distribution
X11()
hist(x2, xlim=c(0,1), freq=FALSE)
# discretize into 10 categories and estimate entropy

#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_Rate_lambda_LL_E.pdf",sep=""))
X11()
hist(drawLL$lambda[,,],col=colourR[2],main="Lambda")
hist(drawNL$lambda[,,],col=colourR[3],add=TRUE)
hist(drawDL$lambda[,,],col=colourR[1],add=TRUE)

X11()
hist(drawLL$gamma[,,],col=colourR[2],main="gamma")
hist(drawNL$gamma[,,],col=colourR[3],add=TRUE)
hist(drawDL$gamma[,,],col=colourR[1],add=TRUE)


X11()
hist(drawLL$phi_0[,,],col=colourR[2],main="Phi 0")
hist(drawNL$phi_0[,,],col=colourR[3],add=TRUE)
hist(drawDL$phi_0[,,],col=colourR[1],add=TRUE)


X11()
hist(drawLL$phi_max[,,],col=colourR[2],main="Phi Max")
hist(drawNL$phi_max[,,],col=colourR[3],add=TRUE)
hist(drawDL$phi_max[,,],col=colourR[1],add=TRUE)
#plot(drawLL$phi_max[3,,])

X11()
hist(drawLL$sigma[,1,,],main="LL")



X11()
#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_StartEnd_u0_NL_E.pdf",sep=""))
hist(drawLL$u1[,,],breaks=50,xlim=c(0,7),col=colourH[2])
hist(drawLL$u0[,,],breaks=50,xlim=c(0,7),add=TRUE,col=colourH[2])

#dev.off()
########################
## NL ###
mNL=jags.model(file="model.tmp",data=dataNL);
drawNL=jags.samples(mNL,steps,thin=thin,variable.names=varnames)


X11()
#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_NL_F.pdf",sep=""))
plotGCRes(drawNL,dataNL,groupID=3)
#dev.off()

## Plot the infered function NL

X11()
#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_Rate_lambda_NL_E.pdf",sep=""))
hist(drawNL$lambda[1,,1],main="NL")
#dev.off()

X11()
#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_StartEnd_u0_NL_E.pdf",sep=""))
hist(drawNL$u1[1,,1],breaks=50,xlim=c(0,7),col="red")
hist(drawNL$u0[1,,1],breaks=50,xlim=c(0,7),add=TRUE,col="red")
#dev.off()





pdf(file= paste(strPlotExportPath,"/stat/stat_densityolinregressionslope.pdf",sep=""))
plot(dDLb,col=colourH[1],xlim=c(0.5,1.2),lwd=3,lty=1,ylim=c(0,20),main="Density Inference of Turn-To-Prey Slope ")
lines(dLLb,col=colourH[2],xlim=c(0.5,1.2),lwd=3,lty=2)

lines(dNLb,col=colourH[3],xlim=c(0.5,1.2),lwd=3,lty=3)
legend("topleft",legend=paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       ,fill=colourL,lty=c(1,2,3))
dev.off()

### PLot Scatter with regression lines with Conf intervals##
#X11()

pdf(file= paste(strPlotExportPath,"/stat/stat_TurnToPrey_LinearRegression.pdf",sep=""))
plot(lFirstBoutPoints[["DL"]][,1], lFirstBoutPoints[["DL"]][,2],
     main=paste("Turn Size Vs Bearing To Prey ", sep=""),
     xlab="Bearing To Prey prior to Bout",ylab="Bearing Change After Bout",xlim=c(-100,100),
     ylim=c(-100,100),
     col=colourP[1] ,pch=pchL[1]) ##boutSeq The order In Which The Occurred Coloured from Dark To Lighter
##Draw 0 Vertical Line
segments(0,-90,0,90); segments(-90,0,90,0); segments(-90,-90,90,90,lwd=1,lty=2);
#text(lFirstBoutPoints[["DL"]][,1]+2,lFirstBoutPoints[["DL"]][,2]+5,labels=lFirstBoutPoints[["DL"]][,3],cex=0.8,col="darkblue")
abline(lm(lFirstBoutPoints[["DL"]][,2] ~ lFirstBoutPoints[["DL"]][,1]),col=colourH[4],lwd=1.0) ##Fit Line / Regression
abline(a=muDLa,b=muDLb,col=colourH[1],lwd=1.5) ##Fit Line / Regression
abline(a=quantile(drawDL$beta[,(steps-ind):steps,1][1,])[2],b=quantile(drawDL$beta[,(steps-ind):steps,1][2,])[2],col=colourR[1],lwd=4.0) ##Fit Line / Regression
abline(a=quantile(drawDL$beta[,(steps-ind):steps,1][1,])[3],b=quantile(drawDL$beta[,(steps-ind):steps,1][2,])[3],col=colourR[1],lwd=4.0) ##Fit Line / Regression

#abline( lsfit(lFirstBoutPoints[["DL"]][,2], lFirstBoutPoints[["DL"]][,1] ) ,col=colourH[1],lwd=2.0)
##LL
points(lFirstBoutPoints[["LL"]][,1], lFirstBoutPoints[["LL"]][,2],pch=pchL[2],col=colourP[2])
#text(lFirstBoutPoints[["LL"]][,1]+2,lFirstBoutPoints[["LL"]][,2]+5,labels=lFirstBoutPoints[["LL"]][,3],cex=0.8,col="darkgreen")
abline(lm(lFirstBoutPoints[["LL"]][,2] ~ lFirstBoutPoints[["LL"]][,1]),col=colourH[4],lwd=1.0)
abline(a=muLLa,b=muLLb,col=colourH[2],lwd=1.5) ##Fit Line / Regression
abline(a=quantile(drawLL$beta[,(steps-ind):steps,1][1,])[2],b=quantile(drawLL$beta[,(steps-ind):steps,1][2,])[2],col=colourR[2],lwd=4.0) ##Fit Line / Regression
abline(a=quantile(drawLL$beta[,(steps-ind):steps,1][1,])[3],b=quantile(drawLL$beta[,(steps-ind):steps,1][2,])[3],col=colourR[2],lwd=4.0) ##Fit Line / Regression

#abline(lsfit(lFirstBoutPoints[["LL"]][,2], lFirstBoutPoints[["LL"]][,1] ) ,col=colourH[2],lwd=2.0)
##NL
points(lFirstBoutPoints[["NL"]][,1], lFirstBoutPoints[["NL"]][,2],pch=pchL[3],col=colourP[3])
#text(lFirstBoutPoints[["NL"]][,1]+2,lFirstBoutPoints[["NL"]][,2]+5,labels=lFirstBoutPoints[["NL"]][,3],cex=0.8,col="darkred")
abline(lm(lFirstBoutPoints[["NL"]][,2] ~ lFirstBoutPoints[["NL"]][,1]),col=colourH[4],lwd=1.0)
abline(a=muNLa,b=muNLb,col=colourH[3],lwd=1.5) ##Fit Line / Regression
abline(a=quantile(drawNL$beta[,(steps-ind):steps,1][1,])[2],b=quantile(drawNL$beta[,(steps-ind):steps,1][2,])[2],col=colourR[3],lwd=4.0) ##Fit Line / Regression
abline(a=quantile(drawNL$beta[,(steps-ind):steps,1][1,])[3],b=quantile(drawNL$beta[,(steps-ind):steps,1][2,])[3],col=colourR[3],lwd=4.0) ##Fit Line / Regression
#abline( lsfit(lFirstBoutPoints[["NL"]][,2], lFirstBoutPoints[["NL"]][,1] ) ,col=colourH[3],lwd=2.0)
legend("topleft",legend=paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       , pch=pchL,col=colourL)

dev.off()





##Plot Densities Summary

sampNL <- coda.samples(mNL,                      variable.names=c("beta","sigma"),                      n.iter=20000, progress.bar="none")
sampDL <- coda.samples(mDL,                      variable.names=c("beta","sigma"),                      n.iter=20000, progress.bar="none")
X11()
plot(sampLL)
X11()
plot(sampNL)
X11()
plot(sampDL,main="DL")


