## Processing Of Tracker Eye TrackerData, Using R 
# Kostasl Nov 2017
#TODO: Add Colour Marker of Hunting On Trajectories
library(tools)
library(RColorBrewer);
library("MASS");
#library(hexbin)
rm("temp","subsetDat","TrackerData","ReyeAll","LeyeAll","frameNAll","ReyeAngleAll","LeyeAngleAll");

nFrWidth = 300;
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')));
r <- rf(30);
#### FUNCTIONS ###

medianf <- function(t,k) {tproc=rep(NA,length(t)); for(i in 1:length(t)) tproc[i]=median(t[max(1,i-k):i]); return(tproc)}
meanf <- function(t,k) {tproc=rep(NA,length(t)); for(i in 1:length(t)) tproc[i]=mean(t[max(1,i-k):i]); return(tproc)}

histj<- function(x,y,x.breaks,y.breaks){
  c1 = as.numeric(cut(x,breaks=x.breaks));
  c2 = as.numeric(cut(y,breaks=y.breaks));
  mat<-matrix(0,ncol=length(y.breaks)-1,nrow=length(x.breaks)-1);
  mat[cbind(c1,c2)] = 1;
  return(mat)
}  

colTraj <- function(x){
  ids = unique(x);
  rcol = rf(length(ids))
  
  #z_scl <- (x - min(x, na.rm=T))/(max(x, na.rm=T) - min(x, na.rm=T))
  #return(r[z_scl*length(r)])
  return(rcol[which(ids %in% x)])
}



### Hunting Episode Analysis ####
source("HuntingEventAnalysis.r")
#################################

### TRAJECTORIES Indicating Hunting  - With distinct colour for each larva ####
source("plotTrackScatterAndDensities.r")
####################


## GLOBAL VARS ##

TrackerData <- list();	
#ttL	<- list();
#ttR	<- list();
#ttV <- list();
LeyeAll <- list();
ReyeAll <- list();
ReyeAngleAll <-vector();
LeyeAngleAll <-vector();
datProcessed <- list();
ttEyeAnglesAll <- list();
frameNAll <- vector();
fileIdxAll <- vector();

lHuntStat <- list();

#### List Of Data files / and result label assuming organized in Directory Structure ###
srcdatList = list()
strCondR  <- "*.csv"; 
srcdatList[["LE"]] <- list(c( list.files(path="./Tracked22-11-17/LiveFed/Empty/", pattern=strCondR,full.names = TRUE),
                              list.files(path="./Tracked30-11-17/LiveFed/Empty/", pattern=strCondR,full.names = TRUE)),"-LiveFed-Empty")
                           
srcdatList[["LL"]] <- list(c(list.files(path="./Tracked22-11-17/LiveFed/Live/", pattern=strCondR,full.names = TRUE),
                             list.files(path="./Tracked30-11-17/LiveFed/Live/", pattern=strCondR,full.names = TRUE) ),"-LiveFed-Live")

srcdatList[["NE"]] <- list(c(list.files(path="./Tracked22-11-17/NotFed/Empty/", pattern=strCondR,full.names = TRUE),
                             list.files(path="./Tracked30-11-17/NotFed/Empty/", pattern=strCondR,full.names = TRUE)),"-NotFed-Empty")

srcdatList[["NL"]] <- list(c(list.files(path="./Tracked22-11-17/NotFed/Live/", pattern=strCondR,full.names = TRUE),
                           list.files(path="./Tracked30-11-17/NotFed/Live/", pattern=strCondR,full.names = TRUE)),"-NotFed-Live")

srcdatList[["DE"]] <- list(c(list.files(path="./Tracked22-11-17/DryFed/Empty/", pattern=strCondR,full.names = TRUE),
                             list.files(path="./Tracked30-11-17/DryFed/Empty/", pattern=strCondR,full.names = TRUE)),"-DryFed-Empty")

srcdatList[["DL"]] <- list(c(list.files(path="./Tracked22-11-17/DryFed/Live/", pattern=strCondR,full.names = TRUE),
                           list.files(path="./Tracked30-11-17/DryFed/Live/", pattern=strCondR,full.names = TRUE)),"-DryFed-Live")

#srcdatList[["LE5"]] <- list(list.files(path="./Tracked30-11-17/LiveFed/Empty/", pattern=strCondR,full.names = TRUE),"-LiveFed-Empty")
#srcdatList[["LL5"]] <- list(list.files(path="./Tracked30-11-17/LiveFed/Live/", pattern=strCondR,full.names = TRUE),"-LiveFed-Live")
#srcdatList[["NE5"]] <- list(list.files(path="./Tracked30-11-17/NotFed/Empty/", pattern=strCondR,full.names = TRUE),"-NotFed-Empty")
#srcdatList[["NL5"]] <- list(list.files(path="./Tracked30-11-17/NotFed/Live/", pattern=strCondR,full.names = TRUE),"-NotFed-Live")
#srcdatList[["DE5"]] <- list(list.files(path="./Tracked30-11-17/DryFed/Empty/", pattern=strCondR,full.names = TRUE),"-DryFed-Empty")
#srcdatList[["DL5"]] <- list(list.files(path="./Tracked30-11-17/DryFed/Live/", pattern=strCondR,full.names = TRUE),"-DryFed-Live")

##CHANGE HASH/ID to select between datasets/groups ##
strCondTags = names(srcdatList);


##Load Files To Tracker Data and Filter Them Out##
procDatIdx = 1;
groupDatIdx = 1;
for (i in strCondTags)
{
  message(paste("#### Load Data Files Of Group ",i," ###############"))

  subsetDat = srcdatList[[i]];
  temp <-  unlist(subsetDat[1])
  TrackerData[[i]] = lapply(temp, read.delim)
  nDat = length(TrackerData[[i]])
  

  groupDatIdx = 1;
  procDatFrames = 0;
  ## FOR EACH DATA FIle IN Group - Filter Data And combine into Single DataFrame For Group ##
  for (j in 1:nDat)
  {
    message(paste(j,". Filtering Data :",  temp[[j]]))
    procDatFrames = procDatFrames + length(TrackerData[[i]][[j]]$frameN);
    message(paste("Found #Rec:",  length(TrackerData[[i]][[j]]$frameN) ))
    
    ##Extract Larva ID
    brokenname = strsplit(temp[[j]],"_")
    larvaID =  as.numeric(brokenname[[1]][length(brokenname[[1]])-3]);
    eventID = as.numeric(brokenname[[1]][length(brokenname[[1]])-2]);
    #Filter Out Empty Files - ones with less than 300 frames ( ~1 sec of data )
    if ( length(TrackerData[[i]][[j]]$frameN) > 300)
    {     
      datProcessed[[procDatIdx]] = data.frame(LEyeAngle= medianf(TrackerData[[i]][[j]]$EyeLDeg,nFrWidth),
                                              REyeAngle= medianf(TrackerData[[i]][[j]]$EyeRDeg,nFrWidth),
                                              posX = TrackerData[[i]][[j]]$Centroid_X,
                                              posY =TrackerData[[i]][[j]]$Centroid_Y,
                                              frameN=TrackerData[[i]][[j]]$frameN,
                                              fileIdx=rep(i,length(TrackerData[[i]][[j]]$EyeLDeg)),
                                              larvaID=rep(larvaID,length(TrackerData[[i]][[j]]$EyeLDeg)),
                                              eventID=rep(eventID,length(TrackerData[[i]][[j]]$EyeLDeg)),
                                              group=rep(i,length(TrackerData[[i]][[j]]$EyeLDeg))
                                              );
      
      procDatIdx = procDatIdx+1; ##INcreased Count Of Processed Files
      groupDatIdx = groupDatIdx + 1;
    }   
  } ##For Each File In Group ##
  datAllFrames = do.call(rbind,datProcessed);
  #xxList = list(LEyeAngle=unlist(ttEyeAnglesAll[,1]),REyeAngle=unlist(ttEyeAnglesAll[,2]),frameN=unlist(ttEyeAnglesAll[,3]),fileIdx=unlist(ttEyeAnglesAll[,4]))
  message(paste("Usable Data files Count :",  groupDatIdx, " total Frames :",procDatFrames));
  
  message(paste("###### Finished Loading Data Files Of Group ",i," ###############"))
  
} ##For Each Group Tag ##
message("#### Loading Complete ###")
message(paste("Total Usable Data files Count :",  procDatIdx, " total Frames :",procDatFrames));


### Process Files #####
source("HuntingEventAnalysis.r")
for (i in strCondTags)
{
  message(paste("#### Process Group ",i," ###############"))
  subsetDat = srcdatList[[i]];
  strCond   <- paste(strCondR,subsetDat[2],collapse=NULL);
  strSerNum <- "25";
  
  ##Filter To Keep Only data in range ###
  filtereddatAllFrames <- datAllFrames[datAllFrames$group == i & datAllFrames$REyeAngle > -35 & datAllFrames$REyeAngle <  35 & datAllFrames$LEyeAngle > -35 & datAllFrames$LEyeAngle <  35,]
  vlarvaID = unique(filtereddatAllFrames$larvaID)
  #
  
  lHuntStat[[i]] = calcHuntStat(filtereddatAllFrames,vlarvaID)
  #plotGroupMotion(filtereddatAllFrames,lHuntStat[[i]],vlarvaID)
}

## Hunt Statistics Summary - Combine Rows ##
datHuntStat = do.call(rbind,lHuntStat)


#
##Bar Plot Mean Hunting Events Per Animal
par(mar = c(5, 6, 4, 5) + 2.5)

datmean <- unlist(datHuntStat[,"meanHuntingEventsPerLarva"],use.names = FALSE)
datse <- unlist(datHuntStat[,"seHuntingEventsPerLarva"],use.names = FALSE)
datlbls <-row.names(datHuntStat)


plotTop <- max(datmean) +
  datse[datmean== max(datmean)] * 3

barCenters <- barplot(height = datmean,
                      names.arg = datlbls,
                      beside = true, las = 2,
                      ylim = c(0, plotTop),
                      cex.names = 0.75, xaxt = "n",
                      main = "Hunting Events Per Larva",
                      ylab = "#",
                      border = "black", axes = TRUE)

# Specify the groupings. We use srt = 45 for a
# 45 degree string rotation
text(x = barCenters, y =  par("usr")[3] -0.5, srt = 45,
     adj = 1, labels = datlbls, xpd = TRUE)

segments(barCenters, datmean - datse * 2, barCenters,
         datmean + datse * 2, lwd = 1.5)

#####################

par(mar = c(5, 6, 4, 5) + 2.5)

datmean <- unlist(datHuntStat[,"meanHuntRatioPerLarva"],use.names = FALSE)
datse <- unlist(datHuntStat[,"seHuntRatioPerLarva"],use.names = FALSE)
datlbls <-row.names(datHuntStat)


plotTop <- max(datmean) +
  datse[datmean== max(datmean)] * 3

barCenters <- barplot(height = datmean,
                      names.arg = datlbls,
                      beside = true, las = 2,
                      ylim = c(0, plotTop),
                      cex.names = 0.75, xaxt = "n",
                      main = "Ratio of Time spent Hunting Over all Frames",
                      ylab = "#",
                      border = "black", axes = TRUE)

# Specify the groupings. We use srt = 45 for a
# 45 degree string rotation
text(x = barCenters, y = par("usr")[3] - 0.01, srt = 45,
     adj = 1, labels = datlbls, xpd = TRUE)

segments(barCenters, datmean - datse * 2, barCenters,
         datmean + datse * 2, lwd = 1.5)



#######################################################################
###  EYE - PLOT Scatter and Eye Densities #####
#source("EyeScatterAndDensities.r")
#####


#############



# #Calculate EyeVergence Index
# #(Vectors) Are (RightEyeAngle,LeftEyeAngle)
# vRestCentre = c(-10,10);
# vVerg = c(-30,30)-vRestCentre
# 
# Vidx <- function(x) {CosA=((x-vRestCentre)%*%vVerg)/sqrt(sum((x-vRestCentre)^2))*sqrt(sum((vVerg)^2)); return(CosA)}
# 
# ##Obtain Change In Eye Vergence 
# vDEyeVergence = meanf(diff(LeyeAngleAll-ReyeAngleAll),30);
# #Pick + VE Changes 
# huntFiles = fileIdxAll[vDEyeVergence > 50]
# huntFrames = fileNAll[vDEyeVergence > 50]

#temp[[huntFiles]]
#X11();
#plot(ReyeAll,LeyeAll,cex=.1,xlim=c(-50,20),ylim=c(-40,60),col='blue')

#
#hL=hist(ttL[[1]],breaks=seq(-200,200,length=500),xlim=c(-50,100),mfg=c(1,2,3, 3))
#X11();
#hR=hist(ttR[[1]],breaks=seq(-200,200,length=500),xlim=c(-50,100),mfg=c(1,3,3, 3))


#Mean Filtered TrackerData
#meanf <- function(t,k) {tproc=rep(NA,length(t)); for(i in 1:length(t)) tproc[i]=mean(t[max(1,i-k):i]); return(tproc)}
#meanttL=meanf(TrackerData[1]$EyeLDeg,200);
#meanttR=meanf(TrackerData[1]$EyeRDeg,200);
#plot(meanttR,meanttL,cex=.1,xlim=c(-50,20),ylim=c(-40,60));



###HISTOGRAM 2D
##Try HexBin ##
# Create hexbin object and plot
#eyeDm <- data.frame(ttL[[1]],ttR[[1]]) 
#h <- hexbin(eyeDm)
#plot(h, colramp=rf)


##### OPTION 5: The Hard Way (DIY) #######
# http://stackoverflow.com/questions/18089752/r-generate-2d-histogram-from-raw-data
#nbins <- 25
#x.bin <- seq(floor(min(df[,1])), ceiling(max(df[,1])), length=nbins)
#y.bin <- seq(floor(min(df[,2])), ceiling(max(df[,2])), length=nbins)

#freq <-  as.data.frame(table(findInterval(df[,1], x.bin),findInterval(df[,2], y.bin)))
#freq[,1] <- as.numeric(freq[,1])
#freq[,2] <- as.numeric(freq[,2])

#freq2D <- diag(nbins)*0
#freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]

# Normal
#image(x.bin, y.bin, freq2D, col=r)

# Log
#image(x.bin, y.bin, log(freq2D), col=r)

##### Addendum: 2D Histogram + 1D on sides (from Computational ActSci w R) #######
#http://books.google.ca/books?id=YWcLBAAAQBAJ&pg=PA60&lpg=PA60&dq=kde2d+log&source=bl&ots=7AB-RAoMqY&sig=gFaHSoQCoGMXrR9BTaLOdCs198U&hl=en&sa=X&ei=8mQDVPqtMsi4ggSRnILQDw&redir_esc=y#v=onepage&q=kde2d%20log&f=false


