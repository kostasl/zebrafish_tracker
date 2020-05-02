### Library Functions For Plotting Hunt Statistics ###

## Adding Legend outside plotting area by overlaying a new plot layer
## SRC: https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


pieChartLabelledEvents <- function(tblRes,GroupID,colourL=NA)
{
  
  
  ## Find Tbl Indexes Indicating Success 
  tblIdxSuccess <- which (grepl("Success",row.names(tblRes) ) ) 
  tblIdxFail <- which (grepl("Fail",row.names(tblRes) ) ) 
  tblidxValidHuntEvents <- c(6,7,tblIdxSuccess,tblIdxFail)
  #
  ##Summarize COmbine Labels ###
  # Success Together, And Fails Together, A No-Target and and Escape separatelly
  #DLRes=c(sum(tblRes[c(3,12),GroupID]) ,sum(tblRes[c(4,10,11),GroupID]),sum(tblRes[c(5),GroupID]),sum(tblRes[c(7),GroupID]))
  
  DLRes=c(sum(tblRes[tblIdxSuccess,GroupID]) ,sum(tblRes[tblIdxFail,GroupID]),sum(tblRes[c(5),GroupID]),sum(tblRes[c(7),GroupID]))
  
  #NLRes=c(sum(tblRes[c(3,12),"NL"]) ,sum(tblRes[c(4,10,11),"NL"]),sum(tblRes[c(5),"NL"]),sum(tblRes[c(7),"NL"]))
  #LLRes=c(sum(tblRes[c(3,12),"LL"]) ,sum(tblRes[c(4,10,11),"LL"]),sum(tblRes[c(5),"LL"]),sum(tblRes[c(7),"LL"]))
  
  #nLabelledDL <- sum(tblRes[c(3,12,4,10,11,5,7),GroupID])
  #nLabelledLL <- sum(tblRes[c(3,12,4,10,11,5,7),"LL"])
  #nLabelledNL <- sum(tblRes[c(3,12,4,10,11,5,7),"NL"])
  nLabelledDL <- sum(tblRes[tblidxValidHuntEvents,GroupID])
  
  ScoreLabels <- c("Success","Fail","No Target","Escape")
  
  #rfc <- colorRampPalette(rev(brewer.pal(8,'Set2')));
  if (is.na(colourL))
    colourL <- c(rfc(NROW(ScoreLabels)),"#FF0000");
  
  
  pie(DLRes , labels = paste(""," ",round((DLRes/nLabelledDL)*100),"%",sep=""),cex=2.8,cex.main=2.8,clockwise = TRUE,
      main=paste(GroupID," #",nLabelledDL,"/",nLabelledDL+sum(tblRes[c(1),GroupID]) ),
      radius=1.0,col=colourL) 
  #pie(NLRes , labels = paste(ScoreLabels," %",round((NLRes/nLabelledNL)*100)/100,sep=""),clockwise = TRUE,main=paste("NL #",nLabelledNL),radius=1.08)
  #pie(LLRes , labels = paste(ScoreLabels," %",round((LLRes/nLabelledLL)*100)/100,sep=""),clockwise = TRUE,main=paste("LL #",nLabelledLL),radius=1.08)
  
}


##Show Results comparing Success only of Identified Hunt Events that tracked Prey
## Return Number of Labelled Events
pieChartLabelledSuccessVsFails_StrikeBreakDown <- function(tblRes,GroupID,colourL=NA)
{
  ## Find Tbl Indexes Indicating Success 
  tblIdxSuccess <- which (grepl("Success",row.names(tblRes) ) ) 
  tblIdxFail <- which (grepl("Fail",row.names(tblRes) ) ) 
  tblidxValidHuntEvents <- c(6,7,tblIdxSuccess,tblIdxFail)
  tblIdxSuccess_Strike <- which (grepl("Success-OnStrike",row.names(tblRes) ) )
  tblIdxSuccess_NoStrike <- which (grepl("Success-OnApproach",row.names(tblRes) ) ) 
  tblIdxFail_Strike <- which (grepl("Fail-With Strike",row.names(tblRes) ) )
  tblIdxFail_NoStrike <- which (grepl("Fail-No Strike",row.names(tblRes) ) ) 
  
  
  ##Summarize COmbine Labels ###
  # Success Together, And Fails Together
  DLRes=c(sum(tblRes[tblIdxSuccess_Strike,GroupID]) ,
          sum(tblRes[tblIdxSuccess_NoStrike,GroupID] ),
          sum(tblRes[tblIdxFail_Strike,GroupID] ),
          sum(tblRes[tblIdxFail_NoStrike,GroupID] )) 
  nLabelledDL <- sum(DLRes)  ##sum(tblRes[c(3,12,4,10,11,5,7),GroupID])
  
  ScoreLabels <- c("Success","Fail")
  
  if (is.na(colourL))
    colourL <- rev(brewer.pal(4,'Paired'))
  
    #colourL <-  c("#66C2A5FF","#66C2A5A1","#A3A3A3FF","#A3A3A3A1") #c(rfc(NROW(ScoreLabels)),"#FF0000");
  
  pie(DLRes , labels = paste("","",prettyNum(  (DLRes/(nLabelledDL) *100),digits=3   ),"%",sep=""),
      cex=1.4,cex.main=1.4,clockwise = TRUE,
      #main=paste(GroupID," #",nLabelledDL,"/",nLabelledDL+sum(tblRes[1,GroupID]) ),
      radius=1.0,col=colourL,border=F )

  return(nLabelledDL)
}



##Show Results comparing Success only of Identified Hunt Events that tracked Prey
## Return Number of Labelled Events
pieChartLabelledStrikeBreakDown <- function(tblRes,GroupID,colourL=NA)
{
  ## Find Tbl Indexes Indicating Success 
  tblIdxSuccess_Strike <- which (grepl("Success-OnStrike",row.names(tblRes) ) )
  tblIdxSuccess_NoStrike <- which (grepl("Success-OnApproach",row.names(tblRes) ) ) 
  tblIdxFail_Strike <- which (grepl("Fail-With Strike",row.names(tblRes) ) )
  tblIdxFail_NoStrike <- which (grepl("Fail-No Strike",row.names(tblRes) ) ) 

  
  ##Summarize COmbine Labels ###
  # Success Together, And Fails Together
  DLRes=c(sum(tblRes[tblIdxSuccess_Strike,GroupID],tblRes[tblIdxFail_Strike,GroupID]) ,
          sum(tblRes[tblIdxSuccess_NoStrike,GroupID],tblRes[tblIdxFail_NoStrike,GroupID])
          ) 
  nLabelledDL <- sum(DLRes)  ##sum(tblRes[c(3,12,4,10,11,5,7),GroupID])
  
  ScoreLabels <- c("Success","Fail")
  
  rfc <- colorRampPalette(rev(brewer.pal(8,'Set2')));
  if (is.na(colourL))
    colourL <-  c("#66C2A5FF","#66C2A5A1","#A3A3A3FF","#A3A3A3A1") #c(rfc(NROW(ScoreLabels)),"#FF0000");
  
  pie(DLRes , labels = paste("","",prettyNum(  (DLRes/(nLabelledDL) *100),digits=3   ),"%",sep=""),
      cex=1.4,cex.main=1.4,clockwise = TRUE,
      #main=paste(GroupID," #",nLabelledDL,"/",nLabelledDL+sum(tblRes[1,GroupID]) ),
      radius=1.0,col=colourL,border=TRUE )
  
  return(nLabelledDL)
}

##Show Results comparing Success only of Identified Hunt Events that tracked Prey
## Return Number of Labelled Events
pieChartLabelledSuccessVsFails <- function(tblRes,GroupID,colourL=NA)
{
  ## Find Tbl Indexes Indicating Success 
  tblIdxSuccess <- which (grepl("Success",row.names(tblRes) ) ) 
  tblIdxFail <- which (grepl("Fail",row.names(tblRes) ) ) 

  ##Summarize COmbine Labels ###
  # Success Together, And Fails Together
  DLRes=c(sum(tblRes[tblIdxSuccess,GroupID]) ,
          sum(tblRes[tblIdxFail,GroupID] )) 
  nLabelledDL <- sum(DLRes)  ##sum(tblRes[c(3,12,4,10,11,5,7),GroupID])

  ScoreLabels <- c("Success","Fail")
  
  #rfc <- colorRampPalette(rev(brewer.pal(8,'Set2')));
  colPaired <- rev(brewer.pal(4,'Paired')) ## chosen for colorblindness
  if (!is.array(colourL))
     #colourL <-  c("#66C2A5FF","#B3B3B3FF","#B3B3B381") #c(rfc(NROW(ScoreLabels)),"#FF0000");
     colourL <- c(colPaired[1],colPaired[3]) 
  
  pie(DLRes , labels = paste("","",round((DLRes/nLabelledDL)*100),"%",sep=""),
      cex=2.8,cex.main=2.8,clockwise = TRUE,
      #main=paste(GroupID," #",nLabelledDL,"/",nLabelledDL+sum(tblRes[1,GroupID]) ),
      radius=1.0,col=colourL,border=F )
  #pie(NLRes , labels = paste(ScoreLabels," %",round((NLRes/nLabelledNL)*100)/100,sep=""),clockwise = TRUE,main=paste("NL #",nLabelledNL),radius=1.08)
  #pie(LLRes , labels = paste(ScoreLabels," %",round((LLRes/nLabelledLL)*100)/100,sep=""),clockwise = TRUE,main=paste("LL #",nLabelledLL),radius=1.08)
  
  return(nLabelledDL)
}


pieChartStrikeVsNonStrike_Success <- function(tblRes,GroupID,colourL)
{
  tblIdxSuccess <- which ("Success" == row.names(tblRes)  ) ##Generic
  tblIdxSuccess_Strike <- which (grepl("Success-OnStrike",row.names(tblRes) ) )
  tblIdxSuccess_NoStrike <- which (grepl("Success-OnApproach",row.names(tblRes) ) ) 

  tblIdxFail <- which ("Fail" == row.names(tblRes) ) 
  tblIdxFail_Strike <- which (grepl("Fail-With Strike",row.names(tblRes) ) )
  tblIdxFail_NoStrike <- which (grepl("Fail-No Strike",row.names(tblRes) ) ) 
  tblidxValidHuntEvents <- c(6,7,tblIdxSuccess,tblIdxFail)
  
  ##Summarize COmbine Labels ###
  # Success Together, And Fails Together
  DLRes=c(
   # sum(tblRes[tblIdxSuccess,GroupID]),
    sum(tblRes[tblIdxSuccess_Strike,GroupID]),
    sum(tblRes[tblIdxSuccess_NoStrike,GroupID])
    #sum(tblRes[tblIdxFail_Strike,GroupID]),
    #sum(tblRes[tblIdxFail_NoStrike,GroupID])
  )
  
  #NLRes=c(sum(tblRes[c(3,12),"NL"]) ,sum(tblRes[c(4,10,11),"NL"]),sum(tblRes[c(5),"NL"]),sum(tblRes[c(7),"NL"]))
  #LLRes=c(sum(tblRes[c(3,12),"LL"]) ,sum(tblRes[c(4,10,11),"LL"]),sum(tblRes[c(5),"LL"]),sum(tblRes[c(7),"LL"]))
  ##Here We Condition on the Fact that these were Hunt Events Tracking Prey / 
  nLabelledDL <- sum(DLRes)  ##sum(tblRes[c(3,12,4,10,11,5,7),GroupID])
  #nLabelledLL <- sum(tblRes[c(3,12,4,10,11,5,7),"LL"])
  #nLabelledNL <- sum(tblRes[c(3,12,4,10,11,5,7),"NL"])
  
  #ScoreLabels <- c("Success Noclass","Success Strike","Success No Strike","Fail Strike","Fail No Strike") ##"Fail Noclass" Removed as it is 0
  ScoreLabels <- c("Success Strike","Success No-Strike") ##"Fail Noclass" Removed as it is 0
  rfc <- colorRampPalette(rev(brewer.pal(8,'Set2')));
  if (is.na(colourL))
    colourL <-  c("#66C2A5","#B3B3B3") #c(rfc(NROW(ScoreLabels)),"#FF0000");
  
  pie(DLRes , labels = c(paste0("",round((DLRes[1]/nLabelledDL)*100),"%"),
                         paste0("",round((DLRes[2]/nLabelledDL)*100),"%") ),
      cex=1.4,cex.main=2.8,clockwise = TRUE,
      radius=1.0,col=colourL) 
  
  
  
  return(nLabelledDL)
  
}  
  

pieChartStrikeVsNonStrike_Fail <- function(tblRes,GroupID,colourL)
{
  tblIdxFail <- which ("Fail" == row.names(tblRes) ) 
  tblIdxFail_Strike <- which (grepl("Fail-With Strike",row.names(tblRes) ) )
  tblIdxFail_NoStrike <- which (grepl("Fail-No Strike",row.names(tblRes) ) ) 
  tblidxValidHuntEvents <- c(6,7,tblIdxSuccess,tblIdxFail)
  
  ##Summarize COmbine Labels ###
  # Success Together, And Fails Together
  DLRes=c(
    # sum(tblRes[tblIdxSuccess,GroupID]),
    #sum(tblRes[tblIdxSuccess_Strike,GroupID]),
    #sum(tblRes[tblIdxSuccess_NoStrike,GroupID])
    sum(tblRes[tblIdxFail_Strike,GroupID]),
    sum(tblRes[tblIdxFail_NoStrike,GroupID])
  )
  
  nLabelledDL <- sum(DLRes)  ##sum(tblRes[c(3,12,4,10,11,5,7),GroupID])

  #ScoreLabels <- c("Success Noclass","Success Strike","Success No Strike","Fail Strike","Fail No Strike") ##"Fail Noclass" Removed as it is 0
  ScoreLabels <- c("Fail Strike","Fail No-Strike") ##"Fail Noclass" Removed as it is 0
  rfc <- colorRampPalette(rev(brewer.pal(8,'Set2')));
  if (is.na(colourL))
    colourL <-  c("#66C2A5","#B3B3B3") #c(rfc(NROW(ScoreLabels)),"#FF0000");
  
  pie(DLRes , labels = c(paste0("",round((DLRes[1]/nLabelledDL)*100),"%"),
                         paste0("",round((DLRes[2]/nLabelledDL)*100),"%") ),
      cex=1.4,cex.main=2.8,clockwise = TRUE,
      #main=paste(GroupID," #",nLabelledDL,"/",nLabelledDL+sum(tblRes[1,GroupID]) ),
      radius=1.0,col=colourL) 
  
  
  
  return(nLabelledDL)
  
}  


##\todo convert Means to BoxPlots
## Box Plots Used to Compare Conditions On Mean Stats - Saves Output As Pdf
barplotPerCondition <- function(datStat,datMean,datSe,strtitle,strsubt,stroutFileName,plotTop)
{
  
  par(mar = c(5, 6, 4, 5) + 2.5)
  
  datN = vector()
  if (is.null(dim(datStat))) ##This is probably A List
  {
    datlbls <-names(datStat) ##Input is LisFt
    for (i in seq(1:NROW(datStat)))
    {
      datN[i] <-NROW( datStat[[i]][!is.na(datStat[[i]]) ]   )
    }
    # datN <- NROW(datStat[[1]])
  }
  else ##A Known Data Frame Struct
  { 
    datN    <- unlist(datStat[,"nLarva"],use.names = FALSE)
    datlbls <-row.names(datStat) ##Multidim Use Row Names
  }
  ##Add N Numbers to Labels
  datlbls <- paste(datlbls,"\nn=",datN,sep="")
  
  if(missing(plotTop)) 
  {
    plotTop <- max(datMean) +
      unique(datSe[datMean== max(datMean)]) * 3
  }else
  {
    plotTop <- 1.1*plotTop ##Increase by 10%
    
  }
  plotBottom <- 2*min(min(datMean),0)
  
  barCenters <- barplot(height = datMean,
                        names.arg = datlbls,
                        beside = true, las = 2,
                        ylim = c(plotBottom, plotTop),
                        cex.names = 0.75, xaxt = "n",
                        main = strtitle,
                        sub = strsubt,
                        ylab = "#",
                        border = "black", axes = TRUE)
  
  # Specify the groupings. We use srt = 45 for a
  # 45 degree string rotation
  #text(x = barCenters, y = par("usr")[3] - 0.01, srt = 45,    adj = 1, labels = datlbls, xpd = TRUE)
  text(x = barCenters+0.2, y = 0, srt = 45,
       adj = 1.8, labels = datlbls, xpd = TRUE)
  
  segments(barCenters, datMean - datSe * 2, barCenters,
           datMean + datSe * 2, lwd = 1.5)
  #title(sub=strsubt)
  #  dev.off()
  
  message(strsubt)
  
  return(barCenters)
}




plotConnectedPointsPairs <- function(vIDTable,vDat,strCondTags,xbarcenters)
{
  for (gIdx in seq(1,NROW(strCondTags),2)  ) ##Iterated Through LF DF And NF Groups
  {
    gE <- strCondTags[gIdx] ##Empty Condution
    gL <- strCondTags[gIdx+1] ##With ROtifers Test Condition 
    
    #stopifnot(NROW(vIDTable[[gE]]$dataSetID) == NROW(vIDTable[[gL]]$dataSetID) )
    
    datSetID <- levels(vIDTable[[gE]]$dataSetID)[ vIDTable[[gE]]$dataSetID[vIDTable[[gE]]$larvaID == vIDTable[[gL]]$larvaID & vIDTable[[gE]]$dataSetID == vIDTable[[gL]]$dataSetID] ]
    
    idsE <- vIDTable[[gE]]$expID[vIDTable[[gE]]$larvaID == vIDTable[[gL]]$larvaID & vIDTable[[gE]]$dataSetID ==vIDTable[[gL]]$dataSetID]
    idsL <- vIDTable[[gL]]$expID[vIDTable[[gL]]$larvaID == vIDTable[[gL]]$larvaID & vIDTable[[gE]]$dataSetID == vIDTable[[gL]]$dataSetID]
    ptSrc  <- vDat[[gE]][levels(idsE)[idsE]]
    ptDest <- vDat[[gL]][levels(idsL)[idsL]]
    
    ##OPTIONAL: Replace NAs With 0 when Plotting So as to To Show Direction For All points
    #ptDest <- replace(ptDest,is.na(ptDest),0)
    ##Otherwise Do not Plot NA connecting line
    
    ##Plot The Lines Connect Each Empty Tested Larva With Itself In THe Live Fed Conditions 
    idxSrc  <- match(gE,strCondTags) ##Bar Center Idx for Each Condition E. Fed
    idxDest <- match(gL ,strCondTags)           
    
    for (p in 1:NROW(ptSrc))
    {
      pcolour <- rDataset[as.numeric( replace( datSetID[p],is.na(datSetID[p]) ,1 ) ) ]
      points(xbarcenters[idxSrc],ptSrc[p],pch=as.numeric(datSetID[p]),
             col=pcolour )
      points(xbarcenters[idxDest],ptDest[p],pch=as.numeric(datSetID[p]),
             col=pcolour )
      
      ccLine <- rbind(c(xbarcenters[idxSrc],ptSrc[p]),c(xbarcenters[idxDest],ptDest[p] ) ) 
      
      lines(ccLine,col=pcolour)
    }
    
  } ## Go Through Pairs Of Conditions ##
}


plotPairedChangeHistogram <- function(vIDTable,vDat,strCondTags,uLim,lLim)
{
  lChangePairs <- list()
  idx = 1
  for (gIdx in seq(1,NROW(strCondTags),2)  ) ##Iterated Through LF DF And NF Groups
  {
    gE <- strCondTags[gIdx] ##Empty Condution
    gL <- strCondTags[gIdx+1] ##With ROtifers Test Condition 
    
    datSetID <- levels(vIDTable[[gE]]$dataSetID) [ vIDTable[[gE]]$dataSetID[vIDTable[[gE]]$larvaID == vIDTable[[gL]]$larvaID & vIDTable[[gE]]$dataSetID == vIDTable[[gL]]$dataSetID] ]
    idsE <- vIDTable[[gE]]$expID[vIDTable[[gE]]$larvaID == vIDTable[[gL]]$larvaID && vIDTable[[gE]]$dataSetID ==vIDTable[[gL]]$dataSetID]
    idsL <- vIDTable[[gL]]$expID[vIDTable[[gL]]$larvaID == vIDTable[[gL]]$larvaID && vIDTable[[gE]]$dataSetID == vIDTable[[gL]]$dataSetID]
    ptSrc <-  vDat[[gE]][levels(idsE)[idsE]]
    ptDest <- vDat[[gL]][levels(idsL)[idsL]]
    
    ##OPTIONAL: Replace NAs With 0 when Plotting So as to To Show Direction For All points
    #ptSrc <- replace(ptSrc,is.na(ptSrc),0)
    #ptDest <- replace(ptDest,is.na(ptDest),0)
    ##Otherwise Do not Plot NA connecting line
    
    ##Plot The Lines Connect Each Empty Tested Larva With Itself In THe Live Fed Conditions 
    idxSrc  <- match(gE,strCondTags) ##Bar Center Idx for Each Condition E. Fed
    idxDest <- match(gL ,strCondTags)           
    
    ##PLot Histogram Of THis Pair
    ptDelta <- ptDest-ptSrc
    ##Saturate Limits 
    ptDelta[ptDelta > uLim] = uLim
    ptDelta[ptDelta < lLim] = lLim
    
    lChangePairs[[idx]] <- ptDelta
    idx = idx+1
    hist(ptDelta,ylim = c(0,NROW(ptDelta)/3),breaks = seq(lLim,uLim,(uLim-lLim)/30),
         main =paste("                                                \t \t Change ",gE,"-",gL,sep="") )
    
  } ## Go Through Pairs Of Conditions ##
  
  names(lChangePairs) <- c("LE-LL","NE-NL","DE-DL")
  
  return (lChangePairs)
}


## Show In Which Prey Count Density Most Hunt Events Occurred - 
## Note: NOT INITIAL #Prey count but rather prey count sample
plotHuntEventPreyCountHist <- function(strCondTags,dataSetsToProcess)
{
  layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE))
  
  for (i in strCondTags)
  {
    strDataFileName <- paste("out/setn",NROW(dataSetsToProcess),"HuntEvents",i,sep="-") ##To Which To Save After Loading
    message(paste(" Loading Hunt Events: ",strDataFileName))
    ##ExPORT 
    load(file=paste(strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
    hist(unlist(datHuntEvent[,"PreyCount"],use.names = FALSE),main=paste(i," Events #",NROW(datHuntEvent[,"PreyCount"]),sub = ""),
         breaks=seq(0,85,5),xlab="# Prey",xlim=c(0,70),ylim=c(0,300) )
  }
  
}



## Show the mean Hunt Rate of a group's Larva vs Prey Count Density - 
## Note: NOT INITIAL #Prey count but rather prey count sample on each Hunt event is used here and then this is divided by 
##' the number of Larvae that did these events in that bin
plotMeanHuntEventPerLarvaVsPreyCountHist <- function(datAllHuntEvent)
{
  
  layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE))
  yl <- 35 
  step <- 10
  strCondTags <- unique(datAllHuntEvent$groupID)
  for (g in strCondTags)
  {
    #strDataFileName <- paste("out/setn",NROW(dataSetsToProcess),"HuntEvents",g,sep="-") ##To Which To Save After Loading
    #message(paste(" Loading Hunt Events: ",strDataFileName))
    ##ExPORT 
    #load(file=paste(strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
    datHuntEventFilt<- datAllHuntEvent[datAllHuntEvent$groupID == g,]
    datHuntEventFilt$expID <- factor(datHuntEventFilt$expID)
    
    histHuntFq <- vector()
    histHuntFqSE <- vector()
    histHuntPrey <- vector()
    histHuntN  <- vector()
    
    nn <- 0
    for (i in seq(0,70,step))
    {
      nn <- nn + 1 
      datExpIdSlice <- datHuntEventFilt$expID[round(datHuntEventFilt$PreyCount) >= i & round(datHuntEventFilt$PreyCount) < i+step]
      datHuntEventSlice <- datHuntEventFilt[round(datHuntEventFilt$PreyCount) >= i & round(datHuntEventFilt$PreyCount) < i+step,]
      
      ##Summarize Number of Hunt Events Per Larvae 
      tblHuntEventCount <- table(datHuntEventSlice$expID)
      
      histHuntPrey[nn]  <- i
      #histHuntFq[nn] <- ifelse(NROW(unique(datExpIdSlice)) > 0,NROW(datExpIdSlice)/NROW(unique(datExpIdSlice)),0) 
      ##There is 1 Event 
      histHuntFq[nn] <- mean(tblHuntEventCount[tblHuntEventCount > 0],na.rm = TRUE)
      histHuntFqSE[nn] <- sd(tblHuntEventCount[tblHuntEventCount > 0],na.rm = TRUE)/sqrt(length(datHuntEventSlice))
      histHuntN[nn] <- NROW(unique(datExpIdSlice))
      
      
    }
    #plot(histHuntPrey,histHuntFq,type='l',xlab="#Prey",ylab="#Hunts/#Larva",main=paste(g,"  ") )
    barCenters <- barplot(histHuntFq,names.arg="",xlab="#Prey",ylab="#Hunts/#Larva",main=paste(g," #",NROW(datHuntEventFilt)),ylim=c(0,yl) )
    
    
    segments(barCenters, histHuntFq - histHuntFqSE * 2, barCenters,
             histHuntFq + histHuntFqSE * 2, lwd = 1.5)
    
    text(x = barCenters-0.1, y = -0.9, srt = 0,
         adj = 1.8, labels = paste(histHuntPrey,sep=""), xpd = TRUE)
    
    text(x = barCenters+0.7, y = -0.05, srt = 45,
         adj = 1.8, labels = paste("\nn=",histHuntN,sep=""), xpd = TRUE)
  }
}



## Show the mean Hunt INTERVALS of a group's Larva vs Prey Count Density - 
## Note: NOT INITIAL #Prey count but rather prey count sample on each Hunt event is used here and then this is divided by 
##' the number of Larvae that did these events in that bin
plotMeanHuntIntervalPerLarvaVsPreyCountHist <- function(datAllHuntEvent)
{
  
  layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE))
  yl <- 10000/G_APPROXFPS
  step <- 20
  
  strCondTags <- unique(datAllHuntEvent$groupID)
  for (g in strCondTags)
  {
    # strDataFileName <- paste("out/setn",NROW(dataSetsToProcess),"HuntEvents",g,sep="-") ##To Which To Save After Loading
    #  message(paste(" Loading Hunt Events: ",strDataFileName))
    ##ExPORT 
    #  load(file=paste(strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
    
    datHuntEventFilt<- datAllHuntEvent[datAllHuntEvent$groupID == g,]
    datHuntEventFilt$expID <- factor(datHuntEventFilt$expID)
    
    
    
    histHuntInterval    <- vector()
    histHuntIntervalSE  <- vector()
    histHuntPrey        <- vector()
    histHuntN           <- vector()
    
    nn <- 0
    for (i in seq(0,70,step))
    {
      nn <- nn + 1 
      datHuntEventSlice <- datHuntEventFilt[round(datHuntEventFilt$PreyCount) >= i & round(datHuntEventFilt$PreyCount) < i+step,]
      #datExpIdSlice <- datHuntEvent$expID[round(datHuntEvent$PreyCount) >= i & round(datHuntEvent$PreyCount) < i+step]
      
      histHuntPrey[nn]  <- i
      #histHuntFq[nn] <- ifelse(NROW(unique(datExpIdSlice)) > 0,NROW(datExpIdSlice)/NROW(unique(datExpIdSlice)),0)
      ##Get Mean Hunt Intervals 
      tblMeanHuntIntervals   <- tapply(datHuntEventSlice$nextHuntFrame-datHuntEventSlice$endFrame, datHuntEventSlice$expID,mean,na.rm=TRUE)
      histHuntInterval[nn]   <- mean(tblMeanHuntIntervals,na.rm=TRUE)/G_APPROXFPS ##Normalize to Seconds
      histHuntIntervalSE[nn] <- (sd(tblMeanHuntIntervals,na.rm=TRUE)/(sqrt(NROW(datHuntEventSlice))) ) / G_APPROXFPS ##Normalize to Seconds
      histHuntN[nn]          <- NROW(unique(datHuntEventSlice$expID)) ##Number of Larvae this belongs to 
      
    }
    #yl <- max(histHuntInterval,na.rm=TRUE)
    #plot(histHuntPrey,histHuntFq,type='l',xlab="#Prey",ylab="#Hunts/#Larva",main=paste(g,"  ") )
    barCenters <- barplot(histHuntInterval,names.arg="",xlab="#Prey",ylab="#Mean Hunts Interval (sec)",main=paste(g," #",NROW(datHuntEventFilt)),ylim=c(0,yl) )
    
    segments(barCenters, histHuntInterval - histHuntIntervalSE * 2, barCenters,
             histHuntInterval + histHuntIntervalSE * 2, lwd = 1.5)
    
    text(x = barCenters-0.1, y = -0.9, srt = 0,
         adj = 1.8, labels = paste(histHuntPrey,sep=""), xpd = TRUE)
    
    text(x = barCenters+0.5, y = -0.04, srt = 45,
         adj = 1.8, labels = paste("\nn=",histHuntN,sep=""), xpd = TRUE) ##Number of Larvae Involved
  }
}



## Show the mean Hunt INTERVALS of a group's Larva vs Prey Count Density -
## Uses the Next HuntFrame data - Which is managed during import processes - and not corrected by labelling
## Note: NOT INITIAL #Prey count but rather prey count sample on each Hunt event is used here and then this is divided by 
##' the number of Larvae that did these events in that bin
plotInterHuntIntervalPerLarvaVsPreyCountHist <- function(datAllHuntEvent)
{
  
  layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE))
  yl <- 10000/G_APPROXFPS
  step <- 20
  
  strCondTags <- unique(datAllHuntEvent$groupID)
  for (g in strCondTags)
  {
    # strDataFileName <- paste("out/setn",NROW(dataSetsToProcess),"HuntEvents",g,sep="-") ##To Which To Save After Loading
    #  message(paste(" Loading Hunt Events: ",strDataFileName))
    ##ExPORT 
    #  load(file=paste(strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
    
    datHuntEventFilt<- datAllHuntEvent[datAllHuntEvent$groupID == g,]
    datHuntEventFilt$expID <- factor(datHuntEventFilt$expID)
    
    histHuntInterval    <- list()
    histHuntIntervalSE  <- vector()
    histHuntPrey        <- vector()
    histHuntN           <- vector()
    
    nn <- 0
    preybreaks <- seq(0,70,step)
    for (i in preybreaks)
    {
      nn <- nn + 1 
      datHuntEventSlice <- datHuntEventFilt[round(datHuntEventFilt$PreyCount) >= i & round(datHuntEventFilt$PreyCount) < i+step,]
      #datExpIdSlice <- datHuntEvent$expID[round(datHuntEvent$PreyCount) >= i & round(datHuntEvent$PreyCount) < i+step]
      
      histHuntPrey[nn]  <- i
      #histHuntFq[nn] <- ifelse(NROW(unique(datExpIdSlice)) > 0,NROW(datExpIdSlice)/NROW(unique(datExpIdSlice)),0)
      ##Get Mean Hunt Intervals - 
      tblMeanHuntIntervals   <- tapply(datHuntEventSlice$nextHuntFrame-datHuntEventSlice$endFrame, datHuntEventSlice$expID,mean,na.rm=TRUE)
      histHuntInterval[[nn]]   <- tblMeanHuntIntervals[!is.na(tblMeanHuntIntervals)]/G_APPROXFPS ##Normalize to Seconds
      histHuntN[nn]          <- NROW(unique(datHuntEventSlice$expID)) ##Number of Larvae this belongs to 
      
    }
    #yl <- max(histHuntInterval,na.rm=TRUE)
    #plot(histHuntPrey,histHuntFq,type='l',xlab="#Prey",ylab="#Hunts/#Larva",main=paste(g,"  ") )
    lbl <- paste(preybreaks,"\n#",histHuntN)
    barCenters <- boxplot(histHuntInterval,names=lbl ,xlab="#Prey",ylab="#Mean IHI (sec)",main=paste(g," #",NROW(datHuntEventFilt)),ylim=c(0,yl) )
    
    #segments(barCenters, histHuntInterval - histHuntIntervalSE * 2, barCenters, histHuntInterval + histHuntIntervalSE * 2, lwd = 1.5)
    
    # text(x = barCenters-0.1, y = -0.9, srt = 0,
    #      adj = 1.8, labels = paste(histHuntPrey,sep=""), xpd = TRUE)
    
    #text(x = preybreaks+0.5, y = -0.04, srt = 45, adj = 1.8, labels = paste("\nn=",histHuntN,sep=""), xpd = TRUE) ##Number of Larvae Involved
  }
}


boxPlotHuntEpisodeDuration <- function(datAllHuntEvent)
{
  
  layout(matrix(c(1,2,3,3,4,5), 3, 2, byrow = TRUE))
  
  ##tblHuntEpisodeDuration <- tapply(datAllHuntEvent$nextHuntFrame-datAllHuntEvent$endFrame,datAllHuntEvent$groupID,mean,na.rm=TRUE)
  datAllHuntEvent_local <- datAllHuntEvent
  datAllHuntEvent_local$huntScore <- convertToScoreLabel(datAllHuntEvent_local$huntScore)
  
  datAllHuntEventSucc <- datAllHuntEvent_local[datAllHuntEvent_local$huntScore == "Success" |
                                                 datAllHuntEvent_local$huntScore == "Success-SpitBackOut",]
  ##Merge Success
  datAllHuntEventSucc[datAllHuntEventSucc$huntScore == "Success-SpitBackOut",]$huntScore <- "Success" 
  
  boxplot((datAllHuntEventSucc$nextHuntFrame-datAllHuntEventSucc$endFrame)/G_APPROXFPS ~  datAllHuntEventSucc$groupID,
          main=paste("Successful Episode T Per Group"," #",NROW(datAllHuntEventSucc)), ylab="(sec)",ylim=c(0,40))
  
  datAllHuntEventFail <- datAllHuntEvent_local[datAllHuntEvent_local$huntScore == "Fail-With Strike" |
                                                 datAllHuntEvent_local$huntScore == "Fail-No Strike" |
                                                 datAllHuntEvent_local$huntScore == "Fail" |
                                                 datAllHuntEvent_local$huntScore == "No_Target",]
  
  ##Merge Fail
  datAllHuntEventFail[datAllHuntEventFail$huntScore == "Fail",]$huntScore <- "Fail-No Strike" 
  
  
  boxplot((datAllHuntEventFail$nextHuntFrame-datAllHuntEventFail$endFrame)/G_APPROXFPS ~  datAllHuntEventFail$groupID,
          main=paste("Failed Episode T Per Group"," #",NROW(datAllHuntEventFail)),ylab="(sec)",ylim=c(0,40))
  
  
  ##Per Label  
  boxplot((datAllHuntEvent_local$nextHuntFrame-datAllHuntEvent_local$endFrame)/G_APPROXFPS ~  datAllHuntEvent_local$huntScore,
          main=paste("T Per Label"," #",NROW(datAllHuntEvent_local)),ylab="(sec)",ylim=c(0,40))
  
  
  lSuccessVsFail <- list()
  
  lSuccessVsFail[["Success"]] <- datAllHuntEventSucc
  lSuccessVsFail[["Fail"]] <- datAllHuntEventFail
  datSuccessVsFail <- do.call(rbind,lSuccessVsFail)
  
  
  ##Merge Relevant Scores
  
  
  datSuccessVsFail$huntScore <- factor(datSuccessVsFail$huntScore) #factor(x=datSuccessVsFail$huntScore,levels=c(0,2,12,4,3,9,10),labels=c("UnLabelled","Success","Success/SpitOut","No_Target","Fail","Fail-No Strike","Fail-With Strike") ) 
  
  
  ##Per Success Fail Label  
  boxplot((datSuccessVsFail$nextHuntFrame-datSuccessVsFail$endFrame)/G_APPROXFPS ~  datSuccessVsFail$huntScore,
          main=paste("Episode T Per Outcome"," #",NROW(datSuccessVsFail)),ylab="(sec)",ylim=c(0,40))
  
  
  
  ##Success Per Group
  boxplot((datAllHuntEventSucc$nextHuntFrame-datAllHuntEventSucc$endFrame)/G_APPROXFPS ~  datAllHuntEventSucc$groupID,
          main=paste("Successful Hunt T Per Group"," #",NROW(datAllHuntEventSucc)),ylab="(sec)",ylim=c(0,40))
  
  
}




## The figure with CDF of hunt efficiency - plotted as power, so as to show differences in actual consumption 
## instead of just ratio of success vs fail
plotHuntPowerDataCDF <- function(datHuntEventAllGroupToLabel)
{
  datFishSuccessRate <- getHuntSuccessPerFish(datHuntEventAllGroupToLabel)
  ## Hunt POwer ##
  vScoreIdx        <- datFishSuccessRate[,"HuntPower"] # ((datFishSuccessRate[,"Success"]*datFishSuccessRate[,"Success"])/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails"]))
  vEfficiencyRatio <- (datFishSuccessRate[,"Success"]/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails"]))
  vEfficiencyRatio_Strike <- (datFishSuccessRate[,"Success"]/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails_WS"]))
  vEfficiencyRatio_NStrike <- (datFishSuccessRate[,"Success"]/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails_NS"]))
  #vScoreIdx[is.nan(vScoreIdx) ] <- 0
  tblEventsTracked <- table(datHuntEventAllGroupToLabel$expID, datHuntEventAllGroupToLabel$markTracked,useNA="always" )
  datFishSuccessRateMerged <- cbind(datFishSuccessRate,
                                    markUnTrackable=data.frame(tblEventsTracked[row.names(datFishSuccessRate),1]),
                                    markTracked=data.frame(tblEventsTracked[row.names(datFishSuccessRate),2]),
                                    notTracked=data.frame(tblEventsTracked[row.names(datFishSuccessRate),3]))
  
  datFishSuccessRateMerged <- cbind(datFishSuccessRateMerged,vScoreIdx,vEfficiencyRatio,vEfficiencyRatio_Strike,vEfficiencyRatio_NStrike)
  
  
  
  ## Subset only the active/Larvae - ones that have hunted 
  datFishSuccessRateActive <- datFishSuccessRateMerged[!is.nan(datFishSuccessRateMerged$vScoreIdx),]
  
  ## Plot Density of Hunting POWER S^2/(S+F)
  densDLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vScoreIdx)
  densNLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vScoreIdx)
  densLLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vScoreIdx)
  
  cdfDLScore <- ecdf(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vScoreIdx)
  cdfNLScore <- ecdf(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vScoreIdx)
  cdfLLScore <- ecdf(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vScoreIdx)
  
  
  #par(mar = c(3.9,4.3,1,1)) 
  plot(cdfNLScore,lty=2,lwd=3,col=colourLegL [1],xlim=c(0,12),pch=pchL[1],ylim=c(0.03,1.01),
       main=NA,ylab=NA,  xlab=NA,cex.main =cex,cex.axis=cex,cex=cex )
  
  plot(cdfLLScore,add=T,lty=1,lwd=3,col=colourLegL[2],pch=pchL[2],ylim=c(0,1.01),cex=cex)
  plot(cdfDLScore,add=T,lty=1,lwd=3,col=colourLegL[3],pch=pchL[3],ylim=c(0,1.01),cex=cex)
  #axis(side=1)
  mtext(side = 1,cex=cex, line = line, expression( "Hunt power " ~ N[S]^2/(N[S]+N[F]) ,paste("") )   )
  mtext(side = 2,cex=cex, line = line, expression("Cumulative function " ))
  
  legend("bottomright",legend=paste(c("NL #","LL #","DL #"),c(densNLScore$n,densLLScore$n,densDLScore$n) ),
         col = colourLegL,pch=pchL,cex=cex+0.2)
}




## The figure with CDF of hunt efficiency - plotted as power, so as to show differences in actual consumption 
## instead of just ratio of success vs fail
plotHuntEfficiencyDataCDF <- function(datHuntEventAllGroupToLabel,showLegend=FALSE)
{
  datFishSuccessRate <- getHuntSuccessPerFish(datHuntEventAllGroupToLabel)
  vScoreIdx        <- ((datFishSuccessRate[,"Success"]*datFishSuccessRate[,"Success"])/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails"]))
  vEfficiencyRatio <- (datFishSuccessRate[,"Success"]/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails"]))
  vEfficiencyRatio_Strike <- (datFishSuccessRate[,"Success"]/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails_WS"]))
  vEfficiencyRatio_NStrike <- (datFishSuccessRate[,"Success"]/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails_NS"]))
  #vScoreIdx[is.nan(vScoreIdx) ] <- 0
  tblEventsTracked <- table(datHuntEventAllGroupToLabel$expID, datHuntEventAllGroupToLabel$markTracked,useNA="always" )
  datFishSuccessRateMerged <- cbind(datFishSuccessRate,
                                    markUnTrackable=data.frame(tblEventsTracked[row.names(datFishSuccessRate),1]),
                                    markTracked=data.frame(tblEventsTracked[row.names(datFishSuccessRate),2]),
                                    notTracked=data.frame(tblEventsTracked[row.names(datFishSuccessRate),3]))
  
  datFishSuccessRateMerged <- cbind(datFishSuccessRateMerged,vScoreIdx,vEfficiencyRatio,vEfficiencyRatio_Strike,vEfficiencyRatio_NStrike)
  
  
  
  ## Subset only the active/Larvae - ones that have hunted 
  datFishSuccessRateActive <- datFishSuccessRateMerged[!is.nan(datFishSuccessRateMerged$vEfficiencyRatio),]
  
  ## Plot Density of Hunting POWER S^2/(S+F)
  densDLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vEfficiencyRatio)
  densNLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vEfficiencyRatio)
  densLLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vEfficiencyRatio)
  
  cdfDLScore <- ecdf(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vEfficiencyRatio)
  cdfNLScore <- ecdf(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vEfficiencyRatio)
  cdfLLScore <- ecdf(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vEfficiencyRatio)
  
  
  message("LL:")
  print(quantile(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vEfficiencyRatio))
  
  message("NL:")
  print(quantile(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vEfficiencyRatio))
  message("DL:")
  print(quantile(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vEfficiencyRatio) )
  
  #par(mar = c(3.9,4.3,1,1)) 
  plot(cdfNLScore,lty=2,lwd=3,col=colourLegL [1],xlim=c(0,1),pch=pchL[1],ylim=c(0.03,1.01),
       main=NA,ylab=NA,  xlab=NA,cex.main =cex,cex.axis=cex,cex=cex )
  
  plot(cdfLLScore,add=T,lty=1,lwd=3,col=colourLegL[2],pch=pchL[2],ylim=c(0,1.01),cex=cex)
  plot(cdfDLScore,add=T,lty=1,lwd=3,col=colourLegL[3],pch=pchL[3],ylim=c(0,1.01),cex=cex)
  #axis(side=1)
  mtext(side = 1,cex=cex, line = line, expression( "Capture Efficiency  "  ,paste("") )   )
  mtext(side = 2,cex=cex, line = line, expression("Cumulative function " ))
  if (showLegend) 
  {
    legend("bottomright",legend=paste(c("NL #","LL #","DL #"),c(densNLScore$n,densLLScore$n,densDLScore$n) ),
           col = colourLegL,pch=pchL,cex=cex+0.2)
  }
}




## The figure with CDF of hunt efficiency - plotted as power, so as to show differences in actual consumption 
## instead of just ratio of success vs fail
plotHuntEfficiencyDataPDF <- function(datHuntEventAllGroupToLabel)
{
  datFishSuccessRate <- getHuntSuccessPerFish(datHuntEventAllGroupToLabel)
  vScoreIdx        <- ((datFishSuccessRate[,"Success"]*datFishSuccessRate[,"Success"])/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails"]))
  vEfficiencyRatio <- (datFishSuccessRate[,"Success"]/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails"]))
  vEfficiencyRatio_Strike <- (datFishSuccessRate[,"Success"]/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails_WS"]))
  vEfficiencyRatio_NStrike <- (datFishSuccessRate[,"Success"]/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails_NS"]))
  #vScoreIdx[is.nan(vScoreIdx) ] <- 0
  tblEventsTracked <- table(datHuntEventAllGroupToLabel$expID, datHuntEventAllGroupToLabel$markTracked,useNA="always" )
  datFishSuccessRateMerged <- cbind(datFishSuccessRate,
                                    markUnTrackable=data.frame(tblEventsTracked[row.names(datFishSuccessRate),1]),
                                    markTracked=data.frame(tblEventsTracked[row.names(datFishSuccessRate),2]),
                                    notTracked=data.frame(tblEventsTracked[row.names(datFishSuccessRate),3]))
  
  datFishSuccessRateMerged <- cbind(datFishSuccessRateMerged,vScoreIdx,vEfficiencyRatio,vEfficiencyRatio_Strike,vEfficiencyRatio_NStrike)
  
  
  
  ## Subset only the active/Larvae - ones that have hunted 
  datFishSuccessRateActive <- datFishSuccessRateMerged[!is.nan(datFishSuccessRateMerged$vEfficiencyRatio),]
  
  ## Plot Density of Hunting POWER S^2/(S+F)
  densDLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vEfficiencyRatio)
  densNLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vEfficiencyRatio)
  densLLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vEfficiencyRatio)
  
  #par(mar = c(3.9,4.3,1,1)) 
  plot(densNLScore,lty=2,lwd=3,col=colourLegL [1],xlim=c(0,1),pch=pchL[1],ylim=c(0.01,3.01),
       main=NA,ylab=NA,  xlab=NA,cex.main =cex,cex.axis=cex,cex=cex )
  
  lines(densLLScore,lty=1,lwd=3,col=colourLegL[2],pch=pchL[2],ylim=c(0,3.01),cex=cex)
  lines(densDLScore,lty=1,lwd=3,col=colourLegL[3],pch=pchL[3],ylim=c(0,3.01),cex=cex)
  #axis(side=1)
  mtext(side = 1,cex=1.5, line = 3.5, expression( "Capture Efficiency  "  ,paste("") )   )
  mtext(side = 2,cex=1.5, line = 3.5, expression("Cumulative function " ))
  
  legend("topright",legend=paste(c("NL #","LL #","DL #"),c(densNLScore$n,densLLScore$n,densDLScore$n) ),
         col = colourLegL,pch=pchL,cex=cex+0.2)
}



