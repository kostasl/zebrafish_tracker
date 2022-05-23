

#######################################################################
###  For Each LARVA in the Group - PLOT Scatter and Eye Densities #####
#if (!exists(filtereddatAllFrames))
#{
#  error("EyeScatteAndDensities Expects filtereddatAllFrames data frame")
#}
# [[Hello Palette]]
library("MASS")
source("TrackerDataFilesImport_lib.r") ##For clip Eye Range

rfHot <- colorRampPalette(rev(brewer.pal(11,'Spectral')));

histj<- function(x,y,x.breaks,y.breaks){
  c1 = as.numeric(cut(x,breaks=x.breaks));
  c2 = as.numeric(cut(y,breaks=y.breaks));
  mat<-matrix(0,ncol=length(y.breaks)-1,nrow=length(x.breaks)-1);
  mat[cbind(c1,c2)] = 1;
  return(mat)
}  


hbinRL = list();
idx = 1;

TotalprocDatFrames <- 0 ##Running Sum of all frames
for (i in vexpID)
{
  print(i)  
  
  datLarvalAllFrames <- datAllGroupFrames[datAllGroupFrames$expID == i 
                                          &datAllGroupFrames$REyeAngle > -G_THRESHCLIPEYEDATA &
                                          datAllGroupFrames$REyeAngle <  G_THRESHCLIPEYEDATA & 
                                          datAllGroupFrames$LEyeAngle > -G_THRESHCLIPEYEDATA &
                                          datAllGroupFrames$LEyeAngle <  G_THRESHCLIPEYEDATA
                                            ,]
  
  procDatFrames <- NROW(datLarvalAllFrames$frameN)
  if (procDatFrames < 10)
    next
  
  lMax <- G_THRESHCLIPEYEDATA ##Min Max Angle Allowed For Each Eye
  lMin <- -G_THRESHCLIPEYEDATA
  
  datLarvalAllFrames$LEyeAngle <-medianf(datLarvalAllFrames$LEyeAngle,nEyeFilterWidth)
  datLarvalAllFrames$LEyeAngle <- clipEyeRange(datLarvalAllFrames$LEyeAngle,lMin,lMax)
  datLarvalAllFrames[is.na(datLarvalAllFrames$LEyeAngle),"LEyeAngle"] <- G_THRESHCLIPEYEDATA
  
  lMax <- G_THRESHCLIPEYEDATA
  lMin <- -G_THRESHCLIPEYEDATA
  

  datLarvalAllFrames$REyeAngle <-medianf(datLarvalAllFrames$REyeAngle,nEyeFilterWidth)
  datLarvalAllFrames$REyeAngle <- clipEyeRange(datLarvalAllFrames$REyeAngle,lMin,lMax)
  datLarvalAllFrames[is.na(datLarvalAllFrames$REyeAngle),"REyeAngle"] <- G_THRESHCLIPEYEDATA
  
  strScatterplotFileName <- paste(strPlotExportPath,"/scatter/EyeAngleScatter-Set-",strCond,"-lID_",i,".pdf",collapse=NULL,sep="");
  strDensityplotFileName <- paste(strPlotExportPath,"/densities/EyeAngleDensity-Set-",strCond,"-lID_",i,".pdf",collapse=NULL,sep="");
  
  ## Eye Trajectory Scatter Plot For all events from This Larva ##
#  pdf(strScatterplotFileName,width=8,height=8)
    sampleSize= length(unique(datLarvalAllFrames$fileIdx));
 
    TotalprocDatFrames = TotalprocDatFrames+procDatFrames
  
    plot(datLarvalAllFrames$REyeAngle,datLarvalAllFrames$LEyeAngle,cex=.1,xlim=c(-50,20),ylim=c(-40,60),asp=1,
         xlab="Right Eye Angle",ylab="Left Eye Angle")
    title(paste(strCond,"R-L Eye Density lID=",i," #n=", sampleSize, " T:",round(procDatFrames/G_APPROXFPS),"sec"),collapse=NULL);
#  dev.off();
  
  hR <- hist(datLarvalAllFrames$REyeAngle, breaks=seq(-G_THRESHCLIPEYEDATA-1,G_THRESHCLIPEYEDATA+1,length=60), plot=F)
  hL <- hist(datLarvalAllFrames$LEyeAngle, breaks=seq(-G_THRESHCLIPEYEDATA-1,G_THRESHCLIPEYEDATA+1,length=60), plot=F)
  ##Do Binarized Histogram Add Each LArva In the group To A list of matrices
  hbinRL[[idx]] <- histj(datLarvalAllFrames$REyeAngle,datLarvalAllFrames$LEyeAngle,(-G_THRESHCLIPEYEDATA:G_THRESHCLIPEYEDATA),(-G_THRESHCLIPEYEDATA:G_THRESHCLIPEYEDATA))
  
  
  # ## DO Eye Density Plot for all events from This Larva ##  
  # if (length(datLarvalAllFrames$REyeAngle) > 10)
  # {
  #   
  #   top <- max(hR$counts, hL$counts)
  #   bw <- bandwidth.nrd(datLarvalAllFrames$REyeAngle)
  #   
  #   bw <- ifelse(is.na(bw),0,bw)
  #   message(paste("kde BWdth:",bw));
  #   if (bw==0)
  #   {    bw <- 1.2
  #   message(paste("**Warning changed kde BWdth to fixed value -> ",bw));
  #   }
  #   
  #   k <- kde2d(datLarvalAllFrames$REyeAngle,datLarvalAllFrames$LEyeAngle,h=bw, n=60, lims=c(range(-35,35),range(-35,35)) )
  #   pdf(strDensityplotFileName,width=8,height=8)
  #   # margins
  #   oldpar <- par()
  #   par(mar=c(3,3,1,1))
  #   layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
  #   image(k, col=r) # axes=TRUE, add=TRUE,ylab="Angle of Left Eye" plot the image dev.off()
  #   #axis(1,outer=TRUE,"Angle of Right Eye")
  #   title(paste(strCond,"lID:",i," #e=", procDatIdx, " #F:",procDatFrames),collapse=NULL);
  #   par(mar=c(0,2,1,0))
  #   #Horizontal Axis #Need To Fix Order From -ve To +ve Angle On Histogramme
  #   barplot(hR$counts, axes=F, ylim=c(0, top), space=0, col='red')
  #   par(mar=c(2,0,0.5,1))
  #   #Vertical Axis
  #   barplot(hL$counts, axes=F, xlim=c(0, top), space=0, col='green', horiz=T)
  #   
  #   dev.off()
  # }
  idx <-idx+1;
    
} ##For Each Exp ID


###### BINARIZED HISTOGRAM PER GROUP ###
## Now Sum All LArva Binarized Response and Display Heat Map
hGroupbinDensity <- Reduce('+', hbinRL)
strDensityplotFileName <- paste(strPlotExportPath,"/binDensity/EyeAngleDensity-BINSet-",strCond,".pdf",collapse=NULL,sep="");
#pdf(strDensityplotFileName,width=8,height=8)
  sampleSize  <- length(vexpID) #Number of Larvae Used 
  hotMap <- c(rfHot(sampleSize),"#FF0000");
  image((-G_THRESHCLIPEYEDATA:G_THRESHCLIPEYEDATA),(-G_THRESHCLIPEYEDATA:G_THRESHCLIPEYEDATA),hGroupbinDensity,axes=TRUE,
      col=hotMap,xlab="Right Eye Angle",ylab="Left Eye Angle",asp=1)
  title(paste(strCond,"R-L Eye Density #n=", sampleSize,  " T:",round(TotalprocDatFrames/G_APPROXFPS),"sec"),collapse=NULL);
#dev.copy(jpeg,filename=paste(strDensityplotFileName,"-plot.jpg"));
#dev.off()
###

#### Eye Density - With Frequencies Over the whole Group #####
 strDensityplotFileName <- paste(strPlotExportPath,"/densities/EyeAngleDensity-Set-",strCond,".pdf",collapse=NULL,sep="");
#pdf(strDensityplotFileName,width=8,height=8)
  
  procDatFrames <- NROW(datAllGroupFrames)
  datAllGroupFrames$LEyeAngle <- clipEyeRange(datAllGroupFrames$LEyeAngle,lMin,lMax)
  datAllGroupFrames$LEyeAngle <-medianf(datAllGroupFrames$LEyeAngle,nEyeFilterWidth)
  
  datAllGroupFrames$REyeAngle <- clipEyeRange(datAllGroupFrames$REyeAngle,lMin,lMax)
  datAllGroupFrames$REyeAngle <-medianf(datAllGroupFrames$REyeAngle,nEyeFilterWidth)
  datAllGroupFrames[is.na(datAllGroupFrames$LEyeAngle),"LEyeAngle"] <- 0
  datAllGroupFrames[is.na(datAllGroupFrames$REyeAngle),"REyeAngle"] <- 0
  
  bw <- bandwidth.nrd(datAllGroupFrames[!is.na(datAllGroupFrames$LEyeAngle),"LEyeAngle"])
     
  bw <- ifelse(is.na(bw),0,bw)
   message(paste("kde BWdth:",bw));
  if (bw==0)
  {    bw <- 1.2
   message(paste("**Warning changed kde BWdth to fixed value -> ",bw));
  }
   eGroupDens <- kde2d(datAllGroupFrames$REyeAngle,datAllGroupFrames$LEyeAngle,h=bw, n=60, 
                       lims=c(range(-G_THRESHCLIPEYEDATA,G_THRESHCLIPEYEDATA),range(-G_THRESHCLIPEYEDATA,G_THRESHCLIPEYEDATA)) ) 
   #(-lMax:lMin),(-lMin:lMax)
  hotMap <- c(rfHot(20),"#FF0000");
  image(eGroupDens,col=hotMap,xlab="Right Eye Angle",ylab="Left Eye Angle",asp=1)
  sampleSize  <- length(vexpID) #Number of Larvae Used 
  title(paste(strCond,"R-L Eye Density #n=", sampleSize,  " T:",round(procDatFrames/G_APPROXFPS),"sec"),collapse=NULL);
#dev.off()
################################

