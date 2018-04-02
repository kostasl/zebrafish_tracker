## Used to Manually Label The hunt Events Stored In datHuntEvent ##
##strVideoFilePath = "/home/kostasl/workspace/build-zebraprey_track-Desktop-Debug"
## Kostas Lagogiagiannis 2018 Jan
## Run The tracker Specifically on video frames isolating the Hunt Events - let the user label if the event was succesful or not

##To Execute The QT tracker application We may need to give the QT library Path - (xcb error)
#Sys.setenv(LD_LIBRARY_PATH="/home/kostasl/Qt/5.9.2/gcc_64/lib/" )
##Check If Qt Is already Added To Exec Path
if (grepl("Qt",Sys.getenv("LD_LIBRARY_PATH") )  == FALSE) 
{
  #Sys.setenv(LD_LIBRARY_PATH=paste(Sys.getenv("LD_LIBRARY_PATH"),"/home/kostasl/Qt/5.9.2/gcc_64/lib/",sep=":" ) )
  Sys.setenv(LD_LIBRARY_PATH="")
  Sys.setenv(LD_LIBRARY_PATH=paste(Sys.getenv("LD_LIBRARY_PATH"),"/opt/Qt/5.9/5.9/gcc_64/lib",sep=":" ) ) ##Home PC
  Sys.setenv(LD_LIBRARY_PATH=paste(Sys.getenv("LD_LIBRARY_PATH"),"/usr/lib/x86_64-linux-gnu/",sep=":" ) )
}

vHuntEventLabels <- c("NA","Success","Fail","No_Target","Not_HuntMode/Delete","Escape","Out_Of_Range","Duplicate/Overlapping","Fail-No Strike","Fail-With Strike","Success-SpitBackOut")
huntLabels <- factor(x=5,levels=c(1,2,3,4,5,6,7,8,9,10,11),labels=vHuntEventLabels )##Set To NoTHuntMode

labelHuntEvents <- function(datHuntEvent,strDataFileName,strVideoFilePath,strTrackerPath,strTrackOutputPath)
{
  message(paste(NROW(datHuntEvent[datHuntEvent$huntScore >0,]),"/",NROW(datHuntEvent), " Data has already been labelled" ) )
  nLabelledSuccess <- NROW(datHuntEvent[datHuntEvent$huntScore == which(levels(huntLabels) == "Success") | datHuntEvent$huntScore == which(levels(huntLabels) == "Success-SpitBackOut"),])
  readline(prompt="-.Begin Data labelling.-")
  
  
  for (i in  (1:NROW(datHuntEvent)) )
  {
    rec <- datHuntEvent[i,] 
    
    if (rec$huntScore != 0 && rec$huntScore != which(levels(huntLabels)=="NA") )
        next ##SKip Record if previously Labelled

    ##For Larva That Did not register any sufficient Hunting Events -  An Empty Record has been added To Acknowledge 
    if (rec$eventID == 0 & rec$huntScore == 0)
    {##Set To Not Hunt Event/Delete - So as to Ignore In Hunt Event Counts
      datHuntEvent[i,]$huntScore = huntLabels
      next
    }
      
    ##Get Respective Video  Filename / Path
    strVideoFile <- list.files(path =strVideoFilePath, pattern = rec$filenames , all.files = FALSE,
               full.names = TRUE, recursive = TRUE,
               ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
    
    message(paste("\n",i,". Examining Hunt Event of Larva:",rec$expID," Event:",rec$eventID, "Video:",rec$filenames, " -s:",max(0,rec$startFrame-1)," -e:",rec$endFrame) )
    ##--
    strArgs = paste("--ModelBG=0 --invideofile=",strVideoFile," --outputdir=",strTrackOutputPath," --startframe=",max(0,rec$startFrame-1)," --stopframe=",rec$endFrame,sep="")
    message(paste(strTrackerPath,"/zebraprey_track",strArgs,sep=""))
    execres <- base::system2(command=paste(strTrackerPath,"/zebraprey_track",sep=""),args =  strArgs,stdout="",stderr=TRUE)
    
    ## execres contains all of the stdout - so cant be used for exit code
    #stopifnot(execres == 0 ) ##Stop If Application Exit Status is not success
    ##Show Labels And Ask Uset input after video is examined
    Keyc = 1000 ##Start With Out Of Range Value
    failInputCount = 0
    flush(con=stdin())
    flush(con=stdout())
    
    ##Log Updates to Separate File ## 
    bColNames = FALSE
  
    while ( Keyc != 'c' & failInputCount < 3 | nchar(Keyc) > 1) #is.na(factor(levels(huntLabels)[as.numeric(Keyc)] ) ) 
    {
      l <- 0
      for (g in levels(huntLabels) )
      {
        l=l+1
        message(paste(l,g,sep="-"))
      }
      l=l+1
      message(paste("f","Fix Frame Range",sep="-"))
      l=l+1
      message(paste("c","End Labelling Process",sep="-"))
      
      failInputCount <- failInputCount + 1
      Keyc <- readline(prompt="###Was this Hunt Succesfull? (# / c to END) :")
      #message(paste(failInputCount,Keyc) )
      
      ##Check for Stop Loop Signal
      print("###")
      if (Keyc == 'c') 
      {
        message(paste(Keyc,"~End Label Process Here") )
        rec$huntScore <- 0
        return(datHuntEvent)
        
      }else {
        rec$huntScore <- as.numeric(Keyc) ##factor(levels(huntLabels)[as.numeric(c)]
      }
      ##Fix Frame Range Manually
      if (Keyc == 'f')
      {
        rec$startFrame <- as.numeric(readline(prompt=paste(" Enter new start frame (",rec$startFrame, "):") ) )
        if (!is.na(rec$startFrame))
          if (rec$startFrame > 0)
            datHuntEvent[i,"startFrame"] <- rec$startFrame
          
        rec$endFrame <- as.numeric( readline(prompt=paste(" Enter new end frame (",rec$endFrame, "):") ) )
        if (!is.na(rec$endFrame))
          if (rec$endFrame > 0)
            datHuntEvent[i,"endFrame"] <- rec$endFrame
        Keyc = 1000
        
        message(paste("*New start:",datHuntEvent[i,"startFrame"]," end frame:", datHuntEvent[i,"endFrame"],"\n") )
      }
      ##User Has selected Label? Then Break From menu loop
      if (!is.na(factor(levels(huntLabels)[as.numeric(Keyc)] ) ) )
        break

    } ##End Menu Loop
    

    datHuntEvent[i,"huntScore"] <-rec$huntScore
    message(datHuntEvent[i,"huntScore"])
     if (rec$huntScore == which(levels(huntLabels)=="Success") )
     {
       message("~Mark Succesfull")
     }
    if (rec$huntScore == 0 || rec$huntScore == which(levels(huntLabels)=="NA"))
     {
       message("~Leave Unlabelled")
     }
    if (rec$huntScore == which(levels(huntLabels)=="Fail"))
     {
       message("~Failed To Capture Prey")
     }
    # 
    if (Keyc == 'c')  ##Stop Event Loop Here if c was pressed
    {
      message(" Stop Labelling Loop Here " )
      return(datHuntEvent)
      #break
    }
    else
       message(paste(levels(huntLabels)[as.numeric(Keyc)] , "-Proceeding to Next Video.") )
     
     #####################################################################################################
    ##### Save With Dataset Idx Identifier On Every Labelling As An Error Could Lose Everything  ########
    save(datHuntEvent,file=paste(strDataFileName,".RData",sep="" ))      
    strOutFileName <- paste(strDataFileName,"-updates.csv",sep="")
    bColNames = FALSE
    if (!file.exists(paste(strDataFileName,"-updates.csv",sep="" )) )
      bColNames = TRUE

    ##Need to Use write.table if you want to append down this list- I Found out after a lot of messing around
     write.table(datHuntEvent[i,],file=strOutFileName, append = TRUE,dec='.',sep=',',col.names = bColNames,quote=FALSE) ##Append Labelled records to a Log CSV File
    # }else ##Work Around Cause we cannot Append To Csvs!
    # {
    #   flogUpdates <- file(strOutFileName,'a',blocking = TRUE) ##Open File
    #   write((rec),file=flogUpdates,sep=",",append=TRUE,ncolumns=length(rec) )
    #   #writeChar('',con=flogUpdates)
    #   close(con=flogUpdates)
    # }
    ###########################################################################################
    
  } ## For Each Hunt Event Detected - 
  
  ##Return Modified Data Frame
  return(datHuntEvent)
  
}

