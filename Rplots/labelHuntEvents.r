## Used to Manually Label The hunt Events Stored In datHuntEvent ##
##strVideoFilePath = "/home/kostasl/workspace/build-zebraprey_track-Desktop-Debug"
## Kostas Lagogiagiannis 2018 Jan
## Run The tracker Specifically on video frames isolating the Hunt Events - let the user label if the event was succesful or not

##To Execute The QT tracker application We may need to give the QT library Path - (xcb error)
#Sys.setenv(LD_LIBRARY_PATH="/home/kostasl/Qt/5.9.2/gcc_64/lib/" )
##Check If Qt Is already Added To Exec Path
if (grepl("Qt",Sys.getenv("LD_LIBRARY_PATH") )  == FALSE) 
{
  Sys.setenv(LD_LIBRARY_PATH=paste(Sys.getenv("LD_LIBRARY_PATH"),"/home/kostasl/Qt/5.9.2/gcc_64/lib/",sep=":" ) )
  Sys.setenv(LD_LIBRARY_PATH=paste(Sys.getenv("LD_LIBRARY_PATH"),"/usr/lib/x86_64-linux-gnu/",sep=":" ) )
}

huntLabels <- factor(x=0,levels=c(1,-1,2,3,5,0),labels=c("Success","Fail","No_Target","Not_HuntMode","Escape","NA"))

labelHuntEvents <- function(datHuntEvent,strVideoFilePath,strTrackerPath,strTrackOutputPath)
{
  
  for (i in  1:NROW(datHuntEvent) )
  {
    rec <- datHuntEvent[i,] 
    
    if (rec$huntScore != 0 && rec$huntScore != which(huntLabels=="NA") )
        next ##SKip Record if previously Labelled
        
    ##Get Respective Video  Filename / Path
    strVideoFile <- list.files(path =strVideoFilePath, pattern = rec$filenames , all.files = FALSE,
               full.names = TRUE, recursive = TRUE,
               ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
    
    message(paste("\n Examining Hunt Event of Larva:",rec$expID," Event:",rec$eventID, "Video:",rec$filenames) )
    ##--
    strArgs = paste("--ModelBG=0 --invideofile=",strVideoFile," --outputdir=",strTrackOutputPath," --startframe=",max(0,rec$startFrame-1)," --stopframe=",rec$endFrame,sep="")
    message(strArgs)
    execres <- base::system2(command=paste(strTrackerPath,"/zebraprey_track",sep=""),args =  strArgs,stdout="",stderr=TRUE)
    #stopifnot(execres == 0 ) ##Stop If Application Exit Status is not success
    ##Show Labels And Ask Uset input after video is examined
    l <- 0
    for (g in levels(huntLabels) )
    {
      l=l+1
      message(paste(l,g,sep="-"))
    }
    
    c = -1
    while (!is.element(c,huntLabels )  || c != 'c')
    {
      c <- readline(prompt="###Was this Hunt Succesfull? Type # or c to end :")
    }
    rec$huntScore <- levels(huntLabels)[c]
    if (!is.numeric(rec$huntScore)) 
    {
      message("~End Label Process Here")
      return(0)
    }
      ## 
    
    datHuntEvent[i,"huntScore"] <-  rec$huntScore
    if (rec$huntScore == which(huntLabels=="Success"))
    {
      message("~Mark Succesfull")
    }
    if (rec$huntScore == 0 ||rec$huntScore == which(huntLabels=="NA"))
    {
      message("~Leave Unlabelled")
    }
    if (rec$huntScore == which(huntLabels=="Fail"))
    {
      message("~Failed To Capture Prey")
    }
    
     message("-Proceeding to Next Video.")
  } ## For Each Hunt Event Detected - 
  
  
}
