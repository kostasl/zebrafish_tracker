## Used to Manually Label The hunt Events Stored In datHuntEvent ##
##strVideoFilePath = "/home/kostasl/workspace/build-zebraprey_track-Desktop-Debug"
## Kostas Lagogiagiannis 2018 Jan
## Run The tracker Specifically on video frames isolating the Hunt Events - let the user label if the event was succesful or not

##To Execute The QT tracker application We may need to give the QT library Path - (xcb error)
#Sys.setenv(LD_LIBRARY_PATH="/home/kostasl/Qt/5.9.2/gcc_64/lib/" )
Sys.setenv(LD_LIBRARY_PATH=paste(Sys.getenv("LD_LIBRARY_PATH"),"/home/kostasl/Qt/5.9.2/gcc_64/lib/",sep=":" ) )
Sys.setenv(LD_LIBRARY_PATH=paste(Sys.getenv("LD_LIBRARY_PATH"),"/usr/lib/x86_64-linux-gnu/",sep=":" ) )


labelHuntEvents <- function(datHuntEvent,strVideoFilePath,strTrackerPath)
{
  
  for (i in  1:NROW(datHuntEvent) )
  {
    rec <- datHuntEvent[i,] 
    
    ##Get Respective Video  Filename / Path
    strVideoFile <- list.files(path =strVideoFilePath, pattern = rec$filenames , all.files = FALSE,
               full.names = TRUE, recursive = TRUE,
               ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
    
    message(paste("\n Examining Hunt Event of Larva:",rec$expID," Event:",rec$eventID, "Video:",rec$filenames) )
    ##--
    strArgs = paste(" --invideofile=",strVideoFile," --outputdir=",strTrackeroutPath," --startframe=",max(0,rec$startFrame-1)," --stopframe=",rec$endFrame,sep="")
    message(strArgs)
    execres <- base::system2(command=paste(strTrackerPath,"/zebraprey_track",sep=""),args =  strArgs,stdout="",stderr=FALSE)
    #stopifnot(execres == 0 ) ##Stop If Application Exit Status is not success
    ##Set Hunt Score - Given by user input after video is examined
    rec$huntScore <- as.numeric(readline(prompt="###Was this Hunt Succesfull? Type 1 (Yes)/-1(No)/0(Leave Undefined) :"))
    stopifnot(is.numeric(rec$huntScore) == TRUE )
    
    if (rec$huntScore == 1)
    {
      datHuntEvent[i,"huntScore"] <-  rec$huntScore
      message("~Mark Succesfull")
    }
    if (rec$huntScore == 0)
    {
      message("~Leave Unlabelled")
      datHuntEvent[i,"huntScore"] <-  rec$huntScore
    }
    if (rec$huntScore == -1)
    {
      message("~Failed To Capture Prey")
      datHuntEvent[i,"huntScore"] <-  rec$huntScore
    }
    
     message("-Proceeding to Next Video.")
  } ## For Each Hunt Event Detected - 
  
  
}