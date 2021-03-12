## Helper Script to Run The Trajectory Analyses
setwd("/home/kostasl/workspace/zebrafishtrack/Rplots")
#setwd(here())
source("config_lib.R")
source("TrajectoryAnalysis.r")
source("HuntingEventAnalysis_lib.r")

setEnvFileLocations("OFFICE") #HOME,OFFICE,#LAPTOP
## To compile to different destination run:
# rmarkdown::render("ForagingStateAnalysis.Rmd",output_dir = paste0(strDataExportDir,'../foragingAnalysis') )

loaddatAllFrames <- function(forceReload = FALSE)
{
  if (!exists("datAllFrames") | forceReload)
  {
    attach(paste(strDatDir,"datAllFramesFix1_Ds-5-19.RData",sep="/"))
    attach(paste(strDatDir,"groupsrcdatListPerDataSet_Ds-5-19.RData",sep="/"))
  }
}


##Re Run All Path Calculation 
generate_DispersionFiles <- function(timePointSequence = seq(2,30,2))
{
  ##Generate Files From Start 
  loaddatAllFrames()
  for (tt in timePointSequence )
    calcTrajectoryDispersions(datAllFrames,tt)
}

generate_DispersionFiles(2)