## Helper Script to Run The Trajectory Analyses
setwd("/home/kostasl/workspace/zebrafishtrack/Rplots")
#setwd(here())
source("config_lib.R")
source("TrajectoryAnalysis.r")
source("HuntingEventAnalysis_lib.r")

setEnvFileLocations("HOME") #HOME,OFFICE,#LAPTOP
## To compile to different destination run:
# rmarkdown::render("ForagingStateAnalysis.Rmd",output_dir = paste0(strDataExportDir,'../foragingAnalysis') )

##Re Run All Path Calculation 
generate_DispersionFiles <- function(timePointSequence = seq(2,30,2))
{
  ##Generate Files From Start 
  loaddatAllFrames()
  for (tt in timePointSequence )
    calcTrajectoryDispersions(datAllFrames,tsec_timeWindow)
}