## Organize Manuscript Figures ### 
#### Kostas Lagogiannis 2019 
## \brief Make a scipt clarifying the script files used to produce each figure Used in the MS 



library(tools)
library(RColorBrewer);
library("MASS");
library(extrafont) ##For F


source("config_lib.R")
setEnvFileLocations("HOME") #OFFICE,#LAPTOP


####################
#source("TrackerDataFilesImport.r")
### Hunting Episode Analysis ####


### Fig 1 ####
## The kinematics was produced by selecting one of the figure produced
## from Hunt event analysis loop in : runHuntEpisodeAnalysis.r
## 
# Time Line manually designed - revised by Martin 


### Fig 2 ####
source("Stats/stat_HuntRateInPreyRange.R")
source("Stats/stat_HuntDuration.R")


### Fig 3 ####
source("Stats/stat_HuntEfficiency.r")

### Fig 4 ####
source("DataLabelling/plotLabelledDataResults.R")
source("Stats/stat_CaptureSpeedVsDistanceToPrey.R")

### Fig 5 ####
source("Stats/stat_LinRegression_TurnVsBearing.R")

## Fig 6  clustering Speed, TurnRatio and Distance to Prey using 3D 2xGaussian mixture method####
source("Stats/stat_ClusterCaptureSpeedVsUndershootAndDistance.r")

### Fig 7 Show Covariance using Gaussian 3D non clustering model aong wit ####
source("Stats/stat_CaptureSpeedVsUndershootAndDistance.r")


