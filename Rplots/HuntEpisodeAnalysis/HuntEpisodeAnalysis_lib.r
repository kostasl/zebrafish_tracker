### Library For Analysis of Hunt Episodes ###

library(signal)
library(MASS)
library(mclust,quietly = TRUE)

#library(sBIC)

citation("mclust")


##Clusters Fish Speed Measurements into Bout And Non Bout
##Use 3 For Better Discrimination When  There Are Exist Bouts Of Different Size
detectMotionBouts <- function(dEventSpeed_smooth)
{
  prior_factor3 <- 0.1 ## Adds a prior shift in the threshold Of Classification
  prior_factor2 <- 0.8 ## Adds a prior shift in the threshold Of Classification 
  colClass <- c("#FF0000","#00FF22","#0000FF")
  
  #t <- datRenderHuntEvent$frameN
  
  #BIC <- mclustBIC(dEventSpeed)
  
  ### INcreased to 3 Clusters TO Include Other Non-Bout Activity
  fit <- Mclust(dEventSpeed_smooth,G=3 ) #modelNames = "V" prior = priorControl(shrinkage = 0) 
  summary(fit)
  
  region <- min(NROW(t),NROW(dEventSpeed_smooth))
  X11()
  plot(fit, what="density", main="", xlab="Velocity (Mm/s)")
  rug(dEventSpeed)
  
  #X11()
  #boutClass <- fit$classification
  #plot(dEventSpeed[1:region],type='l',col=colClass[1])
  #points(which(boutClass == 2), dEventSpeed[boutClass == 2],type='p',col=colClass[2])
  
  #points(which( fit$z[,2]> fit$z[,1]*prior_factor ), dEventSpeed[ fit$z[,2]> fit$z[,1]*prior_factor  ],type='p',col=colClass[3])
  ## Add Prior Bias to Selects from Clusters To The 
  return (which( fit$z[,3]> fit$z[,1]*prior_factor3 | fit$z[,3]> fit$z[,2]*prior_factor2   ))
  
}