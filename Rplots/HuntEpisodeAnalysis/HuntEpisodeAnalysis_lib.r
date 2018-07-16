### Library For Analysis of Hunt Episodes ###

library(signal)
library(MASS)
library(mclust,quietly = TRUE)

library(sBIC)

citation("mclust")


##Clusters Fish Speed Measurements into Bout And Non Bout
detectMotionBouts <- function(dEventSpeed)
{
  prior_factor <- 0.05 ## Adds a prior shift in the threshold Of Classification 
  colClass <- c("#FF0000","#00FF22","#0000FF")
  
  #t <- datRenderHuntEvent$frameN
  
  BIC <- mclustBIC(dEventSpeed)
  
  fit <- Mclust(dEventSpeed,G=2)
  summary(fit)
  
  region <- min(NROW(t),NROW(dEventSpeed))
  X11()
  plot(fit, what="density", main="", xlab="Velocity (Mm/s)")
  rug(dEventSpeed)
  
  X11()
  boutClass <- fit$classification
  plot(dEventSpeed[1:region],type='l',col=colClass[1])
  points(which(boutClass == 2), dEventSpeed[boutClass == 2],type='p',col=colClass[2])
  
  points(which( fit$z[,2]> fit$z[,1]*prior_factor ), dEventSpeed[ fit$z[,2]> fit$z[,1]*prior_factor  ],type='p',col=colClass[3])
  
  return (which( fit$z[,2]> fit$z[,1]*prior_factor  ))
  
}