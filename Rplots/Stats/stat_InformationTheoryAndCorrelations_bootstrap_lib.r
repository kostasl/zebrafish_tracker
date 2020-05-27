
##
## Calc Info In single Sample & Hunt Event
InfoCalc_get2DFreq <- function(datX,datY,XRange,YRange,plot=FALSE)
{
  
  ### Tally 2D data points in the grid
  nXbins <- 16 ##ie 2^3 - x bit encoding of Distance
  nYbins <- 16 ##ie 2^3 - x bit encoding of Distance
  
  x.bin <- seq(floor(min(XRange)), (max(XRange)), length=nXbins)
  y.bin <- seq(floor(min(YRange)), (max(YRange)), length=nYbins)
  
  
  freq <-  as.data.frame(table(X=findInterval(datX, x.bin),Y=findInterval(datY, (y.bin) )))
  freq[,1] <- as.numeric(as.character(freq[,1]))
  freq[,2] <- as.numeric(as.character(freq[,2]))
  
  ##Place frequencies on 2D matrix 
  freq2D <- diag(nrow=nYbins,ncol=nXbins)*0
  ##Reverse the Y, so matrix pos matches speed labels
  freq2D[cbind(as.numeric( as.character(freq[,2]) ), as.numeric(as.character(freq[,1] )  ) ) ] <- freq[,3]
  
  ##Draw Matrix
  if (plot)
  {
    par(mfrow=c(1,3),pty="s")
    plot(datX,datY,xlim=XRange,ylim=YRange)
    image(x.bin, (y.bin), t(freq2D), col=topo.colors(max(freq2D)+1))
    contour(x.bin, y.bin, t(freq2D), add=TRUE,nlevels=5, col=rgb(1,1,1,.7))
    palette(rainbow(max(freq2D)))
    #cols <- (freq2D[-1,-1] + freq2D[-1,-(nYbins-1)] + freq2D[-(nXbins-1),-(nYbins-1)] + freq2D[-(nXbins-1),-1])/4
    persp(freq2D ,theta=45) #col=cols
  }
  
  
  return(freq2D)
}##
##Debug ##
#datX = datCapture_DL$Undershoot
#datY = datCapture_DL$CaptureSpeed


##Shannon Entropy
H_entropy<- function(in_freq) 
{
  in_freq_norm <- in_freq/sum(in_freq)
  in_freq_norm = in_freq_norm[in_freq_norm >0 ]
  
  ## Here we assume uniform distribution among X, which may be incorrect
  H = -sum(in_freq_norm* log2(in_freq_norm))  
  
  return(H)
}
## Normalize Grid To Sum of Points to obtain Joint Prob P(Capt Speed| Distance)

## Calc Empririca Marginal Entropies and Mutual information in 2D frequency/Count matrix - 
##Measures average reduction in uncertainty about X that results from learning Y, or vice-versa - information of X on Y
## I(X;Y) = I(Y;X)
calcMIEntropy <- function(freqM)
{
  H_Y <- H_entropy(rowSums(freqM) )
  H_X <- H_entropy(colSums(freqM))
  H_XY <- H_entropy(freqM)
  ## Calc Mutual Info
  H_X + H_Y - H_XY
  ##Conditional Entropy is related to Joint Entropy
  H_XCY <- H_XY - H_Y
  ## Calc Mutual Info between X and Y
  I_XY <- H_X -  H_XCY
  
  return (list(H_X=H_X, H_Y=H_Y, H_XY=H_XY, H_XGivenY=H_XCY, MutualInf_XY=I_XY) )
}



## Bootstrap Data analysis from 2 chosen columns of datCapture_NL (lFirstBout)  to get stats on correlations and Mutual Information
## Can use Pearson to examine correlations in value, 
## Also Returns a control distribution obtained by suffling the sampled Y columns Breaking pairs against X
## Or spearman to examine correlation in rank (ie value increases correlate and not corr not influenced by the absolute value of each data point)
bootStrap_stat <- function(datCapture_X,datCapture_Y,N,XRange,YRange,corMethod="pearson")
{
  datCapture <- data.frame(cbind(datCapture_X,datCapture_Y))
  
  l_sampleXYAnalysis <- list()
  for (i in 1:N)
  {
    #freqM_DF <- InfoCalc_get2DFreq(datCapture_DL$Undershoot,datCapture_DL$CaptureSpeed,XRange,YRange)
    #freqM_LF <- InfoCalc_get2DFreq(datCapture_LL$Undershoot,datCapture_LL$CaptureSpeed,XRange,YRange)
    idxSample <- sample(1:NROW(datCapture),size=floor(NROW(datCapture)*0.80))
    datSub <- datCapture[idxSample,]

    freqM_NF <- InfoCalc_get2DFreq(datSub[,1],datSub[,2],XRange,YRange)
    infC <- calcMIEntropy(freqM_NF) ##Mutail INformation
    corrXY <- cor(datSub[,1],datSub[,2],method=corMethod) #method="spearman"
    corrXY_suffled <- cor(datSub[,1],sample(datSub[,2]) ,method=corMethod) ## Calculate Control Corr By Suffling the data
    l_sampleXYAnalysis[[i]] <- data.frame(MI = infC$MutualInf_XY,entropy_X = infC$H_X,entropy_Y = infC$H_Y,corr=corrXY,corr_suffled=corrXY_suffled)
    #inf_LF <-calcMIEntropy(freqM_LF)
    #inf_DF <-calcMIEntropy(freqM_DF)
  }
  datXYAnalysis <- do.call(rbind,l_sampleXYAnalysis)
  
  return(datXYAnalysis)
}
