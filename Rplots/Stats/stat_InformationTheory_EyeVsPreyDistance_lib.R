### 24 Nov 2018 - Separated functions used to calculate the information content of eye vergence regressed from hunt events for the distance to prey 
## Assumes uniform distribution for probability of distance from prey - 5bit encoding of distance vector.
#library("entropy")


## Estimate Response Based On Model's mean response ## 
phi_hat <- function(x,Ulist){
  return(Ulist$phi_0  
         +(Ulist$phi_max - Ulist$phi_0)/( 1 + exp( -Ulist$lambda*( Ulist$tau - x)  ) )
         +Ulist$C*exp(Ulist$gamma*( Ulist$tau - x)) )
}

#Returns the FrequenciesAround mean (phihat)
phiDens <- function(phi,x,Ulist)
{
  ### Gaussian around mean response  ###
  return( dnorm(phi,mean=phi_hat(x,Ulist),sd=Ulist$sigma   ) ) #sqrt()
}

##
## Calc Info In single Sample & Hunt Event
InfoCalc <- function(DistRange,minDist,maxDist,Ulist)
{
  PhiRange <- seq(0,100,1) ##We limit The information Obtained To Reasonable Ranges Of Phi (Vergence Angle)
  
  
  #Generates a list with all possible pairs of in the range of Phi And Distance X 
  Grid <- expand.grid(PhiRange,DistRange)
  PVec=rep(0,NROW(Grid)) 
  ##Calc Density for all input Space
  for (i in 1:NROW(Grid) )
  {
    PVec[i] <- phiDens(Grid[i,1],Grid[i,2],Ulist)
  }
  ##Normalize to Probability
  PVec=PVec/sum(PVec)
  
  # Convert Pvec to a matrix / 
  PMatrix=matrix(PVec,nrow=length(PhiRange),ncol=length(DistRange),byrow = FALSE)
  
  ##Image shows a transpose Of the PMatrix, - In reality Rows contain Phi, and Cols are X 
  #image(t(PMatrix),y=PhiRange,x=DistRange)
  MargVec=rowSums(PMatrix) ### Marginalize across X to obtain P(Response/Phi)
  
  Iloc=PMatrix/MargVec*length(DistRange) ##Information On Local x/For Each X - Assume X is unif. and so Prob[X]=1/Length(X)
  
  ###row sum
  sel=PMatrix>0
  #INFO=sums(PMatrix[sel]*log2(Iloc[sel]) )
  ### Return Marginals I_xPhi
  mInfo <- PMatrix*log2(Iloc)
  
  ##Calc Marginals only on relevant Region - up to min distance from Target
  ##Make Binary Array indicating OutOf Bound regions
  vIntRegion<-as.numeric( DistRange >= minDist)
  mIntRegion <- t(matrix(vIntRegion,dim(mInfo)[2],dim(mInfo)[1]))
  
  INFO=colSums(mInfo*mIntRegion,na.rm=TRUE )
  
  return(INFO)
}


## Calculates Mutual Information X to Phi 
##  Fitted Function instances for each Hunt event
##  Assume distance is encoded in Dx steps, in N bits
## Returns a matrix Giving N Inf. Measures for each Hunt event
calcInfoOfHuntEvent <- function(drawS,dataSubset,n=NA,groupID)
{
  bPlot <- TRUE
  DistMin = 0.5
  DistMax =  5
  ##Assume X distance is encoded using 5 bits
  DistRange <- seq(DistMin,DistMax, (DistMax-DistMin )/(2^5-1)  )
  
  
  NHuntEvents <- NROW(unique(dataSubset$hidx) )
  
  vRegIdx <- unique(dataSubset$RegistrarIdx) ##Get Vector Of RegIdx That Associate with the sample Sequence
  vsampleP <- unique(dataSubset$hidx) # Obtain the seq id, so we can link to drawSamples
  
  if (is.na(n))
  {
    n <- NROW(unique(dataSubset$hidx))
    ##Note:Suffled and so the content of vsampleP maintains the hidx from where we can recover the RegIdx
    vsampleP <- sample(unique(dataSubset$hidx),n) 
  }
  
  
  vsub <- which (dataSubset$hidx %in% vsampleP)
  
  
  
  mInfMatrix <- matrix(nrow=NSamples,ncol=NROW(vsampleP) )
  vsampleRegisterIdx <- vector()
  
  hCnt <- 0
  idxChain <- 1
  for (h in vsampleP)
  {
    print(h)
    hCnt <- hCnt + 1
    ##Hold RegIDx For This Hunt Event 
    vsampleRegisterIdx[hCnt] <- vRegIdx[h] ##Save the Registry Idx So we can refer back to which Fish This belongs to 
    
    vphi_0_sub <- tail(drawS$phi_0[h,,idxChain],n=NSamples)
    vphi_max_sub <- tail(drawS$phi_max[h,,idxChain],n=NSamples)
    vgamma_sub <- tail(drawS$gamma[h,,idxChain],n=NSamples)
    vlambda_sub <- tail(drawS$lambda[h,,idxChain],n=NSamples)
    vtau_sub   <- tail(drawS$tau[h,,idxChain],n=NSamples)
    vsigma_sub <- tail(drawS$sigma[h,,idxChain],n=NSamples)
    valpha_sub <- tail(drawS$alpha[h,,idxChain],n=NSamples)
    
    ##Use the minimum Distance to prey, to obtain the effective integration Zone That information should be measured from 
    DistRangeEffective_min <- min(dataSubset$distP[dataSubset$hidx == h])
    DistRangeEffective_max <- max(dataSubset$distP[dataSubset$hidx == h])
    
    #DistRangeEffective_min <- 0.6 ##Set to fixed int interval
    
    vPP <- which (dataSubset$hidx == h)
    #pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_",strGroupID[groupID],"_Sigmoid_",pp,".pdf",sep="")) 
    if (bPlot)
    {
      pdf(file= paste(strPlotExportPath,"/stat/stat_InfMeasure_EyeVsDistance5mm_",strGroupID[groupID],"_SigExp_",vRegIdx[h],"_",strTag,".pdf",sep="")) 
      par(mar = c(5,5,2,5))
      plot(dataSubset$distP[vPP],dataSubset$phi[vPP],pch=19,xlim=c(0,max(DistRange)),ylim=c(0,90),
           main=paste(strGroupID[groupID],h," (",vRegIdx[h],")") , 
           bg=colourP[2],col=colourP[1],
           cex=0.5,
           xlab="Distance From Prey",ylab=expression(paste(Phi)) )  
      #Plot The Fitted Function
      vY  <-  mean(valpha_sub)*exp(mean(vgamma_sub)*( mean(vtau_sub)-DistRange) )+ mean(vphi_0_sub)   +  ( mean(vphi_max_sub) - mean(vphi_0_sub)  )/(1+exp( -( mean(vlambda_sub)   *( mean(vtau_sub) -DistRange )   ) ) ) 
      lines( DistRange ,vY,type="l",col=colourR[3],lwd=3)
    }
    
    lPlist <- list()
    for (i in 1:NSamples)
    {
      
      lPlist[[i]] <- list(phi_0=vphi_0_sub[i],phi_max=vphi_max_sub[i], gamma= vgamma_sub[i],
                          tau=vtau_sub[i],
                          lambda=vlambda_sub[i], 
                          C=valpha_sub[i] #vC_sub[i]
                          , sigma=vsigma_sub[i] )
      
      ##Information Is in The Sum / Marginal Across X ##
      ##Keep Distance Fixed accross so we are comparing against the same amount of information baseline, 
      ## from which mutual inf gives us the reduction in uncertainty
      ## Ex. We need log2(NROW(DistRange)) of bits to encode distance at the res. defined in DistRange
      
      vPPhiPerX <- InfoCalc(DistRange,DistRangeEffective_min,DistRangeEffective_max,Ulist = lPlist[[i]]) ##Get MutInf Per X
      
      ##Plot The Information Content Of The fitted Function ##
      if (bPlot)
      {
        par(new=T) 
        plot(DistRange,vPPhiPerX,ylim=c(0,2.5),xlim=c(0,max(DistRange) ),axes=F,type="p",pch=19,col=colourP[4], xlab=NA,
             ylab=NA,sub=paste("x requires ", round(100*log2(NROW(DistRange)) )/100,"bits"  ) )
        ##Intergate Inf Across Distance
        lines(DistRange,rev(cumsum(rev(vPPhiPerX))),ylim=c(0,2.5),xlim=c(0,max(DistRange) ),type="l",col=colourR[4], xlab=NA, ylab=NA )
      }
      
      mInfMatrix[i,hCnt] <- sum(vPPhiPerX) ## Mutual Inf Of X and Phi
    } 
    
    if (bPlot)
    {
      axis(side=4)
      mtext(side = 4, line = 3, 'Mutual information (bits)')
      dev.off()
    }
  }##For Each Hunt Event   
  
  
  ##Next Is integrate Over sampled points, and Make Ii hidx matrix 
  lInfoMatStruct <- list(infoMatrix=mInfMatrix,vsampleRegisterIdx=vsampleRegisterIdx,vsamplePSeqIdx=vsampleP)
  return(lInfoMatStruct)
}


##Returns Data frame With Mean Information in EyeVergence Along With Undershoot Ratio For each Tracked Event
mergeFirstTurnWithInformation <- function(datFirstBouts,lInfStruct)
{
  
  ##Subset the Event We have a first Turn Bout and Have Measured Eye INf (which are filtered for doing a strike attack)
  datBoutSubset <- datFirstBouts[datFirstBouts$RegistarIdx %in%  lInfStruct$vsampleRegisterIdx,]
  
  ## Retrieve Turn Bouts in right order to match the information Matric Column Order
  lboutSet <- list()
  idx <- 0
  for (idxR in lInfStruct$vsampleRegisterIdx)
  {
    idx <- idx +1
    if (NROW(datFirstBouts[datFirstBouts$RegistarIdx == idxR, ]) == 0 )
      next()
    
    rec <- datFirstBouts[datFirstBouts$RegistarIdx == idxR, ]
    vInfMeasure <- mean(lInfStruct$infoMatrix[,idx])
    
    lboutSet[[idx]] <- list(RegistarIdx=idxR,hidx= lInfStruct$vsamplePSeqIdx[idx], 
                            UnderShootAngle=abs( (rec$OnSetAngleToPrey)-abs(rec$Turn) ),
                            UnderShootRatio = (abs(rec$Turn)/abs(rec$OnSetAngleToPrey) ) ,
                            MInf=vInfMeasure) 
  }
  
  ##Combine Into A dataframe
  datFirstBoutVsInf<- data.frame(do.call(rbind,lboutSet ))
  
  
  return(datFirstBoutVsInf)
}

################# End of Library Functions ##### 
