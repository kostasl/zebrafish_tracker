## 23-10-2018
### Consider Information content in eye vergence on the distance to prey
## We use Model and the sampled parameter values to obtain an estimate of the mean 
##  information content in each hunt episode 
## We calculate the mutual information using P(Phi | X), assuming X is uniform

source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")


#### CalcInformation ##
load(file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RJAgsOUt_4.RData",sep=""))

## Can Load Last Inf Matrix 
load(file=paste(strDataExportDir,"/stat_infoMat_EyeVergenceVsDistance_sigmoidFit5mm-5bit.RData",sep=""))
          

####Select Subset Of Data To Analyse
strRegisterDataFileName <- paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register",".rds",sep="") #Processed Registry on which we add 
message(paste(" Importing Retracked HuntEvents from:",strRegisterDataFileName))
datTrackedEventsRegister <- readRDS(strRegisterDataFileName) ## THis is the Processed Register File On 


##For the 3 Groups 
colourH <- c(rgb(0.01,0.01,0.9,0.8),rgb(0.01,0.7,0.01,0.8),rgb(0.9,0.01,0.01,0.8),rgb(0.00,0.00,0.0,1.0)) ##Legend
colourP <- c(rgb(0.01,0.01,0.8,0.5),rgb(0.01,0.6,0.01,0.5),rgb(0.8,0.01,0.01,0.5),rgb(0.20,0.40,0.5,1.0)) ##points DL,LL,NL
colourR <- c(rgb(0.01,0.01,0.9,0.4),rgb(0.01,0.7,0.01,0.4),rgb(0.9,0.01,0.01,0.4),rgb(0.00,0.00,0.0,0.3)) ##Region (Transparency)
pchL <- c(16,2,4)
ltL <-  c(1,2,3) ##line Types

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
InfoCalc <- function(DistRange,Ulist)
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
  INFO=colSums(PMatrix*log2(Iloc),na.rm=TRUE )
  return(INFO)
}


## Calculates Mutual Information X to Phi 
##  Fitted Function instances for each Hunt event
##  Assume distance is encoded in Dx steps, in N bits
## Returns a matrix Giving N Inf. Measures for each Hunt event
calcInfoOfHuntEvent <- function(drawS,dataSubset,n=NA,groupID)
{
  DistMin = 0.5
  DistMax =  5
  ##Assume X distance is encoded using 5 bits
  DistRange <- seq(DistMin,DistMax, (DistMax-DistMin )/(2^5-1)  )
  
  NSamples <-250
  NHuntEvents <- NROW(unique(dataSubset$hidx) )
  
  vsampleP <- unique(dataSubset$hidx)
  
  if (is.na(n))
  {
    n <- NROW(unique(dataSubset$hidx))
    vsampleP <- sample(unique(dataSubset$hidx),n)
  }
  
  
  vsub <- which (dataSubset$hidx %in% vsampleP)
  
  
  
  mInfMatrix <- matrix(nrow=NSamples,ncol=NROW(vsampleP) )
  vsampleRegisterIdx <- vector()
  hCnt <- 0
  for (h in vsampleP)
  {
    print(h)
    hCnt <- hCnt + 1
    vsampleRegisterIdx[hCnt] <- tail(drawS$RegistrarIdx[h,,],n=1) ##Save the Registry Idx So we can refer back to which Fish This belongs to 
    
    vphi_0_sub <- tail(drawS$phi_0[h,,],n=NSamples)
    vphi_max_sub <- tail(drawS$phi_max[h,,],n=NSamples)
    vgamma_sub <- tail(drawS$gamma[h,,],n=NSamples)
    vlambda_sub <- tail(drawS$lambda[h,,],n=NSamples)
    vtau_sub   <- tail(drawS$tau[h,,],n=NSamples)
    vsigma_sub <- tail(drawS$sigma[h,,],n=NSamples)
    valpha_sub <- tail(drawS$alpha[h,,],n=NSamples)
    
    
    
    vPP <- which (dataSubset$hidx == h)
    #pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_",strGroupID[groupID],"_Sigmoid_",pp,".pdf",sep="")) 
    #pdf(file= paste(strPlotExportPath,"/stat/stat_InfMeasure_EyeVsDistance5mm_",strGroupID[groupID],"_SigExp_",h,".pdf",sep="")) 
    #par(mar = c(5,5,2,5))
    #plot(dataSubset$distP[vPP],dataSubset$phi[vPP],pch=19,xlim=c(0,max(DistRange)),ylim=c(0,90),main=paste(strGroupID[groupID],h), 
    #     bg=colourP[2],col=colourP[1],
    #     cex=0.5,
    #     xlab="Distance From Prey",ylab=expression(paste(Phi)) )  
    ##Plot The Fitted Function
    #vY  <-  mean(valpha_sub)*exp(mean(vgamma_sub)*( mean(vtau_sub)-DistRange) )+ mean(vphi_0_sub)   +  ( mean(vphi_max_sub) - mean(vphi_0_sub)  )/(1+exp( -( mean(vlambda_sub)   *( mean(vtau_sub) -DistRange )   ) ) ) 
    #lines( DistRange ,vY,type="l",col=colourR[3],lwd=3)
    
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
      
      
      vPPhiPerX <- InfoCalc(DistRange,Ulist = lPlist[[i]]) ##Get MutInf Per X

      ##Plot The Information Content Of The fitted Function ##
      #par(new=T) 
      #  plot(DistRange,vPPhiPerX,ylim=c(0,2.5),xlim=c(0,max(DistRange) ),axes=F,type="p",pch=19,col=colourP[4], xlab=NA,
      #       ylab=NA,sub=paste("x requires ", round(100*log2(NROW(DistRange)) )/100,"bits"  ) )
      # lines(DistRange,rev(cumsum(rev(vPPhiPerX))),ylim=c(0,2.5),xlim=c(0,max(DistRange) ),type="l",col=colourR[4], xlab=NA, ylab=NA )
      
      mInfMatrix[i,hCnt] <- sum(vPPhiPerX) ## Mutual Inf Of X and Phi
    } 
    
    #axis(side=4)
    #mtext(side = 4, line = 3, 'Mutual information (bits)')
    #dev.off()
  }##For Each Hunt Event   
  
  
  ##Next Is integrate Over sampled points, and Make Ii hidx matrix 
  lInfoMatStruct <- list(infoMatrix=mInfMatrix,vsampleRegisterIdx=vsampleRegisterIdx,vsamplePSeqIdx=vsampleP)
  return(lInfoMatStruct)
}



## Sample Matrices Of Information / Retursn Struct Containing Mat and Id vectors ##
lInfStructLL <- calcInfoOfHuntEvent(drawLL,dataLL,groupID=2)
lInfStructNL <- calcInfoOfHuntEvent(drawNL,dataNL,groupID=3)
lInfStructDL <- calcInfoOfHuntEvent(drawDL,dataDL,groupID=1)

mInfMatrixLL <- lInfStructLL$infoMatrix
mInfMatrixNL <- lInfStructNL$infoMatrix
mInfMatrixDL <- lInfStructDL$infoMatrix


X11()
hist(mInfMatrixLL,col=colourH[2],xlim=c(0,3),breaks = seq(0,3,1/20))
X11()
hist(mInfMatrixDL,col=colourH[1],xlim=c(0,3),breaks = seq(0,3,1/20))
X11()
hist(mInfMatrixNL,col=colourH[3],xlim=c(0,3),breaks = seq(0,3,1/20))


X11()
hist(colMeans(mInfMatrixLL),col=colourH[2],xlim=c(0,3),breaks = seq(0,3,1/20),main="mean Inf Per Event ")
X11()
hist(colMeans(mInfMatrixDL),col=colourH[1],xlim=c(0,3),breaks = seq(0,3,1/20))
X11()
hist(colMeans(mInfMatrixNL),col=colourH[3],xlim=c(0,3),breaks = seq(0,3,1/20))


### Plot CDF ###
## Match the N 
pdf(file= paste(strPlotExportPath,"/stat/stat_InfSigmoidExp_EyeVsDistance_CDF_4.pdf",sep=""))
subset_mInfMatrixLL <- mInfMatrixLL[,sample(1:58,58)]
plot(ecdf(mInfMatrixDL),col=colourH[1],main="Information In Eye Vergence CDF",
     xlab="Information (bits) ",lty=1,lwd=2,xlim=c(0,2.5))
plot(ecdf(mInfMatrixLL),col=colourH[2],add=T,lty=2,lwd=2)
plot(ecdf(mInfMatrixNL),col=colourH[3],add=T,lty=3,lwd=2)
legend("topleft",legend=paste(c("DL n=","LL n=","NL n="),c(NCOL(mInfMatrixDL),NCOL(mInfMatrixLL) ,NCOL(mInfMatrixNL) ) ) 
       ,col=colourH,lty=c(1,2,3),lwd=2)
dev.off()


##Plot Density 
dLLphi<-density(mInfMatrixLL)
dDLphi<-density(mInfMatrixDL)
dNLphi<-density(mInfMatrixNL)

#X11()
pdf(file= paste(strPlotExportPath,"/stat/stat_InfSigmoidExp_EyeVsDistance_Density_4.pdf",sep=""))
plot(dLLphi,col=colourH[2],type="l",lwd=3,lty=2,
     ylim=c(0,2),main="Mutual Information Distance to Eye Vergence ",
     xlab=expression(paste(" (bits)" ) ) )
lines(dNLphi,col=colourH[3],lwd=3,lty=3)
lines(dDLphi,col=colourH[1],lwd=3,lty=1)
legend("topleft",legend=paste(c("DL n=","LL n=","NL n="),c(NCOL(mInfMatrixDL),NCOL(subset_mInfMatrixLL) ,NCOL(mInfMatrixNL) ) ) 
       ,col=colourH,lty=c(1,2,3),lwd=2)
dev.off()

save(lInfStructLL,lInfStructDL,lInfStructNL,drawLL,drawDL,drawNL,file=paste(strDataExportDir,"/stat_infoMat_EyeVergenceVsDistance_sigmoidFit5mm-5bit_4.RData",sep=""))      






## Entropy ##
#if you use this package please cite: Jean Hausser and Korbinian Strimmer. 2009. Entropy inference
#and the James-Stein estimator, with application to nonlinear gene association networks.  J. Mach.
#Learn.  Res.10:  1469-1484.  Available online fromhttp://jmlr.csail.mit.edu/papers/v10/  hausser09a.html.
library("entropy")

KL.empirical(binNLphi,binLLphi,unit="log2")


mp <- barplot(c(entropy(binDLphi,unit="log2"),entropy(binLLphi,unit="log2"),entropy(binNLphi,unit="log2")),col=colourR,xlab="Group ",
              main=expression(paste("Entropy of eye vergence angle during hunt" ) ),sub=expression(paste("Resolution: ", Phi,": 1 deg" ) ) ,
              ylab="(bits)", ylim=c(0,8) ) #legend.text = strGroupID
text(mp,y=c(-0.3,-0.3,-0.3),labels=strGroupID , xpd = TRUE, col = "black")


## Mutual Information ##
binLLphiVsDist <- discretize2d(dataLL$phi, dataLL$distP, numBins1=80, numBins2=10,r1=c(0,80),r2=c(0.5,2))
binNLphiVsDist <- discretize2d(dataNL$phi, dataNL$distP, numBins1=80, numBins2=10,r1=c(0,80),r2=c(0.5,2))
binDLphiVsDist <- discretize2d(dataDL$phi, dataDL$distP, numBins1=80, numBins2=10,r1=c(0,80),r2=c(0.5,2))

## Joint Entropy
H12_LL = entropy(binLLphiVsDist,unit="log2" )
H12_NL = entropy(binNLphiVsDist ,unit="log2")
H12_DL = entropy(binDLphiVsDist,unit="log2" )
##Plot Mutual Information
mp <- barplot(c(H12_DL,H12_LL,H12_NL),col=colourR,xlab="Group ",
              main=expression(paste("Joint entropy prey distance and eye vergence" ) ),sub=expression(paste("Resolution: ", Phi,": 1 deg, Distance : 0.5mm" ) ) ,
              ylab="(bits)", ylim=c(0,10) ) #legend.text = strGroupID
text(mp,y=c(-0.3,-0.3,-0.3),labels=strGroupID , xpd = TRUE, col = "black")

# mutual information
mi.empirical(binLLphiVsDist,unit="log2" ) # //Mutual Informations
mi.empirical(binNLphiVsDist,unit="log2" ) # approximately zero
mi.empirical(binDLphiVsDist,unit="log2" ) # approximately zero

# another way to compute mutual information
# compute marginal entropies
H12_LL = entropy(binLLphiVsDist,unit="log2" )
H1_LL = entropy(rowSums(binLLphiVsDist),unit="log2")
H2_LL = entropy(colSums(binLLphiVsDist),unit="log2")
H1_LL+H2_LL-H12_LL # mutual entropy

H12_NL = entropy(binNLphiVsDist,unit="log2" )
H1_NL = entropy(rowSums(binNLphiVsDist),unit="log2")
H2_NL = entropy(colSums(binNLphiVsDist),unit="log2")
H1_NL+H2_NL-H12_NL # mutual entropy


H12_DL = entropy(binDLphiVsDist,unit="log2" )
H1_DL = entropy(rowSums(binDLphiVsDist),unit="log2")
H2_DL = entropy(colSums(binDLphiVsDist),unit="log2")
H1_DL+H2_DL-H12_DL # mutual entropy



# sample from continuous uniform distribution
x1 = runif(10000)
X11()
hist(x1, xlim=c(0,1), freq=FALSE)
# discretize into 10 categories
y1 = discretize(x1, numBins=10, r=c(0,1))
# compute entropy from counts
entropy(y1) # empirical estimate near theoretical maximu
log(10)

#discretize2d(x1, x2, numBins1=10, numBins2=10)

# sample from a non-uniform distribution
X11()
hist(x2, xlim=c(0,1), freq=FALSE)
# discretize into 10 categories and estimate entropy

#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_Rate_lambda_LL_E.pdf",sep=""))
X11()
hist(drawLL$lambda[,,],col=colourR[2],main="Lambda")
hist(drawNL$lambda[,,],col=colourR[3],add=TRUE)
hist(drawDL$lambda[,,],col=colourR[1],add=TRUE)

X11()
hist(drawLL$gamma[,,],col=colourR[2],main="gamma")
hist(drawNL$gamma[,,],col=colourR[3],add=TRUE)
hist(drawDL$gamma[,,],col=colourR[1],add=TRUE)


X11()
hist(drawLL$phi_0[,,],col=colourR[2],main="Phi 0")
hist(drawNL$phi_0[,,],col=colourR[3],add=TRUE)
hist(drawDL$phi_0[,,],col=colourR[1],add=TRUE)


X11()
hist(drawLL$phi_max[,,],col=colourR[2],main="Phi Max")
hist(drawNL$phi_max[,,],col=colourR[3],add=TRUE)
hist(drawDL$phi_max[,,],col=colourR[1],add=TRUE)
#plot(drawLL$phi_max[3,,])

X11()
hist(drawLL$sigma[,1,,],main="LL")



X11()
#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_StartEnd_u0_NL_E.pdf",sep=""))
hist(drawLL$u1[,,],breaks=50,xlim=c(0,7),col=colourH[2])
hist(drawLL$u0[,,],breaks=50,xlim=c(0,7),add=TRUE,col=colourH[2])

#dev.off()
########################
## NL ###
mNL=jags.model(file="model.tmp",data=dataNL);
drawNL=jags.samples(mNL,steps,thin=thin,variable.names=varnames)


X11()
#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_NL_F.pdf",sep=""))
plotGCRes(drawNL,dataNL,groupID=3)
#dev.off()

## Plot the infered function NL

X11()
#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_Rate_lambda_NL_E.pdf",sep=""))
hist(drawNL$lambda[1,,1],main="NL")
#dev.off()

X11()
#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_StartEnd_u0_NL_E.pdf",sep=""))
hist(drawNL$u1[1,,1],breaks=50,xlim=c(0,7),col="red")
hist(drawNL$u0[1,,1],breaks=50,xlim=c(0,7),add=TRUE,col="red")
#dev.off()





pdf(file= paste(strPlotExportPath,"/stat/stat_densityolinregressionslope.pdf",sep=""))
plot(dDLb,col=colourH[1],xlim=c(0.5,1.2),lwd=3,lty=1,ylim=c(0,20),main="Density Inference of Turn-To-Prey Slope ")
lines(dLLb,col=colourH[2],xlim=c(0.5,1.2),lwd=3,lty=2)

lines(dNLb,col=colourH[3],xlim=c(0.5,1.2),lwd=3,lty=3)
legend("topleft",legend=paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       ,fill=colourL,lty=c(1,2,3))
dev.off()

### PLot Scatter with regression lines with Conf intervals##
#X11()

pdf(file= paste(strPlotExportPath,"/stat/stat_TurnToPrey_LinearRegression.pdf",sep=""))
plot(lFirstBoutPoints[["DL"]][,1], lFirstBoutPoints[["DL"]][,2],
     main=paste("Turn Size Vs Bearing To Prey ", sep=""),
     xlab="Bearing To Prey prior to Bout",ylab="Bearing Change After Bout",xlim=c(-100,100),
     ylim=c(-100,100),
     col=colourP[1] ,pch=pchL[1]) ##boutSeq The order In Which The Occurred Coloured from Dark To Lighter
##Draw 0 Vertical Line
segments(0,-90,0,90); segments(-90,0,90,0); segments(-90,-90,90,90,lwd=1,lty=2);
#text(lFirstBoutPoints[["DL"]][,1]+2,lFirstBoutPoints[["DL"]][,2]+5,labels=lFirstBoutPoints[["DL"]][,3],cex=0.8,col="darkblue")
abline(lm(lFirstBoutPoints[["DL"]][,2] ~ lFirstBoutPoints[["DL"]][,1]),col=colourH[4],lwd=1.0) ##Fit Line / Regression
abline(a=muDLa,b=muDLb,col=colourH[1],lwd=1.5) ##Fit Line / Regression
abline(a=quantile(drawDL$beta[,(steps-ind):steps,1][1,])[2],b=quantile(drawDL$beta[,(steps-ind):steps,1][2,])[2],col=colourR[1],lwd=4.0) ##Fit Line / Regression
abline(a=quantile(drawDL$beta[,(steps-ind):steps,1][1,])[3],b=quantile(drawDL$beta[,(steps-ind):steps,1][2,])[3],col=colourR[1],lwd=4.0) ##Fit Line / Regression

#abline( lsfit(lFirstBoutPoints[["DL"]][,2], lFirstBoutPoints[["DL"]][,1] ) ,col=colourH[1],lwd=2.0)
##LL
points(lFirstBoutPoints[["LL"]][,1], lFirstBoutPoints[["LL"]][,2],pch=pchL[2],col=colourP[2])
#text(lFirstBoutPoints[["LL"]][,1]+2,lFirstBoutPoints[["LL"]][,2]+5,labels=lFirstBoutPoints[["LL"]][,3],cex=0.8,col="darkgreen")
abline(lm(lFirstBoutPoints[["LL"]][,2] ~ lFirstBoutPoints[["LL"]][,1]),col=colourH[4],lwd=1.0)
abline(a=muLLa,b=muLLb,col=colourH[2],lwd=1.5) ##Fit Line / Regression
abline(a=quantile(drawLL$beta[,(steps-ind):steps,1][1,])[2],b=quantile(drawLL$beta[,(steps-ind):steps,1][2,])[2],col=colourR[2],lwd=4.0) ##Fit Line / Regression
abline(a=quantile(drawLL$beta[,(steps-ind):steps,1][1,])[3],b=quantile(drawLL$beta[,(steps-ind):steps,1][2,])[3],col=colourR[2],lwd=4.0) ##Fit Line / Regression

#abline(lsfit(lFirstBoutPoints[["LL"]][,2], lFirstBoutPoints[["LL"]][,1] ) ,col=colourH[2],lwd=2.0)
##NL
points(lFirstBoutPoints[["NL"]][,1], lFirstBoutPoints[["NL"]][,2],pch=pchL[3],col=colourP[3])
#text(lFirstBoutPoints[["NL"]][,1]+2,lFirstBoutPoints[["NL"]][,2]+5,labels=lFirstBoutPoints[["NL"]][,3],cex=0.8,col="darkred")
abline(lm(lFirstBoutPoints[["NL"]][,2] ~ lFirstBoutPoints[["NL"]][,1]),col=colourH[4],lwd=1.0)
abline(a=muNLa,b=muNLb,col=colourH[3],lwd=1.5) ##Fit Line / Regression
abline(a=quantile(drawNL$beta[,(steps-ind):steps,1][1,])[2],b=quantile(drawNL$beta[,(steps-ind):steps,1][2,])[2],col=colourR[3],lwd=4.0) ##Fit Line / Regression
abline(a=quantile(drawNL$beta[,(steps-ind):steps,1][1,])[3],b=quantile(drawNL$beta[,(steps-ind):steps,1][2,])[3],col=colourR[3],lwd=4.0) ##Fit Line / Regression
#abline( lsfit(lFirstBoutPoints[["NL"]][,2], lFirstBoutPoints[["NL"]][,1] ) ,col=colourH[3],lwd=2.0)
legend("topleft",legend=paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       , pch=pchL,col=colourL)

dev.off()





##Plot Densities Summary

sampNL <- coda.samples(mNL,                      variable.names=c("beta","sigma"),                      n.iter=20000, progress.bar="none")
sampDL <- coda.samples(mDL,                      variable.names=c("beta","sigma"),                      n.iter=20000, progress.bar="none")
X11()
plot(sampLL)
X11()
plot(sampNL)
X11()
plot(sampDL,main="DL")











## N : vector of number of points in Hunt Event
modelExpInd  <- "model{
##Prior

# Prior Sigma On Eye Angle when  In Or Out of hunt region 
for(i in 1:max(hidx)) {
for(j in 1:2){
#inv.var[j] ~ dgamma(0.01, 0.01)  ##Prior for inverse variance
sigma[i,j] ~ dgamma(0.01, 0.01) ##Draw 
}
}

# Likelihood / Separate for each Hunt Event
for(i in 1:N){
phi_0[hidx[i]] ~ dnorm(10,2) # Idle Eye Position
phi_max[hidx[i]] ~ dnorm(15,5) # Max Eye Vergence Angle
lambda[hidx[i]] ~ dgamma(1, 1) # RiseRate of Eye Vs Prey Distance
limDist[hidx[i]] <- max(distMax)
u1[hidx[i]] ~ dunif(0, limDist[hidx[i]]) ## End Hunt Distance - Close to prey
u0[hidx[i]] ~ dunif(u1, limDist[hidx[i]]) ##Start Hunt Distance -Far 


##Make indicator if hunt event is within sampled Range 
#if (u1[hidx[i]] < distP[i]  & distP[i] < u0) 
s[hidx[i],i] <- step( distP[i]-u1[hidx[i]])*step(u0[ hidx[i] ]-distP[i]  ) 

phi_hat[hidx[i],i] <- phi_0[hidx[i]] + s[hidx[i],i] * phi_max[hidx[i]]* (1-exp(-lambda[ hidx[i] ]*(distMax[i] - distP[i] ) )) 
phi[hidx[i],i] ~ dnorm( phi_hat[hidx[i],i], sigma[s[hidx[i],i]+1] ) ##choose sigma 

}"


# 
# 
# # Step 1
# # Reading the the data
# FathersLoveData <- read.csv("C:/Users/zzo1/Dropbox/PBnRTutorial/lovedata.csv")
# # Visual check data of the data
# head(FathersLoveData) 
# # the unnamed column shows the row number: each participant's data corresponds to one row
# # Y1-Y4 are the love scores for each person
# # 'Positivity' shows the positivity categories: 1: low, 2: medium, 3: high
# # X1 takes value 1 for low positivity, otherwise 0
# # X2 takes value 1 for high positivity, otherwise 0
# 
# # Checking the number of rows (i.e., 2nd dimension) of `FathersLoveData' data 
# N <- dim(FathersLoveData)[1] # count number of rows to get number of subjects
# # Creating a matrix with only the love scores: Y1-Y4 (columns 1-4):
# # We `unlist' all N rows of selected columns 1:4 from the data set, 
# # then we transform these values into numeric entries of a matrix
# data <- matrix(as.numeric(unlist(FathersLoveData[,1:4])), nrow = N)
# # Creating a variable that saves the number of time points
# nrT <- 4
# # Saving X1 and X2 as separate variables (same unlisting etc. as explained above)
# grouping <- matrix(as.numeric(unlist(FathersLoveData[,6:7])), nrow = N)
# X1 <- grouping[,1] # 1 when person had low positivity before baby
# X2 <- grouping[,2] # 1 when person had high positivity before baby
# # Creating a time vector for the measurement waves
# time <- c(-3, 3, 9, 36) # time vector (T) based on the time of the measurements
# # Now we have all the data needed to be passed to JAGS
# # Creating a list of all the variables that we created above
# jagsData <- list("Y"=data,"X1"=X1,"X2"=X2,"N"=N,"nrT"=nrT,"time"=time)
# 
# # Step 2
# LinearGrowthCurve = cat("
#                         model {
#                         # Starting loop over participants
#                         for (i in 1:N) { 
#                         # Starting loop over measurement occasions
#                         eps[i,1] <- Y[i, 1] - (betas[i,1] + betas[i,2]*time[1])
#                         Y[i, 1] ~ dnorm(betas[i,1] + betas[i,2]*time[1], precAC)
#                         for (t in 2:nrT) {  
#                         # The likelihood function, corresponding to Equation 1:
#                         
#                         eps[i,t] <- Y[i,t] - (betas[i,1] + betas[i,2]*time[t] + acorr*eps[i,t-1])
#                         Y[i, t] ~ dnorm(betas[i,1] + betas[i,2]*time[t] + acorr*eps[i,t-1], precAC)} 
#                         # end of loop for the observations
#                         
#                         # Describing the level-2 bivariate distribution of intercepts and slopes
#                         betas[i,1:2] ~ dmnorm(Level2MeanVector[i,1:2], interpersonPrecisionMatrix[1:2,1:2])
#                         # The mean of the intercept is modeled as a function of positivity group membership
#                         Level2MeanVector[i,1] <- MedPInt + betaLowPInt*X1[i] + betaHighPInt*X2[i]
#                         Level2MeanVector[i,2] <- MedPSlope + betaLowPSlope*X1[i] + betaHighPSlope*X2[i]
#                         } # end of loop for persons
#                         # Specifying priors  distributions
#                         MedPInt ~ dnorm(0, 0.01)
#                         MedPSlope ~ dnorm(0, 0.01)
#                         betaLowPInt ~ dnorm(0, 0.01)
#                         betaHighPInt ~ dnorm(0, 0.01)
#                         betaLowPSlope ~ dnorm(0, 0.01)
#                         betaHighPSlope ~ dnorm(0, 0.01)
#                         
#                         sd1 ~ dunif(0, 100)
#                         precAC <- 1/pow(sd1,2)*(1-acorr*acorr)
#                         acorr ~ dunif(-1,1)
#                         
#                         sdIntercept  ~ dunif(0, 100)
#                         sdSlope  ~ dunif(0, 100)
#                         corrIntSlope ~ dunif(-1, 1)
#                         # Transforming model parameters
#                         ## Defining the elements of the level-2 covariance matrix
#                         interpersonCovMatrix[1,1] <- sdIntercept * sdIntercept
#                         interpersonCovMatrix[2,2] <- sdSlope * sdSlope
#                         interpersonCovMatrix[1,2] <- corrIntSlope * sdIntercept* sdSlope
#                         interpersonCovMatrix[2,1] <- interpersonCovMatrix[1,2]
#                         ## Taking the inverse of the covariance to get the precision matrix
#                         interpersonPrecisionMatrix <- inverse(interpersonCovMatrix)
#                         ## Creating a variables representing
#                         ### low positivity intercept
#                         LowPInt <- MedPInt + betaLowPInt 
#                         ### high positivity intercept
#                         HighPInt <- MedPInt + betaHighPInt 
#                         ### low positivity slope
#                         LowPSlope <- MedPSlope + betaLowPSlope
#                         ### high positivity slope
#                         HighPSlope <- MedPSlope + betaHighPSlope
#                         ### contrasts terms between high-low, medium-low, high-medium intercepts and slopes
#                         HighLowPInt <- HighPInt - LowPInt
#                         MedLowPInt <- MedPInt - LowPInt
#                         HighMedPInt <- HighPInt - MedPInt
#                         HighLowPSlope <- HighPSlope- LowPSlope
#                         MedLowPSlope <- MedPSlope - LowPSlope
#                         HighMedPSlope <- HighPSlope - MedPSlope
#                         }
#                         ",file = "GCM.txt")
# 
# 
# # Step 3
# # Collecting the model parameters of interest
# parameters  <- c("MedPSlope","betaLowPInt",
#                  "betaHighPInt","betaLowPSlope", 
#                  "betaHighPSlope", "MedPInt", 
#                  "sdIntercept", "sdSlope", 
#                  "corrIntSlope", "betas",
#                  "LowPInt","HighPInt","LowPSlope", "HighPSlope",
#                  "HighLowPInt","HighMedPInt","MedLowPInt",
#                  "HighLowPSlope","HighMedPSlope","MedLowPSlope",
#                  "acorr", "sd1")
# # Sampler settings
# adaptation  <- 2000 # Number of steps to "tune" the samplers
# chains  <- 6    # Re-start the exploration "chains" number of times
# #  with different starting values
# burnin  <- 1000 # Number of steps to get rid of the influence of initial values
# # Define the number of samples drawn from the posterior in each chain
# thinning <- 20
# postSamples <- 60000
# nrOfIter <- ceiling((postSamples * thinning)/chains)
# 
# 
# fixedinits<- list(list(.RNG.seed=5,.RNG.name="base::Mersenne-Twister"),list(.RNG.seed=6,.RNG.name="base::Mersenne-Twister"),list(.RNG.seed=7,.RNG.name="base::Mersenne-Twister"),list(.RNG.seed=8,.RNG.name="base::Mersenne-Twister"),list(.RNG.seed=9,.RNG.name="base::Mersenne-Twister"),list(.RNG.seed=10,.RNG.name="base::Mersenne-Twister"))
# 
# # Step 4
# # loading the rjags package
# library(rjags)            
# # creating JAGS model object
# jagsModel<-jags.model("GCM.txt",data=jagsData,n.chains=chains,n.adapt=adaptation,inits=fixedinits)
# # running burn-in iterations
# update(jagsModel,n.iter=burnin)
# # drawing posterior samples
# codaSamples<-coda.samples(jagsModel,variable.names=parameters,thin = thinning, n.iter=nrOfIter,seed=5)
# 
# source("C:/Users/zzo1/Dropbox/PBnRTutorial/posteriorSummaryStats.R")
# # Part 1: Check convergence
# resulttable <- summarizePost(codaSamples)
# saveNonConverged <- resulttable[resulttable$RHAT>1.1,]
# if (nrow(saveNonConverged) == 0){
#   print("Convergence criterion was met for every parameter.")
# }else{ 
#   print("Not converged parameter(s):")
#   show(saveNonConverged)
# }
# # Part 2: Display summary statistics for selected parameters (regexp)
# show(summarizePost(codaSamples, filters =  c("^Med","^Low","^High","^sd","^corr", "sd1", "acorr"))) 
# 
# save.image(sprintf("GCMPBnR%s.Rdata", Sys.Date()))
# 


