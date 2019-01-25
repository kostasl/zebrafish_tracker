
##THe Growth Model : Carlin and Gelfand (1991) present a nonconjugate Bayesian analysis of the following data set from Ratkowsky (1983):
modelGCSigmoidInd  <- "model
{
  
  for( i in 1 : N ) {
    phi_hat[ hidx[i],i] <-  phi_0[ hidx[i] ] +   (phi_max[hidx[i]] - phi_0[ hidx[i] ])/( 1 + exp( -lambda[ hidx[i] ]*( ( tau[ hidx[i] ] - distP[i]    ) ) ) )
  
    ###OUT Set Region Of Exp Growth Model ##
    # s[hidx[i],i] <- step( distP[i] - u1[hidx[i]] )*step( tau[ hidx[i] ] -distP[i] )     # step( phi_max[hidx[i]] - phi_0[hidx[i]] ) #step(u0[ hidx[i] ] - distP[i]  )  
  
    ## Define Exp Growth Model 
    phi_exp[ hidx[i],i] <- alpha[hidx[i]]*exp( gamma[ hidx[i] ]* ( tau[ hidx[i] ] -  distP[i]))
  
    ### Conditionally Include the exp Model
    phi[i] ~ dnorm(phi_exp[ hidx[i],i]  + phi_hat[ hidx[i],i]  , var_inv[hidx[i]] ) #s[hidx[i],i]+1 
  
  }
  
  
  ## Priors
  limDist <- max(distMax)
  
  
  for(i in 1:max(hidx) ) { 
    phi_max[i] ~ dnorm(65,1e-3)T(0,150) ##I(0,100) # Max Eye Vergence Angle
    phi_0[i] ~ dnorm(0.01, 1e-3)T(0,60)  # Idle Eye Position
    lambda[i] ~ dgamma(100, 1) #dnorm(100.0, 1e-3)T(0,) # RiseRate of Eye Vs Prey Distance Sigmoid
    gamma[i] ~ dgamma(1, 0.5) #dnorm(0.5, 1e-3)I(0,)  # RiseRate of Eye Vs Prey Distance After Sig Rise dunif(0.5, 0.000001)
    alpha[i] ~ dunif(1,3)
    tau[i] ~ dnorm(distMax[i], 1e-1) ##inflexion point, sample from where furthest point of Hunt event is found
    var_inv[i] ~ dgamma(0.001, 0.001) ##Draw   ##Precision
  
    sigma[i] <- 1 / sqrt(var_inv[ i])    
  
  }
  
  
  
  
  
}"
  
  ##Init Vals List -For Run Jags
  # A list of 8 randomly generated starting values for m:
  
  initfunct <- function(nchains,N)
  {
    initlist <- replicate(nchains,list(phi_0=c(rnorm(N,10,5)),
                                       phi_max=rnorm(N,40,5),
                                       lambda=rgamma(N,100,1), ## Sigmoid Rise Rate
                                       gamma=rgamma(N,1,1),
                                       alpha=runif(N,1,3),
                                       tau=rnorm(N,2,0.5)  ),
                          simplify=FALSE)
    
    return(initlist)
  }
  
  ## Implements Regressed Function, evaluated at vX points given the array of parameters
  eyeVregressor <- function(params,vX)
  {
    return(vY  <-  params$ealpha*exp(params$egamma*(params$etau-vX) )+ params$ephi0   +  (params$ephimax -params$ephi0  )/(1+exp( -(params$elambda   *(params$etau -vX )   ) ) ) )
  }
  
  ## Plot the Eye Vs Distance data points and the regression variations ##
  ## Returns Mean Square Error R^2 ## 
  plotEyeGCFit <- function(pp,strGroup,dataSubset,drawS)
  {
    vRegIdx <- unique(dataSubset$RegistrarIdx) ##Get Vector Of RegIdx That Associate with the sample Sequence
    vSqError <- vector()
    
    vX  <- seq(0,5,by=0.01)
    vPP <- which (dataSubset$hidx == pp)
    
    etau    <- (tail(drawS$tau[pp,,],n=100))
    ephimax <- (tail(drawS$phi_max[pp,,],n=100))
    ephi0   <- (tail(drawS$phi_0[pp,,],n=100))
    elambda <- (tail(drawS$lambda[pp,,],n=100))
    egamma  <- (tail(drawS$gamma[pp,,],n=100))
    ealpha  <- (tail(drawS$alpha[pp,,],n=100))
    
    plot(dataSubset$distP[vPP],dataSubset$phi[vPP],pch=19,xlim=c(0,5),ylim=c(0,85),
         main=paste(strGroup,pp," (",vRegIdx[pp],")"), 
         bg=colourP[2],col=colourP[1],cex=0.5)
    
    
    
    ## Draw The 100 Variotons before the fit converged      
    for (k in 1:NROW(etau) )
    {
      params <- list(etau    = etau[k],
                     ephimax = ephimax[k],
                     ephi0   = ephi0[k],
                     elambda = elambda[k],
                     egamma  = egamma[k],
                     ealpha  = ealpha[k]
      )
      #vY  <-  ealpha[k]*exp(egamma[k]*(etau[k]-vX) )+ ephi0[k]   +  (ephimax[k] -ephi0[k]  )/(1+exp( -(elambda[k]   *(etau[k] -vX )   ) ) ) 
      vY          <- eyeVregressor(params,vX)
      ## Obtain Regressor Y at Data points X, and Measure Error To Actual Eye Vergence Data point at X
      vRError     <- sum( (eyeVregressor(params,dataSubset$distP[vPP])-dataSubset$phi[vPP])  ^2 ) / NROW(dataSubset$phi[vPP])
      vSqError[k] <- vRError
      #vY_l  <- quantile(drawS$phi_0[pp,,])[1]   - ( quantile (drawS$lambda[pp])[1] )*((( quantile(drawS$gamma[pp,,])[1] )^( quantile(drawS$u0[pp])[1] - (vX) ) ) ) #
      #vY_u  <- quantile(drawS$phi_0[pp,,])[5]   - (quantile (drawS$lambda[pp,,])[5])*((( quantile(drawS$gamma[pp,,])[5] )^( quantile(drawS$u0[pp,,])[5] - (vX) ) ) ) #
      #      #points(dataSubset$distP[vPP],dataSubset$phi[vPP],pch=19,xlim=c(0,5),ylim=c(-85,85),main=paste("L",pp), bg=colourP[2],col=colourP[1],cex=0.5)
      lines( vX ,vY,type="l",col=colourR[3],lwd=1)
      #        lines( vX ,vY_u,type="l",col=colourR[4],lwd=1)
    }
    
    return(vSqError) 
  }
  
  
  
  ## Plot average regressed function ##
  ## plot( exp(0.1*(-vx+80))+  10 + (90-10)/(1+exp(-100*(60-vx) ))   ,ylim=c(0,400))
  plotGCSig <- function (drawS,dataSubset,n=NA,groupID,nPlotPoints=50){
    
    
    bPlotIndividualEvents <- FALSE
    ## compute 2D kernel density, see MASS book, pp. 130-131
    max_x <- 7
    nlevels <- 12
    
    if (is.na(n))
      n <- NROW(unique(dataSubset$hidx))
    
    vRegIdx <- unique(dataSubset$RegistrarIdx) ##Get Vector Of RegIdx That Associate with the sample Sequence
    
    vsampleP <- sample(unique(dataSubset$hidx),n)
    vsub <- which (dataSubset$hidx %in% vsampleP)
    
    z <- kde2d(dataSubset$distP, dataSubset$phi, n=80)
    
    # X11()
    
    plot(dataSubset$distP[vsub],dataSubset$phi[vsub],pch=21,xlim=c(0,max_x),ylim=c(0,80),
         main=paste("Model Fit : Eye Vergence Vs Distance Data ",strGroupID[groupID]),
         ylab=expression(paste("Eye Vergence ",Phi," (degrees)") ),
         xlab=expression(paste("Distance from Prey (mm)") ),
         bg=colourR[groupID],col="#FFFFFFAA",cex=0.5)
    #points(dataSubset$distToPrey[vsub],dataSubset$vAngle[vsub],pch=21,xlim=c(0,5),ylim=c(0,80),main="LL", bg=colourP[4],col=colourP[1],cex=0.5)
    contour(z, drawlabels=FALSE, nlevels=nlevels,add=TRUE)
    ## Plot The Mean Curve of the selected subset of curves
    vX  <- seq(0,max_x,by=0.01)##max(drawS$u0[vsampleP])
    ##  (phi_max[hidx[i]] - phi_0[hidx[i]] )/(1-exp(-lambda[ hidx[i] ]*(u0[ hidx[i] ]  - distP[i] ) ))
    
    etau    <- (tail(drawS$tau[vsampleP,1,1],n=nPlotPoints))
    ephimax <- (tail(drawS$phi_max[vsampleP,1,1],n=nPlotPoints))
    ephi0   <- (tail(drawS$phi_0[vsampleP,1,1],n=nPlotPoints))
    elambda <- (tail(drawS$lambda[vsampleP,1,1],n=nPlotPoints))
    egamma  <- (tail(drawS$gamma[vsampleP,1,1],n=nPlotPoints))
    ealpha  <- (tail(drawS$alpha[vsampleP,1,1],n=nPlotPoints))
    ## Draw The 100 Variotons before the fit converged      
    for (k in 1:NROW(etau) )
    {
      params <- list(etau    = etau[k],
                     ephimax = ephimax[k],
                     ephi0   = ephi0[k],
                     elambda = elambda[k],
                     egamma  = egamma[k],
                     ealpha  = ealpha[k]
      )
      #vY  <-  ealpha[k]*exp(egamma[k]*(etau[k]-vX) )+ ephi0[k]   +  (ephimax[k] -ephi0[k]  )/(1+exp( -(elambda[k]   *(etau[k] -vX )   ) ) ) 
      vY          <- eyeVregressor(params,vX)
      ## Obtain Regressor Y at Data points X, and Measure Error To Actual Eye Vergence Data point at X
      lines( vX ,vY,type="l",col=colourL[3],lwd=1)
      #        lines( vX ,vY_u,type="l",col=colourR[4],lwd=1)
    }
    dev.off()
    # 
    # 
    # vY  <-    ealpha*exp(egamma*(etau-vX) ) + ephi0 + (ephimax -ephi0  )/(1+exp( -(elambda )  *(etau -(vX)   ) ) ) 
    # 
    # etau <- quantile((drawS$tau[vsampleP]))[2]
    # vY_l  <-  ealpha*exp(egamma*(etau-vX) )+  ephi0   +  (ephimax -ephi0  )/(1+exp( -(elambda )  *(etau -(vX)   ) ) ) 
    # 
    # etau <- quantile((drawS$tau[vsampleP]))[4]
    # vY_u  <-  ealpha*exp(egamma*(etau-vX) )+  ephi0   +  (ephimax -ephi0  )/(1+exp( -(elambda )  *(etau -(vX)   ) ) ) 
    # 
    # 
    # #vY_u <-  quantile(drawS$phi_0[vsampleP])[4]-(quantile(drawS$lambda[vsampleP])[4])*((quantile(drawS$gamma[vsampleP])[4]^( quantile(drawS$u0[vsampleP])[4] - (vX) ) ) )
    # #vY_l <-  quantile(drawS$phi_0[vsampleP])[2]-(quantile(drawS$lambda[vsampleP])[2])*((quantile(drawS$gamma[vsampleP])[2]^( quantile(drawS$u0[vsampleP])[2] - (vX) ) ) )
    # lines( vX ,vY,xlim=c(0,max_x),ylim=c(0,80),type="l",col="black",lwd=3)
    # lines( vX ,vY_u,xlim=c(0,max_x),ylim=c(0,80),type="l",col="red",lwd=0.5)
    # lines( vX ,vY_l,xlim=c(0,max_x),ylim=c(0,80),type="l",col="red",lwd=0.5)
    # 
    
    
    if (!bPlotIndividualEvents)
      return(NA)
    ##plot individual Curve Fits
    for (pp in vsampleP) ##Go through Each hidx 
    {
      #phi_0[hidx[i]] - lambda[ hidx[i] ] * pow(gamma[hidx[i]],distMax[i] - distP[i] )   
      #      #X11()
      
      pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_",strGroupID[groupID],"_Sigmoid_",vRegIdx[pp],".pdf",sep="")) 
      
      plotEyeGCFit(pp,strGroupID[groupID],dataSubset,drawS)
      
      dev.off()
      
    } ##For Each Sampled Hunt Event 
    
  }##END oF Function 
  
  plotConvergenceDiagnostics <- function(strGroupID,drawS,dataS)
  {
    
    lFitScores <- list()
    vRegIdx <- unique(dataS$RegistrarIdx) ##Get Vector Of RegIdx That Associate with the sample Sequence
    N <- NROW(drawS$tau[,1,1])
    vSqError <- 0
    for (idxH in 1:N)
    {
      
      ## Plot multipage for param convergence ##
      pdf(onefile=TRUE,file= paste(strPlotExportPath,"/stat/diag/stat_SigExpFit_",strGroupID,"_",vRegIdx[idxH],".pdf",sep="")) 
      ## plot the regression lines and the data
      vSqError <- plotEyeGCFit(idxH,strGroupID,dataS,drawS) 
      
      plot(drawS$tau[idxH,,1],type='l',ylim=c(0,4),main=paste("Chains of tau",idxH," (",vRegIdx[idxH],")") )
      lines(drawS$tau[idxH,,2],type='l',col="red")
      lines(drawS$tau[idxH,,3],type='l',col="blue")
      #dev.off()
      
      
      ##gelmal rubin diag ##
      print(paste(idxH,". E[R^2] error:",mean(vSqError) )  )
      chains_tau <- mcmc(drawS$tau[idxH,,],thin=thin)
      lmcmc_tau <- mcmc.list(chains_tau[,1],chains_tau[,2],chains_tau[,3])
      taugdiag_psrf <- gelman.diag(lmcmc_tau,autoburnin=TRUE )$psrf
      #pdf(file= paste(strPlotExportPath,"/stat/diag/stat_gelman_SigExpFit_gamma",idxH,".pdf",sep="")) 
      gelman.plot( lmcmc_tau,autoburnin=TRUE,max.bins=100, ylim=c(0.99,1.5),
                   main=paste("tau psrf:", round(taugdiag_psrf[1]*100)/100 ) )
      
      
      
      #pdf(file= paste(strPlotExportPath,"/stat/diag/stat_SigExpFit_gamma",idxH,".pdf",sep="")) 
      plot(drawS$gamma[idxH,,1],type='l',ylim=c(0,4),main=paste("V rise rate gamma ",idxH," (",vRegIdx[idxH],")") )
      lines(drawS$gamma[idxH,,2],type='l',col="red")
      lines(drawS$gamma[idxH,,3],type='l',col="blue")
      #dev.off()
      
      chains_gamma <- mcmc(drawS$gamma[idxH,,],thin=thin)
      lmcmc_gamma <- mcmc.list(chains_gamma[,1],chains_gamma[,2],chains_gamma[,3])
      gammadiag_psrf <- gelman.diag(lmcmc_gamma,autoburnin=TRUE )$psrf
      #pdf(file= paste(strPlotExportPath,"/stat/diag/stat_gelman_SigExpFit_gamma",idxH,".pdf",sep="")) 
      gelman.plot( lmcmc_gamma,autoburnin=TRUE,max.bins=100, ylim=c(0.99,1.5),
                   main=paste("gamma psrf:", round(gammadiag_psrf[1]*100)/100 ) )
      
      ### LAMBDA ###
      plot(drawS$lambda[idxH,,1],type='l',ylim=c(0,max(drawS$lambda[idxH,,])+1),main=paste("Chains of lambda",idxH," (",vRegIdx[idxH],")") )
      lines(drawS$lambda[idxH,,2],type='l',col="red")
      lines(drawS$lambda[idxH,,3],type='l',col="blue")
      
      
      chains_lambda <- mcmc(drawS$lambda[idxH,,],thin=thin)
      lmcmc_lambda <- mcmc.list(chains_lambda[,1],chains_lambda[,2],chains_lambda[,3])
      lambdadiag_psrf <- gelman.diag(lmcmc_lambda,autoburnin=TRUE )$psrf
      #pdf(file= paste(strPlotExportPath,"/stat/diag/stat_gelman_SigExpFit_gamma",idxH,".pdf",sep="")) 
      gelman.plot( lmcmc_lambda,autoburnin=TRUE,max.bins=100, ylim=c(0.99,1.5),
                   main=paste("lambda psrf:", round(lambdadiag_psrf[1]*100)/100 ) )
      
      
      
      lFitScores[[idxH]] <- list(sampleIdx=idxH,RegistryIdx=vRegIdx[idxH],
                                 gammaPsrf=round(gammadiag_psrf[1]*100)/100,
                                 tauPsrf=round(taugdiag_psrf[1]*100)/100,
                                 lambdaPsrf=round(lambdadiag_psrf[1]*100)/100,
                                 meansqFitError = mean(vSqError),
                                 sdsqFitError = sd(vSqError))
      dev.off()
    }
    
    return(lFitScores)
  }
  
  
  
  #This function converts an mcmc list correctly into an mcarray (It may not be needed).
  `as.mcmc.list.mcarray` <-
    function (x, ...) 
    {
      if (is.null(dim(x)) || is.null(names(dim(x)))) {
        NextMethod()
      }
      xdim <- dim(x)
      ndim <- length(xdim)
      dn <- names(xdim)
      which.iter <- which(dn == "iteration")
      if (length(which.iter) != 1) {
        stop("Bad iteration dimension in mcarray")
      }
      which.chain <- which(dn == "chain")
      if (length(which.chain) > 1) {
        stop("Bad chain dimension in mcarray")
      }
      niter <- xdim[which.iter]
      if (length(which.chain) == 0) {
        perm <- c(which.iter,(1:ndim)[-which.iter])
        x <- matrix(aperm(x, perm), nrow = niter)
        ans <- mcmc.list(mcmc(x))
      }
      else {
        nchain <- xdim[which.chain]
        ans <- vector("list", nchain)
        len <- prod(xdim[-which.chain])
        perm <- c(which.iter, (1:ndim)[-c(which.iter, which.chain)], 
                  which.chain)
        x <- aperm(x, perm)
        for (i in 1:nchain) {
          ans[[i]] <- mcmc(matrix(x[1:len + (i - 1) * len], nrow = niter))
        }
        ans <- mcmc.list(ans)
      }
      return(ans)
    }
  environment(as.mcmc.list.mcarray)<-environment(rjags:::as.mcmc.list.mcarray) #coda 
  
  as.mcarray.mcmc.list <- function(x){
    var.index <- dimnames(x[[1]])[[2]]
    vars <- strsplit(var.index,"\\[|\\]|,")
    var.names <- sapply(vars,function(xx)xx[1])
    var.names <- rle(var.names)
    start <- c(1,cumsum(var.names$lengths)+1)
    start <- start[-length(start)] # remove last one
    stopat <- c(cumsum(var.names$lengths))
    var.dims <- sapply(vars,function(xx)length(xx))-1
    var.dims <- var.dims[start]
    chain.length = dim(x[[1]])[1]
    result <- list()
    
    for (i in 1:length(start)){
      tmp <- switch(var.dims[i]+1,
                    array(c(x[[1]][,start[i]],x[[2]][,start[i]],x[[3]][,start[i]]),
                          dim=c(1,chain.length,3)),
                    aperm(array(c(x[[1]][,start[i]:stopat[i]],x[[2]][,start[i]:stopat[3]],x[[1]][,start[i]:stopat[i]]),
                                dim=c(chain.length,var.names$length[i],3)),c(2,1,3)),
                    aperm(array(c(x[[1]][,start[i]:stopat[i]],x[[2]][,start[i]:stopat[3]],x[[1]][,start[i]:stopat[i]]),
                                dim=c(chain.length,as.numeric(vars[[stopat[i]]][2:3]),3)),c(2,3,1,4)))
      class(tmp)<- "mcarray"
      names(dim(tmp)) <- switch(var.dims[i]+1,
                                c("","iteration","chain"),
                                c("","iteration","chain"),
                                c("","","iteration","chain"))
      result[[var.names$values[i]]] <- tmp
    }
    return(result)
  }
  