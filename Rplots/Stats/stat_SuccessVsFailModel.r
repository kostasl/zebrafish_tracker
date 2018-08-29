## stat Model Success Vs Failure - Plots Distribution For HuntRate and Success Rate
# Makes A model Assuming Poisson Rate lambda and Hunt Success Probability q
# Note: We can exclude the fish that produced no events when excluding eventID == 0- 

source("TrackerDataFilesImport_lib.r")
source("DataLabelling/labelHuntEvents_lib.r")
##wiki suggest gamma prior for lambda
model1="model {
         
         for(i in 1:3) {
             lambda[i] ~ dexp(1) 
             q[i] ~ dbeta(1,1)
         }

         for(i in 1:NTOT){
             Events[i] ~ dpois(lambda[ID[i]])
             Success[i] ~ dbinom(q[ID[i]],Events[i])
         }
}"

#strProcDataFileName <- "setn14-HuntEventsFixExpID-SB-Updated"
#strProcDataFileName <-paste("setn-12-HuntEvents-SB-ALL_19-07-18",sep="") ## Latest Updated HuntEvent Labelled data
strProcDataFileName <- "setn15-HuntEvents-SB-Updated-Merged"
message(paste(" Loading Hunt Event List to Analyse... "))
#load(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
datHuntLabelledEventsSB <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))

##We Can Choose To Exclude The Fish That Produced No Hunting Events
datHuntLabelledEventsSB <- datHuntLabelledEventsSB[datHuntLabelledEventsSB$eventID != 0 &
                                                     datHuntLabelledEventsSB$groupID %in% c("LL","NL","DL") ,]
datFishSuccessRate <- getHuntSuccessPerFish(datHuntLabelledEventsSB)
datFishSuccessRate$groupID <- factor(datFishSuccessRate$groupID)
strGroups <-levels(datFishSuccessRate$groupID)


datatest=list(Success=datFishSuccessRate$Success,
              Events=datFishSuccessRate$Fails+datFishSuccessRate$Success,
              ID=as.numeric(datFishSuccessRate$groupID),
              NTOT=nrow(datFishSuccessRate));

varnames1=c("q","lambda")
burn_in=1000;
steps=10000;
thin=10;

library(rjags)
strModelName = "model1.tmp"
fileConn=file(strModelName)
writeLines(model1,fileConn);
close(fileConn)

m=jags.model(file=strModelName,data=datatest);
update(m,burn_in)
draw=jags.samples(m,steps,thin=thin,variable.names=varnames1)


#######PLOT RESULTS
#X11()
colourH <- c("#0303E663","#03B30363","#E6030363")
colourD <- c("#0303E623","#03B30323","#E6030323")
colourL <- c("#0303E6AF","#03B303AF","#E60303AF")
#for(i in 1:3) hist(draw$q[i,,1],breaks=seq(0,0.5,0.01),col=colourH[i],add=!(i==1))
#legend("topright", legend=strGroups,fill=colourL)

#X11()
#for(i in 1:3) hist(draw$lambda[i,,1],breaks=seq(5,20,0.2),col=colourH[i],add=!(i==1))
#legend("topright", legend=strGroups,fill=colourL)

##Density in 2D of Success Vs Fail
#X11()

strPlotName = paste(strPlotExportPath,"/stat/stat_HuntRateAndEfficiencyEstimation.pdf",sep="")
pdf(strPlotName,width=8,height=8,title="Bayesian Inference on distribution of hunt rate parameter and probability of success, based on labelled data set",onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))

  plot(draw$q[1,,1], draw$lambda[1,,1],col=colourD[1],ylim=c(5,16),xlim=c(0,0.5),pch=19,
       main="Bayesian Estimation for Hunt Rate and Efficiency",
       xlab="Probability of Success q",
       ylab=(expression(paste("Hunt Rate ",lambda ) ) )  ) #paste("Hunt Rate", )
  points(draw$q[2,,1], draw$lambda[2,,1],col=colourD[2],ylim=c(5,16),xlim=c(0,0.5),pch=19)
  points(draw$q[3,,1], draw$lambda[3,,1],col=colourD[3],ylim=c(5,16),xlim=c(0,0.5),pch=19)
  legend("topright", legend=strGroups,fill=colourL)
dev.off()
#hist(draw$q[1,,1],breaks=seq(0,30,length=100),col=colourH[1])
#hist(drawNL$q[1,,1],breaks=seq(0,30,length=100),add=T,col=colourH[2])
#hist(drawDL$q[1,,1],breaks=seq(0,30,length=100),add=T,col=colourH[3])
#legend(2,500,legend = c("LL","NL","DL"),fill=c(colourH[1],colourH[2],colourH[3]))

#hist(drawLL$qq[1,,1],breaks=seq(0,50,length=100))
#hist(drawNL$qq[1,,1],breaks=seq(0,50,length=100),add=T,col=rgb(1,0,0,.4))
#hist(drawDL$qq[1,,1],breaks=seq(0,50,length=100),add=T,col=rgb(0,1,0,.4))

