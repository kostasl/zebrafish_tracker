## stat Model Success Vs Failure

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

strProcDataFileName <- "setn14-HuntEventsFixExpID-SB-Updated"
#strProcDataFileName <-paste("setn-12-HuntEvents-SB-ALL_19-07-18",sep="") ## Latest Updated HuntEvent Labelled data
message(paste(" Loading Hunt Event List to Analyse... "))
#load(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
datHuntLabelledEventsSB <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))

datFishSuccessRate <- getHuntSuccessPerFish(getHuntSuccessPerFish)
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

X11()
colourH <- c("#0303E633","#03B30333","#E6030333")
for(i in 1:3) hist(draw$q[i,,1],breaks=seq(0,0.5,0.01),col=colourH[i],add=!(i==1))
legend("topright", legend=strGroups,fill=colourH)

X11()
for(i in 1:3) hist(draw$lambda[i,,1],breaks=seq(5,20,0.2),col=colourH[i],add=!(i==1))
legend("topright", legend=strGroups,fill=colourH)

##Density in 2D of Success Vs Fail
X11()
plot(draw$q[1,,1], draw$lambda[1,,1],col=colourH[1])
plot(draw$q[1,,1], draw$lambda[1,,1],col=colourH[1],ylim=c(5,20),xlim=c(0,0.5),pch=19)
points(draw$q[2,,1], draw$lambda[2,,1],col=colourH[2],ylim=c(5,20),xlim=c(0,0.5),pch=19)
points(draw$q[3,,1], draw$lambda[3,,1],col=colourH[3],ylim=c(5,20),xlim=c(0,0.5),pch=19)
legend("topright", legend=strGroups,fill=colourH)
#hist(draw$q[1,,1],breaks=seq(0,30,length=100),col=colourH[1])
#hist(drawNL$q[1,,1],breaks=seq(0,30,length=100),add=T,col=colourH[2])
#hist(drawDL$q[1,,1],breaks=seq(0,30,length=100),add=T,col=colourH[3])
#legend(2,500,legend = c("LL","NL","DL"),fill=c(colourH[1],colourH[2],colourH[3]))

#hist(drawLL$qq[1,,1],breaks=seq(0,50,length=100))
#hist(drawNL$qq[1,,1],breaks=seq(0,50,length=100),add=T,col=rgb(1,0,0,.4))
#hist(drawDL$qq[1,,1],breaks=seq(0,50,length=100),add=T,col=rgb(0,1,0,.4))

