modelI="model { 
                qq ~ dnorm(10,0.001)T(0,400)
		tt ~ dnorm(0,0.001)T(0,100)

         	for(j in 1:NTOT){
		  q[j] ~ dnorm(qq,tt) 
		  n[j] ~ dpois(q[j])
		 }
}"

modelG="model { 
                q ~ dnorm(10,0.001)T(0,400)
                
                for(j in 1:NTOT){
                  n[j] ~ dpois(q)
                 }
                 

}"


nLL=read.table("testLL.dat")$Freq
nNL=read.table("testNL.dat")$Freq
nDL=read.table("testDL.dat")$Freq
ntest=c(nLL,nNL);

dataLL=list(n=nLL,NTOT=length(nLL));
dataNL=list(n=nNL,NTOT=length(nNL));
dataDL=list(n=nDL,NTOT=length(nDL));
datatest=list(n=ntest,NTOT=length(ntest));
varnames1=c("q")
burn_in=1000;
steps=100000;
thin=10;

library(rjags)
strModelName = "modelG.tmp"
fileConn=file(strModelName)
writeLines(modelG,fileConn);
close(fileConn)

mLL=jags.model(file=strModelName,data=dataLL);
mNL=jags.model(file=strModelName,data=dataNL);
mDL=jags.model(file=strModelName,data=dataDL);
mtest=jags.model(file=strModelName,data=datatest);
update(mLL,burn_in)
update(mNL,burn_in)
update(mDL,burn_in)
update(mtest,burn_in)
drawLL=jags.samples(mLL,steps,thin=thin,variable.names=varnames1)
drawNL=jags.samples(mNL,steps,thin=thin,variable.names=varnames1)
drawDL=jags.samples(mDL,steps,thin=thin,variable.names=varnames1)
drawtest=jags.samples(mtest,steps,thin=thin,variable.names=varnames1)
##q[idx,sampleID,chainID]

X11()
colourH <- c("#03B30333","#E6030333","#0303E633")
hist(drawLL$q[1,,1],breaks=seq(0,30,length=100),col=colourH[1])
hist(drawNL$q[1,,1],breaks=seq(0,30,length=100),add=T,col=colourH[2])
hist(drawDL$q[1,,1],breaks=seq(0,30,length=100),add=T,col=colourH[3])
legend(2,500,legend = c("LL","NL","DL"),fill=c(colourH[1],colourH[2],colourH[3]))

#hist(drawLL$qq[1,,1],breaks=seq(0,50,length=100))
#hist(drawNL$qq[1,,1],breaks=seq(0,50,length=100),add=T,col=rgb(1,0,0,.4))
#hist(drawDL$qq[1,,1],breaks=seq(0,50,length=100),add=T,col=rgb(0,1,0,.4))

