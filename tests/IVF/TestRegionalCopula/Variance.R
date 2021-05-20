N_site=20
N_time=50
burnCoef=0.5
upper=0.95
lower=0.05
NF=45
folderForge="C:/Mydoc/Forge/"
folderTest="hydrostat/trunk/Fortran/test/Tests_IVF/TestRegionalCopula/"
a<-c()
b<-c()
ha<-c()
la<-c()
hb<-c()
lb<-c()
#mcmc<- read.table(paste(folder,"hydrostat/trunk/Fortran/test/Tests_IVF/TestRegionalCopula/MCMC/SimulationCop",1,".txt",sep=""),header=T)
#plot(density(mcmc$par2[(burnCoef*Nite):Nite]))
#mcmc<- read.table(paste(folder,"hydrostat/trunk/Fortran/test/Tests_IVF/TestRegionalCopula/MCMC/Simulation",1,".txt",sep=""),header=T)
#lines(density(mcmc$par2[(burnCoef*Nite):Nite]),col=2)
for (k in 1:NF){
    mcmc<- read.table(paste(folderForge,folderTest,"MCMC/SimulationCop",k,".txt",sep=""),header=T)
    Nite=length(mcmc$par2)
    ha[k]<-quantile(mcmc$par2[(burnCoef*Nite):Nite],upper)
    la[k]<-quantile(mcmc$par2[(burnCoef*Nite):Nite],lower)
    a[k]<-var(mcmc$par2[(burnCoef*Nite):Nite])
#    lines(density(mcmc$par2[(burnCoef*Nite):Nite]),col=k)
    mcmc<- read.table(paste(folderForge,folderTest,"MCMC/Simulation",k,".txt",sep=""),header=T)
    Nite=length(mcmc$par2)
    hb[k]<-quantile(mcmc$par2[(burnCoef*Nite):Nite],upper)
    lb[k]<-quantile(mcmc$par2[(burnCoef*Nite):Nite],lower)
    b[k]<-var(mcmc$par2[(burnCoef*Nite):Nite])
}
#for (k in 51:60){
#    mcmc<- read.table(paste(folderForge,folderTest,"MCMC/SimulationCop",k,".txt",sep=""),header=T)
#    Nite=length(mcmc$par2)
#    ha[k-5]<-quantile(mcmc$par2[(burnCoef*Nite):Nite],upper)
#    la[k-5]<-quantile(mcmc$par2[(burnCoef*Nite):Nite],lower)
#    a[k-5]<-var(mcmc$par2[(burnCoef*Nite):Nite])
##    lines(density(mcmc$par2[(burnCoef*Nite):Nite]),col=k)
#    mcmc<- read.table(paste(folderForge,folderTest,"MCMC/Simulation",k,".txt",sep=""),header=T)
#    Nite=length(mcmc$par2)
#    hb[k-5]<-quantile(mcmc$par2[(burnCoef*Nite):Nite],upper)
#    lb[k-5]<-quantile(mcmc$par2[(burnCoef*Nite):Nite],lower)
#    b[k-5]<-var(mcmc$par2[(burnCoef*Nite):Nite])
#}
phi<- read.table(paste(folderForge,folderTest,"phi.txt",sep=""),header=T)
#phi2<- read.table(paste(folderForge,folderTest,"phi2.txt",sep=""),header=T)
#x11()
#plot(rbind(phi,phi2)$x,a,type='l')
#lines(rbind(phi,phi2)$x,b,col=2)

x11()
par(mar=c(5,5,4,2) + 0.1)
plot(phi$x,a,type='l',main="Gaussian",ylab="Variance",xlab=expression(phi),cex.main=2,cex.lab=2,cex.axis=1.5,lwd=2)
lines(phi$x,b,col=2,lwd=2)
legend("topleft",c("Regression with copula","Regression without copula"),lty=c(1,1),pch=c(NA,NA),col=c(1,2),cex=1.5,lwd=2)


x11(20,10)
par(mar=c(5,5,4,2) + 0.1)
plot(phi$x,ha,ylim=c(0.03,0.17),main="Gaussian",ylab="Value",xlab=expression(phi),cex.main=2,cex.lab=2,cex.axis=1.5,lwd=2)
points(phi$x,la)
points(phi$x,hb,col=3)
points(phi$x,lb,col=3)
for (i in 1:NF){
    lines(c(phi$x[i],phi$x[i]),c(ha[i],la[i]),lwd=8,col=5)
    lines(c(phi$x[i],phi$x[i]),c(hb[i],lb[i]),col=6,lwd=5)
}
abline(h=0.1,col=2)
legend("topleft",c("True value","90% Credibility Interval (with copula)","90% Credibility Interval (without copula)"),lty=c(1,1,1),pch=c(NA,NA,NA),col=c(2,5,6),cex=1.5,lwd=c(1,5,5))


phibox1=matrix(0,5000*20,2)
for (k in 1:20){
    mcmc<- read.table(paste(folderForge,folderTest,"MCMC/SimulationCop",k,".txt",sep=""),header=T)
    phibox1[(5000*(k-1)+1):(5000*k),1]=k
    phibox1[(5000*(k-1)+1):(5000*k),2]=mcmc$par4[15001:20000]
}
phibox2=matrix(0,5000*15,2)
for (k in 1:15){
    mcmc<- read.table(paste(folderForge,folderTest,"MCMC/SimulationCop",(k+20),".txt",sep=""),header=T)
    phibox2[(5000*(k-1)+1):(5000*k),1]=k+20
    phibox2[(5000*(k-1)+1):(5000*k),2]=mcmc$par4[15001:20000]
}
phibox3=matrix(0,5000*10,2)
for (k in 1:10){
    mcmc<- read.table(paste(folderForge,folderTest,"MCMC/SimulationCop",(k+35),".txt",sep=""),header=T)
    phibox3[(5000*(k-1)+1):(5000*k),1]=k+35
    phibox3[(5000*(k-1)+1):(5000*k),2]=mcmc$par4[15001:20000]
}
phibox1b=data.frame(phibox1)
phibox2b=data.frame(phibox2)
phibox3b=data.frame(phibox3)
names(phibox1b)=c("phi","value")
names(phibox2b)=c("phi","value")
names(phibox3b)=c("phi","value")
phival=read.table(paste(folderForge,folderTest,"phi.txt",sep=""),header=T)

x11()
boxplot(value ~ phi,phibox1b,outline=F,names=phival[1:20,1],cex.lab=1.75,cex.axis=1.5,ylim=c(0,1)) #,main=Tmain1,xlab=Txlab,ylab=Tylab
x11()
boxplot(value ~ phi,phibox2b,outline=F,names=phival[21:35,1],cex.lab=1.75,cex.axis=1.5,ylim=c(0,1)) #,main=Tmain1,xlab=Txlab,ylab=Tylab
x11()
boxplot(value ~ phi,phibox3b,outline=F,names=phival[36:45,1],cex.lab=1.75,cex.axis=1.5,ylim=c(0,1)) #,main=Tmain1,xlab=Txlab,ylab=Tylab
    