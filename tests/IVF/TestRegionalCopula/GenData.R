N_site=20
N_time=50
folder="C:/Mydoc/Forge/hydrostat/trunk/Fortran/test/Tests_IVF/TestRegionalCopula/"
#phi=c(seq(0.99,0.8,-0.01),seq(0.78,0.5,-0.02),seq(0.45,0,-0.05))
phi=c(seq(0.999,0.99,-0.001))

l=length(phi)
library(mvtnorm)
simuData<-matrix(0,N_time, N_site)

for (k in 1:l){
m<-rep(0,N_site)
sig<-matrix(0,N_site,N_site)
for (j in 1:N_site){
    sig[j,j]<-25
}
for (i in 1:N_site){
    for (j in 1:(i-1)){
            sig[i,j]<-phi[k]^abs(i-j)*sqrt(sig[i,i])*sqrt(sig[j,j])
            sig[j,i]<-sig[i,j]
    }
}                           
for (i in 1:N_time){
    for (j in 1:N_site){
        m[j]<-20+0.1*i
    }
    temp<- rmvnorm(n=1, mean=m,sigma=sig)
    simuData[i,]<-temp[1,]
}
write.table(simuData, paste(folder,"Data/SimuData_",(k+50),".txt",sep=""),quote=F,row.names=F)
#write.table(simuData, paste(folder,"SimuData_",k,".txt",sep=""),quote=F,row.names=F)
}

write.table(phi,paste(folder,"phi2.txt",sep=""),quote=F,row.names=F)
#write.table(phi,paste(folder,"phi.txt",sep=""),quote=F,row.names=F)

