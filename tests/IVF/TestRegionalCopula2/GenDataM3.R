N_site=20
N_time=50
folder=""
phi=c(seq(0.99,0.8,-0.01),seq(0.78,0.5,-0.02),seq(0.45,0,-0.05))
l=length(phi)

library(mvtnorm)
NsimuData<-matrix(0,N_time, N_site)
simuData<-matrix(0,N_time, N_site)

for (k in 1:l){
m<-matrix(0,N_time,N_site)
sig<-matrix(0,N_site,N_site)
for (j in 1:N_site){
    sig[j,j]<-1
}
for (i in 1:N_site){
    for (j in 1:(i-1)){
        sig[i,j]<-phi[k]^abs(i-j)
        sig[j,i]<-sig[i,j]
    }
}    

NsimuData<-rmvnorm(n=N_time, mean=rep(0,N_site),sigma=sig) 
cdfN<-pnorm(NsimuData)

library(evd)                   
for (i in 1:N_time){
    for (j in 1:N_site){
        m[i,j]<-25+0.2*i
        #simuData[i,j]<-rgev(1,m[i,j],10,0.1)
        simuData[i,j]<-qgev(p=cdfN[i,j], loc=m[i,j], scale=10, shape=0.1, lower.tail = TRUE)
    }
}                                       
write.table(simuData, paste(folder,"SimuData_",k,".txt",sep=""),quote=F,row.names=F)
}

write.table(phi,paste(folder,"phi.txt",sep=""),quote=F,row.names=F)     
