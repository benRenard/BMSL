# View help
system('./distribution -h')

# Options usage
system('./distribution -name GEV -par 100,50,-0.2 -act d -x 20,300,100 -rf myResultFile.txt')
res=read.table('myResultFile.txt')
plot(res[,1],res[,2])

# Check against R computation
library(HydroPortailStats)
checkAgainstR <- function(dist,pars,distR=dist){
  actions=c('d','p','q','r')
  actionsR=c('GetPdf','GetCdf','GetQuantile','Generate')
  cmd0=paste0('./distribution -name ',dist,' -par ',paste0(pars,collapse=','),' -act ')
  par(mfrow=c(1,4))
  for(i in 1:4){
    action=actions[i]
    system(paste0(cmd0,action))
    res=read.table('distribution_result.txt')
    if(action=='r'){
      y=res[,1]
      yR=sapply(distR,actionsR[i],pars,100)
      plot(y);lines(yR,col='red')
    } else {
      x=res[,1];y=res[,2]
      yR=sapply(x,actionsR[i],distR,pars)
      plot(x,y);lines(x,yR,col='red')
    }
  }
}

checkAgainstR('Gaussian',c(100,25),distR='Normal')
checkAgainstR('GEV',c(100,50,-0.2))
checkAgainstR('Uniform',c(5,10))
checkAgainstR('LogNormal',c(0,0.5))
checkAgainstR('GPD',c(100,50,-0.2),distR='GPD3')
