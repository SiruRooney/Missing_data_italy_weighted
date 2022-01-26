setwd("C:/Users/darun/OneDrive - The University of Melbourne/Codes/Missing_data_italy_weighted")
rm(list=ls())
load("C:/Users/darun/OneDrive - The University of Melbourne/Codes/Missing_data_italy_weighted/step3(1)_GibMisinRstudio9.2.RData")
bicQ8820=result1[[1]]
bicQ8820=lapply(bicQ8820,function(x){x[1:8820,]})
load("C:/Users/darun/OneDrive - The University of Melbourne/Codes/Missing_data_italy_weighted/step3(1)_GibMisinRstudio9.3.RData")
bicobs8820=result2[[1]]
bicobs8820=lapply(bicobs8820,function(x){x[1:8820,]})

###################################################################################################################################
plot(bicQ8820[[1]][,1],type="l",col="blue",ylim=c(9300,10220),ylab="BIC_{Q}",xlab="iteration")
abline(h=min(bicQ8820[[1]][,1]),col="red")

axis(2,seq(9300,10220,300),seq(9300,10220,300))
title(main="I-chart for the generated BIC sequences in Italydata")
par(new=T)
plot(bicobs8820[[1]][,1],type="l",col="darkgreen",ylim=c(5000,5800),ylab="",xlab="",xaxt="n",yaxt="n")#

abline(h=min(bicobs8820[[1]][,1]),col="orange")
axis(4,seq(5000,5800,150),seq(5000,5800,150))
mtext("BIC_{obs}",4)
legend("topright",lty=c(1,1),col=c("blue","darkgreen"),legend=c("BIC_{Q}","BIC_{obs}"))

sort(table(bicQ8820[[1]][,1]))#377
sort(table(bicobs8820[[1]][,1]))#2533
NROW(table(bicQ8820[[1]][,1]))#322
NROW(table(bicobs8820[[1]][,1]))#224
##############################################################################
library(MASS)
source("MisGibparameterwgh.R")
opmdl.Q=lapply(bicQ8820,function(x){y=x[x[,1]==sort(bicQ8820[[1]][,1])[1],];y})
opmdl.Q=lapply(opmdl.Q,function(x){colSums(x[,-1])})
opmdlQ.list=lapply(opmdl.Q,function(x){x[-1][(x[-1]/263)<0.9]=rep(0,NROW(x[-1][(x[-1]/263)<0.9]));x[-1][(x[-1]/263)>=0.9]=rep(1,NROW(x[-1][(x[-1]/263)>=0.9]));x[1]=1;x})
opmdlQ.esti=EM.BIC(Italy.dummy,BAPE(opmdlQ.list),ind.BAPE(opmdlQ.list),rs.indicator,Weigh)

opmdl.obs=lapply(bicobs8820,function(x){y=x[x[,1]==sort(bicobs8820[[1]][,1])[1],];y})
opmdl.obs=lapply(opmdl.obs,function(x){colSums(x[,-1])})
opmdlobs.list=lapply(opmdl.obs,function(x){x[-1][(x[-1]/2533)<0.9]=rep(0,NROW(x[-1][(x[-1]/2533)<0.9]));x[-1][(x[-1]/2533)>=0.9]=rep(1,NROW(x[-1][(x[-1]/2533)>=0.9]));x[1]=1;x})
opmdlobs.esti=EM.BIC(Italy.dummy,BAPE(opmdlobs.list),ind.BAPE(opmdlobs.list),rs.indicator,Weigh)

source("fisherinformation.R")
opmdlQ.fish=fisherinfo(Italy.dummy,opmdlQ.esti$par.coef,ind.BAPE(opmdlQ.list),rs.indicator,Weigh)
sqrt(diag(solve(opmdlQ.fish)))
2-2*pnorm(abs(unlist(opmdlQ.esti$par.coef))/sqrt(diag(solve(opmdlQ.fish))))
opmdlobs.fish=fisherinfo(Italy.dummy,opmdlobs.esti$par.coef,ind.BAPE(opmdlobs.list),rs.indicator,Weigh)
sqrt(diag(solve(opmdlobs.fish)))
2-2*pnorm(abs(unlist(opmdlobs.esti$par.coef))/sqrt(diag(solve(opmdlobs.fish))))
save.image("AnalysisforGibMis(1).RData")
