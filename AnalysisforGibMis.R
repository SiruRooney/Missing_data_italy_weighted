setwd("C:/Users/darun/OneDrive - The University of Melbourne/Codes/Missing_data_italy_weighted")
load("step3_GibMisinRstudio9.3.RData")
bicobs4550=tab.out[[1]]
bicobs4550=lapply(bicobs4550,function(x){x[1:8190,]})
load("step3_GibMisinRstudio9.2.RData")
bicQ4550=tab.out[[1]]
bicQ4550=lapply(bicQ4550,function(x){x[1:8190,]})


load("C:/Users/darun/OneDrive - The University of Melbourne/Codes/Missing_data_italy_weighted/step3_GibMisinRstudio9.21.RData")
bicQ11830=tab.out2[[1]]
bicQ11830=lapply(bicQ11830,function(x){x[1:3640,]})
bicQ11830[[1]]=rbind(bicQ4550[[1]],bicQ11830[[1]])
bicQ11830[[2]]=rbind(bicQ4550[[2]],bicQ11830[[2]])
bicQ11830[[3]]=rbind(bicQ4550[[3]],bicQ11830[[3]])
bicQ11830[[4]]=rbind(bicQ4550[[4]],bicQ11830[[4]])
bicQ11830[[5]]=rbind(bicQ4550[[5]],bicQ11830[[5]])
bicQ11830[[6]]=rbind(bicQ4550[[6]],bicQ11830[[6]])

plot(bicQ11830[[1]][,1],type="l",col="blue",ylim=c(11200,11950),ylab="BIC_{Q}",xlab="iteration")
abline(h=min(bicQ11830[[1]][,1]),col="red")
sort(table(bicQ11830[[1]][,1]))#415
NROW(table(bicQ11830[[1]][,1]))#342


###################################################################################################################################
plot(bicQ4550[[1]][,1],type="l",col="blue",ylim=c(11200,11950),ylab="BIC_{Q}",xlab="iteration")
abline(h=min(bicQ4550[[1]][,1]),col="red")

axis(2,seq(11200,11950,300),seq(11200,11950,300))
title(main="I-chart for the generated BIC sequences in Italydata")
par(new=T)
plot(bicobs4550[[1]][,1],type="l",col="darkgreen",ylim=c(6000,6650),ylab="",xlab="",xaxt="n",yaxt="n")
abline(h=min(bicobs4550[[1]][,1]),col="orange")
axis(4,seq(6000,6650,150),seq(6000,6650,150))
mtext("BIC_{obs}",4)
legend("topright",lty=c(1,1),col=c("blue","darkgreen"),legend=c("BIC_{Q}","BIC_{obs}"))

sort(table(bicQ4550[[1]][,1]))#3
sort(table(bicobs4550[[1]][,1]))#503
NROW(table(bicQ4550[[1]][,1]))#238
NROW(table(bicobs4550[[1]][,1]))#259
##############################################################################

source("MisGibparameterwgh.R")
#################################################################
opmdl.Q2=lapply(bicQ11830,function(x){y=x[x[,1]==sort(bicQ11830[[1]][,1])[1],];y})
opmdl.Q2=lapply(opmdl.Q2,function(x){colSums(x[,-1])})

opmdlQ2.list=lapply(opmdl.Q2,function(x){x[-1][(x[-1]/415)<0.9]=rep(0,NROW(x[-1][(x[-1]/415)<0.9]));x[-1][(x[-1]/415)>=0.9]=rep(1,NROW(x[-1][(x[-1]/415)>=0.9]));x[1]=1;x})
library(MASS)
opmdlQ2.esti=EM.BIC(Italy.dummy,BAPE(opmdlQ2.list),ind.BAPE(opmdlQ2.list),rs.indicator,Weigh)

##################################################################
opmdl.Q=lapply(bicQ4550,function(x){y=x[x[,1]==min(bicQ4550[[1]][,1]),];y})
opmdl.Q=lapply(opmdl.Q,function(x){colSums(x[,-1])})
opmdlQ.list=lapply(opmdl.Q,function(x){x[-1][(x[-1]/3)<0.9]=rep(0,NROW(x[-1][(x[-1]/3)<0.9]));x[-1][(x[-1]/3)>=0.9]=rep(1,NROW(x[-1][(x[-1]/3)>=0.9]));x[1]=1;x})
opmdlQ.esti=EM.BIC(Italy.dummy,BAPE(opmdlQ.list),ind.BAPE(opmdlQ.list),rs.indicator,Weigh)

opmdl.obs=lapply(bicobs4550,function(x){y=x[x[,1]==min(bicobs4550[[1]][,1]),];y})
opmdl.obs=lapply(opmdl.obs,function(x){colSums(x[,-1])})
opmdlobs.list=lapply(opmdl.obs,function(x){x[-1][(x[-1]/503)<0.9]=rep(0,NROW(x[-1][(x[-1]/503)<0.9]));x[-1][(x[-1]/503)>=0.9]=rep(1,NROW(x[-1][(x[-1]/503)>=0.9]));x[1]=1;x})
opmdlobs.esti=EM.BIC(Italy.dummy,BAPE(opmdlobs.list),ind.BAPE(opmdlobs.list),rs.indicator,Weigh)

source("fisherinformation.R")

opmdlQ2.fish=fisherinfo(Italy.dummy,opmdlQ2.esti$par.coef,ind.BAPE(opmdlQ2.list),rs.indicator,Weigh)
sqrt(diag(solve(opmdlQ2.fish)))
2-2*pnorm(abs(unlist(opmdlQ2.esti$par.coef))/sqrt(diag(solve(opmdlQ2.fish))))


opmdlQ.fish=fisherinfo(Italy.dummy,opmdlQ.esti$par.coef,ind.BAPE(opmdlQ.list),rs.indicator,Weigh)
sqrt(diag(solve(opmdlQ.fish)))
2-2*pnorm(abs(unlist(opmdlQ.esti$par.coef))/sqrt(diag(solve(opmdlQ.fish))))
opmdlobs.fish=fisherinfo(Italy.dummy,opmdlobs.esti$par.coef,ind.BAPE(opmdlobs.list),rs.indicator,Weigh)
sqrt(diag(solve(opmdlobs.fish)))
2-2*pnorm(abs(unlist(opmdlobs.esti$par.coef))/sqrt(diag(solve(opmdlobs.fish))))
save.image("AnalysisforGibMis.RData")
