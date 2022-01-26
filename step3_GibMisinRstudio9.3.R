#step3.2.1_GibbsModelselectUnix
#when the value of assigned indicator is 1, the BIC of specified model 


#R session 1 for Gibbs
setwd("F://Missingdatacodes//Missing_data_italy_weighted")
rm(list=ls())
load("F://Missingdatacodes//Missing_data_italy_weighted//step1_datageneration.RData")
source("F://Missingdatacodes//Missing_data_italy_weighted//Gibbsampler3.R")

ptm <- proc.time();ptm
rs.indicator=cbind(r.sex,r.fecund,r.contra)
Italy.dummy=as.matrix(Italy.dummy)
indicator.italy=list(rep(1,NCOL(Italy.dummy)),rep(1,NCOL(Italy.dummy)-2),rep(1,NCOL(Italy.dummy)-1),rep(1,NCOL(Italy.dummy)+1),rep(1,NCOL(Italy.dummy)+2),rep(1,NCOL(Italy.dummy)+3))
indicator.italy[[3]][c(13:14)+1]=0
indicator.italy[[4]][c(15,16)+1]=0
indicator.italy[[5]][c(18)+1]=0
indicator.italy[[6]][c(12,14,17)+1]=0
#40
#15
library(MASS)#introduce function ginv()
tab.out=Gibbs.sampler(Italy.dummy,rs.indicator,indicator.italy,90,Weigh)
save.image("F://Missingdatacodes//Missing_data_italy_weighted//step3_GibMisinRstudio9.3.RData")
proc.time()-ptm
