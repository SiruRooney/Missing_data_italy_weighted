
#Analyze the initial full model with weighted
setwd("F://Missingdatacodes//Missing_data_italy_weighted//step3_GibbsforModelselect")
load("F://Missingdatacodes//Missing_data_italy_weighted//step1_datageneration.RData")
source("F://Missingdatacodes//Missing_data_italy_weighted//MisGibparameterwgh.R")

source("F://Missingdatacodes//Missing_data_italy_weighted//fisherinformation.R")
ptm <- proc.time();ptm
rs.indicator=cbind(r.sex,r.fecund,r.contra)
Italy.dummy=as.matrix(Italy.dummy)
indicator.italy=list(rep(1,NCOL(Italy.dummy)),rep(1,NCOL(Italy.dummy)-2),rep(1,NCOL(Italy.dummy)-1),rep(1,NCOL(Italy.dummy)+1),rep(1,NCOL(Italy.dummy)+2),rep(1,NCOL(Italy.dummy)+3))


library(MASS)#introduce function ginv()
initial=EM.BIC(Italy.dummy,BAPE(indicator.italy),ind.BAPE(indicator.italy),rs.indicator,Weigh)
initial.fisher=fisherinfo(Italy.dummy,initial$par.coef,ind.BAPE(indicator.italy),rs.indicator,Weigh)



indicator.italy[[3]][c(13:14)+1]=0
indicator.italy[[4]][c(15,16)+1]=0
indicator.italy[[5]][c(18)+1]=0
indicator.italy[[6]][c(12,14,17)+1]=0
initial1=EM.BIC(Italy.dummy,BAPE(indicator.italy),ind.BAPE(indicator.italy),rs.indicator,Weigh)
################################################################
initial1.fisher=fisherinfo(Italy.dummy,initial1$par.coef,ind.BAPE(indicator.italy),rs.indicator,Weigh)





#test the performance of the initial full model without weights
initial.unwgh=EM.BIC(Italy.dummy,BAPE(indicator.italy),ind.BAPE(indicator.italy),rs.indicator,rep(1,3371))
indicator.italy=list(rep(1,NCOL(Italy.dummy)),rep(1,NCOL(Italy.dummy)-2),rep(1,NCOL(Italy.dummy)-1),rep(1,NCOL(Italy.dummy)+1),rep(1,NCOL(Italy.dummy)+2),rep(1,NCOL(Italy.dummy)+3))
initial.unwgh1=EM.BIC(Italy.dummy,BAPE(indicator.italy),ind.BAPE(indicator.italy),rs.indicator,rep(1,3371))
library(nnet)

Italydata$agegrp=relevel(as.factor(Italydata$agegrp),ref="35")
Italydata$partners=relevel(as.factor(Italydata$partners),ref="0")
Italydata$child=relevel(as.factor(Italydata$child),ref="0")
Italydata$mainact=relevel(as.factor(Italydata$mainact),ref="6")



initial.glmw=glm(contra~agegrp+relat+partners+child+morechild+mainact+sex+fecund,data=Italydata,family=binomial,weights=Weigh)
initial.glm=glm(contra~agegrp+relat+partners+child+morechild+mainact+sex+fecund,data=Italydata,family=binomial)
summary(initial.glmw)

