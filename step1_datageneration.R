#data preprocess--transform dummy variables

#install.packages("nnet") #to generate dummy variables
#install.packages("plyr")
setwd("F://Missingdatacodes//Missing_data_italy_weighted")
rm(list=ls())

Italydata=read.csv("F://Missingdatacodes//Italydata//italy.out",head=TRUE)
dim(Italydata)#3374*11


Italydata[(Italydata==-99)]=NA
Italydata=Italydata[(is.na(Italydata$mainact)!=TRUE),]
Weigh=Italydata$WEIGHT
Italydata=Italydata[,-c(1,2)]
dim(Italydata)#3374*9
dim(Italydata)#3371*9
#analyze the Italydata
cat("colnames(Italydata)",colnames(Italydata),"\n")
#contra\agegrp\relat\fecund\parners\children\sex\mainact\morechild
table(Italydata$contra)#0-1428,1-1730
table(Italydata$agegrp)#20-890,25-877,30-840,35-764
table(Italydata$partners)#0-1390,1-1921,2-60
table(Italydata$children)#0-1748,1-710,2-913
table(Italydata$mainact)#1-1593,2-495,3-857,4-412,5-7,6-7
############################################################################
Italydata$fecund[Italydata$fecund==4]=0

agegrp.dummy=apply(Italydata,1,function(x)
{dummy=switch(as.character(x[2]),"20"=c(1,0,0),"25"=c(0,1,0),"30"=c(0,0,1),"35"=c(0,0,0));return(dummy)})
agegrp.dummy=t(agegrp.dummy)
colnames(agegrp.dummy)=c("agegrp1","agegrp2","agegrp3");agegrp.dummy=as.data.frame(agegrp.dummy)
#############################################################################
partners.dummy=apply(Italydata,1,function(x)
{dummy=switch(as.character(x[5]),"0"=c(0,0),"1"=c(1,0),"2"=c(0,1));return(dummy)})
partners.dummy=t(partners.dummy)
colnames(partners.dummy)=c("partner1","partner2");partners.dummy=as.data.frame(partners.dummy)
#############################################################################
children.dummy=apply(Italydata,1,function(x)
{dummy=switch(as.character(x[6]),"0"=c(0,0),"1"=c(1,0),"2"=c(0,1));return(dummy)})
children.dummy=t(children.dummy)
colnames(children.dummy)=c("children1","children2");children.dummy=as.data.frame(children.dummy)
#############################################################################
mainact.dummy=apply(Italydata,1,function(x)
{dummy=switch(as.character(x[8]),"1"=c(1,0,0,0,0),"2"=c(0,1,0,0,0),"3"=c(0,0,1,0,0),"4"=c(0,0,0,1,0),"5"=c(0,0,0,0,1),"6"=c(0,0,0,0,0))})
mainact.dummy=t(unlist(mainact.dummy))
colnames(mainact.dummy)=c("mainact1","mainact2","mainact3","mainact4","mainact5");mainact.dummy=as.data.frame((mainact.dummy))
#############################################################################


Italy.dummy=data.frame(agegrp.dummy,Italydata$relat,partners.dummy,children.dummy,Italydata$morechild,mainact.dummy,Italydata$sex,Italydata$fecund,Italydata$contra)
colnames(Italy.dummy)=c("agegrp1","agegrp2","agegrp3","relat","partner1","partner2","children1","children2",
                        "mchild","mainact1","mainact2","mainact3","mainact4","mainact5","sex","fecund","contra")
##############################################################################
#Missingness indicator 
r.contra=Italy.dummy$contra
r.contra[(is.na(Italy.dummy$contra)==TRUE)]=0;r.contra[(is.na(Italy.dummy$contra)==FALSE)]=1
r.fecund=Italy.dummy$fecund
r.fecund[(is.na(Italy.dummy$fecund)==TRUE)]=0;r.fecund[(is.na(Italy.dummy$fecund)==FALSE)]=1
r.sex=Italy.dummy$sex
r.sex[(is.na(Italy.dummy$sex)==TRUE)]=0;r.sex[(is.na(Italy.dummy$sex)==FALSE)]=1
#r=cbind(r.sex,r.fecund)
Italy.dummy=as.matrix(Italy.dummy)






save.image("F://Missingdatacodes//Missing_data_italy_weighted//step1_datageneration.RData")
