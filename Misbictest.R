load("F:/Missingdatacodes/Missing_data_italy_weighted/MisGibparameterwghtest.RData")

italy.obs=as.data.frame(italy.obs)



#BIC definition
BICwgh.def=function(z.cop,y,bape,wgh){
  p=NROW(unlist(bape))
  print(p)
  def.val=logwgh.BIC(z.cop,y,bape,wgh)+p*log(NROW(z.cop))#NROW(z.cop)
  return(def.val)
}

#log information of BIC
logwgh.BIC=function(z.cop,y,bape,wgh)
{# the indicator of which sample
  #comp.ind=complete.cases(cbind(z.cop,y))
  y.pro=exp(cbind(rep(1,NROW(z.cop)),z.cop)%*%bape)/(1+exp(cbind(rep(1,NROW(z.cop)),z.cop)%*%bape))
  print(sum(log(dbinom(round(y*wgh),round(wgh),y.pro))))#-1060.964
  -2*sum(log(dbinom(round(y*wgh),round(wgh),y.pro)))
}

#relevel
BIC(glm(contra~.,data=italy.obs[,-16],family=binomial,weights=wgh.obs))#2246.326

result=glm(contra~.,data=italy.obs[,-16],family=binomial,weights=wgh.obs)

BICwgh.def(as.matrix(italy.obs[,-c(16:17)]),italy.obs$contra,result$coefficients,wgh.obs)#2246.326






