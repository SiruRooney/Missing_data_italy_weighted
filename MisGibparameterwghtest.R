
#############################################################

EM.BIC1=function(predata,bape,wgh){

  Dy=function(predata,bape,wgh){
    predata.inpt=matrix(0.0,NROW(predata),NCOL(predata)+1);
    predata.inpt[,1]=rep(1,NROW(predata));predata.inpt[,-1]=predata
    #observed data
    print(dim(predata.inpt[,-NCOL(predata.inpt)]))
    print(NROW(bape))
    
    p.y=exp(predata.inpt[,-NCOL(predata.inpt)]%*%bape)/(1+exp(predata.inpt[,-NCOL(predata.inpt)]%*%bape))
    p.y=as.vector(p.y)
    
    D1y=colSums(wgh*(predata.inpt[,NROW(bape)+1]-p.y)*predata.inpt[,-NCOL(predata.inpt)])
    D2y.k=t(apply(predata.inpt[,-NCOL(predata.inpt)],1,function(x){x%o%x}))*(-wgh*p.y*(1-p.y))
    D2y.k=colSums(D2y.k)
    D2y=matrix(D2y.k,NROW(bape),NROW(bape))

    list(D1y,D2y)
  }

  i.EM=1
  bape0=bape
  while((i.EM==1)|(max(abs(bape-bape0))>=10e-5)){#
    bape0=bape
    cat("i.EM",i.EM,"\n")
    D12y=Dy(predata,bape,wgh)#;print(D12y[[2]])
    bape=bape-ginv(D12y[[2]])%*%D12y[[1]]
    i.EM=i.EM+1
  }
 bape
}

italy.obs=Italy.dummy[complete.cases(Italy.dummy[,-16])==T,]
inil.par=rep(0,16)
wgh.obs=Weigh[complete.cases(Italy.dummy[,-16])==T]

EM.BIC1(italy.obs[,-16],inil.par,wgh.obs)
#[1,] -3.347360386
#[2,]  0.988720329
#[3,]  0.525729655
#[4,]  0.300936561
#[5,] -0.081971295
#[6,] -1.252592320
#[7,] -1.234562112
#[8,]  0.293040382
#[9,] -0.009363829
#[10,] -1.224077670
#[11,]  2.169380384
#[12,]  1.775264698
#[13,]  1.914474966
#[14,]  2.547714579
#[15,] 17.772616540
#[16,]  4.096547922
italy.obs=as.data.frame(italy.obs)
glm(contra~.,data=italy.obs[,-16],family=binomial,weights=wgh.obs)
#glm(contra/wgh.obs~.,data=italy.obs[,-16],family=binomial,weights=wgh.obs)
logLik(glm(contra~.,data=italy.obs[,-16],family=binomial,weights=wgh.obs))#-1060.964





