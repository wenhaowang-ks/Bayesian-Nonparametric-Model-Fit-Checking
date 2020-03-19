#package--------------------
library(dplyr)
library(reshape2)
library(mirt)
library(nimble)
library(coda)

#word directory----------------
path= "C:\\Users\\wenhaowang\\Documents\\Publication\\BaysianModelFit\\Revision2019.12\\Simulation"
setwd(path)

#---Conditon-----------------------
N=5000
n=60
Replication=100
MisfitConditon="Null"

#------PPMC value---------------------
#---------theta node---------------
tn<-as.data.frame(seq(-2.5,2.5,0.1))
colnames(tn)<-"thetanode"
h=1.1*N^(-0.2)

#-----NIMBLE related code-----------------
nimblerank<-nimbleFunction(
  run=function(x=double(1),N=double(0)){
    returnType(double(1))
    for (k in 2:N){
      for (kk in 1:N){
        if (x[kk]>x[k]){
          t=x[kk]
          x[kk]=x[k]
          x[k]=t
        }
      }
    }
    return(x)
  }
)

#b<-runif(4)
#nimblerank(b,4)

irtCode<-nimbleCode({
  for (i in 1:N){
    for (j in 1:ndct){
      y[i,j]~dbern(prob[i,j])
    }
    for (j in (1+ndct):n){
      y[i,j]~dcat(probk[i,(j-ndct),1:K])
    }
    
    theta[i]~dnorm(0,1)
    
    for (j in 1:ndct){
      logit(prob[i,j])<-a[j]*(theta[i]-b[j])
      P[i,j,1]<-1
    }
    for (j in (1+ndct):n){
      for (k in 1:(K-1)){
        logit(P[i,(j-ndct),(k+1)])<-a[j]*(theta[i]-bk[(j-ndct),k])
      }
      for (k in 1:K-1){
        probk[i,(j-ndct),k]<-P[i,(j-ndct),k]-P[i,(j-ndct),k+1]
      }
      probk[i,(j-ndct),K]<-P[i,(j-ndct),K]
    }
    
  }
  
  for (j in 1:n){
    a[j]~dlnorm(0,1)
  }
  for (j in 1:ndct){
    b[j]~dnorm(0,1)
  }
  for (j in 1:(n-ndct)){
    for (k in 1:(K-1)){
      bk.star[j,k]~dnorm(0,1)
    }
    bk[j,1:(K-1)]<-nimblerank(bk.star[j,1:(K-1)],(K-1))
  }
  
})



#--------read in data------------------
truetheta<-read.csv(paste0('truetheta_',N,'.csv'),header=T)
#-----read in item paramter----------------
pool<-read.csv(paste0('pool_',n,'.csv'),header=T)

#------random seed-----------------------------
set.seed(123)

#------fit value summary-------------
irtpx2<-NULL
irtpsx2<-NULL
bsp<-NULL
pppv<-NULL



for (rpl in 1:Replication){
  #print (rpl)
  
  #----step 1: Generate Response-----------
  if (MisfitConditon=='Null'){
    p<-merge(truetheta,pool)%>%
      mutate(p0=1,
             p1=1 / (1 + exp(-a*truetheta + a*b1)),
             p2=ifelse(is.na(b2),0,1 / (1 + exp(-a*truetheta + a*b2))),
             p3=ifelse(is.na(b3),0,1 / (1 + exp(-a*truetheta + a*b3))),
             p4=ifelse(is.na(b4),0,1 / (1 + exp(-a*truetheta + a*b4))),
             p5=0)
  }else if (MisfitConditon=='Cubic'){
    p<-merge(truetheta,pool)%>%
      mutate(p0=1,
             p1=ifelse(itemid%in%c(1:(n*0.15),(n*0.8+1):(n*0.8+n*0.05)),1 / (1 + exp(-a*(truetheta**3) + a*b1)),1 / (1 + exp(-a*truetheta + a*b1))),
             p2=ifelse(is.na(b2),0,ifelse(itemid%in%c(1:(n*0.15),(n*0.8+1):(n*0.8+n*0.05)),1 / (1 + exp(-a*(truetheta**3) + a*b2)),1 / (1 + exp(-a*truetheta + a*b2)))),
             p3=ifelse(is.na(b3),0,ifelse(itemid%in%c(1:(n*0.15),(n*0.8+1):(n*0.8+n*0.05)),1 / (1 + exp(-a*(truetheta**3) + a*b3)),1 / (1 + exp(-a*truetheta + a*b3)))),
             p4=ifelse(is.na(b4),0,ifelse(itemid%in%c(1:(n*0.15),(n*0.8+1):(n*0.8+n*0.05)),1 / (1 + exp(-a*(truetheta**3) + a*b4)),1 / (1 + exp(-a*truetheta + a*b4)))),
             p5=0)

  }else if (MisfitConditon=='NonMonotonic'){
    p<-merge(truetheta,pool)%>%
      mutate(p0=1,
             p1=ifelse(itemid%in%c(1:(n*0.15),(n*0.8+1):(n*0.8+n*0.05)),
                       1 / (1 + exp(-0.5*a*((truetheta+1.5)**3-(truetheta+2)- b1))),
                       1 / (1 + exp(-a*truetheta + a*b1))),
             p2=ifelse(is.na(b2),0,ifelse(itemid%in%c(1:(n*0.15),(n*0.8+1):(n*0.8+n*0.05)),
                                          1 / (1 + exp(-0.5*a*((truetheta+1.5)**3-(truetheta+2)- b2))),
                                          1 / (1 + exp(-a*truetheta + a*b2)))),
             p3=ifelse(is.na(b3),0,ifelse(itemid%in%c(1:(n*0.15),(n*0.8+1):(n*0.8+n*0.05)),
                                          1 / (1 + exp(-0.5*a*((truetheta+1.5)**3-(truetheta+2)- b3))),
                                          1 / (1 + exp(-a*truetheta + a*b3)))),
             p4=ifelse(is.na(b4),0,ifelse(itemid%in%c(1:(n*0.15),(n*0.8+1):(n*0.8+n*0.05)),
                                          1 / (1 + exp(-0.5*a*((truetheta+1.5)**3-(truetheta+2)- b4))),
                                          1 / (1 + exp(-a*truetheta + a*b4)))),
             p5=0)
  }
  p$mcp<-runif(dim(p)[1])
  p<-p%>%mutate(score=ifelse(mcp>p1,0,ifelse(mcp>p2,1,ifelse(mcp>p3,2,ifelse(mcp>p4,3,4)))))
  #sc<-p%>%select(itemid,score)%>%distinct()
  #print(dim(sc)[1])
  data<-dcast(p,sid~itemid,value.var="score")%>%arrange(sid)
  #write.table(data,paste0(MisfitConditon,"_T",N,"_P",n,"_rpl",rpl,".txt"),sep="\t",row.names=F,col.names=F)
  
  
  #---step 2 : IRT Estimation Using MIRT-------------------------------------
  #----this steps 1. get item parameters for bootstrapping-------------
  #---------------2. get traditional IRT fit statistics-------------------------
  ncat<- apply(data[,-1], 2, function(x) length(unique(x[!is.na(x)])))
  item.list <- unlist(lapply(ncat, function(x){if(x==2) item.list="2PL" else item.list="graded"}))
  mod <- mirt(data=data[,-1],model=1, method = "EM", SE=TRUE, quadpts = 60, technical = list(NCYCLES=c(1000,200)), itemtype = item.list)
  estparm<-data.frame(cbind(MDISC(mod), MDIFF(mod)))
  colnames(estparm) <- c("a",paste0("b", 1:(ncol(estparm)-1)))
  irtpx2<-cbind(irtpx2,itemfit(mod, fit_stats = "X2")$p.X2)
  irtpsx2<-cbind(irtpsx2,itemfit(mod, fit_stats = "S_X2")$p.S_X2)
  estthe<-fscores(mod, method='MAP')
  

  #---------original data nonparametric estimation---------------
  #--------this will be used in step 3 and step 4-------------
  data2<-melt(data,id='sid')%>%
    mutate(itemid=as.numeric(as.character(variable)))%>%
    select(sid,itemid,score=value)%>%group_by(sid)%>%mutate(total=sum(score))%>%
    ungroup()%>%mutate(totalni=total-score)
  etheta<-data2%>%
    group_by(itemid,totalni)%>%summarize(counts=n())%>%ungroup()%>%
    arrange(itemid,totalni)%>%group_by(itemid)%>%mutate(p=cumsum(counts)/N)%>%
    mutate(etheta=ifelse(qnorm(p)==Inf,4,qnorm(p)))%>%select(itemid,totalni,etheta)%>%
    ungroup()
  data3<-left_join(data2,etheta,by=c('itemid','totalni'))
  data4<-do.call("rbind", replicate(dim(tn)[1], data3, simplify = FALSE))
  tnd4<-do.call("rbind", replicate(dim(data3)[1], tn, simplify = FALSE))%>%
    arrange(thetanode)
  data3<-cbind(data4,tnd4)%>%mutate(w=dnorm((etheta-thetanode)/h))%>%
    mutate(wr1=ifelse(score>=1,w,0),
           wr2=ifelse(score>=2,w,0),
           wr3=ifelse(score>=3,w,0),
           wr4=ifelse(score>=4,w,0))%>%
    group_by(itemid,thetanode)%>%
    summarize(ks1=sum(wr1)/sum(w),
              ks2=sum(wr2)/sum(w),
              ks3=sum(wr3)/sum(w),
              ks4=sum(wr4)/sum(w))%>%ungroup()
  rm(tnd4,data4,data2)
  
  
  #----step 3 : Bootstrapping to get bootstrapping p value------------------------------
  estthe<-as.data.frame(estthe)%>%mutate(sid=row_number())
  names(estthe)<-c('estthe','sid')
  estparm<-estparm%>%mutate(itemid=row_number())
  #------------------Bootstrapping realized RISE-------------------
  estparm2<-do.call("rbind", replicate(dim(tn)[1], estparm, simplify = FALSE))
  tnep2<-do.call("rbind", replicate(dim(estparm)[1], tn, simplify = FALSE))%>%
    arrange(thetanode)
  estmodelp<-cbind(estparm2,tnep2)%>%mutate(itemid=as.numeric(itemid))%>%
    mutate(p1=1 / (1 + exp(-a*thetanode + a*b1)),
           p2=ifelse(is.na(b2),0,1 / (1 + exp(-a*thetanode + a*b2))),
           p3=ifelse(is.na(b3),0,1 / (1 + exp(-a*thetanode + a*b3))),
           p4=ifelse(is.na(b4),0,1 / (1 + exp(-a*thetanode + a*b4))))%>%
    select(itemid,thetanode,p1,p2,p3,p4)
  rm(tnep2,estparm2)
  bsrRISE<-left_join(estmodelp,data3,by=c('itemid','thetanode'))%>%
    mutate(RISEp1=(p1-ks1)^2*dnorm(thetanode),
           RISEp2=(p2-ks2)^2*dnorm(thetanode),
           RISEp3=(p3-ks3)^2*dnorm(thetanode),
           RISEp4=(p4-ks4)^2*dnorm(thetanode))%>%
    group_by(itemid)%>%
    summarize(rRISE=sqrt(sum(RISEp1))+sqrt(sum(RISEp2))+sqrt(sum(RISEp3))+sqrt(sum(RISEp4)))%>%ungroup()
  #-------------------Bootstrapping predictive RISE and calcualte p value-------------------
  bsdatap<-merge(estthe,estparm)%>%
    mutate(p0=1,
           p1=1 / (1 + exp(-a*estthe + a*b1)),
           p2=ifelse(is.na(b2),0,1 / (1 + exp(-a*estthe + a*b2))),
           p3=ifelse(is.na(b3),0,1 / (1 + exp(-a*estthe + a*b3))),
           p4=ifelse(is.na(b4),0,1 / (1 + exp(-a*estthe + a*b4))),
           p5=0)
  bsallp<-NULL
  #starttime=Sys.time()
  for (i in 1:100){
    bsdata<-bsdatap
    bsdata$mcp<-runif(dim(bsdata)[1])
    bsdata<-bsdata%>%
      mutate(score=ifelse(mcp>p1,0,ifelse(mcp>p2,1,ifelse(mcp>p3,2,ifelse(mcp>p4,3,4)))))%>%
      select(sid,itemid,score)%>%
      group_by(sid)%>%mutate(total=sum(score))%>%
      ungroup()%>%mutate(totalni=total-score)
    bsetheta<-bsdata%>%
      group_by(itemid,totalni)%>%summarize(counts=n())%>%ungroup()%>%
      arrange(itemid,totalni)%>%group_by(itemid)%>%mutate(p=cumsum(counts)/N)%>%
      mutate(etheta=ifelse(qnorm(p)==Inf,4,qnorm(p)))%>%select(itemid,totalni,etheta)%>%
      ungroup()
    bsdata2<-left_join(bsdata,bsetheta,by=c('itemid','totalni'))
    bsdata3<-do.call("rbind", replicate(dim(tn)[1], bsdata2, simplify = FALSE))
    bstn<-do.call("rbind", replicate(dim(bsdata2)[1], tn, simplify = FALSE))%>%
      arrange(thetanode)
    bsdata4<-cbind(bsdata3,bstn)%>%mutate(w=dnorm((etheta-thetanode)/h))%>%
      mutate(wr1=ifelse(score>=1,w,0),
             wr2=ifelse(score>=2,w,0),
             wr3=ifelse(score>=3,w,0),
             wr4=ifelse(score>=4,w,0))%>%
      group_by(itemid,thetanode)%>%
      summarize(ks1=sum(wr1)/sum(w),
                ks2=sum(wr2)/sum(w),
                ks3=sum(wr3)/sum(w),
                ks4=sum(wr4)/sum(w))%>%ungroup()
    bspRISE<-left_join(estmodelp,bsdata4,by=c('itemid','thetanode'))%>%
      mutate(RISEp1=(p1-ks1)^2*dnorm(thetanode),
             RISEp2=(p2-ks2)^2*dnorm(thetanode),
             RISEp3=(p3-ks3)^2*dnorm(thetanode),
             RISEp4=(p4-ks4)^2*dnorm(thetanode))%>%
      group_by(itemid)%>%
      summarize(pRISE=sqrt(sum(RISEp1))+sqrt(sum(RISEp2))+sqrt(sum(RISEp3))+sqrt(sum(RISEp4)))%>%ungroup()
    bssinglep<-left_join(bsrRISE,bspRISE,by=c('itemid'))%>%
      mutate(count=ifelse(rRISE<pRISE,1,0))%>%select(itemid,count)%>%
      mutate(index=i)
    bsallp<-rbind(bsallp,bssinglep)
    
  }
  #endtime=Sys.time()
  #runtime: endtime-starttime Time difference of 5.687642 mins
  bsallp<-bsallp%>%group_by(itemid)%>%summarise(bsp=sum(count)/100)
  bsp<-cbind(bsp,bsallp$bsp)
  
  
  #---step 4 : PPMC to get ppp-value---------------------------------------
  x<-as.matrix(cbind(data[,2:(n*0.8+1)],data[,(2+n*0.8):(n+1)]+1))
  #---------constant, combinde data------------------
  irtConsts<-list(N=N,n=n,ndct=n*0.8,K=5,y=x)
  #-----------initial--------------
  irtInits<-list(b=rep(0,n*0.8),a=rep(1,n),theta=rep(0,N), bk.star=t(replicate(n*0.2,c(-2,-1,1,2))))
  #-----------start calibration-----------------------
  irt.mcmcsamples<-nimbleMCMC(code=irtCode,constants = irtConsts,inits=irtInits,
                              niter =1000,nburnin =500,thin = 2,nchains=2,summary=TRUE,WAIC=TRUE)
  #-------judge convergence and set parameters---------------
  
  #-----------get estimation----------------------------
  samples1<-as.matrix(irt.mcmcsamples$samples[[1]])[,1:(n+1.6*n)]
  mcmcsample1<-mcmc(samples1)
  samples2<-as.matrix(irt.mcmcsamples$samples[[2]])[,1:(n+1.6*n)]
  mcmcsample2<-mcmc(samples2)
  sample<-as.data.frame(rbind(as.matrix(mcmcsample1),as.matrix(mcmcsample2)))
  mcmcsample<-mcmc.list(mcmcsample1,mcmcsample2)
  
  gelman.diag(mcmcsample, confidence = 0.95, transform=TRUE, autoburnin=FALSE,multivariate=F)
  #---------calcualte PPP-Value-----------------------
  #----------getting item pool parameter----------------
  a<-sample[,1:n]%>%mutate(index=row_number())
  colnames(a)<-c(paste0("a",1:n),'index')
  a2<-melt(a,id='index')%>%
    mutate(itemid=substr(as.character(variable),2,nchar(as.character(variable))))%>%
    select(index,itemid,a=value)
  b1<-sample[,(1+n):(n+n)]%>%mutate(index=row_number())
  colnames(b1)<-c(paste0("b1",1:n),'index')
  b12<-melt(b1,id='index')%>%
    mutate(itemid=substr(as.character(variable),3,nchar(as.character(variable))))%>%
    select(index,itemid,b1=value)
  b2<-sample[,(1+2*n):(2*n+n*0.2)]%>%mutate(index=row_number())
  colnames(b2)<-c(paste0("b2",(n*0.8+1):n),'index')
  b22<-melt(b2,id='index')%>%
    mutate(itemid=substr(as.character(variable),3,nchar(as.character(variable))))%>%
    select(index,itemid,b2=value)
  b3<-sample[,(1+2*n+n*0.2):(2*n+n*0.2*2)]%>%mutate(index=row_number())
  colnames(b3)<-c(paste0("b3",(n*0.8+1):n),'index')
  b32<-melt(b3,id='index')%>%
    mutate(itemid=substr(as.character(variable),3,nchar(as.character(variable))))%>%
    select(index,itemid,b3=value)
  b4<-sample[,(1+2*n+n*0.2*2):(2*n+n*0.2*3)]%>%mutate(index=row_number())
  colnames(b4)<-c(paste0("b4",(n*0.8+1):n),'index')
  b42<-melt(b4,id='index')%>%
    mutate(itemid=substr(as.character(variable),3,nchar(as.character(variable))))%>%
    select(index,itemid,b4=value)
  itempool<-list(a2,b12,b22,b32,b42) %>%
    Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("index","itemid")), .)
  rm(a,a2,b1,b12,b2,b22,b3,b32,b4,b42)
  #---------clean pool------------
  b12<-itempool%>%filter(b1>b2)%>%mutate(b1=b2,b2=b1)
  other<-itempool%>%filter(!(b1>b2)|is.na(b2))
  itempool<-rbind(b12,other)
  b23<-itempool%>%filter(b2>b3)%>%mutate(b2=b3,b3=b2)
  other<-itempool%>%filter(!(b2>b3)|is.na(b2))
  itempool<-rbind(b23,other)
  b34<-itempool%>%filter(b3>b4)%>%mutate(b3=b4,b4=b3)
  other<-itempool%>%filter(!(b3>b4)|is.na(b2))
  itempool<-rbind(b34,other)
  rm(b12,b23,b34)
  #----------get theta parameter-------------------------
  theta<-sample[,(1+2*n+n*0.6):(2*n+n*0.6+N)]%>%mutate(index=row_number())
  colnames(theta)<-c(paste0("theta",1:N),'index')
  theta2<-melt(theta,id='index')%>%
    mutate(sid=substr(as.character(variable),6,nchar(as.character(variable))))%>%
    select(index,sid,theta=value)
  all<-left_join(theta2,itempool,by='index')%>%mutate(sid=as.numeric(sid))%>%
    mutate(itemid=as.numeric(itemid))
  rm(theta,theta2)
  #-------------------PPMC realized RISE-------------------
  itempool2<-do.call("rbind", replicate(dim(tn)[1], itempool, simplify = FALSE))
  tnip2<-do.call("rbind", replicate(dim(itempool)[1], tn, simplify = FALSE))%>%
    arrange(thetanode)
  modelp<-cbind(itempool2,tnip2)%>%mutate(itemid=as.numeric(itemid))%>%
    mutate(p1=1 / (1 + exp(-a*thetanode + a*b1)),
           p2=ifelse(is.na(b2),0,1 / (1 + exp(-a*thetanode + a*b2))),
           p3=ifelse(is.na(b3),0,1 / (1 + exp(-a*thetanode + a*b3))),
           p4=ifelse(is.na(b4),0,1 / (1 + exp(-a*thetanode + a*b4))))%>%
    select(index,itemid,thetanode,p1,p2,p3,p4)
  #%>%
  #mutate(p1=ifelse(itemid<=n/2,P1,P1-P2),
  #       p2=ifelse(itemid<=n/2,0,P2-P3),
  #       p3=ifelse(itemid<=n/2,0,P3-P4),
  #       p4=ifelse(itemid<=n/2,0,P4))%>%select(index,itemid,thetanode,p1,p2,p3,p4)
  rm(tnip2,itempool2)
  rRISE<-left_join(modelp,data3,by=c('itemid','thetanode'))%>%
    mutate(RISEp1=(p1-ks1)^2*dnorm(thetanode),
           RISEp2=(p2-ks2)^2*dnorm(thetanode),
           RISEp3=(p3-ks3)^2*dnorm(thetanode),
           RISEp4=(p4-ks4)^2*dnorm(thetanode))%>%
    group_by(index,itemid)%>%
    summarize(rRISE=sqrt(sum(RISEp1))+sqrt(sum(RISEp2))+sqrt(sum(RISEp3))+sqrt(sum(RISEp4)))%>%ungroup()
  rm(etheta,other,data3)
  #-------------------PPMC predictive RISE-------------------
  alldata<-all%>%
    mutate(p0=1,
           p1=1 / (1 + exp(-a*theta + a*b1)),
           p2=ifelse(is.na(b2),0,1 / (1 + exp(-a*theta + a*b2))),
           p3=ifelse(is.na(b3),0,1 / (1 + exp(-a*theta + a*b3))),
           p4=ifelse(is.na(b4),0,1 / (1 + exp(-a*theta + a*b4))),
           p5=0)
  alldata$mcp<-runif(dim(alldata)[1])
  alldata<-alldata%>%
    mutate(score=ifelse(mcp>p1,0,ifelse(mcp>p2,1,ifelse(mcp>p3,2,ifelse(mcp>p4,3,4)))))%>%
    select(index,sid,itemid,score)%>%
    group_by(index,sid)%>%mutate(total=sum(score))%>%
    ungroup()%>%mutate(totalni=total-score)
  alletheta<-alldata%>%
    group_by(index,itemid,totalni)%>%summarize(counts=n())%>%ungroup()%>%
    arrange(index,itemid,totalni)%>%group_by(index,itemid)%>%mutate(p=cumsum(counts)/N)%>%
    mutate(etheta=ifelse(qnorm(p)==Inf,4,qnorm(p)))%>%select(index,itemid,totalni,etheta)%>%
    ungroup()
  alldata3<-left_join(alldata,alletheta,by=c('index','itemid','totalni'))
  rm(all,alldata,alletheta)
  
  #starttime=Sys.time()
  alldata4<-NULL
  for (i in 1:dim(tn)[1]){
    thetanode=tn$thetanode[i]
    tnalldata3<-alldata3%>%
          mutate(w=dnorm((etheta-thetanode)/h))%>%
          mutate(wr1=ifelse(score>=1,w,0),
                 wr2=ifelse(score>=2,w,0),
                 wr3=ifelse(score>=3,w,0),
                 wr4=ifelse(score>=4,w,0))%>%
          group_by(index,itemid)%>%
          summarize(ks1=sum(wr1)/sum(w),
                    ks2=sum(wr2)/sum(w),
                    ks3=sum(wr3)/sum(w),
                    ks4=sum(wr4)/sum(w))%>%
      mutate(thetanode=thetanode)%>%ungroup()
    alldata4<-rbind(alldata4,tnalldata3)
  }
  #endtime=Sys.time()
  #runtime: endtime-starttime Time difference of 2.830743 mins
  
 
  
  pRISE<-left_join(modelp,alldata4,by=c('index','itemid','thetanode'))%>%
    mutate(RISEp1=(p1-ks1)^2*dnorm(thetanode),
           RISEp2=(p2-ks2)^2*dnorm(thetanode),
           RISEp3=(p3-ks3)^2*dnorm(thetanode),
           RISEp4=(p4-ks4)^2*dnorm(thetanode))%>%
    group_by(index,itemid)%>%
    summarize(pRISE=sqrt(sum(RISEp1))+sqrt(sum(RISEp2))+sqrt(sum(RISEp3))+sqrt(sum(RISEp4)))%>%ungroup()
  
  #---------calcualte ppp-value---------------------
  pppvalue<-left_join(rRISE,pRISE,by=c('index','itemid'))%>%
    group_by(itemid)%>%summarize(pppv=sum(rRISE<pRISE)/100)%>%ungroup()%>%arrange(itemid)
  
  pppv<-cbind(pppv,pppvalue$pppv)
  
}

#-------judge convergence and set parameters---------------
irt.mcmcsamples<-nimbleMCMC(code=irtCode,constants = irtConsts,inits=irtInits,
                            niter = 15000,nburnin = 2000,thin = 4,nchains=2,summary=TRUE,WAIC=TRUE)

samples1<-as.matrix(irt.mcmcsamples$samples[[1]])[,1:(n+1.6*n)]
mcmcsample1<-mcmc(samples1,start=2751,end=3250)

samples2<-as.matrix(irt.mcmcsamples$samples[[2]])[,1:(n+1.6*n)]
mcmcsample2<-mcmc(samples2,start=2751,end=3250)

sampleforplot<- mcmc.list(mcmcsample1,mcmcsample2)

pdf(file=paste0('Data\\T',N,'P',n,"\\Cubic\\Cubic_T",N,"_P",n,"_rep",rep,'_traceplot.pdf'),width=6,height=8)
plot(sampleforplot)
dev.off()

pdf(file=paste0('Data\\T',N,'P',n,"\\Cubic\\Cubic_T",N,"_P",n,"_rep",rep,'_autocorrelationplot.pdf'),width=50,height=40)
acfplot(sampleforplot,lag.max=50)
dev.off()

#samples<-rbind(as.matrix(irt.mcmcsamples$samples[[1]]),as.matrix(irt.mcmcsamples$samples[[2]]))

samples1<-as.matrix(irt.mcmcsamples$samples[[1]])
mcmcsample1<-mcmc(samples1,start=2751,end=3250)

samples2<-as.matrix(irt.mcmcsamples$samples[[2]])
mcmcsample2<-mcmc(samples2,start=2751,end=3250)

samples<-rbind(as.matrix(mcmcsample1),as.matrix(mcmcsample2))

write.table(samples,paste0('Data\\T',N,'P',n,"\\Cubic\\Cubic_T",N,"_P",n,"_rep",rep,"_MCMCSamples.txt"),sep="\t",row.names=F,col.names=F)




write.table(pppvalue,paste0('Data\\T',N,'P',n,"\\NonMonotonic\\NonMonotonic_T",N,"_P",n,"_rep",rep,"_pppv.txt"),sep="\t",row.names=F,col.names=F)






