library("glmmTMB")
library(TMB) # http://github.com/kaskr/adcomp

# data from Supporting Information 11
load("experimental_data.RData")

# Add treatment (frozen or not)
db$tmt<-ifelse(db$Bucket_ID%in%c("5","6","7","8"), 1, 2)
db$tmtcombo<-factor(factor(db$Description):as.factor(db$tmt))

# recode time to reflect that "no fish" is after 120h
db$time<-ifelse(db$Description=="Blood_water_no_fish",paste0(as.integer(gsub("h","",db$Time_start))+120,"h"),db$Time_start)

# fractions instead of percent 
db$ObsFrac100 <- as.numeric(db$ObsFrac)/100
db$TrueFracW100 <- as.numeric(db$TrueFracW)/100

db$ctime<-as.integer(gsub("h","",db$time))
db$p<-as.integer(as.factor(db$plate_effect))

## Start model search

fit0 <- glmmTMB(ObsFrac100 ~ I(qlogis(TrueFracW100))+ctime*factor(tmtcombo)+(1|Bucket_ID), family=beta_family(link="logit"),data=db) # M

fit1 <- glmmTMB(ObsFrac100 ~ I(qlogis(TrueFracW100))+ctime+factor(tmtcombo)+(1|Bucket_ID), family=beta_family(link="logit"),data=db) # H1

1-pchisq(2*(logLik(fit0)-logLik(fit1)), 2) # cannot reduce pval 0.83%

# write model in a way where it is simpler to get estimates
fitM <- glmmTMB(ObsFrac100 ~ -1+I(qlogis(TrueFracW100))+factor(tmtcombo):ctime+factor(tmtcombo)+(1|Bucket_ID), family=beta_family(link="logit"),data=db) # M
summary(fitM)


# mackerel biomass prediction using model 1 (fitM)
#the raw eDNA quantities 

load("raw_genetic_data_landings.RData")

pdataBW72 <- data.frame(TrueFracW100=seq(0.001, 0.05, length=100), ctime=72, tmtcombo="Blood_water:2", Bucket_ID=NA)
pBW72<-predict(fitM, newdata=pdataBW72, se.fit=TRUE)

x=pdataBW72$TrueFracW100
y=plogis(pBW72$fit)
f<-splinefun(y,x)

sub$predictedWeight<-f(sub$obsFracW/100) 



#model 2
load("modeling_2_data.RData")

dat<-combi_wl
dat<-dat[!is.na(dat$catch_weight),]
dat$type <- factor(dat$type)
dat$tank <- factor(dat$tank)
dat$haul<- factor(c("1.1"=2, "1.2"=1, "1.3"=2, "2.1"=2, "2.2"=1, "3.1"=1, "3.2"=1, "4.2"=3, "4.3"=3)[dat$tank])

compile("self.cpp") # can be extended as much as we like 
dyn.load(dynlib("self"))

samp<-match(paste0(dat$catch_weight,"x",dat$type), unique(paste0(dat$catch_weight, "x", dat$type)))
samp[dat$type%in%c(3,4)]<-0
samp<-as.integer(as.factor(samp))-1


# script for no mixing
  sub<-dat[dat$type==2, ]
  idx<-rep(NA, nrow(sub))
  idx[1]<-1
  for(i in 2:length(idx))if(sub$haul[i]==sub$haul[i-1]){idx[i]<-idx[i-1]}else{idx[i]<-idx[i-1]+1}
  frac<-rep(NA, nrow(sub))
  frac[1]<- 1
  Cshift<- -101
  oldType<- 1
  prev<- 1
  weSee<-matrix(0, nrow=nrow(sub), ncol=length(unique(sub$haul)))
  weSee[1,1]<-1
  for(i in 2:length(frac)){
      if(idx[i]!=idx[i-1]){
          Cshift<-sub$catch_weight[i-1];
          prev<-sub$haul[i-1]
      }
      frac[i]<-min(1, (sub$catch_weight[i-1]-Cshift)/100)
      weSee[i,prev]<-(1-frac[i])
      weSee[i,sub$haul[i]]<-frac[i]
  }

weSeeExp<-matrix(0,nrow=nrow(dat),ncol=length(unique(dat$haul)))
weSeeExp[cbind(1:nrow(dat),dat$haul)]<-1
#weSeeExp[dat$type==2,]<-weSee      

  
mydat<-list(obs=dat$obs, type=as.integer(dat$type)-1, haul=as.integer(dat$haul)-1, samp=as.integer(samp)-1, weSee=weSeeExp)

ntype=length(unique(dat$type))
nhaul=length(unique(dat$haul))

mypar<-list(logitP=rep(0, nhaul),
            logitMu=matrix(0, nrow=ntype, ncol=nhaul),
            logitPhi=rep(0,ntype), logSigma=0, U=rep(0,max(samp)))

obj<-MakeADFun(mydat, mypar, DLL="self", random="U", map=list(logitPhi=factor(c(1,2,3))))

opt<-nlminb(obj$par, obj$fn, obj$gr)

sdr<-sdreport(obj)

frac<-as.list(sdr, "Est", report=TRUE)$fracTab
fracSD<-as.list(sdr, "Std", report=TRUE)$fracTab


# get estimates of the models for the total catch by taking into account the shares of the hauls to the total outcome

A <- c(0.5135135 ,0.3243243 ,0.1621622) #shares of the 3 hauls to the total catch
idx<-which(names(sdr$value)=="fracTab")
est<-A%*%sdr$value[idx[c(1,3,5)]] #genetic_ship
ship_sd<-sqrt(t(A)%*%sdr$cov[idx[c(1,3,5)],idx[c(1,3,5)]]%*%(A))
ship_ci<-est[1,1]+c(-2,2)*ship_sd[1,1]

est<-A%*%sdr$value[idx[c(2,4,6)]] #genetic_factory
factory_sd<-sqrt(t(A)%*%sdr$cov[idx[c(2,4,6)],idx[c(2,4,6)]]%*%(A))
factory_ci<-est[1,1]+c(-2,2)*factory_sd[1,1]
  
  
# script for no mixing
sub<-dat[dat$type==2, ]
idx<-rep(NA, nrow(sub))
idx[1]<-1
for(i in 2:length(idx))if(sub$haul[i]==sub$haul[i-1]){idx[i]<-idx[i-1]}else{idx[i]<-idx[i-1]+1}
frac<-rep(NA, nrow(sub))
frac[1]<- 1
Cshift<- -101
oldType<- 1
prev<- 1
weSee<-matrix(0, nrow=nrow(sub), ncol=length(unique(sub$haul)))
weSee[1,1]<-1
for(i in 2:length(frac)){
  if(idx[i]!=idx[i-1]){
    Cshift<-sub$catch_weight[i-1];
    prev<-sub$haul[i-1]
  }
  frac[i]<-min(1, (sub$catch_weight[i-1]-Cshift)/100)
  weSee[i,prev]<-(1-frac[i])
  weSee[i,sub$haul[i]]<-frac[i]
}

weSeeExp<-matrix(0,nrow=nrow(dat),ncol=length(unique(dat$haul)))
weSeeExp[cbind(1:nrow(dat),dat$haul)]<-1
weSeeExp[dat$type==2,]<-weSee      


mydat<-list(obs=dat$obs, type=as.integer(dat$type)-1, haul=as.integer(dat$haul)-1, samp=as.integer(samp)-1, weSee=weSeeExp)

ntype=length(unique(dat$type))
nhaul=length(unique(dat$haul))

mypar<-list(logitP=rep(0, nhaul),
            logitMu=matrix(0, nrow=ntype, ncol=nhaul),
            logitPhi=rep(0,ntype), logSigma=0, U=rep(0,max(samp)))

obj<-MakeADFun(mydat, mypar, DLL="self", random="U", map=list(logitPhi=factor(c(1,2,3))))

opt<-nlminb(obj$par, obj$fn, obj$gr)

sdr<-sdreport(obj)

frac<-as.list(sdr, "Est", report=TRUE)$fracTab
fracSD<-as.list(sdr, "Std", report=TRUE)$fracTab


# get estimates of the models for the total catch by taking into account the shares of the hauls to the total outcome

A <- c(0.5135135 ,0.3243243 ,0.1621622) #shares of the 3 hauls to the total catch
idx<-which(names(sdr$value)=="fracTab")
est<-A%*%sdr$value[idx[c(1,3,5)]] #genetic_ship
ship_sd<-sqrt(t(A)%*%sdr$cov[idx[c(1,3,5)],idx[c(1,3,5)]]%*%(A))
ship_ci<-est[1,1]+c(-2,2)*ship_sd[1,1]

est<-A%*%sdr$value[idx[c(2,4,6)]] #genetic_factory
factory_sd<-sqrt(t(A)%*%sdr$cov[idx[c(2,4,6)],idx[c(2,4,6)]]%*%(A))
factory_ci<-est[1,1]+c(-2,2)*factory_sd[1,1]


# script for total mixing
sub<-dat[dat$type==2, ]
weSee<-matrix(0, nrow=nrow(sub), ncol=length(unique(sub$haul)))
weSee[1,sub$haul[1]]<-sub$catch_weight[1]
for(i in 2:nrow(sub)){
  weSee[i,]<-weSee[i-1,]    
  weSee[i,sub$haul[i]]<-weSee[i,sub$haul[i]]+(sub$catch_weight[i]-sub$catch_weight[i-1])
}
weSee<-weSee/rowSums(weSee)

weSeeExp<-matrix(0,nrow=nrow(dat),ncol=length(unique(dat$haul)))
weSeeExp[cbind(1:nrow(dat),dat$haul)]<-1
weSeeExp[dat$type==2,]<-weSee        


mydat<-list(obs=dat$obs, type=as.integer(dat$type)-1, haul=as.integer(dat$haul)-1, samp=as.integer(samp)-1, weSee=weSeeExp)

ntype=length(unique(dat$type))
nhaul=length(unique(dat$haul))

mypar<-list(logitP=rep(0, nhaul),
            logitMu=matrix(0, nrow=ntype, ncol=nhaul),
            logitPhi=rep(0,ntype), logSigma=0, U=rep(0,max(samp)))

obj<-MakeADFun(mydat, mypar, DLL="self", random="U", map=list(logitPhi=factor(c(1,2,3,4))))

opt<-nlminb(obj$par, obj$fn, obj$gr)

sdr<-sdreport(obj)

frac<-as.list(sdr, "Est", report=TRUE)$fracTab
fracSD<-as.list(sdr, "Std", report=TRUE)$fracTab

# get estimates of the models for the total catch by taking into account the shares of the hauls to the total outcome

A <- c(0.5135135 ,0.3243243 ,0.1621622) #shares of the 3 hauls to the total catch
idx<-which(names(sdr$value)=="fracTab")
est<-A%*%sdr$value[idx[c(1,3,5)]] #genetic_ship
ship_sd<-sqrt(t(A)%*%sdr$cov[idx[c(1,3,5)],idx[c(1,3,5)]]%*%(A))
ship_ci<-est[1,1]+c(-2,2)*ship_sd[1,1]

est<-A%*%sdr$value[idx[c(2,4,6)]] #genetic_factory
factory_sd<-sqrt(t(A)%*%sdr$cov[idx[c(2,4,6)],idx[c(2,4,6)]]%*%(A))
factory_ci<-est[1,1]+c(-2,2)*factory_sd[1,1]
