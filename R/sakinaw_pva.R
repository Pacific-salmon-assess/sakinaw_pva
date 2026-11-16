#Sakinaw Sockeye - PVA
source('./R/funcs.R')

#datasets####
#smolt-spawner brood table, same dataset as the abundances, but where smolts in +2 years are aligned with each brood year of returning spawners
smolt_bt=read.csv('./data/smolt_broodtable.csv')
#main dataset - adult returns, fry/smolt releases, fry/smolt lake survival, marine survival, and big qualicum release/return data
data=read.csv('./data/sakinaw_abundances.csv')

#Natural dynamics####

#Ricker spawner-smolt curve#
#format data for samEst
smolt_bt=smolt_bt[,1:3]
colnames(smolt_bt)=c('by','S','R')
smolt_bt=smolt_bt[complete.cases(smolt_bt$by),]
smolt_bt$logRS=log(smolt_bt$R/smolt_bt$S)

#estimate autocorrelated Ricker model, extract key parameters
smolt_ric=samEst::ricker_TMB(smolt_bt,ac=TRUE,silent=T)
sak.alpha=c(smolt_ric$sd_report$value[1],smolt_ric$sd_report$sd[1])
sak.beta=c(smolt_ric$sd_report$value[2],smolt_ric$sd_report$sd[2])
sak.rho=c(smolt_ric$sd_report$value[3],smolt_ric$sd_report$sd[3])
sak.sigma=c(smolt_ric$sd_report$value[5],smolt_ric$sd_report$sd[5])

#randomly sample Ricker parameters from estimates and covariances
mu.ric=smolt_ric$sd_report$par.fixed
vcov.ric=smolt_ric$sd_report$cov.fixed

#Fry and Smolt releases####

#fry release
total.fry=apply(data[,8:9],1,sum,na.rm = TRUE)
names(total.fry)=data$ï..year

#use fry releases from 2018 - 2025 representing the new hatchery regime 
total.fry=total.fry[18:25]

sd.fry=sd(total.fry)
mu.fry=mean(total.fry)
sd.log.fry=sqrt(log(1+(sd.fry/mu.fry)^2))
mu.log.fry=mean(log(total.fry))-sd.log.fry^2



#steps 1b: smolt releases + survival to emigration
total.smolts=data$smolt_release
total.smolts=total.smolts[complete.cases(total.smolts)]

#natural scale estimates
mu.smolts=mean(total.smolts)
sd.smolts=sd(total.smolts)

#log-adjusted estimates
sd.log.smolts=sqrt(log(1+(sd.smolts/mu.smolts)^2))
mu.log.smolts=mean(log(total.smolts))-sd.log.smolts^2

mean(exp(rnorm(1e6,mu.log.smolts,sd.log.smolts)))
mu.smolts


#ratio of smolt-to-total releases
rels=data[,c(8:10)]
total.rels=apply(rels,1,sum,na.rm=T)
p.sm=rels$smolt_release/total.rels
p.sm=p.sm[complete.cases(p.sm)]

mu.logit.psm=mean(qlogis(p.sm))
sd.logit.psm=sd(qlogis(p.sm))

summary(plogis(rnorm(1e3,mu.logit.psm,sd.logit.psm)))

#fry survival####

#natural scale estimates
mu.fry.surv=mean(data$fry_to_smolt_surv/100,na.rm=T)
sd.fry.surv=sd(data$fry_to_smolt_surv/100,na.rm=T)

#logit-adjusted estimates
mu.logit.fry.surv=mean(qlogis(data$fry_to_smolt_surv/100),na.rm=T)
sd.logit.fry.surv=sd(qlogis(data$fry_to_smolt_surv/100),na.rm=T)

mean(plogis(rnorm(1e6,mu.logit.fry.surv,sd.logit.fry.surv))) #simulated
mu.fry.surv #empirical mean

#smolt survival to smolt emigration

#natural scale estimates
mu.smolt.surv=mean(data$smolt_release_emig_survival/100,na.rm=T)
sd.smolt.surv=sd(data$smolt_release_emig_survival/100,na.rm=T)

#logit-adjusted estimates
mu.logit.smolt.surv=mean(qlogis(data$smolt_release_emig_survival/100),na.rm=T)
sd.logit.smolt.surv=sd(qlogis(data$smolt_release_emig_survival/100),na.rm=T)

mean(plogis(rnorm(1000,mu.logit.smolt.surv,sd.logit.smolt.surv)))


#Marine survival####
#mean, sd, correlations among survival estimates
survivals=data[,c(15,17)] #extract hatchery (clipped) and 
surv.cor=cor(log(survivals+0.001),use='pairwise.complete.obs')
surv.sd=apply(survivals,2,sd,na.rm=T)
surv.mean=apply(survivals,2,mean,na.rm=T)
sd.log.surv=sqrt(log(1+(surv.sd/surv.mean)^2))
mu.log.surv=log(surv.mean)-0.5*sd.log.surv^2

#variance-covariance matrix for wild/hatchery/mixed marine survival
surv.mat=matrix(ncol=2,nrow=2)
surv.mat[1,1]=sd.log.surv[1]^2
surv.mat[2,2]=sd.log.surv[2]^2
surv.mat[1,2]=surv.cor[1,2]*sd.log.surv[1]*sd.log.surv[2]
surv.mat[2,1]=surv.cor[1,2]*sd.log.surv[1]*sd.log.surv[2]

ms.log=MASS::mvrnorm(n=1e3,mu=mu.log.surv,Sigma=surv.mat,empirical=TRUE)
ms.dists=exp(ms.log)
surv.mean
apply(ms.dists,2,mean)

#unknown

#autocorrelation in marine survival
rho.ms.h=acf(log(survivals[,1]+0.001),na.action=na.pass,plot=F)$acf[2] #1-yr autocorrelation - hatchery survival
rho.ms.n=acf(log(survivals[,2]+0.001),na.action=na.pass,plot=F)$acf[2] #1-yr autocorrelation - wild survival
se.rho.ms.h=1/sqrt(length(na.omit(survivals[,1]))) #standard error estimate for rho
se.rho.ms.n=1/sqrt(length(na.omit(survivals[,2]))) #standard error estimate for rho

#scale both quantities to logistic scale
logistic.rho.ms.h=logistic_transf(rho.ms.h) #one_to_minus_one transformation (logistic)
logistic.rho.ms.n=logistic_transf(rho.ms.n)
#SE using delta method
logistic.se.rho.ms.h=2/(1-rho.ms.h^2)*se.rho.ms.h
logistic.se.rho.ms.n=2/(1-rho.ms.n^2)*se.rho.ms.n

#check for consistency - hatchery rho
mean(inv_logistic_transf(rnorm(1e6,logistic.rho.ms.h,logistic.se.rho.ms.h)))
rho.ms.h
sd(inv_logistic_transf(rnorm(1e6,logistic.rho.ms.h,logistic.se.rho.ms.h)))
se.rho.ms.h

mean(inv_logistic_transf(rnorm(1e6,logistic.rho.ms.n,logistic.se.rho.ms.n)))
rho.ms.n
sd(inv_logistic_transf(rnorm(1e6,logistic.rho.ms.n,logistic.se.rho.ms.n)))
se.rho.ms.n




#Sim 1: Baseline Scenario####
sims1=list()
n.sims=1e3
sim.length=24
sim.summary1=data.frame(sim.id=seq(1:n.sims),py.tspwn.100=NA,py.nspwn.100=NA,py.tspwn.500=NA,py.nspwn.500=NA,m.tspwn=NA,m.nspwn=NA,extinct=NA)
for(i in 1:n.sims){
  sim.dat=data.frame(sim.year=seq(1:c(sim.length+2)),spawners.t=0,spawners.u=0,spawners.h=0,spawners.n=c(rep(NA,3),rep(0,sim.length-1)),smolts.u=NA,smolts.h=NA,smolts.n=NA,ms.h=NA,ms.n=NA,ms.u=NA,fry.rel=NA,smolt.rel=NA,fry.surv=NA,smolt.surv=NA,ricker.logrs=NA,ricker.mu=NA,ricker.eps=NA)
  #first 2 years - estimate returns based on known quantities (of e.g. fry and smolt releases, smolt outmigrants, etc) in the last several years
  
  #initial spawners - conditioned on observed data
  #total spawners
  sim.dat[1:2,2]=data$total_spawners[c(nrow(data)-1):nrow(data)]
  #unclipped spawners
  sim.dat[1:2,3]=data$unknown_origin_spawners[c(nrow(data)-1):nrow(data)]
  #clipped spawners
  sim.dat[1:2,4]=data$clipped_spawners[c(nrow(data)-1):nrow(data)]
  
  #initial smolt counts - conditioned on observed data
  #unclipped smolts (unc. fry and natural origin):
  sim.dat[1:2,6]=data$smolts_cnt_unclipped[c(nrow(data)-1):nrow(data)]
  #clipped smolts (both clp. fry and smolt releases):
  sim.dat[1:2,7]=data$smolts_cnt_clipped_fry[c(nrow(data)-1):nrow(data)]+data$smolt_cnt_release_clipped[c(nrow(data)-1):nrow(data)]
  
  #hatchery releases
  #total fry
  sim.dat[1:2,12]=c(data$fry_release_captive_brood_unclpd[c(nrow(data)-1):nrow(data)]+data$fry_release_captive_brood_clpd[c(nrow(data)-1):nrow(data)])
  
  #smolts
  sim.dat[1:2,13]=data$smolt_release[c(nrow(data)-1):nrow(data)]
  
  #hatchery lake survival
  #fry survival to smolt emigration
  sim.dat[1:2,14]=data$fry_to_smolt_surv[c(nrow(data)-1):nrow(data)]/100
  #smolt survival to smolt emigration
  sim.dat[1:2,15]=data$smolt_release_emig_survival[c(nrow(data)-1):nrow(data)]/100
  
  ##Simulate - marine survival for clipped/unclipped smolts - spawners t+2
  #hatchery and unknown origin survival for these smolt cohorts
  
  #simulate for marine survival for each set of smolts
  sim.dat[1:2,9:10]=exp(MASS::mvrnorm(n=2,mu=mu.log.surv,Sigma=surv.mat))
  sim.dat[1:2,11]=exp(rnorm(n=2,mu.log.u.ms,sd.log.u.ms)) #unknown marine survival for first 2 unclipped smolt cohorts
  
  #estimate returning spawners at t+2 for each smolt cohort
  #unclipped spawners - unclipped outmigrating smolts x unknown origin marine survival
  sim.dat[3:4,3]=sim.dat[1:2,6]*(sim.dat[1:2,11]/100)
  #clipped spawners - clipped outmigrating smolts x hatchery origin marine survival
  sim.dat[3:4,4]=sim.dat[1:2,7]*(sim.dat[1:2,9]/100)
  #total spawners
  sim.dat[3:4,2]=sim.dat[3:4,3]+sim.dat[3:4,4]
  
  #set Ricker parameters for the iteration
  sim.ric=MASS::mvrnorm(n=1,mu=mu.ric,Sigma=vcov.ric)
  sim.ric[2:3]=exp(sim.ric[2:3]) #exponentiate from log-scale
  sim.ric[4]=inv_logistic_transf(sim.ric[4]) #function to scale transformed rho to -1 to 1
  
  #for all future states - next 25 years
  for(t in 1:c(sim.length+2)){
    
    sim.dat[t,2]=sum(sim.dat[t,3:5],na.rm=T) #total spawners
    
    #natural origin smolts from total spawners:
    
    #estimate natural origin smolts t+2 based on Ricker model
    if(t==1){
      #no autocorrelation if first smolt cohort (no previous estimate available to condition autocorrelation in smolt recruitment)
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2],sd=sim.ric[3])
    }else{
      #autocorrelation included
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2]+sim.ric[4]*sim.dat[t-1,18],sd=sim.ric[3])
    }
    #if total spawners < 2 = no natural smolt production
    if(sim.dat[t,2]<2){
      sim.dat[t+2,8]=0
    }else{
      sim.dat[t+2,8]=exp(logRS)*sim.dat[t,2] #natural-origin smolts outmigrating in t+2
    }
    sim.dat[t,16]=logRS
    sim.dat[t,17]=sim.ric[1]-sim.dat[t,2]/sim.ric[2]
    sim.dat[t,18]=logRS-(sim.ric[1]-sim.dat[t,2]/sim.ric[2])
    
    #fry & smolt releases - after first 2 burn-in years
    if(t>2){
      #fry releases - total fry releases x allocation to clipped/unclipped status
      total_fry.t=exp(rnorm(1,mu.log.fry,sd.log.fry)) #log-normal draws for fry releases
      total_fry.t=min(total_fry.t,2e6) #set capacity for fry release of 2 million
      sim.dat[t,12]=total_fry.t
      sim.dat[t,13]=exp(rnorm(1,mu.log.smolts,sd.log.smolts)) #log-normal draws for smolt releases
      
      #fry/smolt survival from lake to outmigration
      sim.dat[t,14]=plogis(rnorm(1,mu.logit.fry.surv,sd.logit.fry.surv)) #fry survival - applies to both clipped/unclipped
      sim.dat[t,15]=plogis(rnorm(1,mu.logit.smolt.surv,sd.logit.smolt.surv)) #smolt release survival to outmigration
      sim.dat[t,7]=sim.dat[t-1,12]*sim.dat[t,14]+sim.dat[t-1,13]*sim.dat[t,15] #fry releases x fry survival + smolt releases x smolt survival
      
      #expected returns (t+2) based on marine survival - natural origin & hatchery origin fish
      rho.ms.h.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.h,logistic.se.rho.ms.h)) #draw autocorrelation coefficient for hatchery origin smolts
      rho.ms.n.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.n,logistic.se.rho.ms.n)) #draw autocorrelation coefficient for natural origin smolts
      
      mu.adj=c(rho.ms.h.sim*(log(sim.dat[t-1,9])-log(surv.mean)[1]), #adjustment to mean survival based on 1-y autocorrelation
               rho.ms.n.sim*(log(sim.dat[t-1,10])-log(surv.mean)[2])) #adjustment to mean survival based on 1-y autocorrelation
      
      sim.dat[t,9:10]=exp(MASS::mvrnorm(n=1,mu=mu.log.surv+mu.adj,Sigma=surv.mat)) #marine survival
      
      #partition based on age at maturity:
      p_age=extraDistr::rdirichlet(1,c(3,87,10))
      #3% age 3, 87% age 4, 10% age 5
      #proportion recruits, add to each return year based on randomized maturity
      #age 3
      sim.dat[t+1,4]=sim.dat[t+1,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[1] #hatchery returns
      sim.dat[t+1,5]=sim.dat[t+1,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[1] #natural origin returns
      #age 4
      sim.dat[t+2,4]=sim.dat[t+2,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[2] #hatchery returns
      sim.dat[t+2,5]=sim.dat[t+2,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[2] #natural origin returns
      #age 5
      sim.dat[t+3,4]=sim.dat[t+3,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[3] #hatchery returns
      sim.dat[t+3,5]=sim.dat[t+3,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[3] #natural origin returns
    }
  }
  sim.dat=sim.dat[1:c(sim.length+2),] #trim extra brood years
  sim.summary1[i,2]=sum(sim.dat[5:c(sim.length+2),2]>100)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where total spawners > 500
  sim.summary1[i,3]=sum(sim.dat[5:c(sim.length+2),5]>100)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary1[i,4]=sum(sim.dat[5:c(sim.length+2),2]>500)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where total spawners > 500
  sim.summary1[i,5]=sum(sim.dat[5:c(sim.length+2),5]>500)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary1[i,6]=exp(mean(log(sim.dat[5:c(sim.length+2),2])))
  sim.summary1[i,7]=exp(mean(log(sim.dat[5:c(sim.length+2),5])))
  sim.summary1[i,8]=ifelse(all(sim.dat[c(sim.length+2-4):c(sim.length+2),2]==0),1,0)
  
  sims1[[i]]=sim.dat #store simulation
}

#Sim 2: Doubled Hatchery Output####

#adjust total mean fry
mu.fry2=mean(total.fry*2)
sd.fry2=sd(total.fry*2)
sd.log.fry2=sqrt(log(1+(sd.fry2/mu.fry2)^2))
mu.log.fry2=mean(log(total.fry*2))-sd.log.fry2^2

mu.smolts2=mean(total.smolts*2)
sd.smolts2=sd(total.smolts*2)
sd.log.smolts2=sqrt(log(1+(sd.smolts2/mu.smolts2)^2))
mu.log.smolts2=mean(log(total.smolts*2))-sd.log.smolts2^2


sims2=list()
n.sims=1e3
sim.length=24
sim.summary2=data.frame(sim.id=seq(1:n.sims),py.tspwn.100=NA,py.nspwn.100=NA,py.tspwn.500=NA,py.nspwn.500=NA,m.tspwn=NA,m.nspwn=NA,extinct=NA)
for(i in 1:n.sims){
  sim.dat=data.frame(sim.year=seq(1:c(sim.length+2)),spawners.t=0,spawners.u=0,spawners.h=0,spawners.n=c(rep(NA,3),rep(0,sim.length-1)),smolts.u=NA,smolts.h=NA,smolts.n=NA,ms.h=NA,ms.n=NA,ms.u=NA,fry.rel=NA,smolt.rel=NA,fry.surv=NA,smolt.surv=NA,ricker.logrs=NA,ricker.mu=NA,ricker.eps=NA)
  #first 2 years - estimate returns based on known quantities (of e.g. fry and smolt releases, smolt outmigrants, etc) in the last several years
  
  #initial spawners - conditioned on observed data
  #total spawners
  sim.dat[1:2,2]=data$total_spawners[c(nrow(data)-1):nrow(data)]
  #unclipped spawners
  sim.dat[1:2,3]=data$unknown_origin_spawners[c(nrow(data)-1):nrow(data)]
  #clipped spawners
  sim.dat[1:2,4]=data$clipped_spawners[c(nrow(data)-1):nrow(data)]
  
  #initial smolt counts - conditioned on observed data
  #unclipped smolts (unc. fry and natural origin):
  sim.dat[1:2,6]=data$smolts_cnt_unclipped[c(nrow(data)-1):nrow(data)]
  #clipped smolts (both clp. fry and smolt releases):
  sim.dat[1:2,7]=data$smolts_cnt_clipped_fry[c(nrow(data)-1):nrow(data)]+data$smolt_cnt_release_clipped[c(nrow(data)-1):nrow(data)]
  
  #hatchery releases
  #total fry
  sim.dat[1:2,12]=c(data$fry_release_captive_brood_unclpd[c(nrow(data)-1):nrow(data)]+data$fry_release_captive_brood_clpd[c(nrow(data)-1):nrow(data)])
  
  #smolts
  sim.dat[1:2,13]=data$smolt_release[c(nrow(data)-1):nrow(data)]
  
  #hatchery lake survival
  #fry survival to smolt emigration
  sim.dat[1:2,14]=data$fry_to_smolt_surv[c(nrow(data)-1):nrow(data)]/100
  #smolt survival to smolt emigration
  sim.dat[1:2,15]=data$smolt_release_emig_survival[c(nrow(data)-1):nrow(data)]/100
  
  ##Simulate - marine survival for clipped/unclipped smolts - spawners t+2
  #hatchery and unknown origin survival for these smolt cohorts
  
  #simulate for marine survival for each set of smolts
  sim.dat[1:2,9:10]=exp(MASS::mvrnorm(n=2,mu=mu.log.surv,Sigma=surv.mat))
  sim.dat[1:2,11]=exp(rnorm(n=2,mu.log.u.ms,sd.log.u.ms)) #unknown marine survival for first 2 unclipped smolt cohorts
  
  #estimate returning spawners at t+2 for each smolt cohort
  #unclipped spawners - unclipped outmigrating smolts x unknown origin marine survival
  sim.dat[3:4,3]=sim.dat[1:2,6]*(sim.dat[1:2,11]/100)
  #clipped spawners - clipped outmigrating smolts x hatchery origin marine survival
  sim.dat[3:4,4]=sim.dat[1:2,7]*(sim.dat[1:2,9]/100)
  #total spawners
  sim.dat[3:4,2]=sim.dat[3:4,3]+sim.dat[3:4,4]
  
  #set Ricker parameters for the iteration
  sim.ric=MASS::mvrnorm(n=1,mu=mu.ric,Sigma=vcov.ric)
  sim.ric[2:3]=exp(sim.ric[2:3]) #exponentiate from log-scale
  sim.ric[4]=inv_logistic_transf(sim.ric[4]) #function to scale transformed rho to -1 to 1
  
  #for all future states - next 25 years
  for(t in 1:c(sim.length+2)){
    
    sim.dat[t,2]=sum(sim.dat[t,3:5],na.rm=T) #total spawners
    
    #natural origin smolts from total spawners:
    
    #estimate natural origin smolts t+2 based on Ricker model
    if(t==1){
      #no autocorrelation if first smolt cohort (no previous estimate available to condition autocorrelation in smolt recruitment)
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2],sd=sim.ric[3])
    }else{
      #autocorrelation included
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2]+sim.ric[4]*sim.dat[t-1,18],sd=sim.ric[3])
    }
    #if total spawners < 2 = no natural smolt production
    if(sim.dat[t,2]<2){
      sim.dat[t+2,8]=0
    }else{
      sim.dat[t+2,8]=exp(logRS)*sim.dat[t,2] #natural-origin smolts outmigrating in t+2
    }
    sim.dat[t,16]=logRS
    sim.dat[t,17]=sim.ric[1]-sim.dat[t,2]/sim.ric[2]
    sim.dat[t,18]=logRS-(sim.ric[1]-sim.dat[t,2]/sim.ric[2])
    
    #fry & smolt releases - after first 2 burn-in years
    if(t>2){
      #fry releases - total fry releases x allocation to clipped/unclipped status
      total_fry.t=exp(rnorm(1,mu.log.fry2,sd.log.fry2)) #log-normal draws for fry releases
      total_fry.t=min(total_fry.t,2e6) #set capacity for fry release of 2 million
      sim.dat[t,12]=total_fry.t
      sim.dat[t,13]=exp(rnorm(1,mu.log.smolts2,sd.log.smolts2)) #log-normal draws for smolt releases
      
      #fry/smolt survival from lake to outmigration
      sim.dat[t,14]=plogis(rnorm(1,mu.logit.fry.surv,sd.logit.fry.surv)) #fry survival - applies to both clipped/unclipped
      sim.dat[t,15]=plogis(rnorm(1,mu.logit.smolt.surv,sd.logit.smolt.surv)) #smolt release survival to outmigration
      sim.dat[t,7]=sim.dat[t-1,12]*sim.dat[t,14]+sim.dat[t-1,13]*sim.dat[t,15] #fry releases x fry survival + smolt releases x smolt survival
      
      #expected returns (t+2) based on marine survival - natural origin & hatchery origin fish
      rho.ms.h.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.h,logistic.se.rho.ms.h)) #draw autocorrelation coefficient for hatchery origin smolts
      rho.ms.n.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.n,logistic.se.rho.ms.n)) #draw autocorrelation coefficient for natural origin smolts
      
      mu.adj=c(rho.ms.h.sim*(log(sim.dat[t-1,9])-log(surv.mean)[1]), #adjustment to mean survival based on 1-y autocorrelation
               rho.ms.n.sim*(log(sim.dat[t-1,10])-log(surv.mean)[2])) #adjustment to mean survival based on 1-y autocorrelation
      
      sim.dat[t,9:10]=exp(MASS::mvrnorm(n=1,mu=mu.log.surv+mu.adj,Sigma=surv.mat)) #marine survival
      
      #partition based on age at maturity:
      p_age=extraDistr::rdirichlet(1,c(3,87,10))
      #3% age 3, 87% age 4, 10% age 5
      #proportion recruits, add to each return year based on randomized maturity
      #age 3
      sim.dat[t+1,4]=sim.dat[t+1,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[1] #hatchery returns
      sim.dat[t+1,5]=sim.dat[t+1,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[1] #natural origin returns
      #age 4
      sim.dat[t+2,4]=sim.dat[t+2,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[2] #hatchery returns
      sim.dat[t+2,5]=sim.dat[t+2,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[2] #natural origin returns
      #age 5
      sim.dat[t+3,4]=sim.dat[t+3,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[3] #hatchery returns
      sim.dat[t+3,5]=sim.dat[t+3,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[3] #natural origin returns
    }
  }
  sim.dat=sim.dat[1:c(sim.length+2),] #trim extra brood years
  sim.summary2[i,2]=sum(sim.dat[5:c(sim.length+2),2]>100)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where total spawners > 500
  sim.summary2[i,3]=sum(sim.dat[5:c(sim.length+2),5]>100)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary2[i,4]=sum(sim.dat[5:c(sim.length+2),2]>500)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where total spawners > 500
  sim.summary2[i,5]=sum(sim.dat[5:c(sim.length+2),5]>500)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary2[i,6]=exp(mean(log(sim.dat[5:c(sim.length+2),2])))
  sim.summary2[i,7]=exp(mean(log(sim.dat[5:c(sim.length+2),5])))
  sim.summary2[i,8]=ifelse(all(sim.dat[c(sim.length+2-4):c(sim.length+2),2]==0),1,0)
  
  sims2[[i]]=sim.dat #store simulation
}

#Sim 3: Halved Hatchery Output####
#adjust total mean fry
mu.fry2=mean(total.fry/2)
sd.fry2=sd(total.fry/2)
sd.log.fry2=sqrt(log(1+(sd.fry2/mu.fry2)^2))
mu.log.fry2=mean(log(total.fry/2))-sd.log.fry2^2

mu.smolts2=mean(total.smolts/2)
sd.smolts2=sd(total.smolts/2)
sd.log.smolts2=sqrt(log(1+(sd.smolts2/mu.smolts2)^2))
mu.log.smolts2=mean(log(total.smolts/2))-sd.log.smolts2^2

sims3=list()
n.sims=1e3
sim.length=24
sim.summary3=data.frame(sim.id=seq(1:n.sims),py.tspwn.100=NA,py.nspwn.100=NA,py.tspwn.500=NA,py.nspwn.500=NA,m.tspwn=NA,m.nspwn=NA,extinct=NA)
for(i in 1:n.sims){
  sim.dat=data.frame(sim.year=seq(1:c(sim.length+2)),spawners.t=0,spawners.u=0,spawners.h=0,spawners.n=c(rep(NA,3),rep(0,sim.length-1)),smolts.u=NA,smolts.h=NA,smolts.n=NA,ms.h=NA,ms.n=NA,ms.u=NA,fry.rel=NA,smolt.rel=NA,fry.surv=NA,smolt.surv=NA,ricker.logrs=NA,ricker.mu=NA,ricker.eps=NA)
  #first 2 years - estimate returns based on known quantities (of e.g. fry and smolt releases, smolt outmigrants, etc) in the last several years
  
  #initial spawners - conditioned on observed data
  #total spawners
  sim.dat[1:2,2]=data$total_spawners[c(nrow(data)-1):nrow(data)]
  #unclipped spawners
  sim.dat[1:2,3]=data$unknown_origin_spawners[c(nrow(data)-1):nrow(data)]
  #clipped spawners
  sim.dat[1:2,4]=data$clipped_spawners[c(nrow(data)-1):nrow(data)]
  
  #initial smolt counts - conditioned on observed data
  #unclipped smolts (unc. fry and natural origin):
  sim.dat[1:2,6]=data$smolts_cnt_unclipped[c(nrow(data)-1):nrow(data)]
  #clipped smolts (both clp. fry and smolt releases):
  sim.dat[1:2,7]=data$smolts_cnt_clipped_fry[c(nrow(data)-1):nrow(data)]+data$smolt_cnt_release_clipped[c(nrow(data)-1):nrow(data)]
  
  #hatchery releases
  #total fry
  sim.dat[1:2,12]=c(data$fry_release_captive_brood_unclpd[c(nrow(data)-1):nrow(data)]+data$fry_release_captive_brood_clpd[c(nrow(data)-1):nrow(data)])
  
  #smolts
  sim.dat[1:2,13]=data$smolt_release[c(nrow(data)-1):nrow(data)]
  
  #hatchery lake survival
  #fry survival to smolt emigration
  sim.dat[1:2,14]=data$fry_to_smolt_surv[c(nrow(data)-1):nrow(data)]/100
  #smolt survival to smolt emigration
  sim.dat[1:2,15]=data$smolt_release_emig_survival[c(nrow(data)-1):nrow(data)]/100
  
  ##Simulate - marine survival for clipped/unclipped smolts - spawners t+2
  #hatchery and unknown origin survival for these smolt cohorts
  
  #simulate for marine survival for each set of smolts
  sim.dat[1:2,9:10]=exp(MASS::mvrnorm(n=2,mu=mu.log.surv,Sigma=surv.mat))
  sim.dat[1:2,11]=exp(rnorm(n=2,mu.log.u.ms,sd.log.u.ms)) #unknown marine survival for first 2 unclipped smolt cohorts
  
  #estimate returning spawners at t+2 for each smolt cohort
  #unclipped spawners - unclipped outmigrating smolts x unknown origin marine survival
  sim.dat[3:4,3]=sim.dat[1:2,6]*(sim.dat[1:2,11]/100)
  #clipped spawners - clipped outmigrating smolts x hatchery origin marine survival
  sim.dat[3:4,4]=sim.dat[1:2,7]*(sim.dat[1:2,9]/100)
  #total spawners
  sim.dat[3:4,2]=sim.dat[3:4,3]+sim.dat[3:4,4]
  
  #set Ricker parameters for the iteration
  sim.ric=MASS::mvrnorm(n=1,mu=mu.ric,Sigma=vcov.ric)
  sim.ric[2:3]=exp(sim.ric[2:3]) #exponentiate from log-scale
  sim.ric[4]=inv_logistic_transf(sim.ric[4]) #function to scale transformed rho to -1 to 1
  
  #for all future states - next 25 years
  for(t in 1:c(sim.length+2)){
    
    sim.dat[t,2]=sum(sim.dat[t,3:5],na.rm=T) #total spawners
    
    #natural origin smolts from total spawners:
    
    #estimate natural origin smolts t+2 based on Ricker model
    if(t==1){
      #no autocorrelation if first smolt cohort (no previous estimate available to condition autocorrelation in smolt recruitment)
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2],sd=sim.ric[3])
    }else{
      #autocorrelation included
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2]+sim.ric[4]*sim.dat[t-1,18],sd=sim.ric[3])
    }
    #if total spawners < 2 = no natural smolt production
    if(sim.dat[t,2]<2){
      sim.dat[t+2,8]=0
    }else{
      sim.dat[t+2,8]=exp(logRS)*sim.dat[t,2] #natural-origin smolts outmigrating in t+2
    }
    sim.dat[t,16]=logRS
    sim.dat[t,17]=sim.ric[1]-sim.dat[t,2]/sim.ric[2]
    sim.dat[t,18]=logRS-(sim.ric[1]-sim.dat[t,2]/sim.ric[2])
    
    #fry & smolt releases - after first 2 burn-in years
    if(t>2){
      #fry releases - total fry releases x allocation to clipped/unclipped status
      total_fry.t=exp(rnorm(1,mu.log.fry2,sd.log.fry2)) #log-normal draws for fry releases
      total_fry.t=min(total_fry.t,2e6) #set capacity for fry release of 2 million
      sim.dat[t,12]=total_fry.t
      sim.dat[t,13]=exp(rnorm(1,mu.log.smolts2,sd.log.smolts2)) #log-normal draws for smolt releases
      
      #fry/smolt survival from lake to outmigration
      sim.dat[t,14]=plogis(rnorm(1,mu.logit.fry.surv,sd.logit.fry.surv)) #fry survival - applies to both clipped/unclipped
      sim.dat[t,15]=plogis(rnorm(1,mu.logit.smolt.surv,sd.logit.smolt.surv)) #smolt release survival to outmigration
      sim.dat[t,7]=sim.dat[t-1,12]*sim.dat[t,14]+sim.dat[t-1,13]*sim.dat[t,15] #fry releases x fry survival + smolt releases x smolt survival
      
      #expected returns (t+2) based on marine survival - natural origin & hatchery origin fish
      rho.ms.h.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.h,logistic.se.rho.ms.h)) #draw autocorrelation coefficient for hatchery origin smolts
      rho.ms.n.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.n,logistic.se.rho.ms.n)) #draw autocorrelation coefficient for natural origin smolts
      
      mu.adj=c(rho.ms.h.sim*(log(sim.dat[t-1,9])-log(surv.mean)[1]), #adjustment to mean survival based on 1-y autocorrelation
               rho.ms.n.sim*(log(sim.dat[t-1,10])-log(surv.mean)[2])) #adjustment to mean survival based on 1-y autocorrelation
      
      sim.dat[t,9:10]=exp(MASS::mvrnorm(n=1,mu=mu.log.surv+mu.adj,Sigma=surv.mat)) #marine survival
      
      #partition based on age at maturity:
      p_age=extraDistr::rdirichlet(1,c(3,87,10))
      #3% age 3, 87% age 4, 10% age 5
      #proportion recruits, add to each return year based on randomized maturity
      #age 3
      sim.dat[t+1,4]=sim.dat[t+1,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[1] #hatchery returns
      sim.dat[t+1,5]=sim.dat[t+1,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[1] #natural origin returns
      #age 4
      sim.dat[t+2,4]=sim.dat[t+2,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[2] #hatchery returns
      sim.dat[t+2,5]=sim.dat[t+2,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[2] #natural origin returns
      #age 5
      sim.dat[t+3,4]=sim.dat[t+3,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[3] #hatchery returns
      sim.dat[t+3,5]=sim.dat[t+3,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[3] #natural origin returns
    }
  }
  sim.dat=sim.dat[1:c(sim.length+2),] #trim extra brood years
  sim.summary3[i,2]=sum(sim.dat[5:c(sim.length+2),2]>100)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where total spawners > 500
  sim.summary3[i,3]=sum(sim.dat[5:c(sim.length+2),5]>100)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary3[i,4]=sum(sim.dat[5:c(sim.length+2),2]>500)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where total spawners > 500
  sim.summary3[i,5]=sum(sim.dat[5:c(sim.length+2),5]>500)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary3[i,6]=exp(mean(log(sim.dat[5:c(sim.length+2),2])))
  sim.summary3[i,7]=exp(mean(log(sim.dat[5:c(sim.length+2),5])))
  sim.summary3[i,8]=ifelse(all(sim.dat[c(sim.length+2-4):c(sim.length+2),2]==0),1,0)
  
  sims3[[i]]=sim.dat #store simulation
}

#Sim 4: Quartered hatchery output####
mu.fry2=mean(total.fry/4)
sd.fry2=sd(total.fry/4)
sd.log.fry2=sqrt(log(1+(sd.fry2/mu.fry2)^2))
mu.log.fry2=mean(log(total.fry/4))-sd.log.fry2^2

mu.smolts2=mean(total.smolts/4)
sd.smolts2=sd(total.smolts/4)
sd.log.smolts2=sqrt(log(1+(sd.smolts2/mu.smolts2)^2))
mu.log.smolts2=mean(log(total.smolts/4))-sd.log.smolts2^2

sims4=list()
n.sims=1e3
sim.length=24
sim.summary4=data.frame(sim.id=seq(1:n.sims),py.tspwn.100=NA,py.nspwn.100=NA,py.tspwn.500=NA,py.nspwn.500=NA,m.tspwn=NA,m.nspwn=NA,extinct=NA)
for(i in 1:n.sims){
  sim.dat=data.frame(sim.year=seq(1:c(sim.length+2)),spawners.t=0,spawners.u=0,spawners.h=0,spawners.n=c(rep(NA,3),rep(0,sim.length-1)),smolts.u=NA,smolts.h=NA,smolts.n=NA,ms.h=NA,ms.n=NA,ms.u=NA,fry.rel=NA,smolt.rel=NA,fry.surv=NA,smolt.surv=NA,ricker.logrs=NA,ricker.mu=NA,ricker.eps=NA)
  #first 2 years - estimate returns based on known quantities (of e.g. fry and smolt releases, smolt outmigrants, etc) in the last several years
  
  #initial spawners - conditioned on observed data
  #total spawners
  sim.dat[1:2,2]=data$total_spawners[c(nrow(data)-1):nrow(data)]
  #unclipped spawners
  sim.dat[1:2,3]=data$unknown_origin_spawners[c(nrow(data)-1):nrow(data)]
  #clipped spawners
  sim.dat[1:2,4]=data$clipped_spawners[c(nrow(data)-1):nrow(data)]
  
  #initial smolt counts - conditioned on observed data
  #unclipped smolts (unc. fry and natural origin):
  sim.dat[1:2,6]=data$smolts_cnt_unclipped[c(nrow(data)-1):nrow(data)]
  #clipped smolts (both clp. fry and smolt releases):
  sim.dat[1:2,7]=data$smolts_cnt_clipped_fry[c(nrow(data)-1):nrow(data)]+data$smolt_cnt_release_clipped[c(nrow(data)-1):nrow(data)]
  
  #hatchery releases
  #total fry
  sim.dat[1:2,12]=c(data$fry_release_captive_brood_unclpd[c(nrow(data)-1):nrow(data)]+data$fry_release_captive_brood_clpd[c(nrow(data)-1):nrow(data)])
  
  #smolts
  sim.dat[1:2,13]=data$smolt_release[c(nrow(data)-1):nrow(data)]
  
  #hatchery lake survival
  #fry survival to smolt emigration
  sim.dat[1:2,14]=data$fry_to_smolt_surv[c(nrow(data)-1):nrow(data)]/100
  #smolt survival to smolt emigration
  sim.dat[1:2,15]=data$smolt_release_emig_survival[c(nrow(data)-1):nrow(data)]/100
  
  ##Simulate - marine survival for clipped/unclipped smolts - spawners t+2
  #hatchery and unknown origin survival for these smolt cohorts
  
  #simulate for marine survival for each set of smolts
  sim.dat[1:2,9:10]=exp(MASS::mvrnorm(n=2,mu=mu.log.surv,Sigma=surv.mat))
  sim.dat[1:2,11]=exp(rnorm(n=2,mu.log.u.ms,sd.log.u.ms)) #unknown marine survival for first 2 unclipped smolt cohorts
  
  #estimate returning spawners at t+2 for each smolt cohort
  #unclipped spawners - unclipped outmigrating smolts x unknown origin marine survival
  sim.dat[3:4,3]=sim.dat[1:2,6]*(sim.dat[1:2,11]/100)
  #clipped spawners - clipped outmigrating smolts x hatchery origin marine survival
  sim.dat[3:4,4]=sim.dat[1:2,7]*(sim.dat[1:2,9]/100)
  #total spawners
  sim.dat[3:4,2]=sim.dat[3:4,3]+sim.dat[3:4,4]
  
  #set Ricker parameters for the iteration
  sim.ric=MASS::mvrnorm(n=1,mu=mu.ric,Sigma=vcov.ric)
  sim.ric[2:3]=exp(sim.ric[2:3]) #exponentiate from log-scale
  sim.ric[4]=inv_logistic_transf(sim.ric[4]) #function to scale transformed rho to -1 to 1
  
  #for all future states - next 25 years
  for(t in 1:c(sim.length+2)){
    
    sim.dat[t,2]=sum(sim.dat[t,3:5],na.rm=T) #total spawners
    
    #natural origin smolts from total spawners:
    
    #estimate natural origin smolts t+2 based on Ricker model
    if(t==1){
      #no autocorrelation if first smolt cohort (no previous estimate available to condition autocorrelation in smolt recruitment)
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2],sd=sim.ric[3])
    }else{
      #autocorrelation included
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2]+sim.ric[4]*sim.dat[t-1,18],sd=sim.ric[3])
    }
    #if total spawners < 2 = no natural smolt production
    if(sim.dat[t,2]<2){
      sim.dat[t+2,8]=0
    }else{
      sim.dat[t+2,8]=exp(logRS)*sim.dat[t,2] #natural-origin smolts outmigrating in t+2
    }
    sim.dat[t,16]=logRS
    sim.dat[t,17]=sim.ric[1]-sim.dat[t,2]/sim.ric[2]
    sim.dat[t,18]=logRS-(sim.ric[1]-sim.dat[t,2]/sim.ric[2])
    
    #fry & smolt releases - after first 2 burn-in years
    if(t>2){
      #fry releases - total fry releases x allocation to clipped/unclipped status
      total_fry.t=exp(rnorm(1,mu.log.fry2,sd.log.fry2)) #log-normal draws for fry releases
      total_fry.t=min(total_fry.t,2e6) #set capacity for fry release of 2 million
      sim.dat[t,12]=total_fry.t
      sim.dat[t,13]=exp(rnorm(1,mu.log.smolts2,sd.log.smolts2)) #log-normal draws for smolt releases
      
      #fry/smolt survival from lake to outmigration
      sim.dat[t,14]=plogis(rnorm(1,mu.logit.fry.surv,sd.logit.fry.surv)) #fry survival - applies to both clipped/unclipped
      sim.dat[t,15]=plogis(rnorm(1,mu.logit.smolt.surv,sd.logit.smolt.surv)) #smolt release survival to outmigration
      sim.dat[t,7]=sim.dat[t-1,12]*sim.dat[t,14]+sim.dat[t-1,13]*sim.dat[t,15] #fry releases x fry survival + smolt releases x smolt survival
      
      #expected returns (t+2) based on marine survival - natural origin & hatchery origin fish
      rho.ms.h.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.h,logistic.se.rho.ms.h)) #draw autocorrelation coefficient for hatchery origin smolts
      rho.ms.n.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.n,logistic.se.rho.ms.n)) #draw autocorrelation coefficient for natural origin smolts
      
      mu.adj=c(rho.ms.h.sim*(log(sim.dat[t-1,9])-log(surv.mean)[1]), #adjustment to mean survival based on 1-y autocorrelation
               rho.ms.n.sim*(log(sim.dat[t-1,10])-log(surv.mean)[2])) #adjustment to mean survival based on 1-y autocorrelation
      
      sim.dat[t,9:10]=exp(MASS::mvrnorm(n=1,mu=mu.log.surv+mu.adj,Sigma=surv.mat)) #marine survival
      
      #partition based on age at maturity:
      p_age=extraDistr::rdirichlet(1,c(3,87,10))
      #3% age 3, 87% age 4, 10% age 5
      #proportion recruits, add to each return year based on randomized maturity
      #age 3
      sim.dat[t+1,4]=sim.dat[t+1,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[1] #hatchery returns
      sim.dat[t+1,5]=sim.dat[t+1,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[1] #natural origin returns
      #age 4
      sim.dat[t+2,4]=sim.dat[t+2,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[2] #hatchery returns
      sim.dat[t+2,5]=sim.dat[t+2,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[2] #natural origin returns
      #age 5
      sim.dat[t+3,4]=sim.dat[t+3,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[3] #hatchery returns
      sim.dat[t+3,5]=sim.dat[t+3,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[3] #natural origin returns
    }
  }
  sim.dat=sim.dat[1:c(sim.length+2),] #trim extra brood years
  sim.summary4[i,2]=sum(sim.dat[5:c(sim.length+2),2]>100)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where total spawners > 500
  sim.summary4[i,3]=sum(sim.dat[5:c(sim.length+2),5]>100)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary4[i,4]=sum(sim.dat[5:c(sim.length+2),2]>500)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where total spawners > 500
  sim.summary4[i,5]=sum(sim.dat[5:c(sim.length+2),5]>500)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary4[i,6]=exp(mean(log(sim.dat[5:c(sim.length+2),2])))
  sim.summary4[i,7]=exp(mean(log(sim.dat[5:c(sim.length+2),5])))
  sim.summary4[i,8]=ifelse(all(sim.dat[c(sim.length+2-4):c(sim.length+2),2]==0),1,0)
  
  sims4[[i]]=sim.dat #store simulation
}


#Sim 5: No Hatchery####
sims5=list()
n.sims=1e3
sim.length=24
sim.summary5=data.frame(sim.id=seq(1:n.sims),py.tspwn.100=NA,py.nspwn.100=NA,py.tspwn.500=NA,py.nspwn.500=NA,m.tspwn=NA,m.nspwn=NA,extinct=NA)
for(i in 1:n.sims){
  sim.dat=data.frame(sim.year=seq(1:c(sim.length+2)),spawners.t=0,spawners.u=0,spawners.h=0,spawners.n=c(rep(NA,3),rep(0,sim.length-1)),smolts.u=NA,smolts.h=NA,smolts.n=NA,ms.h=NA,ms.n=NA,ms.u=NA,fry.rel=NA,smolt.rel=NA,fry.surv=NA,smolt.surv=NA,ricker.logrs=NA,ricker.mu=NA,ricker.eps=NA)
  #first 2 years - estimate returns based on known quantities (of e.g. fry and smolt releases, smolt outmigrants, etc) in the last several years
  
  #initial spawners - conditioned on observed data
  #total spawners
  sim.dat[1:2,2]=data$total_spawners[c(nrow(data)-1):nrow(data)]
  #unclipped spawners
  sim.dat[1:2,3]=data$unknown_origin_spawners[c(nrow(data)-1):nrow(data)]
  #clipped spawners
  sim.dat[1:2,4]=data$clipped_spawners[c(nrow(data)-1):nrow(data)]
  
  #initial smolt counts - conditioned on observed data
  #unclipped smolts (unc. fry and natural origin):
  sim.dat[1:2,6]=data$smolts_cnt_unclipped[c(nrow(data)-1):nrow(data)]
  #clipped smolts (both clp. fry and smolt releases):
  sim.dat[1:2,7]=data$smolts_cnt_clipped_fry[c(nrow(data)-1):nrow(data)]+data$smolt_cnt_release_clipped[c(nrow(data)-1):nrow(data)]
  
  #hatchery releases
  #total fry
  sim.dat[1:2,12]=c(data$fry_release_captive_brood_unclpd[c(nrow(data)-1):nrow(data)]+data$fry_release_captive_brood_clpd[c(nrow(data)-1):nrow(data)])
  
  #smolts
  sim.dat[1:2,13]=data$smolt_release[c(nrow(data)-1):nrow(data)]
  
  #hatchery lake survival
  #fry survival to smolt emigration
  sim.dat[1:2,14]=data$fry_to_smolt_surv[c(nrow(data)-1):nrow(data)]/100
  #smolt survival to smolt emigration
  sim.dat[1:2,15]=data$smolt_release_emig_survival[c(nrow(data)-1):nrow(data)]/100
  
  ##Simulate - marine survival for clipped/unclipped smolts - spawners t+2
  #hatchery and unknown origin survival for these smolt cohorts
  
  #simulate for marine survival for each set of smolts
  sim.dat[1:2,9:10]=exp(MASS::mvrnorm(n=2,mu=mu.log.surv,Sigma=surv.mat))
  sim.dat[1:2,11]=exp(rnorm(n=2,mu.log.u.ms,sd.log.u.ms)) #unknown marine survival for first 2 unclipped smolt cohorts
  
  #estimate returning spawners at t+2 for each smolt cohort
  #unclipped spawners - unclipped outmigrating smolts x unknown origin marine survival
  sim.dat[3:4,3]=sim.dat[1:2,6]*(sim.dat[1:2,11]/100)
  #clipped spawners - clipped outmigrating smolts x hatchery origin marine survival
  sim.dat[3:4,4]=sim.dat[1:2,7]*(sim.dat[1:2,9]/100)
  #total spawners
  sim.dat[3:4,2]=sim.dat[3:4,3]+sim.dat[3:4,4]
  
  #set Ricker parameters for the iteration
  sim.ric=MASS::mvrnorm(n=1,mu=mu.ric,Sigma=vcov.ric)
  sim.ric[2:3]=exp(sim.ric[2:3]) #exponentiate from log-scale
  sim.ric[4]=inv_logistic_transf(sim.ric[4]) #function to scale transformed rho to -1 to 1
  
  #for all future states - next 25 years
  for(t in 1:c(sim.length+2)){
    
    sim.dat[t,2]=sum(sim.dat[t,3:5],na.rm=T) #total spawners
    
    #natural origin smolts from total spawners:
    
    #estimate natural origin smolts t+2 based on Ricker model
    if(t==1){
      #no autocorrelation if first smolt cohort (no previous estimate available to condition autocorrelation in smolt recruitment)
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2],sd=sim.ric[3])
    }else{
      #autocorrelation included
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2]+sim.ric[4]*sim.dat[t-1,18],sd=sim.ric[3])
    }
    #if total spawners < 2 = no natural smolt production
    if(sim.dat[t,2]<2){
      sim.dat[t+2,8]=0
    }else{
      sim.dat[t+2,8]=exp(logRS)*sim.dat[t,2] #natural-origin smolts outmigrating in t+2
    }
    sim.dat[t,16]=logRS
    sim.dat[t,17]=sim.ric[1]-sim.dat[t,2]/sim.ric[2]
    sim.dat[t,18]=logRS-(sim.ric[1]-sim.dat[t,2]/sim.ric[2])
    
    #fry & smolt releases - after first 2 burn-in years
    if(t>2){
      sim.dat[t,7]=0 #no hatchery smolts leaving sakinaw lake
      
      #expected returns (t+2) based on marine survival - natural origin & hatchery origin fish
      rho.ms.h.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.h,logistic.se.rho.ms.h)) #draw autocorrelation coefficient for hatchery origin smolts
      rho.ms.n.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.n,logistic.se.rho.ms.n)) #draw autocorrelation coefficient for natural origin smolts
      
      mu.adj=c(rho.ms.h.sim*(log(sim.dat[t-1,9])-log(surv.mean)[1]), #adjustment to mean survival based on 1-y autocorrelation
               rho.ms.n.sim*(log(sim.dat[t-1,10])-log(surv.mean)[2])) #adjustment to mean survival based on 1-y autocorrelation
      
      sim.dat[t,9:10]=exp(MASS::mvrnorm(n=1,mu=mu.log.surv+mu.adj,Sigma=surv.mat)) #marine survival
      
      #partition based on age at maturity:
      p_age=extraDistr::rdirichlet(1,c(3,87,10))
      #3% age 3, 87% age 4, 10% age 5
      #proportion recruits, add to each return year based on randomized maturity
      #age 3
      sim.dat[t+1,4]=sim.dat[t+1,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[1] #hatchery returns
      sim.dat[t+1,5]=sim.dat[t+1,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[1] #natural origin returns
      #age 4
      sim.dat[t+2,4]=sim.dat[t+2,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[2] #hatchery returns
      sim.dat[t+2,5]=sim.dat[t+2,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[2] #natural origin returns
      #age 5
      sim.dat[t+3,4]=sim.dat[t+3,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[3] #hatchery returns
      sim.dat[t+3,5]=sim.dat[t+3,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[3] #natural origin returns
    }
  }
  sim.dat=sim.dat[1:c(sim.length+2),] #trim extra brood years
  sim.summary5[i,2]=sum(sim.dat[5:c(sim.length+2),2]>100)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where total spawners > 500
  sim.summary5[i,3]=sum(sim.dat[5:c(sim.length+2),5]>100)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary5[i,4]=sum(sim.dat[5:c(sim.length+2),2]>500)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where total spawners > 500
  sim.summary5[i,5]=sum(sim.dat[5:c(sim.length+2),5]>500)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary5[i,6]=exp(mean(log(sim.dat[5:c(sim.length+2),2])))
  sim.summary5[i,7]=exp(mean(log(sim.dat[5:c(sim.length+2),5])))
  sim.summary5[i,8]=ifelse(all(sim.dat[c(sim.length+2-4):c(sim.length+2),2]==0),1,0)
  
  sims5[[i]]=sim.dat #store simulation
}

#Sim 6: Brood-take####
sims6=list()
n.sims=1e3
sim.length=24
sim.summary6=data.frame(sim.id=seq(1:n.sims),py.tspwn.100=NA,py.nspwn.100=NA,py.tspwn.500=NA,py.nspwn.500=NA,m.tspwn=NA,m.nspwn=NA,extinct=NA)
for(i in 1:n.sims){
  sim.dat=data.frame(sim.year=seq(1:c(sim.length+2)),spawners.t=0,spawners.u=0,spawners.h=0,spawners.n=c(rep(NA,3),rep(0,sim.length-1)),smolts.u=NA,smolts.h=NA,smolts.n=NA,ms.h=NA,ms.n=NA,ms.u=NA,fry.rel=NA,smolt.rel=NA,fry.surv=NA,smolt.surv=NA,ricker.logrs=NA,ricker.mu=NA,ricker.eps=NA,broodtake=NA,spawners.r=NA,total.eggs=NA,egg.surv=NA)
  #first 2 years - estimate returns based on known quantities (of e.g. fry and smolt releases, smolt outmigrants, etc) in the last several years
  
  #initial spawners - conditioned on observed data
  #total spawners
  sim.dat[1:2,2]=data$total_spawners[c(nrow(data)-1):nrow(data)]
  #unclipped spawners
  sim.dat[1:2,3]=data$unknown_origin_spawners[c(nrow(data)-1):nrow(data)]
  #clipped spawners
  sim.dat[1:2,4]=data$clipped_spawners[c(nrow(data)-1):nrow(data)]
  
  #initial smolt counts - conditioned on observed data
  #unclipped smolts (unc. fry and natural origin):
  sim.dat[1:2,6]=data$smolts_cnt_unclipped[c(nrow(data)-1):nrow(data)]
  #clipped smolts (both clp. fry and smolt releases):
  sim.dat[1:2,7]=data$smolts_cnt_clipped_fry[c(nrow(data)-1):nrow(data)]+data$smolt_cnt_release_clipped[c(nrow(data)-1):nrow(data)]
  
  #hatchery releases
  #total fry
  sim.dat[1:2,12]=c(data$fry_release_captive_brood_unclpd[c(nrow(data)-1):nrow(data)]+data$fry_release_captive_brood_clpd[c(nrow(data)-1):nrow(data)])
  
  #smolts
  sim.dat[1:2,13]=data$smolt_release[c(nrow(data)-1):nrow(data)]
  
  #hatchery lake survival
  #fry survival to smolt emigration
  sim.dat[1:2,14]=data$fry_to_smolt_surv[c(nrow(data)-1):nrow(data)]/100
  #smolt survival to smolt emigration
  sim.dat[1:2,15]=data$smolt_release_emig_survival[c(nrow(data)-1):nrow(data)]/100
  
  ##Simulate - marine survival for clipped/unclipped smolts - spawners t+2
  #hatchery and unknown origin survival for these smolt cohorts
  
  #simulate for marine survival for each set of smolts
  sim.dat[1:2,9:10]=exp(MASS::mvrnorm(n=2,mu=mu.log.surv,Sigma=surv.mat))
  sim.dat[1:2,11]=exp(rnorm(n=2,mu.log.u.ms,sd.log.u.ms)) #unknown marine survival for first 2 unclipped smolt cohorts
  
  #estimate returning spawners at t+2 for each smolt cohort
  #unclipped spawners - unclipped outmigrating smolts x unknown origin marine survival
  sim.dat[3:4,3]=sim.dat[1:2,6]*(sim.dat[1:2,11]/100)
  #clipped spawners - clipped outmigrating smolts x hatchery origin marine survival
  sim.dat[3:4,4]=sim.dat[1:2,7]*(sim.dat[1:2,9]/100)
  #total spawners
  sim.dat[3:4,2]=sim.dat[3:4,3]+sim.dat[3:4,4]
  
  #set Ricker parameters for the iteration
  sim.ric=MASS::mvrnorm(n=1,mu=mu.ric,Sigma=vcov.ric)
  sim.ric[2:3]=exp(sim.ric[2:3]) #exponentiate from log-scale
  sim.ric[4]=inv_logistic_transf(sim.ric[4]) #function to scale transformed rho to -1 to 1
  
  #for all future states - next 25 years
  for(t in 1:c(sim.length+2)){
    
    sim.dat[t,2]=sum(sim.dat[t,3:5],na.rm=T) #total spawners
    #brood-take branch starting at t=3
    #remove up to half of total spawners for hatchery broodtake, up to maximum capacity
    if(t>2){
      sim.dat[t,20]=sim.dat[t,2] #save original spawners at weir before broodtake
      sim.dat[t,19]= min(0.5*sim.dat[t,2], 750) #broodtake - minimum of half of total spawners or max take 
      sim.dat[t,2]=sim.dat[t,20]-sim.dat[t,19] #remaining spawners for natural production
      if(sim.dat[t,2]<2){sim.dat[t,2]=0}
    }
    #natural origin smolts from total spawners:
    #estimate natural origin smolts t+2 based on Ricker model
    if(t==1){
      #no autocorrelation if first smolt cohort (no previous estimate available to condition autocorrelation in smolt recruitment)
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2],sd=sim.ric[3])
    }else{
      #autocorrelation included
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2]+sim.ric[4]*sim.dat[t-1,18],sd=sim.ric[3])
    }
    #if total spawners < 2 = no natural smolt production
    if(sim.dat[t,2]<2){
      sim.dat[t+2,8]=0
    }else{
      sim.dat[t+2,8]=exp(logRS)*sim.dat[t,2] #natural-origin smolts outmigrating in t+2
    }
    sim.dat[t,16]=logRS
    sim.dat[t,17]=sim.ric[1]-sim.dat[t,2]/sim.ric[2]
    sim.dat[t,18]=logRS-(sim.ric[1]-sim.dat[t,2]/sim.ric[2])
    
    #fry & smolt releases - after first 2 burn-in years
    if(t>2){
      
      sim.dat[t,21]=sim.dat[t,19]*0.5*rnorm(1,2050,250) #total eggs
      sim.dat[t,22]=plogis(rnorm(1,qlogis(0.86),0.2)) #egg-to-fry survival (mean = 86%, sd is guess-timated)
      
      #fry releases - total fry releases x allocation to clipped/unclipped status
      total_fry.t=sim.dat[t,21]*sim.dat[t,22] #total-fry release is now total eggs from broodtake x egg survival in previous year
      total_fry.t=min(total_fry.t,2e6) #set capacity for fry release of 2 million
      p.sm=plogis(rnorm(1,mu.logit.psm,0.2)) #proportion of total retained to raise as smolts
       #proportion retained as smolts and adjusted for some fry-to-smolt mortality
      sim.dat[t+1,12]=total_fry.t*(1-p.sm*c(1/plogis(rnorm(1,qlogis(0.8),0.2)))) #unclipped fry releases - #proportion retained as smolts and adjusted for some fry-to-smolt mortality
      sim.dat[t+1,13]=total_fry.t*p.sm #smolt releases
      
      #fry/smolt survival from lake to outmigration
      sim.dat[t,14]=plogis(rnorm(1,mu.logit.fry.surv,sd.logit.fry.surv)) #fry survival - applies to both clipped/unclipped
      sim.dat[t,15]=plogis(rnorm(1,mu.logit.smolt.surv,sd.logit.smolt.surv)) #smolt release survival to outmigration
      sim.dat[t,7]=sim.dat[t-1,12]*sim.dat[t,14]+sim.dat[t-1,13]*sim.dat[t,15] #fry releases x fry survival + smolt releases x smolt survival
      
      #expected returns (t+2) based on marine survival - natural origin & hatchery origin fish
      rho.ms.h.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.h,logistic.se.rho.ms.h)) #draw autocorrelation coefficient for hatchery origin smolts
      rho.ms.n.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.n,logistic.se.rho.ms.n)) #draw autocorrelation coefficient for natural origin smolts
      
      mu.adj=c(rho.ms.h.sim*(log(sim.dat[t-1,9])-log(surv.mean)[1]), #adjustment to mean survival based on 1-y autocorrelation
               rho.ms.n.sim*(log(sim.dat[t-1,10])-log(surv.mean)[2])) #adjustment to mean survival based on 1-y autocorrelation
      
      sim.dat[t,9:10]=exp(MASS::mvrnorm(n=1,mu=mu.log.surv+mu.adj,Sigma=surv.mat)) #marine survival
      
      #partition based on age at maturity:
      p_age=extraDistr::rdirichlet(1,c(3,87,10))
      #3% age 3, 87% age 4, 10% age 5
      #proportion recruits, add to each return year based on randomized maturity
      #age 3
      sim.dat[t+1,4]=sim.dat[t+1,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[1] #hatchery returns
      sim.dat[t+1,5]=sim.dat[t+1,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[1] #natural origin returns
      #age 4
      sim.dat[t+2,4]=sim.dat[t+2,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[2] #hatchery returns
      sim.dat[t+2,5]=sim.dat[t+2,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[2] #natural origin returns
      #age 5
      sim.dat[t+3,4]=sim.dat[t+3,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[3] #hatchery returns
      sim.dat[t+3,5]=sim.dat[t+3,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[3] #natural origin returns
    }
  }
  sim.dat=sim.dat[1:c(sim.length+2),] #trim extra brood years
  sim.summary6[i,2]=sum(sim.dat[5:c(sim.length+2),2]>100)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where total spawners > 500
  sim.summary6[i,3]=sum(sim.dat[5:c(sim.length+2),5]>100)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary6[i,4]=sum(sim.dat[5:c(sim.length+2),2]>500)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where total spawners > 500
  sim.summary6[i,5]=sum(sim.dat[5:c(sim.length+2),5]>500)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary6[i,6]=exp(mean(log(sim.dat[5:c(sim.length+2),2])))
  sim.summary6[i,7]=exp(mean(log(sim.dat[5:c(sim.length+2),5])))
  sim.summary6[i,8]=ifelse(all(sim.dat[c(sim.length+2-4):c(sim.length+2),2]==0),1,0)
  
  sims6[[i]]=sim.dat #store simulation
}


#Sim 7: Big Qualicum Brood-take - 60k releases####
bigQms=data$bigQsurvival[complete.cases(data$bigQsurvival)]
mu.bq.ms=mean(bigQms)
sd.bq.ms=sd(bigQms)
sd.log.bq.ms=sqrt(log(1+(sd.bq.ms/mean(bigQms))^2))
mu.log.bq.ms=log(mu.bq.ms)-0.5*sd.log.bq.ms^2
bq.release.tar=60000

sims7=list()
n.sims=1e3
sim.length=24
sim.summary7=data.frame(sim.id=seq(1:n.sims),py.tspwn.100=NA,py.nspwn.100=NA,py.tspwn.500=NA,py.nspwn.500=NA,m.tspwn=NA,m.nspwn=NA,extinct=NA)
for(i in 1:n.sims){
  sim.dat=data.frame(sim.year=seq(1:c(sim.length+2)),spawners.t=0,spawners.u=0,spawners.h=0,spawners.n=c(rep(NA,3),rep(0,sim.length-1)),smolts.u=NA,smolts.h=NA,smolts.n=NA,ms.h=NA,ms.n=NA,ms.u=NA,fry.rel=NA,smolt.rel=NA,fry.surv=NA,smolt.surv=NA,ricker.logrs=NA,ricker.mu=NA,ricker.eps=NA,bigQRelease=NA,bigQms=NA,bigQreturn=NA,total.eggs=NA,egg.surv=NA,total.fry=NA,bqFry=NA,sakFry=NA)
  #first 2 years - estimate returns based on known quantities (of e.g. fry and smolt releases, smolt outmigrants, etc) in the last several years
  
  #initial spawners - conditioned on observed data
  #total spawners
  sim.dat[1:2,2]=data$total_spawners[c(nrow(data)-1):nrow(data)]
  #unclipped spawners
  sim.dat[1:2,3]=data$unknown_origin_spawners[c(nrow(data)-1):nrow(data)]
  #clipped spawners
  sim.dat[1:2,4]=data$clipped_spawners[c(nrow(data)-1):nrow(data)]
  
  #initial smolt counts - conditioned on observed data
  #unclipped smolts (unc. fry and natural origin):
  sim.dat[1:2,6]=data$smolts_cnt_unclipped[c(nrow(data)-1):nrow(data)]
  #clipped smolts (both clp. fry and smolt releases):
  sim.dat[1:2,7]=data$smolts_cnt_clipped_fry[c(nrow(data)-1):nrow(data)]+data$smolt_cnt_release_clipped[c(nrow(data)-1):nrow(data)]
  
  #hatchery releases
  #total fry
  sim.dat[1:2,12]=c(data$fry_release_captive_brood_unclpd[c(nrow(data)-1):nrow(data)]+data$fry_release_captive_brood_clpd[c(nrow(data)-1):nrow(data)])
  
  #smolts
  sim.dat[1:2,13]=data$smolt_release[c(nrow(data)-1):nrow(data)]
  
  #hatchery lake survival
  #fry survival to smolt emigration
  sim.dat[1:2,14]=data$fry_to_smolt_surv[c(nrow(data)-1):nrow(data)]/100
  #smolt survival to smolt emigration
  sim.dat[1:2,15]=data$smolt_release_emig_survival[c(nrow(data)-1):nrow(data)]/100
  
  ##Simulate - marine survival for clipped/unclipped smolts - spawners t+2
  #hatchery and unknown origin survival for these smolt cohorts
  
  #simulate for marine survival for each set of smolts
  sim.dat[1:2,9:10]=exp(MASS::mvrnorm(n=2,mu=mu.log.surv,Sigma=surv.mat))
  sim.dat[1:2,11]=exp(rnorm(n=2,mu.log.u.ms,sd.log.u.ms)) #unknown marine survival for first 2 unclipped smolt cohorts
  
  #estimate returning spawners at t+2 for each smolt cohort
  #unclipped spawners - unclipped outmigrating smolts x unknown origin marine survival
  sim.dat[3:4,3]=sim.dat[1:2,6]*(sim.dat[1:2,11]/100)
  #clipped spawners - clipped outmigrating smolts x hatchery origin marine survival
  sim.dat[3:4,4]=sim.dat[1:2,7]*(sim.dat[1:2,9]/100)
  #total spawners
  sim.dat[3:4,2]=sim.dat[3:4,3]+sim.dat[3:4,4]
  
  #set Ricker parameters for the iteration
  sim.ric=MASS::mvrnorm(n=1,mu=mu.ric,Sigma=vcov.ric)
  sim.ric[2:3]=exp(sim.ric[2:3]) #exponentiate from log-scale
  sim.ric[4]=inv_logistic_transf(sim.ric[4]) #function to scale transformed rho to -1 to 1
  
  #bigQ releases 
  sim.dat[1:3,19]=rnorm(3,bq.release.tar,1e3) #hypothetical 50k smolt releases at big q
  sim.dat[1:3,20]=exp(rnorm(3,mu.log.bq.ms,sd.log.bq.ms)) #hypothetical 50k smolt releases at big q
  
  #add in fry release under old regime for burn-in years
  sim.dat[3:4,12]=exp(rnorm(2,mu.log.fry,sd.log.fry)) #initiate 2 years of captive brood fry cohorts
  sim.dat[3:4,13]=exp(rnorm(2,mu.log.smolts,sd.log.smolts)) #initiate 2 years of captive brood smolt cohorts
  
  
  #for all future states - next 25 years
  for(t in 1:c(sim.length+2)){
    
    
    sim.dat[t,2]=sum(sim.dat[t,3:5],na.rm=T) #total spawners
    
    #natural origin smolts from total spawners:
    #estimate natural origin smolts t+2 based on Ricker model
    if(t==1){
      #no autocorrelation if first smolt cohort (no previous estimate available to condition autocorrelation in smolt recruitment)
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2],sd=sim.ric[3])
    }else{
      #autocorrelation included
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2]+sim.ric[4]*sim.dat[t-1,18],sd=sim.ric[3])
    }
    #if total spawners < 2 = no natural smolt production
    if(sim.dat[t,2]<2){
      sim.dat[t+2,8]=0
    }else{
      sim.dat[t+2,8]=exp(logRS)*sim.dat[t,2] #natural-origin smolts outmigrating in t+2
    }
    sim.dat[t,16]=logRS
    sim.dat[t,17]=sim.ric[1]-sim.dat[t,2]/sim.ric[2]
    sim.dat[t,18]=logRS-(sim.ric[1]-sim.dat[t,2]/sim.ric[2])
    
    
    #brood-take branch starting at t=3
    if(t>3){
      #big Q loop
      
      #calculate returns from smolt release 3-years ago
      sim.dat[t,21]=sim.dat[t-3,19]*(sim.dat[t-3,20]/100)
      sim.dat[t,22]=sim.dat[t,21]*0.5*rnorm(1,2050,250) #total eggs from all big q returns
      sim.dat[t,23]=plogis(rnorm(1,qlogis(0.86),0.2)) #egg-to-fry survival
      
      sim.dat[t,24]=sim.dat[t,22]*sim.dat[t,23] #total fry produced from big Q broodtake
      
      fry_sm_surv=plogis(rnorm(1,qlogis(0.8),0.2))
      #big Q releases & marine survival - these are functionally independent draws
      sm.targ=rnorm(1,bq.release.tar,0.05*bq.release.tar) #annual smolt release targe
      
      if(c(sm.targ*(1/fry_sm_surv))<sim.dat[t,24]){  
        sim.dat[t,19]=sm.targ
        sim.dat[t,25]=sim.dat[t,19]*c(1/fry_sm_surv) #fry kept for Big Q smolt releases - adjust based on assumed survival of fry to smolts of 80%
        sim.dat[t,26]=sim.dat[t,24]-sim.dat[t,25] #total - bigQ fry
      }else{
        sim.dat[t,19]=sim.dat[t,24]*c(plogis(rnorm(1,qlogis(0.8),0.2))) #realized smolt releases - adjusted for fry-to-smolt survival
        sim.dat[t,25]=sim.dat[t,24] #fry kept for Big Q smolt releases
        sim.dat[t,26]=0
      }
      #remaining fry released to sakinaw in next year and a percent for smolt releases
      p.sm=plogis(rnorm(1,mu.logit.psm,0.2)) #proportion of total retained to raise as smolts
      #proportion retained as smolts and adjusted for some fry-to-smolt mortality
      sim.dat[t+1,12]=sim.dat[t,26]*(1-p.sm*c(1/plogis(rnorm(1,qlogis(0.8),0.2)))) #unclipped fry releases - #proportion retained as smolts and adjusted for some fry-to-smolt mortality
      sim.dat[t+1,13]=sim.dat[t,26]*p.sm #smolt releases
    }
    if(t>2){
      #marine survival for release cohort
      sim.dat[t,20]=exp(rnorm(1,mu.log.bq.ms,sd.log.bq.ms)) #survival for this release cohort
    
      #fry/smolt survival from lake to outmigration
      sim.dat[t,14]=plogis(rnorm(1,mu.logit.fry.surv,sd.logit.fry.surv)) #fry survival - applies to both clipped/unclipped
      sim.dat[t,15]=plogis(rnorm(1,mu.logit.smolt.surv,sd.logit.smolt.surv)) #smolt release survival to outmigration
      sim.dat[t,7]=sim.dat[t-1,12]*sim.dat[t,14]+sim.dat[t-1,13]*sim.dat[t,15] #fry releases x fry survival + smolt releases x smolt survival
      
      #expected returns (t+2) based on marine survival - natural origin & hatchery origin fish
      rho.ms.h.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.h,logistic.se.rho.ms.h)) #draw autocorrelation coefficient for hatchery origin smolts
      rho.ms.n.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.n,logistic.se.rho.ms.n)) #draw autocorrelation coefficient for natural origin smolts
      
      mu.adj=c(rho.ms.h.sim*(log(sim.dat[t-1,9])-log(surv.mean)[1]), #adjustment to mean survival based on 1-y autocorrelation
               rho.ms.n.sim*(log(sim.dat[t-1,10])-log(surv.mean)[2])) #adjustment to mean survival based on 1-y autocorrelation
      
      sim.dat[t,9:10]=exp(MASS::mvrnorm(n=1,mu=mu.log.surv+mu.adj,Sigma=surv.mat)) #marine survival
      
      #partition based on age at maturity:
      p_age=extraDistr::rdirichlet(1,c(3,87,10))
      #3% age 3, 87% age 4, 10% age 5
      #proportion recruits, add to each return year based on randomized maturity
      #age 3
      sim.dat[t+1,4]=sim.dat[t+1,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[1] #hatchery returns
      sim.dat[t+1,5]=sim.dat[t+1,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[1] #natural origin returns
      #age 4
      sim.dat[t+2,4]=sim.dat[t+2,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[2] #hatchery returns
      sim.dat[t+2,5]=sim.dat[t+2,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[2] #natural origin returns
      #age 5
      sim.dat[t+3,4]=sim.dat[t+3,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[3] #hatchery returns
      sim.dat[t+3,5]=sim.dat[t+3,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[3] #natural origin returns
      }
  }
  sim.dat=sim.dat[1:c(sim.length+2),] #trim extra brood years
  sim.summary7[i,2]=sum(sim.dat[5:c(sim.length+2),2]>100)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where total spawners > 500
  sim.summary7[i,3]=sum(sim.dat[5:c(sim.length+2),5]>100)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary7[i,4]=sum(sim.dat[5:c(sim.length+2),2]>500)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where total spawners > 500
  sim.summary7[i,5]=sum(sim.dat[5:c(sim.length+2),5]>500)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary7[i,6]=exp(mean(log(sim.dat[5:c(sim.length+2),2])))
  sim.summary7[i,7]=exp(mean(log(sim.dat[5:c(sim.length+2),5])))
  sim.summary7[i,8]=ifelse(all(sim.dat[c(sim.length+2-4):c(sim.length+2),2]==0),1,0)
  
  sims7[[i]]=sim.dat #store simulation
}


#Sim 8: Big Qualicum Brood-take - 120k releases####
bigQms=data$bigQsurvival[complete.cases(data$bigQsurvival)]
mu.bq.ms=mean(bigQms)
sd.bq.ms=sd(bigQms)
sd.log.bq.ms=sqrt(log(1+(sd.bq.ms/mean(bigQms))^2))
mu.log.bq.ms=log(mu.bq.ms)-0.5*sd.log.bq.ms^2
bq.release.tar=120000

sims8=list()
n.sims=1e3
sim.length=24
sim.summary8=data.frame(sim.id=seq(1:n.sims),py.tspwn.100=NA,py.nspwn.100=NA,py.tspwn.500=NA,py.nspwn.500=NA,m.tspwn=NA,m.nspwn=NA,extinct=NA)
for(i in 1:n.sims){
  sim.dat=data.frame(sim.year=seq(1:c(sim.length+2)),spawners.t=0,spawners.u=0,spawners.h=0,spawners.n=c(rep(NA,3),rep(0,sim.length-1)),smolts.u=NA,smolts.h=NA,smolts.n=NA,ms.h=NA,ms.n=NA,ms.u=NA,fry.rel=NA,smolt.rel=NA,fry.surv=NA,smolt.surv=NA,ricker.logrs=NA,ricker.mu=NA,ricker.eps=NA,bigQRelease=NA,bigQms=NA,bigQreturn=NA,total.eggs=NA,egg.surv=NA,total.fry=NA,bqFry=NA,sakFry=NA)
  #first 2 years - estimate returns based on known quantities (of e.g. fry and smolt releases, smolt outmigrants, etc) in the last several years
  
  #initial spawners - conditioned on observed data
  #total spawners
  sim.dat[1:2,2]=data$total_spawners[c(nrow(data)-1):nrow(data)]
  #unclipped spawners
  sim.dat[1:2,3]=data$unknown_origin_spawners[c(nrow(data)-1):nrow(data)]
  #clipped spawners
  sim.dat[1:2,4]=data$clipped_spawners[c(nrow(data)-1):nrow(data)]
  
  #initial smolt counts - conditioned on observed data
  #unclipped smolts (unc. fry and natural origin):
  sim.dat[1:2,6]=data$smolts_cnt_unclipped[c(nrow(data)-1):nrow(data)]
  #clipped smolts (both clp. fry and smolt releases):
  sim.dat[1:2,7]=data$smolts_cnt_clipped_fry[c(nrow(data)-1):nrow(data)]+data$smolt_cnt_release_clipped[c(nrow(data)-1):nrow(data)]
  
  #hatchery releases
  #total fry
  sim.dat[1:2,12]=c(data$fry_release_captive_brood_unclpd[c(nrow(data)-1):nrow(data)]+data$fry_release_captive_brood_clpd[c(nrow(data)-1):nrow(data)])
  
  #smolts
  sim.dat[1:2,13]=data$smolt_release[c(nrow(data)-1):nrow(data)]
  
  #hatchery lake survival
  #fry survival to smolt emigration
  sim.dat[1:2,14]=data$fry_to_smolt_surv[c(nrow(data)-1):nrow(data)]/100
  #smolt survival to smolt emigration
  sim.dat[1:2,15]=data$smolt_release_emig_survival[c(nrow(data)-1):nrow(data)]/100
  
  ##Simulate - marine survival for clipped/unclipped smolts - spawners t+2
  #hatchery and unknown origin survival for these smolt cohorts
  
  #simulate for marine survival for each set of smolts
  sim.dat[1:2,9:10]=exp(MASS::mvrnorm(n=2,mu=mu.log.surv,Sigma=surv.mat))
  sim.dat[1:2,11]=exp(rnorm(n=2,mu.log.u.ms,sd.log.u.ms)) #unknown marine survival for first 2 unclipped smolt cohorts
  
  #estimate returning spawners at t+2 for each smolt cohort
  #unclipped spawners - unclipped outmigrating smolts x unknown origin marine survival
  sim.dat[3:4,3]=sim.dat[1:2,6]*(sim.dat[1:2,11]/100)
  #clipped spawners - clipped outmigrating smolts x hatchery origin marine survival
  sim.dat[3:4,4]=sim.dat[1:2,7]*(sim.dat[1:2,9]/100)
  #total spawners
  sim.dat[3:4,2]=sim.dat[3:4,3]+sim.dat[3:4,4]
  
  #set Ricker parameters for the iteration
  sim.ric=MASS::mvrnorm(n=1,mu=mu.ric,Sigma=vcov.ric)
  sim.ric[2:3]=exp(sim.ric[2:3]) #exponentiate from log-scale
  sim.ric[4]=inv_logistic_transf(sim.ric[4]) #function to scale transformed rho to -1 to 1
  
  #bigQ releases 
  sim.dat[1:3,19]=rnorm(3,bq.release.tar,1e3) #hypothetical 50k smolt releases at big q
  sim.dat[1:3,20]=exp(rnorm(3,mu.log.bq.ms,sd.log.bq.ms)) #hypothetical 50k smolt releases at big q
  
  #add in fry release under old regime for burn-in years
  sim.dat[3:4,12]=exp(rnorm(2,mu.log.fry,sd.log.fry)) #initiate 2 years of captive brood fry cohorts
  sim.dat[3:4,13]=exp(rnorm(2,mu.log.smolts,sd.log.smolts)) #initiate 2 years of captive brood smolt cohorts
  
  
  #for all future states - next 25 years
  for(t in 1:c(sim.length+2)){
    
    
    sim.dat[t,2]=sum(sim.dat[t,3:5],na.rm=T) #total spawners
    
    #natural origin smolts from total spawners:
    #estimate natural origin smolts t+2 based on Ricker model
    if(t==1){
      #no autocorrelation if first smolt cohort (no previous estimate available to condition autocorrelation in smolt recruitment)
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2],sd=sim.ric[3])
    }else{
      #autocorrelation included
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2]+sim.ric[4]*sim.dat[t-1,18],sd=sim.ric[3])
    }
    #if total spawners < 2 = no natural smolt production
    if(sim.dat[t,2]<2){
      sim.dat[t+2,8]=0
    }else{
      sim.dat[t+2,8]=exp(logRS)*sim.dat[t,2] #natural-origin smolts outmigrating in t+2
    }
    sim.dat[t,16]=logRS
    sim.dat[t,17]=sim.ric[1]-sim.dat[t,2]/sim.ric[2]
    sim.dat[t,18]=logRS-(sim.ric[1]-sim.dat[t,2]/sim.ric[2])
    
    
    #brood-take branch starting at t=3
    if(t>3){
      #big Q loop
      
      #calculate returns from smolt release 3-years ago
      sim.dat[t,21]=sim.dat[t-3,19]*(sim.dat[t-3,20]/100)
      sim.dat[t,22]=sim.dat[t,21]*0.5*rnorm(1,2050,250) #total eggs from all big q returns
      sim.dat[t,23]=plogis(rnorm(1,qlogis(0.86),0.2)) #egg-to-fry survival
      
      sim.dat[t,24]=sim.dat[t,22]*sim.dat[t,23] #total fry produced from big Q broodtake
      
      fry_sm_surv=plogis(rnorm(1,qlogis(0.8),0.2))
      #big Q releases & marine survival - these are functionally independent draws
      sm.targ=rnorm(1,bq.release.tar,0.05*bq.release.tar) #annual smolt release targe
      
      if(c(sm.targ*(1/fry_sm_surv))<sim.dat[t,24]){  
        sim.dat[t,19]=sm.targ
        sim.dat[t,25]=sim.dat[t,19]*c(1/fry_sm_surv) #fry kept for Big Q smolt releases - adjust based on assumed survival of fry to smolts of 80%
        sim.dat[t,26]=sim.dat[t,24]-sim.dat[t,25] #total - bigQ fry
      }else{
        sim.dat[t,19]=sim.dat[t,24]*c(plogis(rnorm(1,qlogis(0.8),0.2))) #realized smolt releases - adjusted for fry-to-smolt survival
        sim.dat[t,25]=sim.dat[t,24] #fry kept for Big Q smolt releases
        sim.dat[t,26]=0
      }
      #remaining fry released to sakinaw in next year and a percent for smolt releases
      p.sm=plogis(rnorm(1,mu.logit.psm,0.2)) #proportion of total retained to raise as smolts
      #proportion retained as smolts and adjusted for some fry-to-smolt mortality
      sim.dat[t+1,12]=sim.dat[t,26]*(1-p.sm*c(1/plogis(rnorm(1,qlogis(0.8),0.2)))) #unclipped fry releases - #proportion retained as smolts and adjusted for some fry-to-smolt mortality
      sim.dat[t+1,13]=sim.dat[t,26]*p.sm #smolt releases
    }
    if(t>2){
      #marine survival for release cohort
      sim.dat[t,20]=exp(rnorm(1,mu.log.bq.ms,sd.log.bq.ms)) #survival for this release cohort
      
      #fry/smolt survival from lake to outmigration
      sim.dat[t,14]=plogis(rnorm(1,mu.logit.fry.surv,sd.logit.fry.surv)) #fry survival - applies to both clipped/unclipped
      sim.dat[t,15]=plogis(rnorm(1,mu.logit.smolt.surv,sd.logit.smolt.surv)) #smolt release survival to outmigration
      sim.dat[t,7]=sim.dat[t-1,12]*sim.dat[t,14]+sim.dat[t-1,13]*sim.dat[t,15] #fry releases x fry survival + smolt releases x smolt survival
      
      #expected returns (t+2) based on marine survival - natural origin & hatchery origin fish
      rho.ms.h.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.h,logistic.se.rho.ms.h)) #draw autocorrelation coefficient for hatchery origin smolts
      rho.ms.n.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.n,logistic.se.rho.ms.n)) #draw autocorrelation coefficient for natural origin smolts
      
      mu.adj=c(rho.ms.h.sim*(log(sim.dat[t-1,9])-log(surv.mean)[1]), #adjustment to mean survival based on 1-y autocorrelation
               rho.ms.n.sim*(log(sim.dat[t-1,10])-log(surv.mean)[2])) #adjustment to mean survival based on 1-y autocorrelation
      
      sim.dat[t,9:10]=exp(MASS::mvrnorm(n=1,mu=mu.log.surv+mu.adj,Sigma=surv.mat)) #marine survival
      
      #partition based on age at maturity:
      p_age=extraDistr::rdirichlet(1,c(3,87,10))
      #3% age 3, 87% age 4, 10% age 5
      #proportion recruits, add to each return year based on randomized maturity
      #age 3
      sim.dat[t+1,4]=sim.dat[t+1,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[1] #hatchery returns
      sim.dat[t+1,5]=sim.dat[t+1,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[1] #natural origin returns
      #age 4
      sim.dat[t+2,4]=sim.dat[t+2,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[2] #hatchery returns
      sim.dat[t+2,5]=sim.dat[t+2,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[2] #natural origin returns
      #age 5
      sim.dat[t+3,4]=sim.dat[t+3,4]+sim.dat[t,7]*(sim.dat[t,9]/100)*p_age[3] #hatchery returns
      sim.dat[t+3,5]=sim.dat[t+3,5]+sim.dat[t,8]*(sim.dat[t,10]/100)*p_age[3] #natural origin returns
    }
  }
  sim.dat=sim.dat[1:c(sim.length+2),] #trim extra brood years
  sim.summary8[i,2]=sum(sim.dat[5:c(sim.length+2),2]>100)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where total spawners > 500
  sim.summary8[i,3]=sum(sim.dat[5:c(sim.length+2),5]>100)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary8[i,4]=sum(sim.dat[5:c(sim.length+2),2]>500)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where total spawners > 500
  sim.summary8[i,5]=sum(sim.dat[5:c(sim.length+2),5]>500)/length(5:c(sim.length+2)) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary8[i,6]=exp(mean(log(sim.dat[5:c(sim.length+2),2])))
  sim.summary8[i,7]=exp(mean(log(sim.dat[5:c(sim.length+2),5])))
  sim.summary8[i,8]=ifelse(all(sim.dat[c(sim.length+2-4):c(sim.length+2),2]==0),1,0)
  
  sims8[[i]]=sim.dat #store simulation
}


#Table 1####
summary(sim.summary1) #baseline
summary(sim.summary2) #double hatchery
summary(sim.summary3) #half hatchery
summary(sim.summary4) #quarter hatchery
summary(sim.summary5) #all natural
summary(sim.summary6) #broodtake
summary(sim.summary7) #Big Q broodtake - 100k releases
summary(sim.summary8) #Big Q broodtake - 100k releases

#Fig 2####
jpeg('./outputs/Figure_2.jpeg',height=11,width=8.5,res=500,units='in',pointsize=11)
par(mfrow=c(4,2))
hist(sim.summary1$py.tspwn.100,main = 'S1: Status quo',xlab='Prop. of years > 100 spawners',xlim=c(0,1),breaks=20)
hist(sim.summary2$py.tspwn.100,main = 'S2: Doubled hatchery releases',xlab='Prop. of years > 100 spawners',xlim=c(0,1))
hist(sim.summary3$py.tspwn.100,main = 'S3: Halved hatchery releases',xlab='Prop. of years > 100 spawners',xlim=c(0,1),breaks=20)
hist(sim.summary4$py.tspwn.100,main = 'S4: Quartered hatchery releases',xlab='Prop. of years > 100 spawners',xlim=c(0,1),breaks=20)
hist(sim.summary5$py.tspwn.100,main = 'S5: No hatchery',xlab='Prop. of years > 100 spawners',xlim=c(0,1))
hist(sim.summary6$py.tspwn.100,main = 'S6: In situ broodstock',xlab='Prop. of years > 100 spawners',xlim=c(0,1))
hist(sim.summary7$py.tspwn.100,main = 'S7: Ex situ broodstock (60k smolt target)',xlab='Prop. of years > 100 spawners',xlim=c(0,1),breaks=20)
hist(sim.summary8$py.tspwn.100,main = 'S8: Ex situ broodstock (120k smolt target)',xlab='Prop. of years > 100 spawners',xlim=c(0,1),breaks=20)
dev.off()


#Fig 3####

sims1.df=dplyr::bind_rows(sims1)
sims2.df=dplyr::bind_rows(sims2)
sims3.df=dplyr::bind_rows(sims3)
sims4.df=dplyr::bind_rows(sims4)
sims5.df=dplyr::bind_rows(sims5)
sims6.df=dplyr::bind_rows(sims6)
sims7.df=dplyr::bind_rows(sims7)
sims8.df=dplyr::bind_rows(sims8)


library(dplyr)
sims1.avg=sims1.df %>% group_by(sim.year) %>% summarize(m.spawners=median(spawners.t),l90.spawners=quantile(spawners.t,0.05),u90.spawners=quantile(spawners.t,0.95),l80.spawners=quantile(spawners.t,0.1),u80.spawners=quantile(spawners.t,0.9),l50.spawners=quantile(spawners.t,0.25),u50.spawners=quantile(spawners.t,0.75),main='Status quo')
sims1.avg$brood.year=seq(2024,c(2024+nrow(sims1.avg)-1))
sims2.avg=sims2.df %>% group_by(sim.year) %>% summarize(m.spawners=median(spawners.t),l90.spawners=quantile(spawners.t,0.05),u90.spawners=quantile(spawners.t,0.95),l80.spawners=quantile(spawners.t,0.1),u80.spawners=quantile(spawners.t,0.9),l50.spawners=quantile(spawners.t,0.25),u50.spawners=quantile(spawners.t,0.75),main='Status quo')
sims2.avg$brood.year=seq(2024,c(2024+nrow(sims2.avg)-1))
sims3.avg=sims3.df %>% group_by(sim.year) %>% summarize(m.spawners=median(spawners.t),l90.spawners=quantile(spawners.t,0.05),u90.spawners=quantile(spawners.t,0.95),l80.spawners=quantile(spawners.t,0.1),u80.spawners=quantile(spawners.t,0.9),l50.spawners=quantile(spawners.t,0.25),u50.spawners=quantile(spawners.t,0.75),main='Status quo')
sims3.avg$brood.year=seq(2024,c(2024+nrow(sims3.avg)-1))
sims4.avg=sims4.df %>% group_by(sim.year) %>% summarize(m.spawners=median(spawners.t),l90.spawners=quantile(spawners.t,0.05),u90.spawners=quantile(spawners.t,0.95),l80.spawners=quantile(spawners.t,0.1),u80.spawners=quantile(spawners.t,0.9),l50.spawners=quantile(spawners.t,0.25),u50.spawners=quantile(spawners.t,0.75),main='Status quo')
sims4.avg$brood.year=seq(2024,c(2024+nrow(sims4.avg)-1))
sims5.avg=sims5.df %>% group_by(sim.year) %>% summarize(m.spawners=median(spawners.t),l90.spawners=quantile(spawners.t,0.05),u90.spawners=quantile(spawners.t,0.95),l80.spawners=quantile(spawners.t,0.1),u80.spawners=quantile(spawners.t,0.9),l50.spawners=quantile(spawners.t,0.25),u50.spawners=quantile(spawners.t,0.75),main='Status quo')
sims5.avg$brood.year=seq(2024,c(2024+nrow(sims5.avg)-1))
sims6.avg=sims6.df %>% group_by(sim.year) %>% summarize(m.spawners=median(spawners.t),l90.spawners=quantile(spawners.t,0.05),u90.spawners=quantile(spawners.t,0.95),l80.spawners=quantile(spawners.t,0.1),u80.spawners=quantile(spawners.t,0.9),l50.spawners=quantile(spawners.t,0.25),u50.spawners=quantile(spawners.t,0.75),main='Status quo')
sims6.avg$brood.year=seq(2024,c(2024+nrow(sims7.avg)-1))
sims7.avg=sims7.df %>% group_by(sim.year) %>% summarize(m.spawners=median(spawners.t),l90.spawners=quantile(spawners.t,0.05),u90.spawners=quantile(spawners.t,0.95),l80.spawners=quantile(spawners.t,0.1),u80.spawners=quantile(spawners.t,0.9),l50.spawners=quantile(spawners.t,0.25),u50.spawners=quantile(spawners.t,0.75),main='Status quo')
sims7.avg$brood.year=seq(2024,c(2024+nrow(sims7.avg)-1))
sims8.avg=sims8.df %>% group_by(sim.year) %>% summarize(m.spawners=median(spawners.t),l90.spawners=quantile(spawners.t,0.05),u90.spawners=quantile(spawners.t,0.95),l80.spawners=quantile(spawners.t,0.1),u80.spawners=quantile(spawners.t,0.9),l50.spawners=quantile(spawners.t,0.25),u50.spawners=quantile(spawners.t,0.75),main='Status quo')
sims8.avg$brood.year=seq(2024,c(2024+nrow(sims8.avg)-1))


jpeg('./outputs/Figure_3.jpeg',height=11,width=8.5,res=500,units='in',
     pointsize = 11)
par(mfrow=c(4,2))
plot(sims1.avg$m.spawners~sims1.avg$brood.year,type='l',lwd=3,ylim=c(0,600),ylab='Returning total spawners',xlab='Simulation year',main = 'S1: Status quo')
polygon(
  c(sims1.avg$brood.year, rev(sims1.avg$brood.year)),           # x coordinates forward and back
  c(sims1.avg$u90.spawners, rev(sims1.avg$l90.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(sims1.avg$brood.year, rev(sims1.avg$brood.year)),           # x coordinates forward and back
  c(sims1.avg$u80.spawners, rev(sims1.avg$l80.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.3),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(sims1.avg$brood.year, rev(sims1.avg$brood.year)),           # x coordinates forward and back
  c(sims1.avg$u50.spawners, rev(sims1.avg$l50.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.5),  # transparent fill
  border = NA             # no border line
)
text(x=sims1.avg$brood.year[nrow(sims1.avg)-3],y=sims1.avg$u90.spawners[nrow(sims1.avg)-3],'90% CI',cex=1)
text(x=sims1.avg$brood.year[nrow(sims1.avg)-3],y=sims1.avg$u80.spawners[nrow(sims1.avg)-3],'80% CI',cex=1)
text(x=sims1.avg$brood.year[nrow(sims1.avg)-3],y=sims1.avg$u50.spawners[nrow(sims1.avg)-3],'50% CI',cex=1)

plot(sims2.avg$m.spawners~sims2.avg$brood.year,type='l',lwd=3,ylim=c(0,600),ylab='Returning total spawners',xlab='Simulation Year',main = 'S2: Doubled hatchery releases')
polygon(
  c(sims2.avg$brood.year, rev(sims2.avg$brood.year)),           # x coordinates forward and back
  c(sims2.avg$u90.spawners, rev(sims2.avg$l90.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(sims2.avg$brood.year, rev(sims2.avg$brood.year)),           # x coordinates forward and back
  c(sims2.avg$u80.spawners, rev(sims2.avg$l80.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.3),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(sims2.avg$brood.year, rev(sims2.avg$brood.year)),           # x coordinates forward and back
  c(sims2.avg$u50.spawners, rev(sims2.avg$l50.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.5),  # transparent fill
  border = NA             # no border line
)

plot(sims3.avg$m.spawners~sims3.avg$brood.year,type='l',lwd=3,ylim=c(0,600),ylab='Returning total spawners',xlab='Simulation Year',main = 'S3: Halved hatchery releases')
polygon(
  c(sims3.avg$brood.year, rev(sims3.avg$brood.year)),           # x coordinates forward and back
  c(sims3.avg$u90.spawners, rev(sims3.avg$l90.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(sims3.avg$brood.year, rev(sims3.avg$brood.year)),           # x coordinates forward and back
  c(sims3.avg$u80.spawners, rev(sims3.avg$l80.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.3),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(sims3.avg$brood.year, rev(sims3.avg$brood.year)),           # x coordinates forward and back
  c(sims3.avg$u50.spawners, rev(sims3.avg$l50.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.5),  # transparent fill
  border = NA             # no border line
)

plot(sims4.avg$m.spawners~sims4.avg$brood.year,type='l',lwd=3,ylim=c(0,600),ylab='Returning total spawners',xlab='Simulation Year',main = 'S4: Quartered hatchery releases')
polygon(
  c(sims4.avg$brood.year, rev(sims4.avg$brood.year)),           # x coordinates forward and back
  c(sims4.avg$u90.spawners, rev(sims4.avg$l90.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(sims4.avg$brood.year, rev(sims4.avg$brood.year)),           # x coordinates forward and back
  c(sims4.avg$u80.spawners, rev(sims4.avg$l80.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.3),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(sims4.avg$brood.year, rev(sims4.avg$brood.year)),           # x coordinates forward and back
  c(sims4.avg$u50.spawners, rev(sims4.avg$l50.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.5),  # transparent fill
  border = NA             # no border line
)

plot(sims5.avg$m.spawners~sims5.avg$brood.year,type='l',lwd=3,ylim=c(0,600),ylab='Returning total spawners',xlab='Simulation Year',main = 'S5: No hatchery releases')
polygon(
  c(sims5.avg$brood.year, rev(sims5.avg$brood.year)),           # x coordinates forward and back
  c(sims5.avg$u90.spawners, rev(sims5.avg$l90.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(sims5.avg$brood.year, rev(sims5.avg$brood.year)),           # x coordinates forward and back
  c(sims5.avg$u80.spawners, rev(sims5.avg$l80.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.3),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(sims5.avg$brood.year, rev(sims5.avg$brood.year)),           # x coordinates forward and back
  c(sims5.avg$u50.spawners, rev(sims5.avg$l50.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.5),  # transparent fill
  border = NA             # no border line
)

plot(sims6.avg$m.spawners~sims6.avg$brood.year,type='l',lwd=3,ylim=c(0,600),ylab='Returning total spawners',xlab='Simulation Year',main = 'S6: In situ Broodtake')
polygon(
  c(sims6.avg$brood.year, rev(sims6.avg$brood.year)),           # x coordinates forward and back
  c(sims6.avg$u90.spawners, rev(sims6.avg$l90.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(sims6.avg$brood.year, rev(sims6.avg$brood.year)),           # x coordinates forward and back
  c(sims6.avg$u80.spawners, rev(sims6.avg$l80.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.3),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(sims6.avg$brood.year, rev(sims6.avg$brood.year)),           # x coordinates forward and back
  c(sims6.avg$u50.spawners, rev(sims6.avg$l50.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.5),  # transparent fill
  border = NA             # no border line
)


plot(sims7.avg$m.spawners~sims7.avg$brood.year,type='l',lwd=3,ylim=c(0,600),ylab='Returning total spawners',xlab='Simulation Year',main = 'S7: Ex Situ Broodtake (60k smolt target)')
polygon(
  c(sims7.avg$brood.year, rev(sims7.avg$brood.year)),           # x coordinates forward and back
  c(sims7.avg$u90.spawners, rev(sims7.avg$l90.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(sims7.avg$brood.year, rev(sims7.avg$brood.year)),           # x coordinates forward and back
  c(sims7.avg$u80.spawners, rev(sims7.avg$l80.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.3),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(sims7.avg$brood.year, rev(sims7.avg$brood.year)),           # x coordinates forward and back
  c(sims7.avg$u50.spawners, rev(sims7.avg$l50.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.5),  # transparent fill
  border = NA             # no border line
)

plot(sims8.avg$m.spawners~sims8.avg$brood.year,type='l',lwd=3,ylim=c(0,600),ylab='Returning total spawners',xlab='Simulation Year',main = 'S8: Ex Situ Broodtake (120k smolt target)')
polygon(
  c(sims8.avg$brood.year, rev(sims8.avg$brood.year)),           # x coordinates forward and back
  c(sims8.avg$u90.spawners, rev(sims8.avg$l90.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(sims8.avg$brood.year, rev(sims8.avg$brood.year)),           # x coordinates forward and back
  c(sims8.avg$u80.spawners, rev(sims8.avg$l80.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.3),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(sims8.avg$brood.year, rev(sims8.avg$brood.year)),           # x coordinates forward and back
  c(sims8.avg$u50.spawners, rev(sims8.avg$l50.spawners)),   # y upper then lower reversed
  col = adjustcolor('darkgray',alpha.f=0.5),  # transparent fill
  border = NA             # no border line
)
dev.off()

#Supplemental Figures####
#timeseries
##S1 adult returns####
jpeg('./outputs/Figure_S1_returns.jpeg',height=5,width=8,res=500,units='in')
par(mfrow=c(1,1))
plot(data$total_spawners~data$ï..year,type='l',ylab='Adult Returns',bty='l',xlab='Year',lwd=3,col=adjustcolor('darkgray',alpha.f=0.9))
points(data$total_spawners~data$ï..year,pch=21,cex=1.5,col=adjustcolor('darkgray',alpha.f=0.9))
lines(data$unknown_origin_spawners~data$ï..year,lwd=2,col=adjustcolor('purple4'))
points(data$unknown_origin_spawners~data$ï..year,pch=21,cex=1.5,col=adjustcolor('purple4'))
lines(data$natural_origin_spawners~data$ï..year,lwd=2,col=adjustcolor('turquoise4'))
points(data$natural_origin_spawners~data$ï..year,pch=21,cex=1.5,col=adjustcolor('turquoise4'))
lines(data$clipped_spawners~data$ï..year,lwd=2,col=adjustcolor('darkred',alpha.f=0.9))
points(data$clipped_spawners~data$ï..year,pch=21,cex=1.5,col=adjustcolor('darkred'))
text(x=par('usr')[1]+2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.05,'Total',col='darkgray',adj=0)
text(x=par('usr')[1]+2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.1,'Hatchery-origin',col='darkred',adj=0)
text(x=par('usr')[1]+2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.15,'Natural-origin',col='turquoise4',adj=0)
text(x=par('usr')[1]+2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.2,'Unknown-origin',col='purple4',adj=0)
dev.off()

##S2 hatchery releases####
data$total_fry=apply(data[,8:9],1,sum,na.rm=T)

jpeg('./outputs/Figure_S2_releases.jpeg',height=10,width=8,res=500,units='in')
par(mfrow=c(2,1))
plot(data$total_fry~data$ï..year,type='l',ylab='Fry Releases',bty='l',xlab='Year',lwd=3,col=adjustcolor('darkgray',alpha.f=0.9))
points(data$total_fry~data$ï..year,pch=21,cex=1.5,col=adjustcolor('darkgray',alpha.f=0.9))
lines(data$fry_release_captive_brood_unclpd~data$ï..year,lwd=2,col=adjustcolor('purple4'))
points(data$fry_release_captive_brood_unclpd~data$ï..year,pch=21,cex=1.5,col=adjustcolor('purple4'))
lines(data$fry_release_captive_brood_clpd~data$ï..year,lwd=2,col=adjustcolor('darkred',alpha.f=0.9))
points(data$fry_release_captive_brood_clpd~data$ï..year,pch=21,cex=1.5,col=adjustcolor('darkred'))
text(x=par('usr')[1]+2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.05,'Total',col='darkgray',adj=0)
text(x=par('usr')[1]+2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.1,'Clipped',col='darkred',adj=0)
text(x=par('usr')[1]+2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.15,'Unclipped',col='purple4',adj=0)
mtext(side=3,adj=0,'A')

plot(data$smolt_release~data$ï..year,type='l',ylab='Smolt Releases',bty='l',xlab='Year',lwd=3)
points(data$smolt_release~data$ï..year,pch=21,cex=1.5)
mtext(side=3,adj=0,'B')
dev.off()

             

##S3 smolt counts####
jpeg('./outputs/Figure_S3_smolts.jpeg',height=5,width=8,res=500,units='in')
plot(data$smolts_cnt_clipped_fry~data$ï..year,type='l',ylab='Smolt Counts',bty='l',xlab='Year',lwd=2,col='darkred')
points(data$smolts_cnt_clipped_fry~data$ï..year,pch=21,cex=1.2,col='darkred')
lines(data$smolts_cnt_naturalorigin~data$ï..year,lwd=2,col='turquoise4')
points(data$smolts_cnt_naturalorigin~data$ï..year,pch=21,cex=1.2,col='turquoise4')
lines(data$smolts_cnt_unclipped~data$ï..year,type='l',col='purple4',lwd=2)
points(data$smolts_cnt_unclipped~data$ï..year,pch=21,cex=1.2,col='purple4')
text(x=par('usr')[2]-5,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.05,'Hatchery-origin',col='darkred',adj=0)
text(x=par('usr')[2]-5,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.1,'Natural-origin',col='turquoise4',adj=0)
text(x=par('usr')[2]-5,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.15,'Unknown-origin',col='purple4',adj=0)
dev.off()


#S4. Ricker fit
sim.ric=MASS::mvrnorm(n=10000,mu=mu.ric,Sigma=vcov.ric)
sim.ric[,2:3]=exp(sim.ric[,2:3]) #exponentiate from log-scale


#predicted mean and bootstrapped 90% intervals for predicted smolts by spawner abundance
s_pred=seq(0,max(smolt_bt$S))
sm_pred=exp(sak.alpha[1]-sak.beta[1]*s_pred)*s_pred
s.pred.l=matrix(nrow=length(s_pred),ncol=10000)
for(i in 1:10000){
  s.pred.l[,i]=exp(sim.ric[i,1]-s_pred/sim.ric[i,2])*s_pred
}
s.pred.lower=apply(s.pred.l,1,quantile,0.05)
s.pred.upper=apply(s.pred.l,1,quantile,0.95)

##S4: Natural Origin Smolt production####
jpeg('./outputs/Figure_S4_smolt_production.jpeg',height=5,width=10,res=500,units='in')
par(mfrow=c(1,2))
plot(R~S,smolt_bt,type='n',xlab='Total Spawners',ylab='Natural Origin Smolts',ylim=c(0,4e4))
lines(R~S,lwd=0.5,col=adjustcolor('darkgray',alpha.f = 0.5),data=smolt_bt)
lines(sm_pred~s_pred,lwd=2)
lines(s.pred.lower~s_pred,lty=5)
lines(s.pred.upper~s_pred,lty=5)

points(R~S,pch=21,cex=1.5,data=smolt_bt)
lines(c(5000,5000)~c(-50,100))
lines(c(-1200,5000)~c(100,100))

text(x=smolt_bt$S-40,y=smolt_bt$R+700,smolt_bt$by,cex=0.7)

plot(R~S,smolt_bt,type='n',xlab='Total Spawners',ylab='Natural Origin Smolts',xlim=c(0,100),ylim=c(0,5000))
lines(R~S,lwd=0.5,col=adjustcolor('darkgray',alpha.f = 0.5),data=smolt_bt)
lines(sm_pred~s_pred,lwd=2)
points(R~S,pch=21,cex=1.5,data=smolt_bt)
text(x=smolt_bt$S-4,y=smolt_bt$R+200,smolt_bt$by,cex=0.7)
dev.off()

##S5 Lake survival####
jpeg('./outputs/Figure_S5_Lake_survival.jpeg',height=5,width=8,res=500,units='in')
plot(data$fry_to_smolt_surv~data$ï..year,type='l',ylab='Release Survival to Out-migration (%)',bty='l',xlab='Year',lwd=2,col='navy',ylim=c(0,100))
points(data$fry_to_smolt_surv~data$ï..year,pch=21,cex=1.2,col='navy')
lines(data$smolt_release_emig_survival~data$ï..year,lwd=2,col='deepskyblue')
points(data$smolt_release_emig_survival~data$ï..year,pch=21,cex=1.2,col='deepskyblue')
text(x=par('usr')[2]-5,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.05,'Fry',col='navy',adj=0)
text(x=par('usr')[2]-5,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.1,'Smolt',col='deepskyblue',adj=0)
dev.off()

##S6 Marine survival time-series####
jpeg('./outputs/Figure_S6_Marine_survival.jpeg',height=5,width=8,res=500,units='in')
par(mfrow=c(1,1))
plot(data$marine_surv_clipped~data$ï..year,xlab='Outmigration Year',ylab='Marine Survival (%)',type='l',lwd=2,col='darkred',ylim=c(0,1),bty='l')
lines(data$marine_surv_wild~data$ï..year,lwd=2,col='turquoise4')
lines(data$marine_surv_unclipped~data$ï..year,lwd=2,col='purple4')
lines(data$bigQsurvival~data$ï..year,lwd=2,col='goldenrod')
points(data$marine_surv_clipped~data$ï..year,col='darkred',cex=1.2)
points(data$marine_surv_wild~data$ï..year,col='turquoise4',cex=1.2)
points(data$marine_surv_unclipped~data$ï..year,col='purple4',cex=1.2)
points(data$bigQsurvival~data$ï..year,col='goldenrod',cex=2)

text(x=par('usr')[2]-10,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.05,'Hatchery-origin',col='darkred',adj=0)
text(x=par('usr')[2]-10,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.1,'Natural-origin',col='turquoise4',adj=0)
text(x=par('usr')[2]-10,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.15,'Unknown-origin',col='purple4',adj=0)
text(x=par('usr')[2]-10,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.2,'BQ releases',col='goldenrod',adj=0)
dev.off()

##S7 MS distributions####
u.ms=data[complete.cases(data[,16]),16]
mu.u.ms=mean(u.ms)
sd.u.ms=sd(u.ms)
sd.log.u.ms=sqrt(log(1+(sd.u.ms/mean(u.ms))^2))
mu.log.u.ms=log(mu.u.ms)-0.5*sd.log.u.ms^2
u.ms.sim=exp(rnorm(1e3,mu.log.u.ms,sd.log.u.ms))


bigQms=data$bigQsurvival[complete.cases(data$bigQsurvival)]
mu.bq.ms=mean(bigQms)
sd.bq.ms=sd(bigQms)
sd.log.bq.ms=sqrt(log(1+(sd.bq.ms/mean(bigQms))^2))
mu.log.bq.ms=log(mu.bq.ms)-0.5*sd.log.bq.ms^2
bq.ms=exp(rnorm(1e3,mu.log.bq.ms,sd.log.bq.ms))

jpeg('./outputs/Figure_S7_MS_hists.jpeg',height=5,width=10,res=500,units='in')
par(mfrow=c(1,2))
d.h=density(ms.dists[,1],bw=0.05)
d.u=density(u.ms.sim,bw=0.05)
d.w=density(ms.dists[,2],bw=0.05)
d.bq=density(bq.ms,bw=0.05)

plot(d.h$y~d.h$x,type='l',col='darkred',lwd=2,yaxt='n',ylab='',xlim=c(0,1),xlab='Simulated Distributions of Marine Survival (%)')
lines(d.w$y~d.w$x,lwd=2,col='turquoise4')
lines(d.u$y~d.u$x,lwd=2,col='purple4')
lines(d.bq$y~d.bq$x,lwd=2,col='goldenrod')
text(x=par('usr')[2]-0.3,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.05,'Hatchery-origin',col='darkred',adj=0)
text(x=par('usr')[2]-0.3,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.1,'Natural-origin',col='turquoise4',adj=0)
text(x=par('usr')[2]-0.3,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.15,'Unknown-origin',col='purple4',adj=0)
text(x=par('usr')[2]-0.3,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.2,'BQ releases',col='goldenrod',adj=0)

#hatchery vs wild survival
plot(log10(ms.dists[,1])~log10(ms.dists[,2]),ylab='',xlab='',xlim=c(-2,0.8),ylim=c(-2,0.8),yaxt='n',xaxt='n')
minor.ticks=as.vector(1:9 %o% 10^(-2:0))
axis(side=1,at=c(-2,-1,0),labels=c(10^-2,10^-1,10^0))
axis(side=1,at=log10(minor.ticks),labels=NA,tcl=-0.2)
axis(side=2,at=c(-2,-1,0),labels=c(10^-2,10^-1,10^0))
axis(side=2,at=log10(minor.ticks),labels=NA,tcl=-0.2)

mtext(side=1,'Simulated Natural-origin Marine Survival (%)',col='turquoise4',line=2.5)
mtext(side=2,'Simulated Hatchery-origin Marine Survival (%)',col='darkred',line=2.5)
dev.off()

#fixed parameter distributions####
#plot distribution
x <- seq(0, qlnorm(0.999, meanlog = mu.log.fry, sdlog = sd.log.fry),
         length.out = 1000)
plot(
  x,
  dlnorm(x, meanlog = mu.log.fry, sdlog = sd.log.fry),
  type = "l",
  lwd = 2,
  col = "blue",
  xlab = "Fry releases",
  ylab = "Density",
  yaxt='n',
  main = ""
)

x <- seq(0, qlnorm(0.999, meanlog = mu.log.smolts, sdlog = sd.log.smolts),
         length.out = 1000)

plot(
  x,
  dlnorm(x, meanlog = mu.log.smolts, sdlog = sd.log.smolts),
  type = "l",
  lwd = 2,
  col = "blue",
  xlab = "Smolt releases",
  ylab = "Density",
  yaxt='n',
  main = ""
)

x <- seq(-1.4, qnorm(0.999, mean = mu.logit.smolt.surv, sd = sd.logit.smolt.surv),
         length.out = 1000)

plot(
  plogis(x)*100,
  dnorm(x, mean = mu.logit.smolt.surv, sd = sd.logit.smolt.surv),
  type = "l",
  lwd = 2,
  col = "blue",
  xlab = "Smolt Lake survival (%)",
  ylab = "Density",
  yaxt='n',
  main = ""
)

x <- seq(-5, qnorm(0.999, mean = mu.logit.fry.surv, sd = sd.logit.fry.surv),
         length.out = 1000)

plot(
  plogis(x)*100,
  dnorm(x, mean = mu.logit.fry.surv, sd = sd.logit.fry.surv),
  type = "l",
  lwd = 2,
  col = "blue",
  xlab = "Fry Lake survival (%)",
  ylab = "Density",
  yaxt='n',
  main = ""
)

#fecundity and egg-to-fry survival
par(mfrow=c(2,2))
x <- seq(1e3, 3e3,
         length.out = 1000)

plot(
  x,
  dnorm(x, mean = 2050, sd = 250),
  type = "l",
  lwd = 2,
  col = "blue",
  xlab = "Mean Fecundity",
  ylab = "Density",
  yaxt='n',
  main = ""
)

x <- seq(0.5, 2,
         length.out = 1000)

plot(
  plogis(x)*100,
  dnorm(x, mean = qlogis(0.78), sd = 0.2),
  type = "l",
  lwd = 2,
  col = "blue",
  xlab = "Egg-to-Fry Survival (%)",
  ylab = "Density",
  yaxt='n',
  main = ""
)

x <- seq(0.5, 2,
         length.out = 1000)

plot(
  plogis(x)*100,
  dnorm(x, mean = qlogis(0.8), sd = 0.2),
  type = "l",
  lwd = 2,
  col = "blue",
  xlab = "Fry-to-Smolt Survival (%)",
  ylab = "Density",
  yaxt='n',
  main = ""
)

