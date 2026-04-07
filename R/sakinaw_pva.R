#Sakinaw Sockeye - PVA

#sections 1-4: simulating conditions check
#smolt-spawner curve
#marine survial
#fry/smolt releases & survivals
#simulations: baseline, 2x captive brood production, alteration to broodtake


#Ricker spawner-smolt curve####
smolt_bt=read.csv('./data/smolt_broodtable.csv')
smolt_bt=smolt_bt[,1:3]
colnames(smolt_bt)=c('by','S','R')
smolt_bt=smolt_bt[complete.cases(smolt_bt$by),]
smolt_bt$logRS=log(smolt_bt$R/smolt_bt$S)

smolt_ric=samEst::ricker_TMB(smolt_bt,ac=TRUE,silent=T)
sak.alpha=c(smolt_ric$sd_report$value[1],smolt_ric$sd_report$sd[1])
sak.beta=c(smolt_ric$sd_report$value[2],smolt_ric$sd_report$sd[2])
sak.rho=c(smolt_ric$sd_report$value[3],smolt_ric$sd_report$sd[3])
sak.sigma=c(smolt_ric$sd_report$value[5],smolt_ric$sd_report$sd[5])

s_pred=seq(0,max(smolt_bt$S))
sm_pred=exp(sak.alpha[1]-sak.beta[1]*s_pred)*s_pred

plot(R~S,smolt_bt,type='n',xlab='Total Spawners',ylab='Natural Origin Smolts')
lines(R~S,lwd=0.5,col=adjustcolor('darkgray',alpha.f = 0.5),data=smolt_bt)
lines(sm_pred~s_pred,lwd=2)
points(R~S,pch=21,cex=1.5,data=smolt_bt)
#text(x=smolt_bt$S*0.98,y=smolt_bt$R*1.08,smolt_bt$by,cex=0.7)

#Marine survival####
data=read.csv('./data/sakinaw_abundances.csv')

#mean, sd, correlations among survival estimates
survivals=data[,14:16]
surv.cor=cor(survivals,use='pairwise.complete.obs')
surv.sd=apply(survivals,2,sd,na.rm=T)
surv.mean=apply(survivals,2,mean,na.rm=T)
sd.log.surv=sqrt(log(1+(surv.sd/surv.mean)^2))
mu.log.surv=log(surv.mean)-0.5*sd.log.surv^2

#variance-covariance matrix for wild/hatchery/mixed marine survival
surv.mat=matrix(ncol=3,nrow=3)
surv.mat[1,1]=log.surv.sd[1]
surv.mat[2,2]=log.surv.sd[2]
surv.mat[3,3]=log.surv.sd[3]
surv.mat[1,2]=surv.cor[1,2]*log.surv.sd[1]*log.surv.sd[2]
surv.mat[2,1]=surv.cor[1,2]*log.surv.sd[1]*log.surv.sd[2]
surv.mat[1,3]=surv.cor[1,3]*log.surv.sd[1]*log.surv.sd[3]
surv.mat[3,1]=surv.cor[1,3]*log.surv.sd[1]*log.surv.sd[3]
surv.mat[2,3]=0
surv.mat[3,2]=0

#specify 2x2 covariance matrices for each part of the simulation
surv.mat.h.u=matrix(data=c(surv.mat[1,1],surv.mat[1,2],surv.mat[2,1],surv.mat[2,2]),ncol=2,nrow=2)
surv.mat.h.n=matrix(data=c(surv.mat[1,1],surv.mat[1,3],surv.mat[3,1],surv.mat[3,3]),ncol=2,nrow=2)

ms.log=MASS::mvrnorm(n=1e3,mu=c(log(surv.mean[1])-0.5*log.surv.sd[1]^2,log(surv.mean[2])-0.5*log.surv.sd[2]^2,log(surv.mean[3])-0.5*log.surv.sd[3]^2),Sigma=surv.mat,empirical=TRUE)
ms.dists=exp(ms.log)

#check empirical vs simulated means
surv.mean
apply(ms.dists,2,mean)

#autocorrelation in marine survival
rho.ms.h=acf(log(survivals[,1]+0.01),na.action=na.pass,plot=F)$acf[2] #1-yr autocorrelation - hatchery survival
rho.ms.n=acf(log(survivals[,3]+0.01),na.action=na.pass,plot=F)$acf[2] #1-yr autocorrelation - wild survival
se.rho.ms.h=1/sqrt(length(na.omit(survivals[,1]))) #standard error estimate for rho
se.rho.ms.n=1/sqrt(length(na.omit(survivals[,3]))) #standard error estimate for rho

#scale both quantities to logistic scale
logistic.rho.ms.h=logistic_transf(rho.ms.h) #one_to_minus_one transformation (logistic)
logistic.rho.ms.n=logistic_transf(rho.ms.n)
#SE using delta method
logistic.se.rho.ms.h=2/(1-rho.ms.h^2)*se.rho.ms.h
logistic.se.rho.ms.n=2/(1-rho.ms.n^2)*se.rho.ms.n

#check for consistency - hatchery rho
mean(inv_logistic_transf(rnorm(1e6,logistic.rho.ms.h,logistic.se.rho.ms.h)))
sd(inv_logistic_transf(rnorm(1e6,logistic.rho.ms.h,logistic.se.rho.ms.h)))

#plot 1: Marine survival time-series
par(mfrow=c(1,1))
plot(data$marine_surv_clipped~data$ï..year,xlab='Year',ylab='Marine Survival (%)',type='l',lwd=2,col='darkred',ylim=c(0,1))
lines(data$marine_surv_wild~data$ï..year,lwd=2,col='turquoise4')
lines(data$marine_surv_unclipped~data$ï..year,lwd=2,col='purple4')
points(data$marine_surv_clipped~data$ï..year,col='darkred',cex=2)
points(data$marine_surv_wild~data$ï..year,col='turquoise4',cex=2)
points(data$marine_surv_unclipped~data$ï..year,col='purple4',cex=2)
text(x=par('usr')[2]-4,y=par('usr')[4]-0.05,'Origin:',cex=1.4)
text(x=par('usr')[2]-4,y=par('usr')[4]-0.1,'Hatchery',col='darkred',cex=1.4)
text(x=par('usr')[2]-4,y=par('usr')[4]-0.15,'Wild',col='turquoise4',cex=1.4)
text(x=par('usr')[2]-4,y=par('usr')[4]-0.2,'Unknown',col='purple4',cex=1.4)


#Marine survival simulation conditions####
par(mfrow=c(1,2))

d.h=density(ms.dists[,1],bw=0.05)
d.u=density(ms.dists[,2],bw=0.05)
d.w=density(ms.dists[,3],bw=0.05)

plot(d.h$y~d.h$x,type='l',col='darkred',lwd=2,yaxt='n',ylab='',xlim=c(0,1),xlab='Simulated Distribution of Marine Survival (%)')
lines(d.w$y~d.w$x,lwd=2,col='turquoise4')
lines(d.u$y~d.m$x,lwd=2,col='purple4')

#hatchery vs wild survival
plot(log10(ms.dists[,1])~log10(ms.dists[,3]),ylab='',xlab='',xlim=c(-2,0.8),ylim=c(-2,0.8),yaxt='n',xaxt='n')
axis(side=1,at=c(-2,-1,0,log10(3)),labels=c(10^-2,10^-1,10^0,3))
axis(side=2,at=c(-2,-1,0,log10(3)),labels=c(10^-2,10^-1,10^0,3))
mtext(side=1,'Wild-origin Marine Survival (%)',col='turquoise4',line=2.5)
mtext(side=2,'Hatchery-origin Marine Survival (%)',col='darkred',line=2.5)

#Fry and Smolt releases####
#fry release
total.fry=apply(data[,7:8],1,sum,na.rm = TRUE)
total.fry=total_fry[total_fry>0]

sd.fry=sd(total.fry)
mu.fry=mean(total.fry)
sd.log.fry=sqrt(log(1+(sd.fry/mu.fry)^2))
mu.log.fry=mean(log(total.fry))

fry.frac=data[c(nrow(data)-2):nrow(data),7:8] #look at last 3 years of clipped/unclipped fry releases
fry.frac[,3]=apply(fry.frac[,1:2],1,sum)
fry.frac[,4]=fry.frac[,1]/fry.frac[,3] #proportion clipped

mu.logit.fry.frac=mean(qlogis(fry.frac[,4]))
sd.logit.fry.frac=sd(qlogis(fry.frac[,4]))

mean(plogis(rnorm(1e6,mu.logit.fry.frac,sd.logit.fry.frac))) #simulated
mean(fry.frac[,4]) #empirical mean 

#fry survival to smolt emigration
mu.fry.surv=mean(data$fry_to_smolt_surv/100,na.rm=T)
sd.fry.surv=sd(data$fry_to_smolt_surv/100,na.rm=T)

mu.logit.fry.surv=mean(qlogis(data$fry_to_smolt_surv/100),na.rm=T)
sd.logit.fry.surv=sd(qlogis(data$fry_to_smolt_surv/100),na.rm=T)

mean(plogis(rnorm(1e6,mu.logit.fry.surv,sd.logit.fry.surv))) #simulated
mu.fry.surv #empirical mean

#steps 1b: smolt releases + survival to emigration
total.smolts=data$smolt_cnt_release_clipped
total.smolts=total_smolts[complete.cases(total_smolts)]

mu.smolts=mean(total.smolts)
sd.smolts=sd(total.smolts)
sd.log.smolts=sqrt(log(1+(sd.smolts/mu.smolts)^2))
mu.log.smolts=mean(log(total.smolts))

sd(exp(rnorm(1e6,mu.log.smolts,sd.log.smolts)))

#smolt survival to smolt emigration
mu.smolt.surv=mean(data$smolt_release_emig_survival/100,na.rm=T)
sd.smolt.surv=sd(data$smolt_release_emig_survival/100,na.rm=T)

mu.logit.smolt.surv=mean(qlogis(data$smolt_release_emig_survival/100),na.rm=T)
sd.logit.smolt.surv=sd(qlogis(data$smolt_release_emig_survival/100),na.rm=T)

mean(plogis(rnorm(1000,logit.mu.smolt.surv,logit.sd.smolt.surv)))




#Sim 1: Baseline Scenario####
sims1=list()
n.sims=500
sim.length=24
sim.summary1=data.frame(sim.id=seq(1:n.sims),py.tspwn.500=NA,py.nspwn.500=NA,py.tspwn.100=NA,py.nspwn.100=NA,m.tspwn=NA,m.nspwn=NA)
for(i in 1:n.sims){
  sim.dat=data.frame(sim.year=seq(1:c(sim.length+2)),spawners.t=0,spawners.u=0,spawners.c=0,spawners.n=c(rep(NA,3),rep(0,sim.length-1)),smolts.u=NA,smolts.c=NA,smolts.n=NA,ms.h=NA,ms.u=NA,ms.n=NA,fry.rel.u=NA,fry.rel.c=NA,smolt.rel.c=NA,fry.surv=NA,smolt.surv=NA,ricker.logrs=NA,ricker.mu=NA,ricker.eps=NA)
  #first X years - estimate returns based on known quantities in the last several years
  
  #initial spawners - conditioned on observed data
  #total spawners
  sim.dat[1:2,2]=data$total_spawners[c(nrow(data)-1):nrow(data)]
  #unclipped spawners
  sim.dat[1:2,3]=data$unclipped_spawners[c(nrow(data)-1):nrow(data)]
  #clipped spawners
  sim.dat[1:2,4]=data$clipped_spawners[c(nrow(data)-1):nrow(data)]
  
  #initial smolt counts - conditioned on observed data
  #unclipped smolts (unc. fry and natural origin):
  sim.dat[1:2,6]=data$smolts_cnt_unclipped[c(nrow(data)-1):nrow(data)]
  #clipped smolts (both clp. fry and smolt releases):
  sim.dat[1:2,7]=data$smolts_cnt_clipped_fry[c(nrow(data)-1):nrow(data)]+data$smolt_cnt_release_clipped[c(nrow(data)-1):nrow(data)]
  
  #hatchery releases
  #unclipped fry
  sim.dat[1:2,12]=data$fry_release_captive_brood_unclpd[c(nrow(data)-1):nrow(data)]
  #clipped fry
  sim.dat[1:2,13]=data$fry_release_captive_brood_unclpd[c(nrow(data)-1):nrow(data)]
  #clipped smolts
  sim.dat[1:2,14]=data$smolt_release[c(nrow(data)-1):nrow(data)]
  
  #hatchery lake survival
  #fry survival to smolt emigration
  sim.dat[1:2,15]=data$fry_to_smolt_surv[c(nrow(data)-1):nrow(data)]/100
  #smolt survival to smolt emigration
  sim.dat[1:2,16]=data$smolt_release_emig_survival[c(nrow(data)-1):nrow(data)]/100
  
  ##Simulate - marine survival for clipped/unclipped smolts - spawners t+2
  #hatchery and unknown origin survival for these smolt cohorts
  
  #simulate for marine survival for each set of smolts
  #include autocorrelated previous years marine survival deviate for each state
  rho.ms.h.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.h,logistic.se.rho.ms.h))
  mu.adj.h=rho.ms.h.sim*(log(data$marine_surv_clipped[nrow(data)-2])-log(surv.mean)[1])
  mu.adj.u=rho.ms.h.sim*(log(data$marine_surv_unclipped[nrow(data)-2])-log(surv.mean)[2])
  
  sim.dat[1:2,9:11]=exp(MASS::mvrnorm(n=2,mu=mu.log.surv,Sigma=surv.mat))
  
  #estimate returning spawners at t+2 for each smolt cohort
  #unclipped spawners - unclipped outmigrating smolts x unknown origin marine survival
  sim.dat[3:4,3]=sim.dat[1:2,6]*(sim.dat[1:2,11]/100)
  #clipped spawners - clipped outmigrating smolts x hatchery origin marine survival
  sim.dat[3:4,4]=sim.dat[1:2,7]*(sim.dat[1:2,9]/100)
  #total spawners
  sim.dat[3:4,2]=sim.dat[3:4,3]+sim.dat[3:4,4]
  
  #need to add unclipped natural smolts for sim year 3 - these will be stored as 'smolts.n'
  
  #for all future states - next 25 years
  for(t in 1:c(sim.length+2)){
    
    sim.dat[t,2]=sum(sim.dat[t,3:5],na.rm=T) #total spawners
    
    #natural origin smolts from total spawners:
    
    #randomly sample Ricker parameters from estimates and covariances
    sim.ric=MASS::mvrnorm(n=1,mu=mu.ric,Sigma=vcov.ric)
    sim.ric[2:3]=exp(sim.ric[2:3]) #exponentiate from log-scale
    sim.ric[4]=inv_logistic_transf(sim.ric[4]) #function to scale transformed rho to -1 to 1
    #estimate natural origin smolts t+2 based on Ricker model
    if(t==1){
      #no autocorrelation if first smolt cohort (no previous estimate available to condition autocorrelation in smolt recruitment)
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2],sd=sim.ric[3])
    }else{
      #autocorrelation included
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2]+sim.ric[4]*sim.dat[t-1,19],sd=sim.ric[3])
    }
    #if total spawners < 1 = no natural smolt production
    if(sim.dat[t,2]<2){
      sim.dat[t+2,8]=0
    }else{
      sim.dat[t+2,8]=exp(logRS)*sim.dat[t,2]
    }
    sim.dat[t,17]=logRS
    sim.dat[t,18]=sim.ric[1]-sim.dat[t,2]/sim.ric[2]
    sim.dat[t,19]=logRS-(sim.ric[1]-sim.dat[t,2]/sim.ric[2])
    
    #fry & smolt releases
    if(t>2){
      #fry releases - total fry releases x allocation to clipped/unclipped status
      total_fry.t=exp(rnorm(1,mu.log.fry,sd.log.fry)) #log-normal draws for fry releases
      total_fry.t=min(total_fry.t,2e6) #set capacity for fry release of 2 million
      frac.clipped=plogis(rnorm(1,mu.logit.fry.frac,sd.logit.fry.frac)) #proportion clipped/unclipped fry releases
      sim.dat[t,12]=total_fry.t*(1-frac.clipped) #unclipped fry releases
      sim.dat[t,13]=total_fry.t*frac.clipped #clipped fry releases
      sim.dat[t,14]=exp(rnorm(1,mu.log.smolts,sd.log.smolts)) #log-normal draws for smolt releases
      
      #survival of natural, unclipped/clipped fry, clipped smolts from lake to outmigration
      sim.dat[t,15]=plogis(rnorm(1,mu.logit.fry.surv,sd.logit.fry.surv)) #fry survival - applies to both clipped/unclipped
      sim.dat[t,16]=plogis(rnorm(1,mu.logit.smolt.surv,sd.logit.smolt.surv)) #smolt release survival to outmigration
      sim.dat[t,6]=sim.dat[t-1,12]*sim.dat[t,15] #unclipped fry releases x fry survival for 3rd (1st) sim year
      sim.dat[t,7]=sim.dat[t-1,13]*sim.dat[t,15]+sim.dat[t-1,14]*sim.dat[t,16] #clipped fry releases x (same) fry survival + clipped smolt releases x smolt emigration survival
      
      #expected returns (t+2) based on marine survival - natural origin & hatchery origin fish
      rho.ms.h.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.h,logistic.se.rho.ms.h)) #draw autocorrelation coefficient for hatchery origin smolts
      rho.ms.n.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.n,logistic.se.rho.ms.n)) #draw autocorrelation coefficient for natural origin smolts
      
      mu.adj=c(rho.ms.h.sim*(log(sim.dat[t-1,9])-log(surv.mean)[1]), #adjustment to mean survival based on 1-y autocorrelation
               0, #unknown origin marine survival not used beyond simulation year 2, set to 0 for convenience
               rho.ms.n.sim*(log(sim.dat[t-1,10])-log(surv.mean)[3])) #adjustment to mean survival based on 1-y autocorrelation
      
      sim.dat[t,9:11]=exp(MASS::mvrnorm(n=1,mu=mu.log.surv+mu.adj,Sigma=surv.mat)) #marine survival
      
      #partition based on age at maturity:
      p_age=extraDistr::rdirichlet(1,c(3,87,10))
      #3% age 3, 87% age 4, 10% age 5
      #proportion recruits, add to each return year based on randomized maturity
      #age 3
      sim.dat[t+1,3:4]=sim.dat[t+1,3:4]+sim.dat[t,6:7]*(sim.dat[t,9]/100)*p_age[1] #unclipped/clipped hatchery returns
      sim.dat[t+1,5]=sim.dat[t+1,5]+sim.dat[t,8]*(sim.dat[t,11]/100)*p_age[1] #natural origin returns
      #age 4
      sim.dat[t+2,3:4]=sim.dat[t+2,3:4]+sim.dat[t,6:7]*(sim.dat[t,9]/100)*p_age[2] #unclipped/clipped hatchery returns
      sim.dat[t+2,5]=sim.dat[t+2,5]+sim.dat[t,8]*(sim.dat[t,11]/100)*p_age[2] #natural origin returns
      #age 5
      sim.dat[t+3,3:4]=sim.dat[t+3,3:4]+sim.dat[t,6:7]*(sim.dat[t,9]/100)*p_age[3] #unclipped/clipped hatchery returns
      sim.dat[t+3,5]=sim.dat[t+3,5]+sim.dat[t,8]*(sim.dat[t,11]/100)*p_age[3] #natural origin returns
    }
  }
  sim.dat=sim.dat[1:c(sim.length+2),] #trim extra brood years
  sim.summary1[i,2]=sum(sim.dat[3:c(sim.length+2),2]>500)/sim.length #for simulated part of series - proportion of years where total spawners > 500
  sim.summary1[i,3]=sum(sim.dat[5:c(sim.length+2),5]>500)/c(sim.length-2) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary1[i,4]=sum(sim.dat[3:c(sim.length+2),2]>100)/sim.length #for simulated part of series - proportion of years where total spawners > 500
  sim.summary1[i,5]=sum(sim.dat[5:c(sim.length+2),5]>100)/c(sim.length-2) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary1[i,6]=exp(mean(log(sim.dat[3:c(sim.length+2),2])))
  sim.summary1[i,7]=exp(mean(log(sim.dat[5:c(sim.length+2),5])))
  
  sims1[[i]]=sim.dat #store simulation
}

#Sim 2: No Hatchery####
sims2=list()
n.sims=500
sim.length=24
sim.summary2=data.frame(sim.id=seq(1:n.sims),py.tspwn.500=NA,py.nspwn.500=NA,py.tspwn.100=NA,py.nspwn.100=NA,m.tspwn=NA,m.nspwn=NA)
for(i in 1:n.sims){
  sim.dat=data.frame(sim.year=seq(1:c(sim.length+2)),spawners.t=0,spawners.u=0,spawners.c=0,spawners.n=c(rep(NA,3),rep(0,sim.length-1)),smolts.u=NA,smolts.c=NA,smolts.n=NA,ms.h=NA,ms.u=NA,ms.n=NA,fry.rel.u=NA,fry.rel.c=NA,smolt.rel.c=NA,fry.surv=NA,smolt.surv=NA,ricker.logrs=NA,ricker.mu=NA,ricker.eps=NA)
  #first X years - estimate returns based on known quantities in the last several years
  
  #initial spawners - conditioned on observed data
  #total spawners
  sim.dat[1:2,2]=data$total_spawners[c(nrow(data)-1):nrow(data)]
  #unclipped spawners
  sim.dat[1:2,3]=data$unclipped_spawners[c(nrow(data)-1):nrow(data)]
  #clipped spawners
  sim.dat[1:2,4]=data$clipped_spawners[c(nrow(data)-1):nrow(data)]
  
  #initial smolt counts - conditioned on observed data
  #unclipped smolts (unc. fry and natural origin):
  sim.dat[1:2,6]=data$smolts_cnt_unclipped[c(nrow(data)-1):nrow(data)]
  #clipped smolts (both clp. fry and smolt releases):
  sim.dat[1:2,7]=data$smolts_cnt_clipped_fry[c(nrow(data)-1):nrow(data)]+data$smolt_cnt_release_clipped[c(nrow(data)-1):nrow(data)]
  
  #hatchery releases
  #unclipped fry
  sim.dat[1:2,12]=data$fry_release_captive_brood_unclpd[c(nrow(data)-1):nrow(data)]
  #clipped fry
  sim.dat[1:2,13]=data$fry_release_captive_brood_unclpd[c(nrow(data)-1):nrow(data)]
  #clipped smolts
  sim.dat[1:2,14]=data$smolt_release[c(nrow(data)-1):nrow(data)]
  
  #hatchery lake survival
  #fry survival to smolt emigration
  sim.dat[1:2,15]=data$fry_to_smolt_surv[c(nrow(data)-1):nrow(data)]/100
  #smolt survival to smolt emigration
  sim.dat[1:2,16]=data$smolt_release_emig_survival[c(nrow(data)-1):nrow(data)]/100
  
  ##Simulate - marine survival for clipped/unclipped smolts - spawners t+2
  #hatchery and unknown origin survival for these smolt cohorts
  
  #simulate for marine survival for each set of smolts
  #include autocorrelated previous years marine survival deviate for each state
  rho.ms.h.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.h,logistic.se.rho.ms.h))
  mu.adj.h=rho.ms.h.sim*(log(data$marine_surv_clipped[nrow(data)-2])-log(surv.mean)[1])
  mu.adj.u=rho.ms.h.sim*(log(data$marine_surv_unclipped[nrow(data)-2])-log(surv.mean)[2])
  
  sim.dat[1:2,9:11]=exp(MASS::mvrnorm(n=2,mu=mu.log.surv,Sigma=surv.mat))
  
  #estimate returning spawners at t+2 for each smolt cohort
  #unclipped spawners - unclipped outmigrating smolts x unknown origin marine survival
  sim.dat[3:4,3]=sim.dat[1:2,6]*(sim.dat[1:2,11]/100)
  #clipped spawners - clipped outmigrating smolts x hatchery origin marine survival
  sim.dat[3:4,4]=sim.dat[1:2,7]*(sim.dat[1:2,9]/100)
  #total spawners
  sim.dat[3:4,2]=sim.dat[3:4,3]+sim.dat[3:4,4]
  
  #need to add unclipped natural smolts for sim year 3 - these will be stored as 'smolts.n'
  
  #for all future states - next 25 years
  for(t in 1:c(sim.length+2)){
    
    sim.dat[t,2]=sum(sim.dat[t,3:5],na.rm=T) #total spawners
    
    #natural origin smolts from total spawners:
    
    #randomly sample Ricker parameters from estimates and covariances
    sim.ric=MASS::mvrnorm(n=1,mu=mu.ric,Sigma=vcov.ric)
    sim.ric[2:3]=exp(sim.ric[2:3]) #exponentiate from log-scale
    sim.ric[4]=inv_logistic_transf(sim.ric[4]) #function to scale transformed rho to -1 to 1
    #estimate natural origin smolts t+2 based on Ricker model
    if(t==1){
      #no autocorrelation if first smolt cohort (no previous estimate available to condition autocorrelation in smolt recruitment)
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2],sd=sim.ric[3])
    }else{
      #autocorrelation included
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2]+sim.ric[4]*sim.dat[t-1,19],sd=sim.ric[3])
    }
    #if total spawners < 1 = no natural smolt production
    if(sim.dat[t,2]<2){
      sim.dat[t+2,8]=0
    }else{
      sim.dat[t+2,8]=exp(logRS)*sim.dat[t,2]
    }
    sim.dat[t,17]=logRS
    sim.dat[t,18]=sim.ric[1]-sim.dat[t,2]/sim.ric[2]
    sim.dat[t,19]=logRS-(sim.ric[1]-sim.dat[t,2]/sim.ric[2])
    
     if(t>2){
      #fry releases - total fry releases x allocation to clipped/unclipped status
  #    total_fry.t=exp(rnorm(1,mu.log.fry,sd.log.fry)) #log-normal draws for fry releases
   #   total_fry.t=min(total_fry.t,2e6) #set capacity for fry release of 2 million
    #  frac.clipped=plogis(rnorm(1,mu.logit.fry.frac,sd.logit.fry.frac)) #proportion clipped/unclipped fry releases
     # sim.dat[t,12]=total_fry.t*(1-frac.clipped) #unclipped fry releases
      #sim.dat[t,13]=total_fry.t*frac.clipped #clipped fry releases
      #sim.dat[t,14]=exp(rnorm(1,mu.log.smolts,sd.log.smolts)) #log-normal draws for smolt releases
      
      #survival of unclipped/clipped fry, clipped smolts from lake to outmigration
    #  sim.dat[t,15]=plogis(rnorm(1,mu.logit.fry.surv,sd.logit.fry.surv)) #fry survival - applies to both clipped/unclipped
    #  sim.dat[t,16]=plogis(rnorm(1,mu.logit.smolt.surv,sd.logit.smolt.surv)) #smolt release survival to outmigration
      sim.dat[t,6]=0 #unclipped fry releases x fry survival for 3rd (1st) sim year
      sim.dat[t,7]=0 #clipped fry releases x (same) fry survival + clipped smolt releases x smolt emigration survival
      
      #expected returns (t+2) based on marine survival - natural origin & hatchery origin fish
      rho.ms.h.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.h,logistic.se.rho.ms.h)) #draw autocorrelation coefficient for hatchery origin smolts
      rho.ms.n.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.n,logistic.se.rho.ms.n)) #draw autocorrelation coefficient for natural origin smolts
      
      mu.adj=c(rho.ms.h.sim*(log(sim.dat[t-1,9])-log(surv.mean)[1]), #adjustment to mean survival based on 1-y autocorrelation
               0, #unknown origin marine survival not used beyond simulation year 2, set to 0 for convenience
               rho.ms.n.sim*(log(sim.dat[t-1,10])-log(surv.mean)[3])) #adjustment to mean survival based on 1-y autocorrelation
      
      sim.dat[t,9:11]=exp(MASS::mvrnorm(n=1,mu=mu.log.surv+mu.adj,Sigma=surv.mat)) #marine survival
      
      #partition based on age at maturity:
      p_age=extraDistr::rdirichlet(1,c(3,87,10))
      #3% age 3, 87% age 4, 10% age 5
      #proportion recruits, add to each return year based on randomized maturity
      #age 3
      sim.dat[t+1,3:4]=sim.dat[t+1,3:4]+sim.dat[t,6:7]*(sim.dat[t,9]/100)*p_age[1] #unclipped/clipped hatchery returns
      sim.dat[t+1,5]=sim.dat[t+1,5]+sim.dat[t,8]*(sim.dat[t,11]/100)*p_age[1] #natural origin returns
      #age 4
      sim.dat[t+2,3:4]=sim.dat[t+2,3:4]+sim.dat[t,6:7]*(sim.dat[t,9]/100)*p_age[2] #unclipped/clipped hatchery returns
      sim.dat[t+2,5]=sim.dat[t+2,5]+sim.dat[t,8]*(sim.dat[t,11]/100)*p_age[2] #natural origin returns
      #age 5
      sim.dat[t+3,3:4]=sim.dat[t+3,3:4]+sim.dat[t,6:7]*(sim.dat[t,9]/100)*p_age[3] #unclipped/clipped hatchery returns
      sim.dat[t+3,5]=sim.dat[t+3,5]+sim.dat[t,8]*(sim.dat[t,11]/100)*p_age[3] #natural origin returns
    }
  }
  sim.dat=sim.dat[1:c(sim.length+2),] #trim extra brood years
  sim.summary2[i,2]=sum(sim.dat[3:c(sim.length+2),2]>500)/sim.length #for simulated part of series - proportion of years where total spawners > 500
  sim.summary2[i,3]=sum(sim.dat[5:c(sim.length+2),5]>500)/c(sim.length-2) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary2[i,4]=sum(sim.dat[3:c(sim.length+2),2]>100)/sim.length #for simulated part of series - proportion of years where total spawners > 500
  sim.summary2[i,5]=sum(sim.dat[5:c(sim.length+2),5]>100)/c(sim.length-2) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary2[i,6]=exp(mean(log(sim.dat[3:c(sim.length+2),2])))
  sim.summary2[i,7]=exp(mean(log(sim.dat[5:c(sim.length+2),5])))
  
  sims2[[i]]=sim.dat #store simulation
}

#Sim 3: Brood-take####
sims3=list()
n.sims=500
sim.length=24
sim.summary3=data.frame(sim.id=seq(1:n.sims),py.tspwn.500=NA,py.nspwn.500=NA,py.tspwn.100=NA,py.nspwn.100=NA,m.tspwn=NA,m.nspwn=NA)
for(i in 1:n.sims){
  sim.dat=data.frame(sim.year=seq(1:c(sim.length+2)),spawners.t=0,spawners.u=0,spawners.c=0,spawners.n=c(rep(NA,3),rep(0,sim.length-1)),smolts.u=NA,smolts.c=NA,smolts.n=NA,ms.h=NA,ms.u=NA,ms.n=NA,fry.rel.u=NA,fry.rel.c=NA,smolt.rel.c=NA,fry.surv=NA,smolt.surv=NA,ricker.logrs=NA,ricker.mu=NA,ricker.eps=NA,broodtake=NA,spawners.r=NA,total.eggs=NA,egg.surv=NA)
  #first X years - estimate returns based on known quantities in the last several years
  
  #initial spawners - conditioned on observed data
  #total spawners
  sim.dat[1:2,2]=data$total_spawners[c(nrow(data)-1):nrow(data)]
  #unclipped spawners
  sim.dat[1:2,3]=data$unclipped_spawners[c(nrow(data)-1):nrow(data)]
  #clipped spawners
  sim.dat[1:2,4]=data$clipped_spawners[c(nrow(data)-1):nrow(data)]
  
  #initial smolt counts - conditioned on observed data
  #unclipped smolts (unc. fry and natural origin):
  sim.dat[1:2,6]=data$smolts_cnt_unclipped[c(nrow(data)-1):nrow(data)]
  #clipped smolts (both clp. fry and smolt releases):
  sim.dat[1:2,7]=data$smolts_cnt_clipped_fry[c(nrow(data)-1):nrow(data)]+data$smolt_cnt_release_clipped[c(nrow(data)-1):nrow(data)]
  
  #hatchery releases
  #unclipped fry
  sim.dat[1:2,12]=data$fry_release_captive_brood_unclpd[c(nrow(data)-1):nrow(data)]
  #clipped fry
  sim.dat[1:2,13]=data$fry_release_captive_brood_unclpd[c(nrow(data)-1):nrow(data)]
  #clipped smolts
  sim.dat[1:2,14]=data$smolt_release[c(nrow(data)-1):nrow(data)]
  
  #hatchery lake survival
  #fry survival to smolt emigration
  sim.dat[1:2,15]=data$fry_to_smolt_surv[c(nrow(data)-1):nrow(data)]/100
  #smolt survival to smolt emigration
  sim.dat[1:2,16]=data$smolt_release_emig_survival[c(nrow(data)-1):nrow(data)]/100
  
  ##Simulate - marine survival for clipped/unclipped smolts - spawners t+2
  #hatchery and unknown origin survival for these smolt cohorts
  
  #simulate for marine survival for each set of smolts
  #include autocorrelated previous years marine survival deviate for each state
  rho.ms.h.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.h,logistic.se.rho.ms.h))
  mu.adj.h=rho.ms.h.sim*(log(data$marine_surv_clipped[nrow(data)-2])-log(surv.mean)[1])
  mu.adj.u=rho.ms.h.sim*(log(data$marine_surv_unclipped[nrow(data)-2])-log(surv.mean)[2])
  
  sim.dat[1:2,9:11]=exp(MASS::mvrnorm(n=2,mu=mu.log.surv,Sigma=surv.mat))
  
  #estimate returning spawners at t+2 for each smolt cohort
  #unclipped spawners - unclipped outmigrating smolts x unknown origin marine survival
  sim.dat[3:4,3]=sim.dat[1:2,6]*(sim.dat[1:2,11]/100)
  #clipped spawners - clipped outmigrating smolts x hatchery origin marine survival
  sim.dat[3:4,4]=sim.dat[1:2,7]*(sim.dat[1:2,9]/100)
  #total spawners
  sim.dat[3:4,2]=sim.dat[3:4,3]+sim.dat[3:4,4]
  
  #need to add unclipped natural smolts for sim year 3 - these will be stored as 'smolts.n'
  
  #for all future states - next 25 years
  for(t in 1:c(sim.length+2)){
    
    sim.dat[t,2]=sum(sim.dat[t,3:5],na.rm=T) #total spawners
    #brood-take branch starting at t=3
    #remove up to half of total spawners for hatchery broodtake, up to maximum capacity
    if(t>2){
      sim.dat[t,21]=sim.dat[t,2] #save original spawners at wier before broodtake
      sim.dat[t,20]= min(0.5*sim.dat[t,2], 750) #broodtake - minimum of half of total spawners or max take 
      sim.dat[t,2]=sim.dat[t,21]-sim.dat[t,20] #remaining spawners for natural production
    }
    #natural origin smolts from total spawners:
    
    #randomly sample Ricker parameters from estimates and covariances
    sim.ric=MASS::mvrnorm(n=1,mu=mu.ric,Sigma=vcov.ric)
    sim.ric[2:3]=exp(sim.ric[2:3]) #exponentiate from log-scale
    sim.ric[4]=inv_logistic_transf(sim.ric[4]) #function to scale transformed rho to -1 to 1
    #estimate natural origin smolts t+2 based on Ricker model
    if(t==1){
      #no autocorrelation if first smolt cohort (no previous estimate available to condition autocorrelation in smolt recruitment)
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2],sd=sim.ric[3])
    }else{
      #autocorrelation included
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2]+sim.ric[4]*sim.dat[t-1,19],sd=sim.ric[3])
    }
    #if total spawners < 1 = no natural smolt production
    if(sim.dat[t,2]<2){
      sim.dat[t+2,8]=0
    }else{
      sim.dat[t+2,8]=exp(logRS)*sim.dat[t,2]
    }
    sim.dat[t,17]=logRS
    sim.dat[t,18]=sim.ric[1]-sim.dat[t,2]/sim.ric[2]
    sim.dat[t,19]=logRS-(sim.ric[1]-sim.dat[t,2]/sim.ric[2])
    
    if(t>2){
      #fry & smolt releases
      #broodtake - estimate eggs from broodtake and survival to fry
      # broodtake x fixed sex ratio of 50% x randomized eggs per female (mean = 2050, sd = 250)
      sim.dat[t,22]=sim.dat[t,20]*0.5*rnorm(1,2050,250) #total eggs
      sim.dat[t,23]=plogis(rnorm(1,1.26,0.5)) #egg-to-fry survival (mean = 78% logit(1.26), sd is guess-timated)
      
      #fry releases - total fry releases x allocation to clipped/unclipped status
      total_fry.t=sim.dat[t,22]*sim.dat[t,23] #total-fry is now total eggs from broodtake x egg survival
      total_fry.t=min(total_fry.t,2e6) #set capacity for fry release of 2 million
      frac.clipped=plogis(rnorm(1,mu.logit.fry.frac,sd.logit.fry.frac)) #proportion clipped/unclipped fry releases
      sim.dat[t,12]=total_fry.t*(1-frac.clipped) #unclipped fry releases
      sim.dat[t,13]=total_fry.t*frac.clipped #clipped fry releases
      sim.dat[t,14]=exp(rnorm(1,mu.log.smolts,sd.log.smolts)) #log-normal draws for smolt releases
      
      #survival of natural, unclipped/clipped fry, clipped smolts from lake to outmigration
      sim.dat[t,15]=plogis(rnorm(1,mu.logit.fry.surv,sd.logit.fry.surv)) #fry survival - applies to both clipped/unclipped
      sim.dat[t,16]=plogis(rnorm(1,mu.logit.smolt.surv,sd.logit.smolt.surv)) #smolt release survival to outmigration
      sim.dat[t,6]=sim.dat[t-1,12]*sim.dat[t,15] #unclipped fry releases x fry survival for 3rd (1st) sim year
      sim.dat[t,7]=sim.dat[t-1,13]*sim.dat[t,15]+sim.dat[t-1,14]*sim.dat[t,16] #clipped fry releases x (same) fry survival + clipped smolt releases x smolt emigration survival
      
      #expected returns (t+2) based on marine survival - natural origin & hatchery origin fish
      rho.ms.h.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.h,logistic.se.rho.ms.h)) #draw autocorrelation coefficient for hatchery origin smolts
      rho.ms.n.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.n,logistic.se.rho.ms.n)) #draw autocorrelation coefficient for natural origin smolts
      
      mu.adj=c(rho.ms.h.sim*(log(sim.dat[t-1,9])-log(surv.mean)[1]), #adjustment to mean survival based on 1-y autocorrelation
               0, #unknown origin marine survival not used beyond simulation year 2, set to 0 for convenience
               rho.ms.n.sim*(log(sim.dat[t-1,10])-log(surv.mean)[3])) #adjustment to mean survival based on 1-y autocorrelation
      
      sim.dat[t,9:11]=exp(MASS::mvrnorm(n=1,mu=mu.log.surv+mu.adj,Sigma=surv.mat)) #marine survival
      
      #partition based on age at maturity:
      p_age=extraDistr::rdirichlet(1,c(3,87,10))
      #3% age 3, 87% age 4, 10% age 5
      #proportion recruits, add to each return year based on randomized maturity
      #age 3
      sim.dat[t+1,3:4]=sim.dat[t+1,3:4]+sim.dat[t,6:7]*(sim.dat[t,9]/100)*p_age[1] #unclipped/clipped hatchery returns
      sim.dat[t+1,5]=sim.dat[t+1,5]+sim.dat[t,8]*(sim.dat[t,11]/100)*p_age[1] #natural origin returns
      #age 4
                   sim.dat[t+2,3:4]=sim.dat[t+2,3:4]+sim.dat[t,6:7]*(sim.dat[t,9]/100)*p_age[2] #unclipped/clipped hatchery returns
      sim.dat[t+2,5]=sim.dat[t+2,5]+sim.dat[t,8]*(sim.dat[t,11]/100)*p_age[2] #natural origin returns
      #age 5
      sim.dat[t+3,3:4]=sim.dat[t+3,3:4]+sim.dat[t,6:7]*(sim.dat[t,9]/100)*p_age[3] #unclipped/clipped hatchery returns
      sim.dat[t+3,5]=sim.dat[t+3,5]+sim.dat[t,8]*(sim.dat[t,11]/100)*p_age[3] #natural origin returns
    }
  }
  sim.dat=sim.dat[1:c(sim.length+2),] #trim extra brood years
  sim.summary3[i,2]=sum(sim.dat[3:c(sim.length+2),2]>500)/sim.length #for simulated part of series - proportion of years where total spawners > 500
  sim.summary3[i,3]=sum(sim.dat[5:c(sim.length+2),5]>500)/c(sim.length-2) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary3[i,4]=sum(sim.dat[3:c(sim.length+2),2]>100)/sim.length #for simulated part of series - proportion of years where total spawners > 500
  sim.summary3[i,5]=sum(sim.dat[5:c(sim.length+2),5]>100)/c(sim.length-2) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary3[i,6]=exp(mean(log(sim.dat[3:c(sim.length+2),2])))
  sim.summary3[i,7]=exp(mean(log(sim.dat[5:c(sim.length+2),5])))
  
  sims3[[i]]=sim.dat #store simulation
}

#Sim 4: Doubled Hatchery Output####
#adjust total mean fry

sd.log.fry=sqrt(log(1+(sd.fry/mu.fry)^2))
mu.log.fry=mean(log(total.fry*2))
exp(mu.log.fry)

sims4=list()
n.sims=500
sim.length=24
sim.summary4=data.frame(sim.id=seq(1:n.sims),py.tspwn.500=NA,py.nspwn.500=NA,py.tspwn.100=NA,py.nspwn.100=NA,m.tspwn=NA,m.nspwn=NA)
for(i in 1:n.sims){
  sim.dat=data.frame(sim.year=seq(1:c(sim.length+2)),spawners.t=0,spawners.u=0,spawners.c=0,spawners.n=c(rep(NA,3),rep(0,sim.length-1)),smolts.u=NA,smolts.c=NA,smolts.n=NA,ms.h=NA,ms.u=NA,ms.n=NA,fry.rel.u=NA,fry.rel.c=NA,smolt.rel.c=NA,fry.surv=NA,smolt.surv=NA,ricker.logrs=NA,ricker.mu=NA,ricker.eps=NA)
  #first X years - estimate returns based on known quantities in the last several years
  
  #initial spawners - conditioned on observed data
  #total spawners
  sim.dat[1:2,2]=data$total_spawners[c(nrow(data)-1):nrow(data)]
  #unclipped spawners
  sim.dat[1:2,3]=data$unclipped_spawners[c(nrow(data)-1):nrow(data)]
  #clipped spawners
  sim.dat[1:2,4]=data$clipped_spawners[c(nrow(data)-1):nrow(data)]
  
  #initial smolt counts - conditioned on observed data
  #unclipped smolts (unc. fry and natural origin):
  sim.dat[1:2,6]=data$smolts_cnt_unclipped[c(nrow(data)-1):nrow(data)]
  #clipped smolts (both clp. fry and smolt releases):
  sim.dat[1:2,7]=data$smolts_cnt_clipped_fry[c(nrow(data)-1):nrow(data)]+data$smolt_cnt_release_clipped[c(nrow(data)-1):nrow(data)]
  
  #hatchery releases
  #unclipped fry
  sim.dat[1:2,12]=data$fry_release_captive_brood_unclpd[c(nrow(data)-1):nrow(data)]
  #clipped fry
  sim.dat[1:2,13]=data$fry_release_captive_brood_unclpd[c(nrow(data)-1):nrow(data)]
  #clipped smolts
  sim.dat[1:2,14]=data$smolt_release[c(nrow(data)-1):nrow(data)]
  
  #hatchery lake survival
  #fry survival to smolt emigration
  sim.dat[1:2,15]=data$fry_to_smolt_surv[c(nrow(data)-1):nrow(data)]/100
  #smolt survival to smolt emigration
  sim.dat[1:2,16]=data$smolt_release_emig_survival[c(nrow(data)-1):nrow(data)]/100
  
  ##Simulate - marine survival for clipped/unclipped smolts - spawners t+2
  #hatchery and unknown origin survival for these smolt cohorts
  
  #simulate for marine survival for each set of smolts
  #include autocorrelated previous years marine survival deviate for each state
  rho.ms.h.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.h,logistic.se.rho.ms.h))
  mu.adj.h=rho.ms.h.sim*(log(data$marine_surv_clipped[nrow(data)-2])-log(surv.mean)[1])
  mu.adj.u=rho.ms.h.sim*(log(data$marine_surv_unclipped[nrow(data)-2])-log(surv.mean)[2])
  
  sim.dat[1:2,9:11]=exp(MASS::mvrnorm(n=2,mu=mu.log.surv,Sigma=surv.mat))
  
  #estimate returning spawners at t+2 for each smolt cohort
  #unclipped spawners - unclipped outmigrating smolts x unknown origin marine survival
  sim.dat[3:4,3]=sim.dat[1:2,6]*(sim.dat[1:2,11]/100)
  #clipped spawners - clipped outmigrating smolts x hatchery origin marine survival
  sim.dat[3:4,4]=sim.dat[1:2,7]*(sim.dat[1:2,9]/100)
  #total spawners
  sim.dat[3:4,2]=sim.dat[3:4,3]+sim.dat[3:4,4]
  
  #need to add unclipped natural smolts for sim year 3 - these will be stored as 'smolts.n'
  
  #for all future states - next 25 years
  for(t in 1:c(sim.length+2)){
    
    sim.dat[t,2]=sum(sim.dat[t,3:5],na.rm=T) #total spawners
    
    #natural origin smolts from total spawners:
    
    #randomly sample Ricker parameters from estimates and covariances
    sim.ric=MASS::mvrnorm(n=1,mu=mu.ric,Sigma=vcov.ric)
    sim.ric[2:3]=exp(sim.ric[2:3]) #exponentiate from log-scale
    sim.ric[4]=inv_logistic_transf(sim.ric[4]) #function to scale transformed rho to -1 to 1
    #estimate natural origin smolts t+2 based on Ricker model
    if(t==1){
      #no autocorrelation if first smolt cohort (no previous estimate available to condition autocorrelation in smolt recruitment)
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2],sd=sim.ric[3])
    }else{
      #autocorrelation included
      logRS=rnorm(1,mean=sim.ric[1]-sim.dat[t,2]/sim.ric[2]+sim.ric[4]*sim.dat[t-1,19],sd=sim.ric[3])
    }
    #if total spawners < 1 = no natural smolt production
    if(sim.dat[t,2]<2){
      sim.dat[t+2,8]=0
    }else{
      sim.dat[t+2,8]=exp(logRS)*sim.dat[t,2]
    }
    sim.dat[t,17]=logRS
    sim.dat[t,18]=sim.ric[1]-sim.dat[t,2]/sim.ric[2]
    sim.dat[t,19]=logRS-(sim.ric[1]-sim.dat[t,2]/sim.ric[2])
    
    #fry & smolt releases
    if(t>2){
      #fry releases - total fry releases x allocation to clipped/unclipped status
      total_fry.t=exp(rnorm(1,mu.log.fry,sd.log.fry)) #log-normal draws for fry releases
      total_fry.t=min(total_fry.t,4e6) #set capacity for fry release of 2 million
      frac.clipped=plogis(rnorm(1,mu.logit.fry.frac,sd.logit.fry.frac)) #proportion clipped/unclipped fry releases
      sim.dat[t,12]=total_fry.t*(1-frac.clipped) #unclipped fry releases
      sim.dat[t,13]=total_fry.t*frac.clipped #clipped fry releases
      sim.dat[t,14]=exp(rnorm(1,mu.log.smolts,sd.log.smolts)) #log-normal draws for smolt releases
      
      #survival of natural, unclipped/clipped fry, clipped smolts from lake to outmigration
      sim.dat[t,15]=plogis(rnorm(1,mu.logit.fry.surv,sd.logit.fry.surv)) #fry survival - applies to both clipped/unclipped
      sim.dat[t,16]=plogis(rnorm(1,mu.logit.smolt.surv,sd.logit.smolt.surv)) #smolt release survival to outmigration
      sim.dat[t,6]=sim.dat[t-1,12]*sim.dat[t,15] #unclipped fry releases x fry survival for 3rd (1st) sim year
      sim.dat[t,7]=sim.dat[t-1,13]*sim.dat[t,15]+sim.dat[t-1,14]*sim.dat[t,16] #clipped fry releases x (same) fry survival + clipped smolt releases x smolt emigration survival
      
      #expected returns (t+2) based on marine survival - natural origin & hatchery origin fish
      rho.ms.h.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.h,logistic.se.rho.ms.h)) #draw autocorrelation coefficient for hatchery origin smolts
      rho.ms.n.sim=inv_logistic_transf(rnorm(1,logistic.rho.ms.n,logistic.se.rho.ms.n)) #draw autocorrelation coefficient for natural origin smolts
      
      mu.adj=c(rho.ms.h.sim*(log(sim.dat[t-1,9])-log(surv.mean)[1]), #adjustment to mean survival based on 1-y autocorrelation
               0, #unknown origin marine survival not used beyond simulation year 2, set to 0 for convenience
               rho.ms.n.sim*(log(sim.dat[t-1,10])-log(surv.mean)[3])) #adjustment to mean survival based on 1-y autocorrelation
      
      sim.dat[t,9:11]=exp(MASS::mvrnorm(n=1,mu=mu.log.surv+mu.adj,Sigma=surv.mat)) #marine survival
      
      #partition based on age at maturity:
      p_age=extraDistr::rdirichlet(1,c(3,87,10))
      #3% age 3, 87% age 4, 10% age 5
      #proportion recruits, add to each return year based on randomized maturity
      #age 3
      sim.dat[t+1,3:4]=sim.dat[t+1,3:4]+sim.dat[t,6:7]*(sim.dat[t,9]/100)*p_age[1] #unclipped/clipped hatchery returns
      sim.dat[t+1,5]=sim.dat[t+1,5]+sim.dat[t,8]*(sim.dat[t,11]/100)*p_age[1] #natural origin returns
      #age 4
      sim.dat[t+2,3:4]=sim.dat[t+2,3:4]+sim.dat[t,6:7]*(sim.dat[t,9]/100)*p_age[2] #unclipped/clipped hatchery returns
      sim.dat[t+2,5]=sim.dat[t+2,5]+sim.dat[t,8]*(sim.dat[t,11]/100)*p_age[2] #natural origin returns
      #age 5
      sim.dat[t+3,3:4]=sim.dat[t+3,3:4]+sim.dat[t,6:7]*(sim.dat[t,9]/100)*p_age[3] #unclipped/clipped hatchery returns
      sim.dat[t+3,5]=sim.dat[t+3,5]+sim.dat[t,8]*(sim.dat[t,11]/100)*p_age[3] #natural origin returns
    }
  }
  sim.dat=sim.dat[1:c(sim.length+2),] #trim extra brood years
  sim.summary4[i,2]=sum(sim.dat[3:c(sim.length+2),2]>500)/sim.length #for simulated part of series - proportion of years where total spawners > 500
  sim.summary4[i,3]=sum(sim.dat[5:c(sim.length+2),5]>500)/c(sim.length-2) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary4[i,4]=sum(sim.dat[3:c(sim.length+2),2]>100)/sim.length #for simulated part of series - proportion of years where total spawners > 500
  sim.summary4[i,5]=sum(sim.dat[5:c(sim.length+2),5]>100)/c(sim.length-2) #for simulated part of series - proportion of years where natural origin spawners > 500
  sim.summary4[i,6]=exp(mean(log(sim.dat[3:c(sim.length+2),2])))
  sim.summary4[i,7]=exp(mean(log(sim.dat[5:c(sim.length+2),5])))
  
  sims4[[i]]=sim.dat #store simulation
}

#summary of each scenario
summary(sim.summary1) #baseline
summary(sim.summary2) #all natural
summary(sim.summary3) #broodtake
summary(sim.summary4) #double hatchery

