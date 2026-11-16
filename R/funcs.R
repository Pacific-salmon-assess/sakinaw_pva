inv_logistic_transf <- function(x) {
  return(2 * (exp(x)/(1+exp(x))) - 1)
}

logistic_transf<- function(x) {
  log((1+x)/(1-x))
}

rollmean_4yr <- function(x) {
  stats::filter(x, rep(1/4, 4), sides = 1)
}


cbp_sim=function(mean.fry,mean.smolt,sd.fry,sd.smolt,n.sims,length){
  
  mu.log.fry=log(mean.fry)
  sd.log.fry
  
  mu.log.smolts=log(mean.smolt)
  sd.log.smolts=
  
  
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
  
  return(sim.summary1,sims1)
}