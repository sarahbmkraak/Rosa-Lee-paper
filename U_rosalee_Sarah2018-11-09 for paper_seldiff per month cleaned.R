##++++++++++++++++++++++++++++++++
## Preliminary selectivity modelling
## based on CM, SK: 2012
## many modifications 2017 and 2018 CM, SK, JS, and SH 
## parametrisation based on western baltic cod SH
##++++++++++++++++++++++++++++++++++++
rm(list=ls())


library(plyr)
library(truncnorm)
library(ggplot2)
library(readr)
library(mvtnorm)


# Set working directory ---------------------------------------------------

##working directory Sarah
  setwd("C:/Users/kraak/Dropbox/Master Thesis/R_rosalee")
  source("C:/Users/kraak/Dropbox/Master Thesis/R_rosalee/get_det_growth_mat_betamon_cod_20190108.R")
  source("C:/Users/kraak/Dropbox/Master Thesis/R_rosalee/get_survivorship_mat_cod_Sarah.R")
  source("C:/Users/kraak/Dropbox/Master Thesis/R_rosalee/get_catch_mat_cod_Sarah.R")
  source("C:/Users/kraak/Dropbox/Master Thesis/R_rosalee/get_length_dist_age_betamon_cod_Sarah.R")
  source("C:/Users/kraak/Dropbox/Master Thesis/R_rosalee/get_weight_cod_using1seasononly.R")
  source("C:/Users/kraak/Dropbox/Master Thesis/R_rosalee/get_maturity_cod.R")
  source("C:/Users/kraak/Dropbox/Master Thesis/R_rosalee/get_n_c_length_betamon_cod_Sarah.R")
  source("C:/Users/kraak/Dropbox/Master Thesis/R_rosalee/get_B_cod.R")
  source("C:/Users/kraak/Dropbox/Master Thesis/R_rosalee/get_SSB_cod.R")
  source("C:/Users/kraak/Dropbox/Master Thesis/R_rosalee/get_Cw_cod.R")
  source("C:/Users/kraak/Dropbox/Master Thesis/R_rosalee/get_Lw_cod_Sarah.R")
  source("C:/Users/kraak/Dropbox/Master Thesis/R_rosalee/get_Dw_cod_Sarah.R")


##-------
## SETUP
##-------


## These are partly from Kate's tagging study (McQueen et al. 2018) and partly from Steffie's minimization. We can play around with rho as an experimental treatment
linf.val <- 154.56 #Kate
sd.linf <- 16.75 #Kate
k.val <- 0.11 #Kate
sd.k <- 0.025 #we chose it arbitrarily based on visual inspection of 100 growth curves
t0.val <- -0.13 #Kate
sd.l.age <- 0.18 #we chose it arbitrarily; only used in backcalculating recruits-at-length giving rise to N-at-age in year 1; burn-in makes actual scenarios independent of starting N
rho <- 0.7 #we chose it arbitrarily; results did not differ much with values of 0.99, -0.99 and 0.01

### set number of ages and years
## We don't have a plusgroup so we need many ages, also to accommodate scenarions with low F where fish survive to old age; let's take 20.
n.ages<-20 
age.vec<-1:n.ages

## number of years to project
## 'burn-in' to reach equilibrium and independence of the starting numbers:
## number of years for burn-in =n.ages; then after 1 year the gear or F may change according to scenarios
## then 10 years for predictions; this 10 can be changed if desired.
n.years<-n.ages+11

## numbers at age in first year (because of burn-in, results will be independent of these)
##SH200618: Estimated stock numbers from year 2002 from the WGBFAS Report 2017 
N.age<-c(40135, 17445, 11119, 2043, 1542, 342, 35,rep(0.1,n.ages-7))


## number of seasons
## if per quarter or month, set to 4 or 12
# we tried with weeks (52 seasons)
n.seasons <- 12

# ## Selectivity parameters
L50.1 <- 30
SR <- 6
L50.2 <- 40 #this is the new L50 in the default gear-change scenario; it can be set to 20 or stay at 30



b1 <- 2.197/SR
a1 <- -L50.1*b1  
a2 <- -L50.2*b1  
b2 <- b1 # in scenarios where SR is changed to 3 or 9, this is set to 2.197/3 or 2.197/9 respectively

## Fully selected fishing mortality
Ffull <- 1  
## the fully selected F for the F-change scenario
Fdeclined <- 0.75*Ffull # this can be set to 0.9, 0.5 or 1.25 respectively for the different F-change scenarios

## Natural mortality at age
M.age<-c(0.242,rep(0.2,n.ages-1)) 

## length-class break-points
## with many length classes we set by=1, else by=0.5
l.breaks<-seq(2,200, by=1)

## get length class mid-points
get.mids<-function(x){
  (x[2:length(x)]+x[1:(length(x)-1)])/2
}
l.mids<-get.mids(l.breaks)
n.bins<-length(l.mids)

## F - used to generate the numbers at length for first year (took the T90 120 mm selection pattern and N.standard.a from 1st year) 
## WHERE DO THESE COME FROM?
##SH200618: need to check this. Is it used later anyway?
### SK020718 they are used for the back calculation to get the length distributions at age in year 1. But owing to burn-in the simulation results are independent of this
F.age<- c(0.103, 0.637, 1.299, 1.299, 1.299, 1.299, 1.299, rep(1.299, n.ages-7))#c(0.11, 0.368, 1.354, 1.354, 1.354, 1.354, 1.354)

## The following lines have to be used once, to create G; after that we can use the G that has been stored.
 #B number of pseudo-VB curves used in the montecarlo simulation
  G<-llply(1:n.ages,function(z){
   get.det.growth.mat.mon(age=age.vec[z],n.seasons,l.breaks,l.mids,n.bins, linf=linf.val,sd.linf,k=k.val,t0=t0.val,B=10000, sd.k,rho)
                 })
  #save(G, file="G_rosa_a20l200s52.RData")
# one should not just load any G that is saved from previous runs, because it may have different numbers of seasons, ages or length bins
#load("C:/Users/kraak/Dropbox/Master Thesis/R_rosalee/G_rosa_a20l500.RData")

## Main calculations


N.C.statusquo <- get.N.C.length(N.age=N.age, age.vec=age.vec, l.breaks=l.breaks,G=G, l.mids=l.mids,linf.val=linf.val, k.val=k.val, t0.val=t0.val, sd.l.age=sd.l.age, n.years=n.years, n.seasons=n.seasons, a1=a1, a2=a1, b1=b1, b2=b1, M.age=M.age, F.age=F.age, Ffull=Ffull, scenario="statusquo", Fbar,Fdeclined) 
N.C.gearchange <- get.N.C.length(N.age=N.age, age.vec=age.vec, l.breaks=l.breaks,G=G, l.mids=l.mids, linf.val=linf.val, k.val=k.val, t0.val=t0.val, sd.l.age=sd.l.age, n.years=n.years, n.seasons=n.seasons, a1=a1, a2=a2, b1=b1, b2=b2, M.age=M.age, F.age=F.age, Ffull=Ffull, scenario="gearchange", Fbar,Fdeclined) 
N.C.Fchange <- get.N.C.length(N.age=N.age, age.vec=age.vec, l.breaks=l.breaks,G=G, l.mids=l.mids, linf.val=linf.val, k.val=k.val, t0.val=t0.val, sd.l.age=sd.l.age, n.years=n.years, n.seasons=n.seasons, a1=a1, a2=a2, b1=b1, b2=b2, M.age=M.age, F.age=F.age, Ffull=Ffull, scenario="Fchange", Fbar,Fdeclined) 

##---------------------------
## STATUS QUO SCENARIO EXTRACTIONS
##---------------------------
N.statusquo<-N.C.statusquo[["N"]]
C.statusquo<-N.C.statusquo[["C"]]

#Matrix with the numbers of fish for every prediction year in the first season
Numbers.statusquo <-matrix(0,n.years, n.ages)
for(i in 1:n.years){
  Numbers.statusquo[i,1:n.ages] <- apply(N.statusquo[,,1,i],2 ,sum)
}

#Matrix with the numbers of fish caught for every prediction year
Catch.statusquo <-matrix(0,n.years, n.ages)
for(i in 1:n.years){
  Catch.statusquo[i,1:n.ages] <- apply(C.statusquo[,,,i],2 ,sum)
}


## Total Biomass at the end of the last season each year
B.statusquo<-get.B(N.array=N.statusquo,n.years,n.ages,n.seasons, l.mids=l.mids)
B.statusquo.sum <- apply(B.statusquo[,,n.seasons,], 3, sum)

## SSB at the end of the last season each year
SSB.statusquo<-get.SSB(B.array=B.statusquo,n.years,n.ages,n.seasons, l.mids=l.mids)
SSB.statusquo.sum<-apply(SSB.statusquo[,,n.seasons,],3,sum)

## Total catch in weight each year
Cw.statusquo<-get.Cw(C.array=C.statusquo,n.years,n.ages,n.seasons, l.mids=l.mids)
Cw.statusquo.sum<-apply(Cw.statusquo[,,,],4,sum)

## Landings above MCRS each year
Lw.statusquo<-get.Lw(Cw.array=Cw.statusquo,n.years,n.ages,n.seasons, l.mids=l.mids)
Lw.statusquo.sum<-apply(Lw.statusquo[,,,],4,sum)
#LW.statusquo.a.sum_age <-apply(Lw.statusquo,2,sum)

## Catches below MCRS each year
Dw.statusquo<-get.Dw(Cw.array=Cw.statusquo,n.years,n.ages,n.seasons, l.mids=l.mids)
Dw.statusquo.sum<-apply(Dw.statusquo[,,,],4,sum)


##---------------------------
## GEAR CHANGE EXTRACTIONS
##---------------------------
N.gearchange<-N.C.gearchange[["N"]]
C.gearchange<-N.C.gearchange[["C"]]

Numbers.gearchange <-matrix(0,n.years, n.ages)
for(i in 1:n.years){
  Numbers.gearchange[i,1:n.ages] <- apply(N.gearchange[,,1,i],2 ,sum)
}

#Matrix with the numbers of fish caught for every prediction year
Catch.gearchange <-matrix(0,n.years, n.ages)
for(i in 1:n.years){
#  Catch.gearchange[i,1] <- i
  #Catch.gearchange[i,1:7] <- apply(C.gearchange[,,1,i],2 ,sum)
  Catch.gearchange[i,1:n.ages] <- apply(C.gearchange[,,,i],2 ,sum)
}
## Total Biomass
B.gearchange<-get.B(N.array=N.gearchange,n.years,n.ages,n.seasons, l.mids=l.mids)
B.gearchange.sum <- apply(B.gearchange[,,n.seasons,],3, sum)

## SSB
SSB.gearchange<-get.SSB(B.array=B.gearchange,n.years,n.ages,n.seasons, l.mids=l.mids)
SSB.gearchange.sum<-apply(SSB.gearchange[,,n.seasons,],3,sum)

## Total catch in weight
Cw.gearchange<-get.Cw(C.array=C.gearchange,n.years,n.ages,n.seasons, l.mids=l.mids)
Cw.gearchange.sum<-apply(Cw.gearchange[,,,],4,sum)

## Landings above MCRS
Lw.gearchange<-get.Lw(Cw.array=Cw.gearchange,n.years,n.ages,n.seasons, l.mids=l.mids)
Lw.gearchange.sum<-apply(Lw.gearchange[,,,],4,sum)

## Catches below MCRS
Dw.gearchange<-get.Dw(Cw.array=Cw.gearchange,n.years,n.ages,n.seasons, l.mids=l.mids)
Dw.gearchange.sum<-apply(Dw.gearchange[,,,],4,sum)


##---------------------------
## F CHANGE SCENARIO EXTRACTIONS
##---------------------------
N.Fchange<-N.C.Fchange[["N"]]
C.Fchange<-N.C.Fchange[["C"]]

#Matrix with the numbers of fish for every prediction year
#attention if i have more than 1 season because then it will ad up and I will have x4 or x12 the population size
Numbers.Fchange <-matrix(0,n.years, n.ages)
for(i in 1:n.years){
  Numbers.Fchange[i,1:n.ages] <- apply(N.Fchange[,,1,i],2 ,sum)
}

#Matrix with the numbers of fish caught for every prediction year
Catch.Fchange <-matrix(0,n.years, n.ages)
for(i in 1:n.years){
  #Catch.Fchange[i,1:7] <- apply(C.Fchange[,,2,i],2 ,sum)
  Catch.Fchange[i,1:n.ages] <- apply(C.Fchange[,,,i],2 ,sum)
}


## Total Biomass
B.Fchange<-get.B(N.array=N.Fchange,n.years,n.ages,n.seasons, l.mids=l.mids)
B.Fchange.sum <- apply(B.Fchange[,,n.seasons,], 3, sum)

## SSB
SSB.Fchange<-get.SSB(B.array=B.Fchange,n.years,n.ages,n.seasons, l.mids=l.mids)
SSB.Fchange.sum<-apply(SSB.Fchange[,,n.seasons,],3,sum)

## Total catch in weight
Cw.Fchange<-get.Cw(C.array=C.Fchange,n.years,n.ages,n.seasons, l.mids=l.mids)
Cw.Fchange.sum<-apply(Cw.Fchange[,,,],4,sum)

## Landings above MCRS
Lw.Fchange<-get.Lw(Cw.array=Cw.Fchange,n.years,n.ages,n.seasons, l.mids=l.mids)
Lw.Fchange.sum<-apply(Lw.Fchange[,,,],4,sum)
#LW.Fchange.a.sum_age <-apply(Lw.Fchange,2,sum)

## Catches below MCRS
Dw.Fchange<-get.Dw(Cw.array=Cw.Fchange,n.years,n.ages,n.seasons, l.mids=l.mids)
Dw.Fchange.sum<-apply(Dw.Fchange[,,,],4,sum)

##---------------------------
## tracking changes in mean lenhth@age over time
##---------------------------

#tracking changes in mean lenhth@age over time
mean_length_at_age<-array(0,c(n.ages,11,2))
for (j in 1:n.ages){
  for (i in 1:11){
    # for gear change
    mean_length_at_age[j,i,1]<-sum(l.mids*(N.gearchange[,j,n.seasons,n.ages+i]))/sum(N.gearchange[,j,n.seasons,n.ages+i])  
  }
}
for (j in 1:n.ages){
  for (i in 1:11){
    # for F change
    mean_length_at_age[j,i,2]<-sum(l.mids*(N.Fchange[,j,n.seasons,n.ages+i]))/sum(N.Fchange[,j,n.seasons,n.ages+i])  
  }
}
plot(mean_length_at_age[1,,1],type="n", ylim=c(30,140))
abline(h=seq(0,200,by=10), col="grey")
for (j in 1:15){
  lines(mean_length_at_age[j,,1])
}
mtext("mean length at age after gear change")
plot(mean_length_at_age[1,,2],type="n", ylim=c(30,140))
abline(h=seq(0,200,by=10), col="grey")
for (j in 1:15){
  lines(mean_length_at_age[j,,2])
}
mtext("mean length at age after F change")


##----------------------

# CALCULATIONS FOR THE STATIC METHOD --------------------------------------

##----------------------


## projections
start.year <- 2002
pred.years<-start.year:(start.year+(n.years-1))

## Mean length of each age by season
mean_length_at_age_at_season <- sel1<- sel2 <- sel.diff <- array(0,c(n.ages,n.seasons))
  for (j in 1:n.ages){
    for (k in 1:n.seasons){
    mean_length_at_age_at_season[j,k]<- sum(l.mids*(N.statusquo[,j,k,n.ages+1]))/sum(N.statusquo[,j,k,n.ages+1])
    sel1[j,k] <- plogis(a1+b1*(mean_length_at_age_at_season[j,k]))
    sel2[j,k] <- plogis(a2+b2*(mean_length_at_age_at_season[j,k]))
    sel.diff[j,k] <- sel2[j,k]/sel1[j,k]
    }
 }
# this is how people with the static method calculate the new selectivity; they first calculate the ratio at age between the new and the old


# Estimation of F-at-age for the static method -------------------------------------------------

F_jk<- matrix(numeric(n.ages*n.seasons),nrow=n.ages) 
# k is for season
# j is for age
for (j in 1:n.ages){
  for(k in 1:n.seasons){

    if (k < n.seasons){
      F_jk[j,k] <- log(sum(N.statusquo[,j,k,n.ages+1])/sum(N.statusquo[,j,k+1,n.ages+1]))-M.age[j]/n.seasons  #season 1 to 11
    }
    if (k == n.seasons & j< n.ages){
                          F_jk[j,k] <- log(sum(N.statusquo[,j,k,n.ages+1])/sum(N.statusquo[,j+1,1,n.ages+2]))-M.age[j]/n.seasons  #season 12
    }
    
    
    if (j == n.ages & k == n.seasons) {
      #F for last season of last age = previous season.
      F_jk[j,k] <- F_jk[j,n.seasons-1]
      }
    
  }
}


N.at.age.jki <- array(0, c(n.ages,n.seasons,n.years))
C.at.age.jki <- array(0, c(n.ages,n.seasons,n.years))
# i is for year
# the first n.ages+1 years are before the measure change; they are the same between static and dynamic method
# we use a burn-in period to get an equilibrium population; the number of years for that is =n.ages; year=n.ages+1 is the first year after burn-in
for( i in 1:(n.ages+1)){
  
  for (k in 1:n.seasons){
  
    N.at.age.jki[,k,i]<- apply(N.statusquo[,,k,i],2,sum)
    
  }
  
}

# year n.ages+1 (just after "burn in" period) biomasses per age and season to calculate ogives for static method and mean weight@age@season


B_jk<- matrix(numeric(n.ages*n.seasons),nrow=n.ages) 

for (k in 1:n.seasons){
  B_jk[,k] <-apply(B.statusquo[,,k,n.ages+1],2,sum)
  
}


SSB_jk<- matrix(numeric(n.ages*n.seasons),nrow=n.ages) 

for (k in 1:n.seasons){
  SSB_jk[,k] <-apply(SSB.statusquo[,,k,n.ages+1],2,sum)
  
}

Cw_jk<- matrix(numeric(n.ages*n.seasons),nrow=n.ages) 

for (k in 1:n.seasons){
  Cw_jk[,k] <-apply(Cw.statusquo[,,k,n.ages+1],2,sum)
  
}

Lw_jk<- matrix(numeric(n.ages*n.seasons),nrow=n.ages) 

for (k in 1:n.seasons){
  Lw_jk[,k] <-apply(Lw.statusquo[,,k,n.ages+1],2,sum)
  
}

Dw_jk<- matrix(numeric(n.ages*n.seasons),nrow=n.ages) 

for (k in 1:n.seasons){
  Dw_jk[,k] <-apply(Dw.statusquo[,,k,n.ages+1],2,sum)
  
}

dis_prop <-Dw_jk/Cw_jk
mat_prop <-SSB_jk/B_jk

# mean weight in the stock
mean_weight_jk <- matrix(numeric(n.ages*n.seasons),nrow=n.ages)

for (k in 1:n.seasons){
  
  mean_weight_jk[,k] <- apply(B.statusquo[,,k,n.ages+1],2,sum)/apply(N.statusquo[,,k,n.ages+1],2,sum)
  
}

# mean weight in the catch
mean_cweight_jk <- matrix(numeric(n.ages*n.seasons),nrow=n.ages)

for (k in 1:n.seasons){
  
  mean_cweight_jk[,k] <- apply(Cw.statusquo[,,k,n.ages+1],2,sum)/apply(C.statusquo[,,k,n.ages+1],2,sum)
  
}


##----------------------
## STATUS QUO static projection after year n.ages+1 (years n.ages+2 to n.ages+11) for status quo (=no change)
##----------------------

pred.years <- 1:10
Ffull<-1

for (i in pred.years){

   N.at.age.jki[1,1,i+(n.ages+1)] <- 14000  #recruitment

   for (k in 1:n.seasons){

    if (k != 1){

                N.at.age.jki[1,k,i+(n.ages+1)] <- exp(-(F_jk[1,k-1]*Ffull+M.age[1]/n.seasons))*N.at.age.jki[1,k-1,i+(n.ages+1)]  #season 2 to 12, age 1
    }

      for(j in 2:n.ages){


            if (k == 1){
              N.at.age.jki[j,1,i+(n.ages+1)] <- exp(-(F_jk[j-1,n.seasons]*Ffull+M.age[j-1]/n.seasons))*N.at.age.jki[j-1,n.seasons,i+n.ages]  #season 1
          }else{
                N.at.age.jki[j,k,i+(n.ages+1)] <- exp(-(F_jk[j,k-1]*Ffull+M.age[j]/n.seasons))*N.at.age.jki[j,k-1,i+(n.ages+1)]  #season 2 to 12
          }

      }
     for( j in 1: n.ages){

     
     C.at.age.jki[j,k,i+(n.ages+1)] <- F_jk[j,k]*Ffull/(F_jk[j,k]*Ffull+M.age[j]/n.seasons)*(1-exp(-(F_jk[j,k]*Ffull+M.age[j]/n.seasons)))*N.at.age.jki[j,k,i+(n.ages+1)]
   }
}
}

landings<- numeric(10)
discards<- numeric(10)
catches<-numeric(10)
SSB<-numeric(10)

for (i in 1:10){
  SSB[i] <- sum(N.at.age.jki[,n.seasons,(n.ages+1)+i]*mean_weight_jk[,n.seasons]*mat_prop[,n.seasons]) 
  landings[i] <- sum(C.at.age.jki[,,(n.ages+1)+i]*mean_cweight_jk*(1-dis_prop))
  discards [i]<- sum(C.at.age.jki[,,(n.ages+1)+i]*mean_cweight_jk*dis_prop)
  catches[i]<- sum(C.at.age.jki[,,(n.ages+1)+i]*mean_cweight_jk)
}


landings<- c(rep(NA, n.ages), Lw.statusquo.sum[(n.ages+1)], landings)
discards<- c(rep(NA, n.ages), Dw.statusquo.sum[(n.ages+1)], discards)
catches<- c(rep(NA, n.ages), Cw.statusquo.sum[(n.ages+1)], catches)
SSB<- c(rep(NA, n.ages), SSB.statusquo.sum[(n.ages+1)], SSB)

##----------------------
## GEAR CHANGE: static projection after year n.ages+1 (years n.ages+2 to n.ages+11) for the scenario with a gear change
##----------------------
 
pred.years<-1:10
N.at.age.jki.gear <- array(0, c(n.ages,n.seasons,n.years))
C.at.age.jki.gear <- array(0, c(n.ages,n.seasons,n.years))
Ffull <- 1

for( i in 1:(n.ages+1)){
  
  for (k in 1:n.seasons){
    
    N.at.age.jki.gear[,k,i]<- apply(N.statusquo[,,k,i],2,sum)
    
  }
  
}

for (i in pred.years){
  
  N.at.age.jki.gear[1,1,i+(n.ages+1)] <- 14000  #recruitment
  
  for (k in 1:n.seasons){
    
    if (k != 1){
      N.at.age.jki.gear[1,k,i+(n.ages+1)] <- exp(-(F_jk[1,k-1]*Ffull*sel.diff[1,k-1]+M.age[1]/n.seasons))*N.at.age.jki.gear[1,k-1,i+(n.ages+1)]  #season 2 to 12, age 1
    }
    
    for(j in 2:n.ages){
      

      if (k == 1){
        N.at.age.jki.gear[j,1,i+n.ages+1] <- exp(-(F_jk[j-1,n.seasons]*Ffull*sel.diff[j-1,n.seasons]+M.age[j]/n.seasons))*N.at.age.jki.gear[j-1,n.seasons,i+n.ages]  #season 1
      }else{
        N.at.age.jki.gear[j,k,i+n.ages+1] <- exp(-(F_jk[j,k-1]*Ffull*sel.diff[j,k-1]+M.age[j]/n.seasons))*N.at.age.jki.gear[j,k-1,i+n.ages+1]  #season 2 to 12
      }
      
    }
    for( j in 1: n.ages){
      
      C.at.age.jki.gear[j,k,i+n.ages+1] <- F_jk[j,k]*Ffull*sel.diff[j,k]/(F_jk[j,k]*Ffull*sel.diff[j,k]+M.age[j]/n.seasons)*(1-exp(-(F_jk[j,k]*Ffull*sel.diff[j,k]+M.age[j]/n.seasons)))*N.at.age.jki.gear[j,k,i+n.ages+1]
    }
  }
  
}


landings.gear<- numeric(10)
discards.gear<- numeric(10)
catches.gear<- numeric(10)
SSB.gear<- numeric(10)

for (i in 1:10){
  SSB.gear[i] <- sum(N.at.age.jki.gear[,n.seasons,n.ages+1+i]*mean_weight_jk[,n.seasons]*mat_prop[,n.seasons])
  landings.gear[i] <- sum(C.at.age.jki.gear[,,n.ages+1+i]*mean_cweight_jk*(1-dis_prop))
  discards.gear [i]<- sum(C.at.age.jki.gear[,,n.ages+1+i]*mean_cweight_jk*dis_prop)
  catches.gear [i]<- sum(C.at.age.jki.gear[,,n.ages+1+i]*mean_cweight_jk)
}

landings.gear<- c(rep(NA, n.ages), Lw.statusquo.sum[n.ages+1], landings.gear)
discards.gear<- c(rep(NA, n.ages), Dw.statusquo.sum[n.ages+1], discards.gear)
catches.gear<- c(rep(NA, n.ages), Cw.statusquo.sum[n.ages+1], catches.gear)
SSB.gear<- c(rep(NA, n.ages), SSB.statusquo.sum[n.ages+1], SSB.gear)

##----------------------
## F CHANGE: static projection after year n.ages+1 (years n.ages+2 to n.ages+11) for the scenario with a change in F
##----------------------
pred.years<-1:10
N.at.age.jki.newF <- array(0, c(n.ages,n.seasons,n.years))
C.at.age.jki.newF <- array(0, c(n.ages,n.seasons,n.years))

Ffull <- Fdeclined

for( i in 1:n.ages+1){
  
  for (k in 1:n.seasons){
    
    N.at.age.jki.newF[,k,i]<- apply(N.statusquo[,,k,i],2,sum)
    
  }
  
}

for (i in pred.years){
  
  N.at.age.jki.newF[1,1,i+n.ages+1] <- 14000 #recruitment

  
  for (k in 1:n.seasons){
    
    if (k != 1){
      N.at.age.jki.newF[1,k,i+n.ages+1] <- exp(-(F_jk[1,k-1]*Ffull+M.age[1]/n.seasons))*N.at.age.jki.newF[1,k-1,i+n.ages+1]  #season 2 to 12, age 1
    }
  
    for(j in 2:n.ages){
      

      if (k == 1){
        N.at.age.jki.newF[j,1,i+n.ages+1] <- exp(-(F_jk[j-1,n.seasons]*Ffull+M.age[j]/n.seasons))*N.at.age.jki.newF[j-1,n.seasons,i+n.ages]  #season 1
      }else{
        N.at.age.jki.newF[j,k,i+n.ages+1] <- exp(-(F_jk[j,k-1]*Ffull+M.age[j]/n.seasons))*N.at.age.jki.newF[j,k-1,i+n.ages+1]  #season 2 to 12
      }
      
    }
    
    for(j in 1: n.ages){
      
      C.at.age.jki.newF[j,k,i+n.ages+1] <- F_jk[j,k]*Ffull/(F_jk[j,k]*Ffull+M.age[j]/n.seasons)*(1-exp(-(F_jk[j,k]*Ffull+M.age[j]/n.seasons)))*N.at.age.jki.newF[j,k,i+n.ages+1]
    }
  }
  
}


landings.newF<- numeric(10)
discards.newF<- numeric(10)
catches.newF<- numeric(10)
SSB.newF<- numeric(10)

for (i in 1:10){
  SSB.newF[i] <- sum(N.at.age.jki.newF[,n.seasons,n.ages+1+i]*mean_weight_jk[,n.seasons]*mat_prop[,n.seasons])
  landings.newF[i] <- sum(C.at.age.jki.newF[,,n.ages+1+i]*mean_cweight_jk*(1-dis_prop))
  discards.newF [i]<- sum(C.at.age.jki.newF[,,n.ages+1+i]*mean_cweight_jk*dis_prop)
  catches.newF [i]<- sum(C.at.age.jki.newF[,,n.ages+1+i]*mean_cweight_jk)
}

landings.newF<- c(rep(NA, n.ages), Lw.statusquo.sum[n.ages+1], landings.newF)
discards.newF<- c(rep(NA, n.ages), Dw.statusquo.sum[n.ages+1], discards.newF)
catches.newF<- c(rep(NA, n.ages), Cw.statusquo.sum[n.ages+1], catches.newF)
SSB.newF<- c(rep(NA, n.ages), SSB.statusquo.sum[n.ages+1], SSB.newF)


## WHY IS SSB BASED ON MONTH 12 INSTEAD OF 1?
##SH260618: because SSB at the end of month 12 is the same as in the beginning of month 1 the following year. If you take month 1 it would be the SSB after one month
apply(N.at.age.jki[,,n.ages+2]*mean_weight_jk*mat_prop,2,sum)[n.seasons] # static
SSB.statusquo.sum[n.ages+2] #dynamic


##-----------------
## SSB: status quo versus gear change in dynamic and static
##-----------------
par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(5,5,5,5))
pred.year.date <- (1+n.ages-1):(1+n.ages+10)
start.year <- 1+n.ages-1

#png('SSB status quo vs gearchange.png', units="cm", width=8, height=9.5, res=600)
par(mfrow=c(1,1))
plot(SSB.statusquo.sum[20:31],type="n", col="darkblue", xlim=c(start.year,pred.year.date[length(pred.year.date)]), ylim=c(0,25000), xaxt="n", yaxt="n", ann=F)
axis(side=2, at=seq(0,200000,by=5000), font=1,cex.axis=0.8)
axis(side=1, at=seq(0,40,by=1), font=1,cex.axis=0.8)
abline(v=start.year:pred.year.date[length(pred.year.date)], col="grey")
abline(h=seq(0,100000,by=5000), col="grey")
# static
lines(pred.year.date, SSB[20:31], lwd=1)
lines(pred.year.date, SSB.gear[20:31], lwd=1.5, lty=5)
# dynamic
lines(pred.year.date, SSB.gearchange.sum[20:31], lwd=1.5, lty=3)
lines(pred.year.date,SSB.statusquo.sum[20:31], lwd=1.5, lty=1)
mtext(side=2,line=2.2, "SSB (tonnes)")
mtext(side=1, line=2.5,"simulation year")
#dev.off()
##-----------------
## SSB: status quo versus F change in dynamic and static
##-----------------
#png('SSB status quo vs Fchange.png', units="cm", width=8, height=9.5, res=600)
par(mfrow=c(1,1))
plot(pred.year.date,SSB.statusquo.sum[20:31],type="n", col="darkblue", xlim=c(start.year,pred.year.date[length(pred.year.date)]), ylim=c(0,20000), xaxt="n", yaxt="n", ylab="SSB",ann=F)
axis(side=2, at=seq(0,200000,by=5000), font=1,cex.axis=0.8)
axis(side=1, at=seq(0,40,by=1), font=1,cex.axis=0.8)
abline(v=start.year:pred.year.date[length(pred.year.date)], col="grey")
abline(h=seq(0,100000,by=5000), col="grey")
# static
lines(pred.year.date, SSB[20:31], lwd=1)
lines(pred.year.date, SSB.newF[20:31], lwd=1.5, lty=5)
# dynamic
lines(pred.year.date, SSB.Fchange.sum[20:31], lwd=1.5, lty=3)
lines(pred.year.date,SSB.statusquo.sum[20:31], lwd=1.5, lty=1)
mtext(side=2,line=2.2, "SSB (tonnes)")
mtext(side=1, line=2.5,"simulation year")
#dev.off()
##-----------------
## Total catches: status quo versus gear change in dynamic and static
##-----------------
#png('total catch status quo vs gearchange.png', units="cm", width=8, height=9.5, res=600)
par(mfrow=c(1,1))
plot(pred.year.date,Cw.statusquo.sum[20:31],type="n", col="darkblue", xlim=c(start.year,pred.year.date[length(pred.year.date)]), ylim=c(0,15000), xaxt="n", yaxt="n", ylab="catches",ann=F)
axis(side=2, at=seq(0,200000,by=5000), font=1,cex.axis=0.8)
axis(side=1, at=seq(0,40,by=1), font=1,cex.axis=0.8)
abline(v=start.year:pred.year.date[length(pred.year.date)], col="grey")
abline(h=seq(0,100000,by=5000), col="grey")
lines(pred.year.date, catches[20:31], lwd=1)
lines(pred.year.date, catches.gear[20:31], lwd=1.5, lty=5)
lines(pred.year.date, Cw.gearchange.sum[20:31], lwd=1.5, lty=3)
lines(pred.year.date,Cw.statusquo.sum[20:31], lwd=1.5, lty=1)
mtext(side=2,line=2.2, "total catch (tonnes)")
mtext(side=1, line=2.5,"simulation year")
#dev.off()


##-----------------
## Total catches: status quo versus F change in dynamic and static
##-----------------
#png('total catch status quo vs Fchange.png', units="cm", width=8, height=9.5, res=600)
par(mfrow=c(1,1))
plot(pred.year.date,Cw.statusquo.sum[20:31],type="n", col="darkblue", xlim=c(start.year,pred.year.date[length(pred.year.date)]), ylim=c(0,15000), xaxt="n", yaxt="n", ylab="catches",ann=F)
axis(side=2, at=seq(0,200000,by=5000), font=1,cex.axis=0.8)
axis(side=1, at=seq(0,40,by=1), font=1,cex.axis=0.8)
abline(v=start.year:pred.year.date[length(pred.year.date)], col="grey")
abline(h=seq(0,100000,by=5000), col="grey")
lines(pred.year.date, catches[20:31], lwd=1)
lines(pred.year.date, catches.newF[20:31], lwd=1.5, lty=5)
lines(pred.year.date, Cw.Fchange.sum[20:31], lwd=1.5, lty=3)
lines(pred.year.date,Cw.statusquo.sum[20:31], lwd=1.5, lty=1)
mtext(side=2,line=2.2, "total catch (tonnes)")
mtext(side=1, line=2.5,"simulation year")
#dev.off()


##-----------------
## Landings: status quo versus gear change in dynamic and static
##-----------------
#png('landings status quo vs gearchange.png', units="cm", width=8, height=9.5, res=600)
par(mfrow=c(1,1))
plot(pred.year.date,Lw.statusquo.sum[20:31],type="n", col="darkblue", xlim=c(start.year,pred.year.date[length(pred.year.date)]), ylim=c(0,15000), xaxt="n", yaxt="n", ylab="landings",ann=F)
axis(side=2, at=seq(0,200000,by=5000), font=1,cex.axis=0.8)
axis(side=1, at=seq(0,40,by=1), font=1,cex.axis=0.8)
abline(v=start.year:pred.year.date[length(pred.year.date)], col="grey")
abline(h=seq(0,100000,by=5000), col="grey")
lines(pred.year.date, landings[20:31], lwd=1)
lines(pred.year.date, landings.gear[20:31], lwd=1.5, lty=5)
lines(pred.year.date, Lw.gearchange.sum[20:31], lwd=1.5, lty=3)
lines(pred.year.date,Lw.statusquo.sum[20:31], lwd=1.5, lty=1)
mtext(side=2,line=2.2, "catch above MCRS (tonnes)")
mtext(side=1, line=2.5,"simulation year")
#dev.off()


##-----------------
## Landings: status quo versus F change in dynamic and static
##-----------------
#png('landings status quo vs Fchange.png', units="cm", width=8, height=9.5, res=600)
par(mfrow=c(1,1))
plot(pred.year.date,Lw.statusquo.sum[20:31],type="n", col="darkblue", xlim=c(start.year,pred.year.date[length(pred.year.date)]), ylim=c(0,15000), xaxt="n", yaxt="n", ylab="landings",ann=F)
axis(side=2, at=seq(0,200000,by=5000), font=1,cex.axis=0.8)
axis(side=1, at=seq(0,40,by=1), font=1,cex.axis=0.8)
abline(v=start.year:pred.year.date[length(pred.year.date)], col="grey")
abline(h=seq(0,100000,by=5000), col="grey")
lines(pred.year.date, landings[20:31], lwd=1)
lines(pred.year.date, landings.newF[20:31], lwd=1.5, lty=5)
lines(pred.year.date, Lw.Fchange.sum[20:31], lwd=1.5, lty=3)
lines(pred.year.date,Lw.statusquo.sum[20:31], lwd=1.5, lty=1)
mtext(side=2,line=2.2, "catch above MCRS (tonnes)")
mtext(side=1, line=2.5,"simulation year")
#dev.off()

##-----------------
## discards: status quo versus gear change in dynamic and static
##-----------------
#png('discards status quo vs gearchange.png', units="cm", width=8, height=9.5, res=600)
par(mfrow=c(1,1))
plot(pred.year.date,Dw.statusquo.sum[20:31],type="n", col="darkblue", xlim=c(start.year,pred.year.date[length(pred.year.date)]), ylim=c(0,1500), xaxt="n", yaxt="n", ylab="discards",ann=F)
axis(side=2, at=seq(0,200000,by=500), font=1,cex.axis=0.8)
axis(side=1, at=seq(0,40,by=1), font=1,cex.axis=0.8)
abline(v=start.year:pred.year.date[length(pred.year.date)], col="grey")
abline(h=seq(0,100000,by=500), col="grey")
lines(pred.year.date, discards[20:31], lwd=1)
lines(pred.year.date, discards.gear[20:31], lwd=1.5, lty=5)
lines(pred.year.date, Dw.gearchange.sum[20:31], lwd=1.5, lty=3)
lines(pred.year.date,Dw.statusquo.sum[20:31], lwd=1.5, lty=1)
mtext(side=2,line=2.2, "catch below MCRS (tonnes)")
mtext(side=1, line=2.5,"simulation year")
#dev.off()

##-----------------
## discards: status quo versus F change in dynamic and static
##-----------------

#png('discards status quo vs Fchange.png', units="cm", width=8, height=9.5, res=600)
par(mfrow=c(1,1))
plot(pred.year.date,Dw.statusquo.sum[20:31],type="n", col="darkblue", xlim=c(start.year,pred.year.date[length(pred.year.date)]), ylim=c(0,1500), xaxt="n", yaxt="n", ylab="discards",ann=F)
axis(side=2, at=seq(0,200000,by=500), font=1,cex.axis=0.8)
axis(side=1, at=seq(0,40,by=1), font=1,cex.axis=0.8)
abline(v=start.year:pred.year.date[length(pred.year.date)], col="grey")
abline(h=seq(0,100000,by=500), col="grey")
lines(pred.year.date, discards[20:31], lwd=1)
lines(pred.year.date, discards.newF[20:31], lwd=1.5, lty=5)
lines(pred.year.date, Dw.Fchange.sum[20:31], lwd=1.5, lty=3)
lines(pred.year.date,Dw.statusquo.sum[20:31],lwd=1.5, lty=1)
mtext(side=2,line=2.2, "catch below MCRS (tonnes)")
mtext(side=1, line=2.5,"simulation year")
#dev.off()

##
##-----------------
## Fbar: We have to calculate Fbar
##-----------------
# Fbar will be the mean of F ages 3-5
# F will here be the F over one year; for the dynamic, it is calculated by the log of the ratio of the N's from one year to the next, minus M
# for static it will be the status quo F times the sel.dif or the decline (for gearchange and Fchange respectively)
# I need a container with 5 rows (status quo, gearchange static, gearchange dynamic, Fchange static, Fchange dynamic) and 12 columns for years 20-31 for the graphs 
Fbar.container<-array(0,c(5,12))

#dynamic & static
#status quo
Fbar<-c(rep(0,10))
for (y in 1:10){
  F_jk<- c(rep(0,20)) 
  #this is copied from where we first calculate F for status quo but season is left out
  # j is for age
  for (j in 1:n.ages){
     if (j< n.ages){
      F_jk[j] <- log(sum(N.statusquo[,j,1,n.ages+y])/sum(N.statusquo[,j+1,1,n.ages+y+1]))-M.age[j]
    }
    if (j == n.ages) {
      #F for last age = previous age.
      F_jk[j] <- F_jk[j-1]
    }
  }
  Fbar[y]<-mean(F_jk[3:5])
}
Fbar.container[,1:2]<-Fbar[1]
Fbar.container[1,3:12]<-Fbar[1:10]

#F change
#static
Fbar.container[4,3:12]<-Fbar[1:10]*Fdeclined

#Fchange
#dynamic
Fbar<-c(rep(0,10))
for (y in 1:10){
  F_jk<- c(rep(0,20)) 
  # j is for age
  for (j in 1:n.ages){
    if (j< n.ages){
      F_jk[j] <- log(sum(N.Fchange[,j,1,n.ages+y])/sum(N.Fchange[,j+1,1,n.ages+y+1]))-M.age[j]
    }
    if (j == n.ages) {
      # F for last age = previous age.
      F_jk[j] <- F_jk[j-1]
    }
  }
  Fbar[y]<-mean(F_jk[3:5])
}
Fbar.container[5,2:11]<-Fbar[1:10]
Fbar.container[5,12]<-Fbar[10]

#gear change
#dynamic
Fbar<-c(rep(0,10))
for (y in 1:10){
  F_jk<- c(rep(0,20)) 
  # j is for age
  for (j in 1:n.ages){
    if (j< n.ages){
      F_jk[j] <- log(sum(N.gearchange[,j,1,n.ages+y])/sum(N.gearchange[,j+1,1,n.ages+y+1]))-M.age[j]
    }
    if (j == n.ages) {
      # F for last age = previous age.
      F_jk[j] <- F_jk[j-1]
    }
  }
  Fbar[y]<-mean(F_jk[3:5])
}
Fbar.container[3,2:11]<-Fbar[1:10]
Fbar.container[3,12]<-Fbar[10]

#static
#gear change
Fbar<-c(rep(0,10))
for (y in 1:10){
  F_jk<- c(rep(0,20)) 
  # j is for age
  for (j in 1:n.ages){
   if (j< n.ages){
      F_jk[j] <- log(sum(N.statusquo[,j,1,n.ages+y])/sum(N.statusquo[,j+1,1,n.ages+y+1]))-M.age[j]
    }
    if (j == n.ages) {
      # F for last age = previous age.
      F_jk[j] <- F_jk[j-1]
    }
  }
  Fbar[y]<-mean(sel.diff[3:5,]*F_jk[3:5])
}
Fbar
Fbar.container[2,3:12]<-Fbar[1:10]
##
##-----------------
## Fbar: status quo versus gear change in dynamic and static
##-----------------

#png('Fbar status quo vs gearchange.png', units="cm", width=8, height=9.5, res=600)
par(mfrow=c(1,1))
plot(pred.year.date,c(rep(1,12)),type="n", col="darkblue", xlim=c(start.year,pred.year.date[length(pred.year.date)]), ylim=c(0,1.5), xaxt="n", yaxt="n", ylab="discards",ann=F)
axis(side=2, at=seq(0,2,by=1), font=1,cex.axis=0.8)
axis(side=1, at=seq(0,40,by=1), font=1,cex.axis=0.8)
abline(v=start.year:pred.year.date[length(pred.year.date)], col="grey")
abline(h=seq(0,10,by=0.5), col="grey")
lines(pred.year.date, Fbar.container[2,], lwd=1.5, lty=5)
lines(pred.year.date, Fbar.container[3,], lwd=1.5, lty=3)
lines(pred.year.date,Fbar.container[1,], lwd=1.5, lty=1)
#lines(pred.year.date, c(rep(1.16976,2),rep(1.091517,10)), lwd=1.5, lty=5)
#lines(pred.year.date, c(rep(1.1697597,2), 1.082498, 1.076351, 1.087492, 1.092792, rep(1.093671,6)), lwd=1.5, lty=3)
#lines(pred.year.date,c(rep(1.16976,12)), lwd=1.5, lty=1)
#abline(v=2002+n.ages+1,lty=2)
#legend("topleft", c("F = 1", "F = 0.5","F=1, naive", "Blim"), col=c("darkblue", "magenta","green", "black"), pch="-")
#legend("topleft", c("status quo, dynamic", "higher L50, dynamic","status quo, static","higher L50, static"), col=c("magenta","blue","green", "red"), pch="-")
#legend("topleft", c("status quo", "higher L50, dynamic","higher L50, static"), col=c("green","blue","red"), pch="-")
mtext(side=2,line=2.2, expression("F"[3-5]),cex=1.2)
mtext(side=1, line=2.5,"simulation year")
#dev.off()

##-----------------
## Fbar: status quo versus F change in dynamic and static
##-----------------

#png('Fbar status quo vs Fchange.png', units="cm", width=8, height=9.5, res=600)
par(mfrow=c(1,1))
plot(pred.year.date,c(rep(1,12)),type="n", col="darkblue", xlim=c(start.year,pred.year.date[length(pred.year.date)]), ylim=c(0,1.5), xaxt="n", yaxt="n", ylab="discards", ann=F)
axis(side=2, at=seq(0,2,by=1), font=1,cex.axis=0.8)
axis(side=1, at=seq(0,40,by=1), font=1,cex.axis=0.8)
abline(v=start.year:pred.year.date[length(pred.year.date)], col="grey")
abline(h=seq(0,10,by=0.5), col="grey")
lines(pred.year.date, Fbar.container[4,], lwd=1.5, lty=5)
lines(pred.year.date, Fbar.container[5,], lwd=1.5, lty=3)
lines(pred.year.date,Fbar.container[1,], lwd=1.5, lty=1)
#lines(pred.year.date, c(rep(1.16976,2),rep(0.8773198,10)), lwd=1.5, lty=5)
#lines(pred.year.date, c(rep(1.1697597,2), 0.9235584, 0.9247908, 0.9253744, 0.9254978, rep(0.9255135,6)), lwd=1.5, lty=3)
#lines(pred.year.date,c(rep(1.16976,12)), lwd=1.5, lty=1)
#abline(v=2002+n.ages+1,lty=2)
#legend("topleft", c("F = 1", "F = 0.5","F=1, naive", "Blim"), col=c("darkblue", "magenta","green", "black"), pch="-")
#legend("topleft", c("status quo, dynamic", "lower F, dynamic","status quo, static","lower F, static"), col=c("magenta","blue","green", "red"), pch="-")
#legend("topleft", c("status quo", "lower F, dynamic","lower F, static"), col=c("green","blue","red"), pch="-")
mtext(side=2,line=2.2, expression("F"[3-5]),cex=1.2)
mtext(side=1, line=2.5,"simulation year")
#dev.off()
##

###
# All 4 plots of gear-change vs status quo and of F-change vs status quo

png('status quo vs gear-change.png', units="cm", width=17, height=15, res=600)
par(mfrow=c(2,2), mar=c(1,1,0,0), oma=c(3,3,1,4))

plot(SSB.statusquo.sum[20:31],type="n", col="darkblue", xlim=c(start.year,pred.year.date[length(pred.year.date)]), ylim=c(0,25000), xaxt="n", yaxt="n", ann=F)
axis(side=2, at=seq(0,200000,by=5000), font=1,cex.axis=0.8)
abline(v=start.year:pred.year.date[length(pred.year.date)], col="grey")
abline(h=seq(0,100000,by=5000), col="grey")
# static
lines(pred.year.date, SSB[20:31], lwd=1)
lines(pred.year.date, SSB.gear[20:31], lwd=1.5, lty=5)
# dynamic
lines(pred.year.date, SSB.gearchange.sum[20:31], lwd=1.5, lty=3)
lines(pred.year.date,SSB.statusquo.sum[20:31], lwd=1.5, lty=1)
mtext(side=2,line=2.2, "SSB (t)")


plot(pred.year.date,Lw.statusquo.sum[20:31],type="n", col="darkblue", xlim=c(start.year,pred.year.date[length(pred.year.date)]), ylim=c(0,15000), xaxt="n", yaxt="n", ylab="landings",ann=F)
axis(side=4, at=seq(0,200000,by=5000), font=1,cex.axis=0.8)
abline(v=start.year:pred.year.date[length(pred.year.date)], col="grey")
abline(h=seq(0,100000,by=5000), col="grey")
lines(pred.year.date, landings[20:31], lwd=1)
lines(pred.year.date, landings.gear[20:31], lwd=1.5, lty=5)
lines(pred.year.date, Lw.gearchange.sum[20:31], lwd=1.5, lty=3)
lines(pred.year.date,Lw.statusquo.sum[20:31], lwd=1.5, lty=1)
mtext(side=4,line=2.2, "catch above MCRS (t)")


plot(pred.year.date,Dw.statusquo.sum[20:31],type="n", col="darkblue", xlim=c(start.year,pred.year.date[length(pred.year.date)]), ylim=c(0,1500), xaxt="n", yaxt="n", ylab="discards",ann=F)
axis(side=2, at=seq(0,200000,by=500), font=1,cex.axis=0.8)
axis(side=1, at=seq(0,40,by=1), font=1,cex.axis=0.8)
abline(v=start.year:pred.year.date[length(pred.year.date)], col="grey")
abline(h=seq(0,100000,by=500), col="grey")
lines(pred.year.date, discards[20:31], lwd=1)
lines(pred.year.date, discards.gear[20:31], lwd=1.5, lty=5)
lines(pred.year.date, Dw.gearchange.sum[20:31], lwd=1.5, lty=3)
lines(pred.year.date,Dw.statusquo.sum[20:31], lwd=1.5, lty=1)
mtext(side=2,line=2.2, "catch below MCRS (t)")

#tracking changes in mean length@age over time
mean_length_at_age<-array(0,c(n.ages,n.ages+length(pred.year.date)-1,2))
for (j in 1:n.ages){
  for (i in 1:length(pred.year.date)){
    # for gear change
    mean_length_at_age[j,n.ages-1+i,1]<-sum(l.mids*(N.gearchange[,j,n.seasons,n.ages-1+i]))/sum(N.gearchange[,j,n.seasons,n.ages-1+i])  
  }
}
plot(mean_length_at_age[1,,1],type="n", xlim=c(start.year,pred.year.date[length(pred.year.date)]), ylim=c(30,140),xaxt="n", yaxt="n")
axis(side=4, at=seq(0,200,by=20), font=1,cex.axis=0.8)
axis(side=1, at=seq(0,40,by=1), font=1,cex.axis=0.8)
abline(h=seq(0,200,by=10), col="grey")
abline(v=start.year:pred.year.date[length(pred.year.date)], col="grey")
for (j in 1:15){
  lines(pred.year.date,mean_length_at_age[j,20:31,1])
}
mtext(side = 4, line = 2.2, "mean length at age (cm)")
mtext(side=1, line=1.5, outer = TRUE, "simulation year")


dev.off()

png('status quo vs F-change.png', units="cm", width=17, height=15, res=600)
par(mfrow=c(2,2), mar=c(1,1,0,0), oma=c(3,3,1,4))

plot(SSB.statusquo.sum[20:31],type="n", col="darkblue", xlim=c(start.year,pred.year.date[length(pred.year.date)]), ylim=c(0,25000), xaxt="n", yaxt="n", ann=F)
axis(side=2, at=seq(0,200000,by=5000), font=1,cex.axis=0.8)
abline(v=start.year:pred.year.date[length(pred.year.date)], col="grey")
abline(h=seq(0,100000,by=5000), col="grey")
# static
lines(pred.year.date, SSB[20:31], lwd=1)
lines(pred.year.date, SSB.newF[20:31], lwd=1.5, lty=5)
# dynamic
lines(pred.year.date, SSB.Fchange.sum[20:31], lwd=1.5, lty=3)
lines(pred.year.date,SSB.statusquo.sum[20:31], lwd=1.5, lty=1)
mtext(side=2,line=2.2, "SSB (t)")


plot(pred.year.date,Lw.statusquo.sum[20:31],type="n", col="darkblue", xlim=c(start.year,pred.year.date[length(pred.year.date)]), ylim=c(0,15000), xaxt="n", yaxt="n", ylab="landings",ann=F)
axis(side=4, at=seq(0,200000,by=5000), font=1,cex.axis=0.8)
abline(v=start.year:pred.year.date[length(pred.year.date)], col="grey")
abline(h=seq(0,100000,by=5000), col="grey")
lines(pred.year.date, landings[20:31], lwd=1)
lines(pred.year.date, landings.newF[20:31], lwd=1.5, lty=5)
lines(pred.year.date, Lw.Fchange.sum[20:31], lwd=1.5, lty=3)
lines(pred.year.date,Lw.statusquo.sum[20:31], lwd=1.5, lty=1)
mtext(side=4,line=2.2, "catch above MCRS (t)")


plot(pred.year.date,Dw.statusquo.sum[20:31],type="n", col="darkblue", xlim=c(start.year,pred.year.date[length(pred.year.date)]), ylim=c(0,1500), xaxt="n", yaxt="n", ylab="discards",ann=F)
axis(side=2, at=seq(0,200000,by=500), font=1,cex.axis=0.8)
axis(side=1, at=seq(0,40,by=1), font=1,cex.axis=0.8)
abline(v=start.year:pred.year.date[length(pred.year.date)], col="grey")
abline(h=seq(0,100000,by=500), col="grey")
lines(pred.year.date, discards[20:31], lwd=1)
lines(pred.year.date, discards.newF[20:31], lwd=1.5, lty=5)
lines(pred.year.date, Dw.Fchange.sum[20:31], lwd=1.5, lty=3)
lines(pred.year.date,Dw.statusquo.sum[20:31], lwd=1.5, lty=1)
mtext(side=2,line=2.2, "catch below MCRS (t)")

#tracking changes in mean length@age over time
for (j in 1:n.ages){
  for (i in 1:length(pred.year.date)){
    # for F change
    mean_length_at_age[j,n.ages-1+i,2]<-sum(l.mids*(N.Fchange[,j,n.seasons,n.ages-1+i]))/sum(N.Fchange[,j,n.seasons,n.ages-1+i])  
  }
}
plot(mean_length_at_age[1,,2],type="n", xlim=c(start.year,pred.year.date[length(pred.year.date)]), ylim=c(30,140),xaxt="n", yaxt="n")
axis(side=4, at=seq(0,200,by=20), font=1,cex.axis=0.8)
axis(side=1, at=seq(0,40,by=1), font=1,cex.axis=0.8)
abline(h=seq(0,200,by=10), col="grey")
abline(v=start.year:pred.year.date[length(pred.year.date)], col="grey")
for (j in 1:15){
  lines(pred.year.date,mean_length_at_age[j,20:31,2])
}
mtext(side = 4, line = 2.2, "mean length at age (cm)")
mtext(side=1, line=1.5, outer = TRUE, "simulation year")


dev.off()

## calculate the biases in year 31 (= year 10 after the change)
(SSB.gear[31]-SSB.gearchange.sum[31])*100/SSB.gearchange.sum[31]
(landings.gear[31]-Lw.gearchange.sum[31])*100/Lw.gearchange.sum[31]
(discards.gear[31]-Dw.gearchange.sum[31])*100/Dw.gearchange.sum[31]
(SSB.newF[31]-SSB.Fchange.sum[31])*100/SSB.Fchange.sum[31]
(landings.newF[31]-Lw.Fchange.sum[31])*100/Lw.Fchange.sum[31]
(discards.newF[31]-Dw.Fchange.sum[31])*100/Dw.Fchange.sum[31]

#####
# make a plot of three Von Bertalanffy curves, for the 3 growth scenarios
VBages <- c(1:15)
VBlength.slow <- VBlength.fast <- VBlength.default <- VBages
VBlength.default[] <- 154.56 * (1 - exp(-0.11 * (VBages[]-t0.val)))
VBlength.fast[] <- 135 * (1 - exp(-0.165 * (VBages[]-t0.val)))
VBlength.slow[] <- 178.5 * (1 - exp(-0.073 * (VBages[]-t0.val)))
#png('three growth curves.png', units="cm", width=10, height=10, res=600)
plot(VBlength.default,type="n", col="darkblue", xaxt="n", yaxt="n", ann=F)
axis(side=2, at=seq(0,200,by=20), font=1,cex.axis=0.8)
axis(side=1, at=seq(0,40,by=1), font=1,cex.axis=0.8)
abline(h=seq(0,140, by=20), col="grey")
abline(v=seq(0,15,by=1), col="grey")
lines(VBlength.fast, lwd=1.5, lty=5)
lines(VBlength.slow, lwd=1.5, lty=3)
lines(VBlength.default, lwd=1.5, lty=1)
mtext(side=2,line=2.2, "length (cm)")
mtext(side=1,line=2.2, "age (years)")
#dev.off()

#####
# make 100 VB curves, as an illustration for the growth trajectories used for the G matrix for growth increments
VBages <- c(1:15)
VB<-array(0,c(15,100))
vcov <- matrix(c(sd.linf^2, rep(rho * sd.linf * sd.k, 2), sd.k^2), nrow = 2)
## simulate the linf and k draws
linf.k.mat <- rmvnorm(100, mean = c(linf.val, k.val), sigma = vcov)
#png('100 growth curves temp.png', units="cm", width=10, height=10, res=600)
plot(VBlength.default,type="n", col="darkblue", xaxt="n", yaxt="n", ylim=c(0,200), ann=F)
axis(side=2, at=seq(0,350,by=20), font=1,cex.axis=0.8)
axis(side=1, at=seq(0,40,by=1), font=1,cex.axis=0.8)
abline(h=seq(0,200, by=20), col="grey")
abline(v=seq(0,15,by=1), col="grey")
for (i in 1:100){
  VB[,i] <- linf.k.mat[i,1] * (1 - exp(-linf.k.mat[i,2] * (VBages[]-t0.val)))  
  lines(VB[,i], lwd=0.5, lty=1)
}
mtext(side=2,line=2.2, "length (cm)")
mtext(side=1,line=2.2, "age (years)")
#dev.off()
