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
source("C:/Users/kraak/Dropbox/Master Thesis/R_rosalee/get_det_growth_mat_betamon_cod_20190507.R")
source("C:/Users/kraak/Dropbox/Master Thesis/R_rosalee/get_survivorship_mat_cod_Sarah.R")
source("C:/Users/kraak/Dropbox/Master Thesis/R_rosalee/get_catch_mat_cod_Sarah.R")
source("C:/Users/kraak/Dropbox/Master Thesis/R_rosalee/get_length_dist_age_betamon_cod_Sarah_revision.R")
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
n.seasons <- 1

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
l.breaks<-seq(1,100, by=10)

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



linf<-linf.val
k<-k.val
t0=t0.val
B=1000

require(plyr)
require(mvtnorm)

## covariance matrix k and linf
vcov <- matrix(c(sd.linf^2, rep(rho * sd.linf * sd.k, 2), sd.k^2), nrow = 2)

## simulate the linf and k draws
linf.k.mat <- abs(rmvnorm(B, mean = c(linf, k), sigma = vcov))


# Plug-in chunk to avoid negative parameters (from 23 november 2018 version) ------------------------------

for(i in 1:1000){
  
  if(any(linf.k.mat[,1]<0)==T){
    
    linf.k.mat[,1]<- ifelse(linf.k.mat[,1]<0,rmvnorm(1, mean = c(linf, k), sigma = vcov) , linf.k.mat[,1])
    
    
  }
  
  if(any(linf.k.mat[,2]<0)==T){
    
    linf.k.mat[,2]<- ifelse(linf.k.mat[,2]<0,rmvnorm(1, mean = c(linf, k), sigma = vcov) , linf.k.mat[,2])
    
    print(i)
    
  } else cat("\nno negative values found")
  
  
  
}


png('age and size-based length increments.png', units="cm", width=17, height=10, res=600)

par(mfrow=c(2,2), mar=c(1,1,0,0), oma=c(3,3,1,4))

age<-2

ageplus <- age + 1
ages<-c(age,ageplus)
l.mat <- sapply(101:200, FUN = function(z){linf.k.mat[z, 1] * (1 - exp(-linf.k.mat[z, 2] * (ages - t0)))})

plot(l.mat[1,],(l.mat[2,]-l.mat[1,]),xlim=c(0,120),ylim=c(0,25),xaxt="n", yaxt="n", ann=F, pch=20, cex=0.5)
abline(v=seq(0,120,by=20), col="grey")
abline(h=seq(0,25,by=5), col="grey")
axis(side=2, at=seq(0,25,by=5), font=1,cex.axis=1, las=1)
mtext(side=2,line=1.5, outer = TRUE, "Length increment (cm)")
mtext('age = 2', side = 3, line = -1.2, adj = 0.05, cex = 1, col = "black")

age<-3

ageplus <- age + 1
ages<-c(age,ageplus)
l.mat <- sapply(101:200, FUN = function(z){linf.k.mat[z, 1] * (1 - exp(-linf.k.mat[z, 2] * (ages - t0)))})

plot(l.mat[1,],(l.mat[2,]-l.mat[1,]),xlim=c(0,120),ylim=c(0,25),xaxt="n", yaxt="n", ann=F, pch=20, cex=0.5)
abline(v=seq(0,120,by=20), col="grey")
abline(h=seq(0,25,by=5), col="grey")
mtext('age = 3', side = 3, line = -1.2, adj = 0.05, cex = 1, col = "black")

age<-4

ageplus <- age + 1
ages<-c(age,ageplus)
l.mat <- sapply(101:200, FUN = function(z){linf.k.mat[z, 1] * (1 - exp(-linf.k.mat[z, 2] * (ages - t0)))})

plot(l.mat[1,],(l.mat[2,]-l.mat[1,]),xlim=c(0,120),ylim=c(0,25),xaxt="n", yaxt="n", ann=F, pch=20, cex=0.5)
abline(v=seq(0,120,by=20), col="grey")
abline(h=seq(0,25,by=5), col="grey")
axis(side=2, at=seq(0,25,by=5), font=1,cex.axis=1, las=1)
axis(side=1, at=seq(0,120,by=20), font=1,cex.axis=1)
mtext('age = 4', side = 3, line = -1.2, adj = 0.05, cex = 1, col = "black")
mtext(side=1, line=1.5, outer = TRUE, "Length at age (cm)")

age<-5

ageplus <- age + 1
ages<-c(age,ageplus)
l.mat <- sapply(101:200, FUN = function(z){linf.k.mat[z, 1] * (1 - exp(-linf.k.mat[z, 2] * (ages - t0)))})

plot(l.mat[1,],(l.mat[2,]-l.mat[1,]),xlim=c(0,120),ylim=c(0,25),xaxt="n", yaxt="n", ann=F, pch=20, cex=0.5)
abline(v=seq(0,120,by=20), col="grey")
abline(h=seq(0,25,by=5), col="grey")
axis(side=1, at=seq(0,120,by=20), font=1,cex.axis=1)
mtext('age = 5', side = 3, line = -1.2, adj = 0.05, cex = 1, col = "black")
mtext(side=1, line=1.5, outer = TRUE, "Length at age (cm)")

dev.off()
