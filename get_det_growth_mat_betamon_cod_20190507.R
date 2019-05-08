

##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
## get.delta.l.mon: Function for estimating empirical growth matrix using montecarlo 
##
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Arguments:
# age = Fish age evaluated
# seasons = number of seasons evaluated. If n.seasons= 1, the growth will be estimated in a yearly basis
# linf, k, t0 =   Average BV parameters 
# sd.linf, sd.k = Standard deviations of BV parameters
# rho = correlation between parameters
# B number of pseudo-VB curves used in the montecarlo simulation


get.det.growth.mat.mon<-function(age,n.seasons,l.breaks,l.mids,n.bins,linf,sd.linf,k,sd.k,rho,t0,B){
  
  
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
  
  
  
if(n.seasons==1){
  
                   ageplus <- age + 1
                   ages<-c(age,ageplus)
  
                   l.mat <- sapply(1:B, FUN = function(z){linf.k.mat[z, 1] * (1 - exp(-linf.k.mat[z, 2] * (ages - t0)))})
                   ## bin individuals at this age
                   lbin <- cut(l.mat[1, ], breaks = l.breaks)
                   ## bin individuals at next age
                   lbinplus <- cut(l.mat[2, ], breaks = l.breaks)
                   ltab <- table(lbinplus,lbin)
                   G <- unclass(prop.table(ltab, margin = 2))
                   G[is.nan(G)]<-0
                   #G<-t(G)
  
  
                  } else {
                    
                    ages=c(age,age+1:n.seasons/n.seasons)
                    
                    l.mat <- sapply(1:B, FUN = function(z){linf.k.mat[z, 1] * (1 - exp(-linf.k.mat[z, 2] * (ages- t0)))})
                    
                    G<-llply(1:(n.seasons),function(j){
                                                          lbin <- cut(l.mat[j, ], breaks = l.breaks)
                                                          ## bin individuals at next age
                                                          lbinplus <- cut(l.mat[j+1, ], breaks = l.breaks)
                                                          ltab <- table(lbin,lbinplus)
                                                          g <- unclass(prop.table(ltab, margin = 1))
                                                          g[is.nan(g)]<-0
                                                          g<-t(g)
                                                          return(g)
                                                          })
                    
                          }
      
  return(G)}
  
  
  
  
  
  
  
  

