
get.length.dist.age<-function(age, N.vec, F.vec, M.vec, l.breaks, l.mids, G, linf, k, t0, sd.l.age, n.seasons, a, b){
  ##------------------------------------------
  ## calculate the length distribution at age
  ## considering that the cohort was fished
  ## with the status quo gear since recruitment (fish recruit in year 1)
  ##------------------------------------------
  require(truncnorm)
  n.bins<-length(l.mids)
  ## number at age 1, choose from different scenarios
  if(age==1){
    N1 <- N.vec[1]    #Recruitment from year N1
  }else{
    N1<- N.vec[age]*exp(sum(F.vec[1:(age-1)]+M.vec[1:(age-1)]))
  }
    ## spread the N1 out among the size classes
    ## NOTE! fixed lognormal bias correction here
  cdf.l<-plnorm(l.breaks[1:n.bins], meanlog=log(linf*(1-exp(-k*(1-t0))))-sd.l.age^2/2, sdlog=sd.l.age)
  cdf.l1<-plnorm(l.breaks[2:(n.bins+1)], meanlog=log(linf*(1-exp(-k*(1-t0))))-sd.l.age^2/2, sdlog=sd.l.age)
  prop.l<-cdf.l1-cdf.l
  ## number in first season, age 1 
  n.l<-prop.l*N1
  ## age 1
  if(age==1){
    return(n.l)
  }else{
    ## container
    age.vec <- 1:age
    n.ages<-length(age.vec)
    ## check this
    n.l.array<-array(0,c(n.bins,n.seasons,n.ages))
    n.l.array[,1,1]<-n.l
    
    ## grow these forward under F and M and selectivity
    for(i in 1:n.ages){ 
      age.val<-age.vec[i]
      ## If only one season
      if(n.seasons==1){
        S<-get.survivorship.mat(l.mids=l.mids,l.breaks=l.breaks,age.val=age.val,Ffull=Ffull,M=M.vec,a=a,b=b,n.seasons=n.seasons)
        #a.pred and b.pred with selectivity parameters from T0 130, which was mainly used in 2003
        if(i<n.ages){
          n.l.array[,1,i+1]<-G[[i]]%*%S%*%n.l.array[,1,i]
          
        }
      
      }else{
        for(j in 1:(n.seasons-1)){
          age.val<-age.vec[i]

         
          S<-get.survivorship.mat(l.mids=l.mids,l.breaks=l.breaks,age.val=age.val,Ffull=Ffull,M=M.vec,a=a,b=b,n.seasons=n.seasons)
          n.l.array[,j+1,i]<-G[[i]][[j]]%*%S%*%n.l.array[,j,i]
         
          
          
        }
        if(i<n.ages){
          S<-get.survivorship.mat(l.mids=l.mids,l.breaks=l.breaks,age.val=age.val,Ffull=Ffull,M=M.vec,a=a,b=b,n.seasons=n.seasons)
          n.l.array[,1,i+1]<-G[[i]][[j+1]]%*%S%*%n.l.array[,j+1,i]
          ##print(c(linf,k,t0))
        }
      }
    }
    
    return(n.l.array[,1,i])  
  }
  
}
