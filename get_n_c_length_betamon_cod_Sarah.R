get.N.C.length<-function(N.age, age.vec, l.breaks, l.mids,G,linf.val, k.val, t0.val, sd.l.age, n.years, n.seasons, a1, a2, b1, b2, M.age, F.age, Ffull, scenario="statusquo", Fbar, Fdeclined){   
  ##---------------------------------------
  ## numbers and catch at length
  ## Main function to generate numbers at length
  ## arguments
  ## n.age: numbers at age in first year
  ## scenario= "statusquo" or "gearchange" or "Fchange"
  ##---------------------------------------
  #print("Working on projections, message update shortly ...")
  ##-------
  ## SETUP
  ##-------

  ## number of ages
  n.ages<-length(age.vec)
  ## number of length bins
  n.bins<-length(l.mids)
  ## container for population numbers at length

  n.l.array<-array(0,c(n.bins,n.ages,n.seasons,n.years))

  ## container for catch numbers at length
  c.l.array<-n.l.array
  ##container for F, recruitment and ssb (in between calculation)

  rec.array <- rep(0,n.years)
  ssb.l.array <- n.l.array
  ssb.l.sum <- rep(0,n.years)
 

  ##------------
  ## First year
  ##------------
  ## need to distribute numbers at age to numbers at length
  ## according to fishing since birth by old gear
  ## fill out the numbers at length by age
  a.pred<-a1
  b.pred<-b1
  #year 1 calculations
  for(j in 1:n.ages){
    n.l.array[,j,1,1]<-get.length.dist.age(age=age.vec[j], N.vec=N.age, F.vec=F.age, M.vec=M.age, l.breaks=l.breaks, l.mids=l.mids, G=G, linf=linf.val, k=k.val, t0=t0.val, sd.l.age=sd.l.age, n.seasons=n.seasons, a=a.pred, b=b.pred)
  }

  
    # 05.07.2018->Current version: next four lines extracted from the loop below the lines  
  rec.array<- c(rep(14000,n.years))
  cdf.l<-plnorm(l.breaks[1:n.bins], meanlog=log(linf.val*(1-exp(-k.val*(1-t0.val))))-sd.l.age^2/2, sdlog=sd.l.age)
  cdf.l1<-plnorm(l.breaks[2:(n.bins+1)], meanlog=log(linf.val*(1-exp(-k.val*(1-t0.val))))-sd.l.age^2/2, sdlog=sd.l.age)
  prop.l<-cdf.l1-cdf.l
  
  for(i in 1:n.years){
    n.l.array[,1,1,i]<-prop.l*rec.array[i]
                      }

 
  ##------------------
  ## Main projections
  ##------------------
  ## get the correct selectivity parameter

  for(i in 1:n.years){
    
    

    a.pred<-ifelse(scenario=="gearchange" & i > n.ages+1, a2, a1)
    b.pred<-ifelse(scenario=="gearchange" & i > n.ages+1, b2, b1)

    print(paste("Working on projection year", i))

    s.length <- plogis(a.pred + b.pred * l.mids)  #trawl
    
    
## the alternative experiment is that Ffull changes instead of the gear selctivity
    Ffull<- ifelse(scenario=="Fchange" & i > n.ages+1, Fdeclined, Ffull)
    
    for(j in 1:n.ages){
      age.val<-age.vec[j]
      

      
      ## IF ONLY ONE SEASON
      if(n.seasons==1){


        ## survivorship matrix
     
        S<-get.survivorship.mat(l.mids=l.mids,l.breaks=l.breaks,age.val=age.val,Ffull=Ffull,M=M.age,a=a.pred,b=b.pred,n.seasons=n.seasons)
        ## numbers projection
        if(j<n.ages & i<n.years){
            n.l.array[,j+1,1,i+1]<-G[[j]]%*%S%*%n.l.array[,j,1,i]
            

        }
        ## catches
        
        
        C <- get.catch.mat(l.mids=l.mids,l.breaks=l.breaks,age.val=age.val,Ffull=Ffull,M=M.age,a=a.pred,b=b.pred,n.seasons=n.seasons)
        c.l.array[,j,1,i]<-C%*%n.l.array[,j,1,i]
         
        

      }else{
        ## > 1 SEASON
       
        
        for(k in 1:(n.seasons-1)){
          
          cat(" Year=", i,"Month=",k,"\n")

          # survivorship matrix
          S<-get.survivorship.mat(l.mids=l.mids,l.breaks=l.breaks,age.val=age.val,Ffull=Ffull,M=M.age,a=a.pred,b=b.pred,n.seasons=n.seasons)
          # numbers
          n.l.array[,j,k+1,i]<-G[[j]][[k]]%*%S%*%n.l.array[,j,k,i]
          # catches
          
          C <- get.catch.mat(l.mids=l.mids,l.breaks=l.breaks,age.val=age.val,Ffull=Ffull,M=M.age,a=a.pred,b=b.pred,n.seasons=n.seasons)
          c.l.array[,j,k,i]<-C%*%n.l.array[,j,k,i]
        }
        ## numbers projection
        if(j<n.ages & i<n.years){

            S<-get.survivorship.mat(l.mids=l.mids,l.breaks=l.breaks,age.val=age.val,Ffull=Ffull,M=M.age,a=a.pred,b=b.pred,n.seasons=n.seasons)
            n.l.array[,j+1,1,i+1]<-G[[j]][[k+1]]%*%S%*%n.l.array[,j,k+1,i]
        }
        ## catches in last season
        C <- get.catch.mat(l.mids=l.mids,l.breaks=l.breaks,age.val=age.val,Ffull=Ffull,M=M.age,a=a.pred,b=b.pred,n.seasons=n.seasons)
        c.l.array[,j,k+1,i]<-C%*%n.l.array[,j,k+1,i]
      }
      
    }
  }
  
  return(list("N"=n.l.array, "C"=c.l.array, "Rec"=rec.array)) 
}
