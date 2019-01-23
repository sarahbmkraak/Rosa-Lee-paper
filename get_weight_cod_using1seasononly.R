get.weight<-function(length.vec, season, n.seasons){
  ## get predicted cod weight for a given length and quarter
  ## in kgs
#  if(n.seasons%%2==1 & n.seasons!=1){stop("get.weight only implemented for even number of seasons greater than 1")}
#  if(n.seasons==1){
    Intercept <-mean(c(exp(-10.9097),exp(-11.2358),exp(-11.2382),exp(-12.3398)))
    Slope <- mean(c(2.8827,2.9545,2.9425,3.1414))

 ## Intercept and slope from Sven Stötera (commercial data 2016, SD24)   
    
   
#  }else{
#    if(n.seasons==2){stop("no weight data available")
#      if(season==1){
#        Intercept <- mean(c(exp(-10.9097),exp(-11.2358)))
#        Slope <- mean(c(2.8827, 2.9545))
#      }
#      if(season==2)
#        Intercept <- mean(c(exp(-11.2382),exp(-12.3398)))
#      Slope <- mean(c(2.9425,3.1414))
#    }
#    if(n.seasons>=4){
#      ## find out which quarter
#      quarter.cutoff<-cbind(c(0,0.25,0.5,0.75), c(0.25,0.5,0.75,1))
#      quarter<-which(season/n.seasons>quarter.cutoff[,1]&season/n.seasons<=quarter.cutoff[,2])
#      if(quarter==1){
#        Intercept <- exp(-10.9097)
#        Slope <- 2.8827
#      }
#      if(quarter==2){
#        Intercept <- exp(-11.2358)
#        Slope <- 2.9545
#      }
#      if(quarter==3){
#        Intercept <- exp(-11.2382)
#        Slope <- 2.9425
#      }
#      if(quarter==4){
#        Intercept <- exp(-12.3398)
#        Slope <- 3.1414
#      }
#    }
#  }
  ## note should have bias correction here
  #pred.weight<-exp(Intercept+Slope*log(length.vec))
  pred.weight <- Intercept* length.vec ^Slope
  return(pred.weight)    
}

