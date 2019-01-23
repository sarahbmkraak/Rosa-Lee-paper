get.maturity_cod<-function(length.vec, age){
  ## get predicted cod maturity at length for a given age
  prop.mature<-rep(0,length(length.vec))
  
  
  if(age==1){
    L50 <- 19.304
    Slope <- 4.038
    #L50 <- 19.074
    #Slope <- 4.005
     prop.mature<- 1/(1+exp(-Slope*(length.vec-L50)))
    
 
  }
  if(age==2){
    L50 <- 24.523
    Slope <- 4.337
    #L50 <- 24.199
    #Slope <- 4.481
    prop.mature<- 1/(1+exp(-Slope*(length.vec-L50))) 
    
  }
  if(age==3){
    L50 <- 28.583
    Slope <- 7.638
    #L50 <- 28.094
    #Slope <- 7.507
    prop.mature<- 1/(1+exp(-Slope*(length.vec-L50)))  
    
  }
  if(age==4){
    L50 <- 29.594
    Slope <- 12.026
    #L50 <- 31.48
    #Slope <- 6.431
    prop.mature<- 1/(1+exp(-Slope*(length.vec-L50))) 
    
  }
  
  if(age>=5){
    prop.mature<-1   
  }
  return(prop.mature)
}

