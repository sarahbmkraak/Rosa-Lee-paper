
get.catch.mat<-function(l.mids,l.breaks,age.val,Ffull,M,a,b,n.seasons){
  
  cat(" F value is=",Ffull,"\n","Mortality is =",M,"\n","Age is =",age.val,"\n","\n")
  ##---------------------------------------------------------------
  ## function to create a removals matrix
  ## Arguments:
  ## l.mids: vector of midpoints of length classes
  ## l.breaks: vector of breakpoints of length classes
  ## age.val: age
  ## Ffull: fully selected F
  ## M: natural mortality at age (??) vector
  ## a and b: parameters of selectivity model
  ##---------------------------------------------------------------
  n.bins<-length(l.breaks)-1
  C<-matrix(0,n.bins,n.bins)
  for(i in 1: dim(C)[2]){
    ## selectivity at length
    s.length <- 1/(1+exp(-(a+b*l.mids[i])))   #trawl
    ## mortality at age
    M.age <- M[age.val]/n.seasons
    F <- s.length*Ffull/n.seasons #s.length
    Z <- M.age+F

    catch <- F/Z*(1-exp(-Z))
    #C[i,i] <- F/Z*(1-exp(-Z))
    if(is.na(catch)) C[i,i] <- 0 else {C[i,i] <- F/Z*(1-exp(-Z))}
  }
  return(C)
}
