get.survivorship.mat<-function(l.mids,l.breaks,age.val,Ffull,M,a,b,n.seasons){
  ##---------------------------------------------------------------
  ## function to create a survivorship matrix
  ## Arguments:
  ## l.mids: vector of midpoints of length classes
  ## l.breaks: vector of breakpoints of length classes
  ## age.val: age
  ## Ffull: fully selected F
  ## M: natural mortality at age vector
  ## a and b: parameters of selectivity model
  ##---------------------------------------------------------------  
  n.bins<-length(l.breaks)-1
  S<-matrix(0,n.bins,n.bins)
  for(i in 1: dim(S)[2]){
    ## selectivity at length
    s.length<-1/(1+exp(-(a+b*l.mids[i])))  #trawl
   ## mortality at age
    M.age<-M[age.val]/n.seasons
    F <- s.length*Ffull/n.seasons
    Z <- M.age+F
    S[i,i] <- exp(-Z)  
  }
  return(S)
}
