get.SSB<-function(B.array,n.years,n.ages,n.seasons,l.mids){
  SSB<-array(0,c(n.bins,n.ages,n.seasons,n.years))
  for(i in 1:n.years){
    for(j in 1:n.ages){
      for(k in 1:n.seasons){
        SSB[,j,k,i]<-B.array[,j,k,i]*get.maturity_cod(l.mids,j)
      }
    }
  }
  return(SSB)
}
