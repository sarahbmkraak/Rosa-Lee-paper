get.B<-function(N.array,n.years,n.ages,n.seasons,l.mids){
  ## convert numbers at length array to weight
  B<-array(0,c(n.bins,n.ages,n.seasons,n.years))
  for(i in 1:n.years){
    for(j in 1:n.ages){
      for(k in 1:n.seasons){
        B[,j,k,i]<-N.array[,j,k,i]*get.weight(l.mids,season=k, n.seasons=n.seasons)
      }
    }
  }
  return(B)
}
