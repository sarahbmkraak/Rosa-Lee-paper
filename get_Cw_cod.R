get.Cw<-function(C.array,n.years,n.ages,n.seasons,l.mids){
  Cw<-array(0,c(n.bins,n.ages,n.seasons,n.years))
  for(i in 1:n.years){
    print(i)
    for(j in 1:n.ages){
      for(k in 1:n.seasons){
        Cw[,j,k,i]<-C.array[,j,k,i]*get.weight(l.mids,season=k,n.seasons=n.seasons)
      }
    }
  }
  return(Cw)
}
