get.Dw<-function(Cw.array,n.years,n.ages,n.seasons,l.mids){
  Dw<-array(0,c(n.bins,n.ages,n.seasons,n.years))

  for(i in 1:n.years){
    m <- 35
    prop.discard<-ifelse(l.mids<m,1,0)
    #print(i)
    for(j in 1:n.ages){
      for(k in 1:n.seasons){
        Dw[,j,k,i]<-Cw.array[,j,k,i]*prop.discard
      }
    }
  }
  return(Dw)
}
