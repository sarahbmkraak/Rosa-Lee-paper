get.Lw<-function(Cw.array,n.years,n.ages,n.seasons,l.mids){
  Lw<-array(0,c(n.bins,n.ages,n.seasons,n.years))
   #MRLS is 35 cm for cod in the Baltic
  for(i in 1:n.years){
    m <- 35
    prop.land<-ifelse(l.mids>=m,1,0)
    print(i)
    for(j in 1:n.ages){
      for(k in 1:n.seasons){
        Lw[,j,k,i]<-Cw.array[,j,k,i]*prop.land
      }
    }
  }
  return(Lw)
}
