MK.q.byStat<-function (FIs, M){
  kappa <- apply(FIs, 2, function(x) {
    if (all(x == 0)) {
      return(1)
    } else {
      return(which.max(x) - 1)
    }
  })
  tau <- apply(FIs, 2, function(i) max(i) - median(i[-which.max(i)]))
  W <- apply(FIs, 2, function(i) (i[1] - median(i[2:(M+1)])) * ifelse(i[1] >= max(i[2:(M+1)]), 1, 0))
  Rej.Bound=10000
  b<-order(tau,decreasing=T)
  c_0<-kappa[b]==0
  #calculate ratios for top Rej.Bound tau values
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    #if(i==1){temp_0=c_0[i]}
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  #calculate q values for top Rej.Bound values
  q<-rep(1,length(tau))
  for(i in 1:length(b)){
    q[b[i]]<-min(ratio[i:min(length(b),Rej.Bound)])*c_0[i]+1-c_0[i]
    if(i>Rej.Bound){break}
  }
  return(q)
}
