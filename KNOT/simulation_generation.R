library(dplyr)
library(MASS)
sparse.cor <- function(x){
  n <- nrow(x)
  cMeans <- colMeans(x)
  covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
  sdvec <- sqrt(diag(covmat))
  cormat <- covmat/tcrossprod(sdvec)
  list(cov=covmat,cor=cormat)
}

generate_sim_data<-function(effectsize, quan=FALSE, sigwin=5000,n = 10000, p = 1000, para=1,maxld=0.7){
  null=FALSE
  p0 <- 0.01
  if (null) maxld<-1
  CVorRV="cv"
  xchr=FALSE
  library(SKAT)
  data("SKAT.haplotypes")
  X0 = SKAT.haplotypes$Haplotype
  pos0 = SKAT.haplotypes$SNPInfo$CHROM_POS
  maf0<- apply(X0,2,sum)/(nrow(X0))
  X0 = X0[,maf0>0.01]
  pos0 = pos0[maf0>0.01]
  maf0 = maf0[maf0>0.01]
  corX0<-cor(X0)
  totalnsnp<-ncol(X0)
  if (CVorRV=='cvrv') {
    nsnp<-nsnp0<-500
    ncausal<-5
    effect<-0.02
  }
  if (CVorRV=="cv") {
    nsnp<-nsnp0<-p
    if (null) nsnp<-nsnp0<-200
    ncausal<-3 ###3
    effect<-0.015
  }
  if (CVorRV=="rv") {
    nsnp<-nsnp0<-400 #400
    ncausal<-4 #4
    effect<-0.015
  }
  
  Sigma.distance = as.dist(1 - abs(corX0))
  fit = hclust(Sigma.distance, method = "single")
  clusters = cutree(fit, h = 1 - maxld)
  df<-data.frame(ind=1:totalnsnp,clusters)
  rm(corX0)
  set.seed(para) #seed<-2
  tmp<-df %>% group_by(clusters) %>% sample_n(size = 1)
  nsnp<-min(max(clusters),nsnp0)
  snpindex<-sort(sample(tmp$ind,nsnp))
  X<-X0[,snpindex] #summary(as.vector(cor(X))) hist(as.vector(cor(X)))
  maf<-maf0[snpindex]
  pos<-pos0[snpindex]
  rm(X0)
  jth<-sort(sample(1:length(pos),ncausal))
  
  if (quan){
    if (!xchr){
      effect_cv<-effectsize
      effect_rv<-0.38
    } else{
      effect_cv<-0.3
      effect_rv<-0.55 #0.6
    }
  } else{
    if (!xchr){
      effect_cv<-effectsize #0.26
      effect_rv<-0.40 #0.48
    } else{
      effect_cv<-0.35 #0.4
      effect_rv<-0.7
    }
  }
  if (null) effect_cv<-effect_rv<-0
  #SKAT models
  if (CVorRV=="cv"){
    beta<-effect_cv*abs(log10(maf[jth]))
  }
  if (CVorRV=="rv") beta<-effect_rv*abs(log10(maf[jth]))
  if (CVorRV=="cvrv") {
    beta<-effect_rv*abs(log10(maf[jth]))
    ind<-which(maf[jth]>=0.01)
    if (length(ind)>0) beta[ind]<-effect_cv*abs(log10(maf[jth][ind]))
  }
  b0<-log(p0/(1-p0))
  ncase<-1
  dat<-array(dim=c(3*n,nsnp))
  available<-1:nrow(X)
  phasing.dad<-phasing.mom<-rep(1,n)
  ##
  if (!quan){
      dat.hap<-array(dim=c(4*n,nsnp))
      y <- rep(0, 3 * n)
      count_0 <- 0
      count_1 <- 0
      print(maf[jth])
        while (ncase<=n){
          tmp<-sample(available,4)
          hap.dad<-sample(c(1,2),1)
          hap.mom<-sample(c(1,2),1)
          Y<-X[tmp[hap.dad],]+X[tmp[2+hap.mom],]
          P1<-X[tmp[1],]+X[tmp[2],]
          P2<-X[tmp[3],]+X[tmp[4],]
          ##nonlinear
          potential_p1 = (beta %*%P1[jth])**2 + b0
          potential_p2 = (beta %*%P2[jth])**2 + b0
          potential_Y = (beta %*%Y[jth])**2 + b0
          y_p1 <- exp(potential_p1)/(1+exp(potential_p1))
          y_p2 <- exp(potential_p2)/(1+exp(potential_p2))
          y_Y <- exp(potential_Y)/(1+exp(potential_Y))
          if (runif(1)<y_Y & runif(1)>=y_p1 & runif(1)>=y_p2){
            # if (0.5<y_Y & 0.5>=y_p1 & 0.5>=y_p2){
            dat[ncase*3-2,]<-P1
            dat[ncase*3-1,]<-P2
            dat[ncase*3,]<-Y
            y[ncase*3] <- 1
            dat.hap[ncase * 4 - 3, ] <- X[tmp[1], ]
            dat.hap[ncase * 4 - 2, ] <- X[tmp[2], ]
            dat.hap[ncase * 4 - 1, ] <- X[tmp[3], ]
            dat.hap[ncase * 4, ] <- X[tmp[4], ]
            phasing.dad[ncase] <- hap.dad
            phasing.mom[ncase] <- hap.mom
            ncase <- ncase + 1
            count_1 <- count_1 + 1}
        }}else{
    dat.hap<-array(dim=c(4*n,nsnp))
    y<-rep(0,3*n)
    while (ncase<=n){
      tmp<-sample(available,4)
      hap.dad<-sample(c(1,2),1)
      hap.mom<-sample(c(1,2),1)
      P1<-X[tmp[1],]+X[tmp[2],]
      P2<-X[tmp[3],]+X[tmp[4],]
      Y<-X[tmp[hap.dad],]+X[tmp[2+hap.mom],]
      potential_p1 = (beta %*%P1[jth])**2+b0
      potential_p2 = (beta %*%P2[jth])**2+b0
      potential_Y = (beta %*%Y[jth])**2+b0
      y_p1 <- exp(potential_p1)/(1+exp(potential_p1))
      y_p2 <- exp(potential_p2)/(1+exp(potential_p2))
      y_Y <- exp(potential_Y)/(1+exp(potential_Y))
      if (runif(1)<y_Y & runif(1)>=y_p1 & runif(1)>=y_p2){
        y[ncase*3-2]<-potential_p1+rnorm(1,0,2)-b0
        y[ncase*3-1]<-potential_p2+rnorm(1,0,2)-b0
        y[ncase*3]<-potential_Y+rnorm(1,0,2)-b0
        dat[ncase*3-2,]<-P1
        dat[ncase*3-1,]<-P2
        dat[ncase*3,]<-Y
        dat.hap[ncase*4-3,]<-X[tmp[1],]
        dat.hap[ncase*4-2,]<-X[tmp[2],]
        dat.hap[ncase*4-1,]<-X[tmp[3],]
        dat.hap[ncase*4,]<-X[tmp[4],]
        phasing.dad[ncase]<-hap.dad
        phasing.mom[ncase]<-hap.mom
        ncase<-ncase+1
      }
  }}
  
  #remove ld>maxld
  if (!null){
    corX0<-sparse.cor(dat[-seq(3,(3*n),3),])$cor
    corX0[is.na(corX0)] <- 0
    Sigma.distance = as.dist(1 - abs(corX0))
    fit = hclust(Sigma.distance, method = "single")
    clusters = cutree(fit, h = 1 - maxld)
    df<-data.frame(ind=1:nsnp,clusters)
    rm(corX0)
    index<-setdiff(which(df$clusters %in% df$clusters[jth]),jth)
    if (length(index)>0){
      dat<-dat[,-index]
      dat.hap<-dat.hap[,-index]
      pos<-pos[-index]
      jth<-match(jth,(1:nsnp)[-index])
      nsnp<-ncol(dat)
    }
  }
  out<-list()
  out$dat<-dat
  out$dat.hap<-dat.hap
  out$pos<-pos
  if (!null) out$pos_causal<-pos[jth] else out$pos_causal<-NULL
  if (xchr) out$sex<-sex
  out$y<-as.numeric(y)
  out$phasing.dad<-phasing.dad
  out$phasing.mom<-phasing.mom
  return(out)
}