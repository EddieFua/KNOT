library(dplyr)
library(MASS)
ACAT<-function(p){
  p[p>0.99]<-0.99 #p[p>1-1e-2]<-1-1e-2
  is.small<-(p<1e-16) & !is.na(p)
  is.regular<-(p>=1e-16) & !is.na(p)
  temp<-rep(NA,length(p))
  temp[is.small]<-1/p[is.small]/pi
  temp[is.regular]<-as.numeric(tan((0.5-p[is.regular])*pi))
  
  cct.stat<-mean(temp,na.rm=T)
  if(is.na(cct.stat)){return(NA)}
  if(cct.stat>1e+15){return((1/cct.stat)/pi)}else{
    return(1-pcauchy(cct.stat))
  }
}
Compute_expected<-function(dat){
  n <- nrow(dat)/3
  index_dad <- seq(1, (3 * n), 3)
  index_mom <- seq(2, (3 * n), 3)
  index_off <- seq(3, (3 * n), 3)
  parent_matrix <- dat[-index_off, ]
  expected_matrix <- dat[index_off, ]
  for (i in 1:n){
    expected_matrix[i,] = apply(parent_matrix[c(i*2-1,i*2),],2,mean)
  }
  return(expected_matrix)
}
Compute_U<-function(dat){
  n <- nrow(dat)/3
  index_dad <- seq(1, (3 * n), 3)
  index_mom <- seq(2, (3 * n), 3)
  index_off <- seq(3, (3 * n), 3)
  U_matrix <- dat[index_off, ]
  expected_matrix <- dat[-index_off, ]
  for (i in 1:n){
    expected_matrix[i,] = apply(expected_matrix[c(i*2-1,i*2),],2,mean)
    U_matrix[i,] <- U_matrix[i,]-expected_matrix[i,]
  }
  return(U_matrix)
}
MK.statistic<-function (T_0,T_k,method='median'){
  T_0<-as.matrix(T_0);T_k<-as.matrix(T_k)
  T.temp<-cbind(T_0,T_k)
  T.temp[is.na(T.temp)]<-0
  kappa<-apply(T.temp,1,which.max)-1
  if(method=='max'){tau<-apply(T.temp,1,max)-apply(T.temp,1,max.nth,n=2)}
  if(method=='median'){
    Get.OtherMedian<-function(x){median(x[-which.max(x)])}
    tau<-apply(T.temp,1,max)-apply(T.temp,1,Get.OtherMedian)
  }
  return(cbind(kappa,tau))
}

sparse.cor <- function(x){
  n <- nrow(x)
  cMeans <- colMeans(x)
  covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
  sdvec <- sqrt(diag(covmat))
  cormat <- covmat/tcrossprod(sdvec)
  list(cov=covmat,cor=cormat)
}

sparse.cov.cross <- function(x,y){
  n <- nrow(x)
  cMeans.x <- colMeans(x);cMeans.y <- colMeans(y)
  covmat <- (as.matrix(crossprod(x,y)) - n*tcrossprod(cMeans.x,cMeans.y))/(n-1)
  list(cov=covmat)
}

create.MK<-function (X, pos, M = 5, corr_max = 0.75, info,info2, maxN.neighbor = Inf, 
                     maxBP.neighbor = 100000) {
  X <- as.matrix(X)
  sparse.fit <- sparse.cor(X)
  cor.X <- sparse.fit$cor
  cor.X[is.na(cor.X)] <- 0
  cor.X[is.infinite(cor.X)] <- 0
  cov.X <- sparse.fit$cov
  cov.X[is.na(cov.X)] <- 0
  cov.X[is.infinite(cov.X)] <- 0
  Sigma.distance = as.dist(1 - abs(cor.X))
  if (ncol(X) > 1) {
    fit = hclust(Sigma.distance, method = "single")
    corr_max = corr_max
    clusters = cutree(fit, h = 1 - corr_max)
  }
  else {
    clusters <- 1
  }
  X_k <- array(0, dim = c(nrow(X), ncol(X), M))
  index.exist <- c()
  for (k in unique(clusters)) {
    #cluster.fitted <- cluster.residuals <- matrix(NA, nrow(X),sum(clusters == k))
    ind_clusterk<-which(clusters == k)
    for (i in ind_clusterk) {
      index.pos <- which(pos >= max(pos[i] - maxBP.neighbor, 
                                    pos[1]) & pos <= min(pos[i] + maxBP.neighbor, 
                                                         pos[length(pos)]))
      temp <- abs(cor.X[i, ])
      temp[ind_clusterk] <- 0
      temp[-index.pos] <- 0
      index <- order(temp, decreasing = T)
      index <- setdiff(index[1:min(length(index), sum(temp > 
                                                        0.05), floor((nrow(X))^(1/3)), maxN.neighbor)], 
                       i)
      y <- X[, i]
      if (length(index) == 0) {
        fitted.values <- mean(y)
      }
      else {
        x <- X[, index, drop = F]
        temp.xy <- rbind(mean(y), crossprod(x, y)/length(y) - 
                           colMeans(x) * mean(y))
        x.exist <- c()
        for (j in 1:M) {
          x.exist <- cbind(x.exist, X_k[ , intersect(index, 
                                                     index.exist),j])
        }
        temp.xy <- rbind(temp.xy, crossprod(x.exist, 
                                            y)/length(y) - colMeans(x.exist) * mean(y))
        temp.cov.cross <- sparse.cov.cross(x, x.exist)$cov
        temp.cov <- sparse.cor(x.exist)$cov
        temp.xx <- cov.X[index, index]
        temp.xx <- rbind(cbind(temp.xx, temp.cov.cross), 
                         cbind(t(temp.cov.cross), temp.cov))
        temp.xx <- cbind(0, temp.xx)
        temp.xx <- rbind(c(1, rep(0, ncol(temp.xx) - 
                                    1)), temp.xx)
        pca.fit <- princomp(covmat = temp.xx)
        v <- pca.fit$loadings
        cump <- cumsum(pca.fit$sdev^2)/sum(pca.fit$sdev^2)
        n.pc <- which(cump >= 0.999)[1]
        pca.index <- intersect(1:n.pc, which(pca.fit$sdev != 
                                               0))
        temp.inv <- v[, pca.index, drop = F] %*% (pca.fit$sdev[pca.index]^(-2) * 
                                                    t(v[, pca.index, drop = F]))
        temp.beta <- temp.inv %*% temp.xy
        temp.j <- 1
        fitted.values <- temp.beta[1] + crossprod(t(x), 
                                                  temp.beta[(temp.j + 1):(temp.j + ncol(x)), 
                                                            , drop = F]) - sum(colMeans(x) * temp.beta[(temp.j + 
                                                                                                          1):(temp.j + ncol(x)), , drop = F])
        temp.j <- temp.j + ncol(x)
        for (j in 1:M) {
          temp.x <- as.matrix(X_k[, intersect(index, 
                                              index.exist),j])
          if (ncol(temp.x) >= 1) {
            fitted.values <- fitted.values + crossprod(t(temp.x), 
                                                       temp.beta[(temp.j + 1):(temp.j + ncol(temp.x)), 
                                                                 , drop = F]) - sum(colMeans(temp.x) * 
                                                                                      temp.beta[(temp.j + 1):(temp.j + ncol(temp.x)), 
                                                                                                , drop = F])
          }
          temp.j <- temp.j + ncol(temp.x)
        }
      }
      residuals <- y - fitted.values
      X_k[,i,] <- as.vector(fitted.values) + sapply(1:M, 
                                                    function(x) {
                                                      residuals2<-residuals
                                                      if (length(info[[i]])>1) residuals2[info[[i]]]<-residuals2[sample(info[[i]])]
                                                      if (length(info2[[i]])>1) residuals2[info2[[i]]]<-residuals2[sample(info2[[i]])]
                                                      return(residuals2)
                                                    }) #sample(residuals)
      #cluster.fitted[, match(i, which(clusters == k))] <- as.vector(fitted.values)
      #cluster.residuals[, match(i, which(clusters == k))] <- as.vector(residuals)
      index.exist <- c(index.exist, i)
    }
    # cluster.sample.index <- sapply(1:M, function(x) sample(1:nrow(X)))
    # for (j in 1:M) {
    #   X_k[j, , which(clusters == k)] <- cluster.fitted + 
    #     cluster.residuals[cluster.sample.index[, j], 
    #                       , drop = F]
    # }
  }
  return(X_k)
}

create.MK.hap<-function (X, pos, M = 5, corr_max = 0.75, info,info2, info3, maxN.neighbor = Inf, 
                         maxBP.neighbor = 100000) {
  X <- as.matrix(X)
  sparse.fit <- sparse.cor(X)
  cor.X <- sparse.fit$cor
  cor.X[is.na(cor.X)] <- 0
  cor.X[is.infinite(cor.X)] <- 0
  cov.X <- sparse.fit$cov
  cov.X[is.na(cov.X)] <- 0
  cov.X[is.infinite(cov.X)] <- 0
  Sigma.distance = as.dist(1 - abs(cor.X))
  if (ncol(X) > 1) {
    fit = hclust(Sigma.distance, method = "single")
    corr_max = corr_max
    clusters = cutree(fit, h = 1 - corr_max)
  }
  else {
    clusters <- 1
  }
  X_k <- array(0, dim = c(nrow(X), ncol(X), M))
  index.exist <- c()
  for (k in unique(clusters)) {
    #cluster.fitted <- cluster.residuals <- matrix(NA, nrow(X),sum(clusters == k))
    ind_clusterk<-which(clusters == k)
    for (i in ind_clusterk) {
      index.pos <- which(pos >= max(pos[i] - maxBP.neighbor, 
                                    pos[1]) & pos <= min(pos[i] + maxBP.neighbor, 
                                                         pos[length(pos)]))
      temp <- abs(cor.X[i, ])
      temp[ind_clusterk] <- 0
      temp[-index.pos] <- 0
      index <- order(temp, decreasing = T)
      index <- setdiff(index[1:min(length(index), sum(temp > 
                                                        0.05), floor((nrow(X))^(1/3)), maxN.neighbor)], 
                       i)
      y <- X[, i]
      if (length(index) == 0) {
        fitted.values <- mean(y)
      }
      else {
        x <- X[, index, drop = F]
        temp.xy <- rbind(mean(y), crossprod(x, y)/length(y) - 
                           colMeans(x) * mean(y))
        x.exist <- c()
        for (j in 1:M) {
          x.exist <- cbind(x.exist, X_k[ , intersect(index, 
                                                     index.exist),j])
        }
        temp.xy <- rbind(temp.xy, crossprod(x.exist, 
                                            y)/length(y) - colMeans(x.exist) * mean(y))
        temp.cov.cross <- sparse.cov.cross(x, x.exist)$cov
        temp.cov <- sparse.cor(x.exist)$cov
        temp.xx <- cov.X[index, index]
        temp.xx <- rbind(cbind(temp.xx, temp.cov.cross), 
                         cbind(t(temp.cov.cross), temp.cov))
        temp.xx <- cbind(0, temp.xx)
        temp.xx <- rbind(c(1, rep(0, ncol(temp.xx) - 
                                    1)), temp.xx)
        pca.fit <- princomp(covmat = temp.xx)
        v <- pca.fit$loadings
        cump <- cumsum(pca.fit$sdev^2)/sum(pca.fit$sdev^2)
        n.pc <- which(cump >= 0.999)[1]
        pca.index <- intersect(1:n.pc, which(pca.fit$sdev != 
                                               0))
        temp.inv <- v[, pca.index, drop = F] %*% (pca.fit$sdev[pca.index]^(-2) * 
                                                    t(v[, pca.index, drop = F]))
        temp.beta <- temp.inv %*% temp.xy
        temp.j <- 1
        fitted.values <- temp.beta[1] + crossprod(t(x), 
                                                  temp.beta[(temp.j + 1):(temp.j + ncol(x)), 
                                                            , drop = F]) - sum(colMeans(x) * temp.beta[(temp.j + 
                                                                                                          1):(temp.j + ncol(x)), , drop = F])
        temp.j <- temp.j + ncol(x)
        for (j in 1:M) {
          temp.x <- as.matrix(X_k[, intersect(index, 
                                              index.exist),j])
          if (ncol(temp.x) >= 1) {
            fitted.values <- fitted.values + crossprod(t(temp.x), 
                                                       temp.beta[(temp.j + 1):(temp.j + ncol(temp.x)), 
                                                                 , drop = F]) - sum(colMeans(temp.x) * 
                                                                                      temp.beta[(temp.j + 1):(temp.j + ncol(temp.x)), 
                                                                                                , drop = F])
          }
          temp.j <- temp.j + ncol(temp.x)
        }
      }
      residuals <- y - fitted.values
      X_k[,i,] <- as.vector(fitted.values) + sapply(1:M, 
                                                    function(x) {
                                                      residuals2<-residuals
                                                      if (length(info[[i]])>1) residuals2[info[[i]]]<-residuals2[sample(info[[i]])]
                                                      if (length(info2[[i]])>1) residuals2[info2[[i]]]<-residuals2[sample(info2[[i]])]
                                                      if (length(info3[[i]])>1) residuals2[info3[[i]]]<-residuals2[sample(info3[[i]])]
                                                      return(residuals2)
                                                    }) #sample(residuals)
      #cluster.fitted[, match(i, which(clusters == k))] <- as.vector(fitted.values)
      #cluster.residuals[, match(i, which(clusters == k))] <- as.vector(residuals)
      index.exist <- c(index.exist, i)
    }
    # cluster.sample.index <- sapply(1:M, function(x) sample(1:nrow(X)))
    # for (j in 1:M) {
    #   X_k[j, , which(clusters == k)] <- cluster.fitted + 
    #     cluster.residuals[cluster.sample.index[, j], 
    #                       , drop = F]
    # }
  }
  return(X_k)
}

create.MK.original.new <- function(X,pos,M=5,corr_max=0.75,maxN.neighbor=Inf,maxBP.neighbor=100000) {
  
  X <- as.matrix(X)
  sparse.fit <- sparse.cor(X)
  cor.X <- sparse.fit$cor
  cor.X[is.na(cor.X)] <- 0
  cor.X[is.infinite(cor.X)] <- 0
  cov.X <- sparse.fit$cov
  cov.X[is.na(cov.X)] <- 0
  cov.X[is.infinite(cov.X)] <- 0
  Sigma.distance = as.dist(1 - abs(cor.X))
  if(ncol(X)>1){
    fit = hclust(Sigma.distance, method="single")
    corr_max = corr_max
    clusters = cutree(fit, h=1-corr_max)
  }else{clusters<-1}
  
  #temp.M<-matrix(1,M+1,M+1)
  #cov.M<-kronecker(temp.M,cov.X)
  
  #X_k<-matrix(0,nrow(X),ncol(X));index.exist<-c()
  X_k<-array(0,dim=c(nrow(X),ncol(X),M));index.exist<-c()
  for (k in unique(clusters)){
    cluster.fitted<-cluster.residuals<-matrix(NA,nrow(X),sum(clusters==k))
    for(i in which(clusters==k)){
      #print(i)
      index.pos<-which(pos>=max(pos[i]-maxBP.neighbor,pos[1]) & pos<=min(pos[i]+maxBP.neighbor,pos[length(pos)]))
      temp<-abs(cor.X[i,]);
      temp[which(clusters==k)]<-0; #only keep neighboring SNPs that are not in the same cluster (i.e., cor < maxcor)
      temp[-index.pos]<-0
      
      #keep neighboring SNPs that have an absolute correlation between maxcor and 0.05
      index<-order(temp,decreasing=T)
      index<-setdiff(index[1:min(length(index),sum(temp>0.05),floor((nrow(X))^(1/3)),maxN.neighbor)],i)
      print(length(index))
      y<-X[,i]
      if(length(index)==0){fitted.values<-mean(y)}else{
        
        x<-X[,index,drop=F];
        temp.xy<-rbind(mean(y),crossprod(x,y)/length(y)-colMeans(x)*mean(y))
        x.exist<-c()
        for(j in 1:M){
          x.exist<-cbind(x.exist,X_k[,intersect(index,index.exist),j])
        }
        temp.xy<-rbind(temp.xy,crossprod(x.exist,y)/length(y)-colMeans(x.exist)*mean(y))
        
        temp.cov.cross<-sparse.cov.cross(x,x.exist)$cov
        temp.cov<-sparse.cor(x.exist)$cov
        temp.xx<-cov.X[index,index]
        temp.xx<-rbind(cbind(temp.xx,temp.cov.cross),cbind(t(temp.cov.cross),temp.cov))
        
        temp.xx<-cbind(0,temp.xx)
        temp.xx<-rbind(c(1,rep(0,ncol(temp.xx)-1)),temp.xx)
        
        pca.fit<-princomp(covmat=temp.xx)
        v<-pca.fit$loadings
        cump<-cumsum(pca.fit$sdev^2)/sum(pca.fit$sdev^2)
        n.pc<-which(cump>=0.999)[1]#nrow(temp.xx)#nrow(temp.xx)#
        pca.index<-intersect(1:n.pc,which(pca.fit$sdev!=0))#which(cump<=0.99)
        #calculate
        #inverse ZZ matrix
        temp.inv<-v[,pca.index,drop=F]%*%(pca.fit$sdev[pca.index]^(-2)*t(v[,pca.index,drop=F]))
        #beta coefficients
        temp.beta<-temp.inv%*%temp.xy
        
        temp.j<-1
        fitted.values<-temp.beta[1]+crossprod(t(x),temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F])-sum(colMeans(x)*temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F])
        temp.j<-temp.j+ncol(x)
        for(j in 1:M){
          temp.x<-as.matrix(X_k[,intersect(index,index.exist),j])
          if(ncol(temp.x)>=1){
            fitted.values<-fitted.values+crossprod(t(temp.x),temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),,drop=F])-sum(colMeans(temp.x)*temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),,drop=F])
          }
          temp.j<-temp.j+ncol(temp.x)
        }
      }
      residuals<-y-fitted.values
      cluster.fitted[,match(i,which(clusters==k))]<-as.vector(fitted.values)
      cluster.residuals[,match(i,which(clusters==k))]<-as.vector(residuals)
      
      index.exist<-c(index.exist,i)
    }
    #sample mutiple knockoffs
    cluster.sample.index<-sapply(1:M,function(x)sample(1:nrow(X)))
    for(j in 1:M){
      X_k[,which(clusters==k),j]<-cluster.fitted+cluster.residuals[cluster.sample.index[,j],,drop=F]
    }
  }
  return(X_k)
}

getslidingwindow<-function(left,right,size,pos){
  start<-end<-n<-ind3<-ind4<-c()
  step<-round(size/2) #move half of the window size each time
  l<-left-step
  r<-l+size-1
  while (r<=(right+step)){
    start<-c(start,l)
    end<-c(end,r)
    tmp<-which(pos>=l & pos<=r)
    n<-c(n,length(tmp))
    if (length(tmp)==0) tmp<-0
    ind3<-c(ind3,tmp[1])
    ind4<-c(ind4,tmp[length(tmp)])
    l<-l+step
    r<-r+step
  }
  return(data.frame(start,end,ind3,ind4,n))
}

pcontribution <- function(dat,dat1=NA,dosage=FALSE,xchr=FALSE,sex=NA){
  n.row <- nrow(dat)
  dad <- dat[seq.int(1, n.row, 3),, drop=FALSE]
  mom <- dat[seq.int(2, n.row, 3),, drop=FALSE]
  kid <- dat[seq.int(3, n.row, 3),, drop=FALSE]
  if (dosage==FALSE){
    Z.pa<-Z.ma<-array(0,dim=c(n.row/3,ncol(dat)))
    het <- (mom == 1)
    hethom <- het & (dad == 2)
    # ind212<-hethom & (kid == 2)
    # ind211<-hethom & (kid == 1)
    # n212 <- colSums(ind212, na.rm=TRUE)
    # n211 <- colSums(ind211, na.rm=TRUE)
    Z.ma[hethom & (kid == 2)]<-1
    Z.ma[hethom & (kid == 1)]<-(-1)
    hethom <- het & (dad == 0)
    # ind011<-hethom & (kid == 1)
    # ind010<-hethom & (kid == 0)
    # n011 <- colSums(ind011, na.rm=TRUE)
    # n010 <- colSums(ind010, na.rm=TRUE)
    Z.ma[hethom & (kid == 1)]<-1
    Z.ma[hethom & (kid == 0)]<-(-1)
    het <- (dad == 1)
    hethom <- het & (mom == 2)
    # ind122<-hethom & (kid == 2)
    # ind121<-hethom & (kid == 1)
    # n122 <- colSums(ind122, na.rm=TRUE)
    # n121 <- colSums(ind121, na.rm=TRUE)
    Z.pa[hethom & (kid == 2)]<-1
    Z.pa[hethom & (kid == 1)]<-(-1)
    hethom <- het & (mom == 0)
    # ind101<-hethom & (kid == 1)
    # ind100<-hethom & (kid == 0)
    # n101 <- colSums(ind101, na.rm=TRUE)
    # n100 <- colSums(ind100, na.rm=TRUE)
    Z.pa[hethom & (kid == 1)]<-1
    Z.pa[hethom & (kid == 0)]<-(-1)
    
    het <- (mom == 1) & (dad == 1)
    ind112<-het & (kid == 2)
    ind110<-het & (kid == 0)
    # n112 <- colSums(ind112, na.rm=TRUE)
    # n110 <- colSums(ind110, na.rm=TRUE)
    Z.pa[ind112]<-Z.ma[ind112]<-1
    Z.pa[ind110]<-Z.ma[ind110]<-(-1)
    #mat<-cbind(n100,n110,n121,n101,n112,n122,n010,n110,n211,n011,n112,n212)
    #colnames(mat)<-c("p-1","p-1","p-1","p+1","p+1","p+1","m-1","m-1","m-1","m+1","m+1","m+1") #contributions
    
    if (xchr){
      ind <- (dad == 1) & (mom == 1) & (kid==1)
      ind0<-sweep(ind,1,sex==0,'&')
      Z.pa[ind0]<-1
      Z.ma[ind0]<-(-1)
      ind0<-sweep(ind,1,sex==1,'&')
      Z.pa[ind0]<-(-1)
      Z.ma[ind0]<-1
    }
    
    out<-list()
    out$Z.pa<-Z.pa
    out$Z.ma<-Z.ma
    return(out)
  }else{
    dims<-dim(dat1)
    dims[1]<-dims[1]/3
    Z.pa<-Z.ma<-array(0,dim=dims)
    dadk <- dat1[seq.int(1, n.row, 3),,,drop=FALSE]
    momk <- dat1[seq.int(2, n.row, 3),,,drop=FALSE]
    kidk <- dat1[seq.int(3, n.row, 3),,,drop=FALSE]
    
    het <- (mom == 1)
    hethom <- het & (dad == 2)
    ind212<-hethom & (kid == 2)
    ind211<-hethom & (kid == 1)
    # n212 <- colSums(ind212, na.rm=TRUE)
    # n211 <- colSums(ind211, na.rm=TRUE)
    Z.pa[ind212]<-kidk[ind212]-dadk[ind212]
    Z.ma[ind212]<-kidk[ind212]-momk[ind212]
    Z.pa[ind211]<-kidk[ind211]-momk[ind211]
    Z.ma[ind211]<-kidk[ind211]-dadk[ind211]
    
    hethom <- het & (dad == 0)
    ind011<-hethom & (kid == 1)
    ind010<-hethom & (kid == 0)
    # n011 <- colSums(ind011, na.rm=TRUE)
    # n010 <- colSums(ind010, na.rm=TRUE)
    Z.pa[ind011]<-kidk[ind011]-momk[ind011]
    Z.ma[ind011]<-kidk[ind011]-dadk[ind011]
    Z.pa[ind010]<-kidk[ind010]-dadk[ind010]
    Z.ma[ind010]<-kidk[ind010]-momk[ind010]
    
    het <- (dad == 1)
    hethom <- het & (mom == 2)
    ind122<-hethom & (kid == 2)
    ind121<-hethom & (kid == 1)
    # n122 <- colSums(ind122, na.rm=TRUE)
    # n121 <- colSums(ind121, na.rm=TRUE)
    Z.pa[ind122]<-kidk[ind122]-dadk[ind122]
    Z.ma[ind122]<-kidk[ind122]-momk[ind122]
    Z.pa[ind121]<-kidk[ind121]-momk[ind121]
    Z.ma[ind121]<-kidk[ind121]-dadk[ind121]
    
    hethom <- het & (mom == 0)
    ind101<-hethom & (kid == 1)
    ind100<-hethom & (kid == 0)
    # n101 <- colSums(ind101, na.rm=TRUE)
    # n100 <- colSums(ind100, na.rm=TRUE)
    Z.pa[ind101]<-kidk[ind101]-momk[ind101]
    Z.ma[ind101]<-kidk[ind101]-dadk[ind101]
    Z.pa[ind100]<-kidk[ind100]-dadk[ind100]
    Z.ma[ind100]<-kidk[ind100]-momk[ind100]
    
    het <- (mom == 1) & (dad == 1)
    ind112<-het & (kid == 2)
    ind110<-het & (kid == 0)
    # n112 <- colSums(ind112, na.rm=TRUE)
    # n110 <- colSums(ind110, na.rm=TRUE)
    Z.pa[ind112]<-kidk[ind112]-momk[ind112]
    Z.ma[ind112]<-kidk[ind112]-dadk[ind112]
    Z.pa[ind110]<-kidk[ind110]-dadk[ind110]
    Z.ma[ind110]<-kidk[ind110]-momk[ind110]
    
    if (xchr){
      ind <- (dad == 1) & (mom == 1) & (kid==1)
      ind0<-sweep(ind,1,sex==0,'&')
      Z.pa[ind0]<-dadk[ind0]
      Z.ma[ind0]<-(-1)*momk[ind0]
      ind0<-sweep(ind,1,sex==1,'&')
      Z.pa[ind0]<-(-1)*dadk[ind0]
      Z.ma[ind0]<-momk[ind0]
    }
    
    #mat<-cbind(n100,n110,n121,n101,n112,n122,n010,n110,n211,n011,n112,n212)
    #colnames(mat)<-c("p-1","p-1","p-1","p+1","p+1","p+1","m-1","m-1","m-1","m+1","m+1","m+1") #contributions
    out<-list()
    out$Z.pa<-Z.pa
    out$Z.ma<-Z.ma
    return(out)
  }
}

xcontribution <- function(dat,dat1=NA,dosage=FALSE,sex=NA){
  #sex: gender info for the offspring, 0 for females and 1 for males
  #Only knockoffs of mothers are required
  n.row <- nrow(dat)
  dad <- dat[seq.int(1, n.row, 3),, drop=FALSE]
  mom <- dat[seq.int(2, n.row, 3),, drop=FALSE]
  kid <- dat[seq.int(3, n.row, 3),, drop=FALSE]
  if (dosage==FALSE){
    Z<-array(0,dim=c(n.row/3,ncol(dat)))
    het <- (mom == 1)
    hethom <- het & (dad == 0)
    Z[hethom & (kid == 0)]<-(-0.5)
    Z[hethom & (kid == 1)]<-0.5
    hethom <- het & (dad == 1)
    ind111<-hethom & (kid == 1)
    Z[hethom & (kid == 0)]<-(-0.5) #110 trios on chrX can only have male offspring
    Z[sweep(ind111,1,sex==0,'&')]<-(-0.5)
    Z[sweep(ind111,1,sex==1,'&')]<-0.5
    Z[hethom & (kid == 2)]<-0.5 #112 trios on chrX can only have female offspring
  }else{
    dims<-dim(dat1)
    dims[1]<-dims[1]/3
    Z<-array(0,dim=dims)
    dadk <- dat1[seq.int(1, n.row, 3),,,drop=FALSE]
    momk <- dat1[seq.int(2, n.row, 3),,,drop=FALSE]
    kidk <- dat1[seq.int(3, n.row, 3),,,drop=FALSE]
    
    het<-dad==0 & mom==1
    Z[het]<-kidk[het]-0.5*(dadk[het]+momk[het])
    #Z[het]<-kidk[het]-0.5*pmax(dadk[het],momk[het])-pmin(dadk[het],momk[het])
    #Z[het]<-kidk[het]-0.5*pmax(dadk[het],momk[het])
    #Z[het]<-kidk[het]-0.5*momk[het]
    het<-dad==1 & mom==1
    het0<-sweep(het,1,sex==0,'&')
    Z[het0]<-kidk[het0]-0.75*(dadk[het0]+momk[het0])
    #Z[het0]<-kidk[het0]-0.5*(pmin(dadk[het0],momk[het0])+dadk[het0]+momk[het0])
    #Z[het0]<-kidk[het0]-pmin(dadk[het0],momk[het0])-0.5*pmax(dadk[het0],momk[het0])
    het0<-sweep(het,1,sex==1,'&')
    Z[het0]<-kidk[het0]-0.25*(dadk[het0]+momk[het0])
    #Z[het0]<-kidk[het0]-0.5*pmax(dadk[het0],momk[het0])
  }
  return(Z)
}

fbat<-function(dat,adjust_for_covariates=FALSE,y=NA,dosage=FALSE,dat1=NA,xchr=FALSE,sex=NA){
  #if (xchr & is.na(sex)) stop("Gender information is required if FBAT is applied to the X chromosome")
  if (is.null(ncol(dat))) dat<-as.matrix(dat)
  n<-nrow(dat)/3
  index_dad<-seq(1,(3*n),3) #index for parent1
  index_mom<-seq(2,(3*n),3) #index for parent2
  index_off<-seq(3,(3*n),3) #index for offspring
  if (!xchr){
    if (dosage==FALSE) Z<-dat[index_off,,drop=F]-(dat[index_dad,,drop=F]+dat[index_mom,,drop=F])/2 else Z<-dat1[index_off,,,drop=FALSE]-(dat1[index_dad,,,drop=FALSE]+dat1[index_mom,,,drop=FALSE])/2
  } else Z<-xcontribution(dat,dat1,dosage,sex)
  tat<-pcontribution(dat,dat1,dosage,xchr,sex)
  Z.pa<-tat$Z.pa
  Z.ma<-tat$Z.ma
  if (adjust_for_covariates){
    Z<-sweep(Z,1,y,'*')
    Z.pa<-sweep(Z.pa,1,y,'*')
    Z.ma<-sweep(Z.ma,1,y,'*')
  }
  additive<-colSums(Z)/sqrt(colSums(Z^2)) #if dosage=T, dim=c(nsnp,M)
  paternal<-colSums(Z.pa)/sqrt(colSums(Z.pa^2)) #if dosage=T, dim=c(nsnp,M)
  maternal<-colSums(Z.ma)/sqrt(colSums(Z.ma^2)) #if dosage=T, dim=c(nsnp,M)
  additive[is.na(additive)]<-0
  paternal[is.na(paternal)]<-0
  maternal[is.na(maternal)]<-0
  out<-list()
  out$additive<-additive
  out$paternal<-paternal
  out$maternal<-maternal
  out$Z<-Z
  out$Z.pa<-Z.pa
  out$Z.ma<-Z.ma
  return (out)
}

# res_fact = fbat(sim$dat,adjust_for_covariates=FALSE,y=NA,dosage=FALSE,dat1=NA,xchr=FALSE,sex=NA)
# abs(res_fbat$additive)[order(abs(res_fbat$additive),decreasing = T)]
# 2*pnorm(-abs(res_fbat$additive))[c(102,103,104)]

fbat_set<-function(W,V,ind){
  if (length(ind)>0){
    s1<-sum(W[ind])
    s2<-sum(V[ind,ind])
    if (s2==0) s2<-1
    return (unname(s1/sqrt(s2)))
  } else return(NA)
}

calculate_w_kappatau<-function(q1,q2){
  out<-list()
  t1<--log10(q1)
  t2<--log10(q2)
  t2_med<-apply(t2,2,median)
  t2_max<-apply(t2,2,max)
  out$w<-(t1-t2_med)*(t1>=t2_max)
  out$w.raw<-t1-t2_med
  out$kappatau<-MK.statistic(t1,t(t2),method="median")
  #out$q<-MK.q.byStat(out$kappatau[,1],out$kappatau[,2],M=M)
  return(out)
}

get_residuals<-function(dat,sex){
  n<-nrow(dat)/3 #number of trios
  y<-rep(c(1,0),each=n)
  gender<-c(sex,rep(1,n/2),rep(0,n/2))
  df_gender<-data.frame(y,gender)
  lmout<-lm(y~gender,df_gender)
  slm<-summary(lmout)
  y<-slm$residuals[1:n] #residuals after regressing out gender, 0.366 for males, 0.787 for females
  return(y)
}

quality_control<-function(g0,g0.hap,error_number=1,missing_rate=0.05,xchr=FALSE,sex=NA,hwe=FALSE,hwe.threshold=1e-7,gender_dependent=TRUE,gender_dependent_thr=0.05/789504,maf.thr=0){
  #remove variants that have more than 1 trio with Mendelian error or have >5% missing trios
  #threshold: remove snps that have >threshold proportion of trios that are incomplete or have mendelian errors
  #replace the remaining incomplete or erroneous trios with 0 0 0
  #gender_dependent_thr=0.05/789504=0.05/# of non-chrX variants
  #g0 must be a matrix type
  g0[!g0 %in% c(0,1,2)] <- NA
  n<-nrow(g0)/2 #number of trios
  index_dad<-seq(1,(2*n),2)
  index_mom<-seq(2,(2*n),2)
  
  #recode 2 to 1 for males genotpye
  if (xchr){
    index_son<-which(sex==1)*3 #index for male offspring
    #index_dau<-which(sex==0)*3
    ind<-which(g0==2,arr.ind=TRUE)
    ind<-ind[which(ind[,1]%%3==1 | ind[,1] %in% index_son),] #index for males who are coded 2 for # of minor alleles on X
    g0[ind]<-1
    # ind<-which(g0==1,arr.ind=TRUE)
    # ind<-ind[which(ind[,1]%%3==1 | ind[,1]%in% index_son),] #index for males who are coded 1 for # of minor alleles on X
    # if (nrow(ind)>0) g0[ind]<-NA
  }
  
  #For autosomes, replace genotypes with 2-genotypes when maf>0.5 
  #For X chromosome, replace genotypes with 2-genotypes for females and 1-genotypes for males when maf>0.5
  if (!xchr) maf<-apply(g0,2,function(x){return(mean(x,na.rm=T)/2)}) else maf<-(apply(g0[index_dad,,drop=F],2,function(x){return(mean(x,na.rm=T))})+apply(g0[index_mom,,drop=F],2,function(x){return(mean(x,na.rm=T)/2)}))/2
  ind<-which(maf>0.5)
  if (length(ind)>0) {
    if (!xchr) g0[,ind]<-2-g0[,ind] else{
      index_dau<-which(sex==0)*3
      g0[index_dad,ind]<-1-g0[index_dad,ind]
      g0[index_mom,ind]<-2-g0[index_mom,ind]
      g0[index_son,ind]<-1-g0[index_son,ind]
      g0[index_dau,ind]<-2-g0[index_dau,ind]
    }
    g0.hap[,ind]<-1-g0.hap[,ind]
  }
  
  removed<-c()
  n_denovo<-n_error<-missing.rate<-rep(0,ncol(g0))
  
  if (hwe){
    if (!xchr){
      tab2x2.pa = data.frame(AA = colSums(g0[-index_off,,drop=F] == 0,na.rm = T), AB = colSums(g0[-index_off,,drop=F] == 1,na.rm = T), BB = colSums(g0[-index_off,,drop=F] == 2,na.rm = T))
      hwe.pa = apply(X = tab2x2.pa, MARGIN = 1, function(x) HardyWeinberg::HWExact(unlist(x),verbose = F)$pval)
      tab2x2.off = data.frame(AA = colSums(g0[index_off,,drop=F] == 0,na.rm = T), AB = colSums(g0[index_off,,drop=F] == 1,na.rm = T), BB = colSums(g0[index_off,,drop=F] == 2,na.rm = T))
      hwe.off = apply(X = tab2x2.off, MARGIN = 1, function(x) HardyWeinberg::HWExact(unlist(x),verbose = F)$pval)
    } else{
      tab2x2.pa = data.frame(A=colSums(g0[index_dad,,drop=F] == 0,na.rm = T),B=colSums(g0[index_dad,,drop=F] == 1,na.rm = T),AA = colSums(g0[index_mom,,drop=F] == 0,na.rm = T), AB = colSums(g0[index_mom,,drop=F] == 1,na.rm = T), BB = colSums(g0[index_mom,,drop=F] == 2,na.rm = T))
      hwe.pa = apply(X = tab2x2.pa, MARGIN = 1, function(x) HardyWeinberg::HWExact(unlist(x), x.linked=T,verbose = F)$pval)
      male<-which(sex==1)
      female<-which(sex==0)
      tab2x2.off = data.frame(A=colSums(g0[index_off,,drop=F][male,,drop=F] == 0,na.rm = T),B=colSums(g0[index_off,,drop=F][male,,drop=F] == 1,na.rm = T),AA = colSums(g0[index_off,,drop=F][female,,drop=F] == 0,na.rm = T), AB = colSums(g0[index_off,,drop=F][female,,drop=F] == 1,na.rm = T), BB = colSums(g0[index_off,,drop=F][female,,drop=F] == 2,na.rm = T))
      hwe.off = apply(X = tab2x2.off, MARGIN = 1, function(x) HardyWeinberg::HWExact(unlist(x), x.linked=T,verbose = F)$pval)
    }
    removed<-union(removed,which(hwe.pa<hwe.threshold | hwe.off<hwe.threshold))
  }
  
  if (!xchr & gender_dependent){
    index_fs<-which(sex==1)*3-2
    index_ms<-which(sex==1)*3-1
    index_son<-which(sex==1)*3
    index_fd<-which(sex==0)*3-2
    index_md<-which(sex==0)*3-1
    index_dau<-which(sex==0)*3
    fs<-g0[index_fs,,drop=F]
    ms<-g0[index_ms,,drop=F]
    son<-g0[index_son,,drop=F]
    fd<-g0[index_fd,,drop=F]
    md<-g0[index_md,,drop=F]
    dau<-g0[index_dau,,drop=F]
    Z_fs<-Z_ms<-array(FALSE,dim=c(length(index_fs),ncol(g0)))
    Z_fd<-Z_md<-array(FALSE,dim=c(length(index_fd),ncol(g0))) #fs:father-son, ms:mother-son, fd:father-daughter, md:mother-daughter
    hetd<-fs==1
    Z_fs[hetd & ((ms==0 & son==1)|(ms==2 & son==2))]<-TRUE
    hetm<-ms==1
    Z_ms[hetm & ((fs==0 & son==1)|(fs==2 & son==2))]<-TRUE
    ind112<-hetd & hetm & son==2
    Z_fs[ind112]<-Z_ms[ind112]<-TRUE
    
    hetd<-fd==1
    Z_fd[hetd & ((md==0 & dau==1)|(md==2 & dau==2))]<-TRUE
    hetm<-md==1
    Z_md[hetm & ((fd==0 & dau==1)|(fd==2 & dau==2))]<-TRUE
    ind112<-hetd & hetm & dau==2
    Z_fd[ind112]<-Z_md[ind112]<-TRUE
    
    t1<-colSums(Z_fs)
    t2<-colSums(Z_ms)
    t3<-colSums(Z_fd)
    t4<-colSums(Z_md)
    p.fisher<-rep(1,length(t1))
    for (i in 1:length(t1)){
      tab<-matrix(c(t1[i],t2[i],t3[i],t4[i]),nrow=2,ncol=2)
      p.fisher[i]<-fisher.test(tab)$p.value
    }
    removed<-union(removed,which(p.fisher<gender_dependent_thr))
  }
  
  if (!xchr) maf<-apply(g0,2,function(x){return(mean(x,na.rm=T)/2)}) else maf<-(apply(g0[index_dad,,drop=F],2,function(x){return(mean(x,na.rm=T))})+apply(g0[index_mom,,drop=F],2,function(x){return(mean(x,na.rm=T)/2)}))/2
  removed<-union(removed,which(maf<maf.thr))
  
  out<-list()
  if (length(removed)>0) {
    out$g0<-g0[,-removed] 
    out$g0.hap<-g0.hap[,-removed]
  }else {
    out$g0<-g0
    out$g0.hap<-g0.hap
  }
  out$removed<-removed
  out$n_denovo<-n_denovo
  out$n_error<-n_error
  out$missing.rate<-missing.rate
  out$maf<-maf
  if (!xchr & gender_dependent) out$p.fisher<-p.fisher
  return(out)
}

#new version
knockofftrio_create_knockoff<-function(dat,pos,M=10,hap=FALSE,dat.hap=NA,xchr=FALSE,sex=NA,phasing.dad=NA,phasing.mom=NA){
  if (nrow(dat) %% 3!=0) stop("The number of rows of the input matrix must be a multiple of three.")
  #if (!all(dat %in% c(0,1,2))) stop("The input matrix can only contain 0, 1, or 2.")
  n<-nrow(dat)/3 #number of trios
  nsnp<-ncol(dat) #number of variants
  
  dad<-dat[seq(1,(3*n),3),,drop=FALSE]
  mom<-dat[seq(2,(3*n),3),,drop=FALSE]
  kid<-dat[seq(3,(3*n),3),,drop=FALSE]
  P<-(dad+mom)/2 #average of parents
  Z<-kid-P #offspring minus average of parents
  info<-info2<-list() #index for informative and non-informative trios
  for (i in 1:nsnp) {
    info[[i]]<-which(Z[,i]!=0) #informative trios
    info2[[i]]<-which(!1:n %in% info[[i]]) #non-informative trios
  }
  if (!hap){
    #Create knockoffs for parents
    dadk<-create.MK(dad,pos,M,corr_max=0.85,info,info2) #knockoff for parent1, dim=(n,nsnp,M). #0.75
    momk<-create.MK(mom,pos,M,corr_max=0.85,info,info2) #knockoff for parent2, dim=(n,nsnp,M). #0.75
    # momk<-create.MK.original(mom,pos,M,corr_max=0.75)
    # dadk<-create.MK.original(dad,pos,M,corr_max=0.75)
    #Create knockoffs for offspring
    kidk<-array(0,dim=c(n,nsnp,M)) #knockoff for offspring, dim=(n,nsnp,M)
    ind_noninfo<-Z==0
    kidk[ind_noninfo]<-(dadk+momk)[ind_noninfo]/2
    #het<-(dad==0 & mom==1) | (dad==1 & mom==0)
    ind<-P==0.5 & kid==0
    kidk[ind]<-pmin(dadk[ind],momk[ind])
    ind<-P==0.5 & kid==1
    kidk[ind]<-pmax(dadk[ind],momk[ind])
    ind<-dad==1 & mom==1 & kid==2
    kidk[ind]<-dadk[ind]+momk[ind]
    ind<-P==1.5 & kid==1
    kidk[ind]<-pmax(dadk[ind],momk[ind])/2
    ind<-P==1.5 & kid==2
    kidk[ind]<-pmax(dadk[ind],momk[ind])/2+pmin(dadk[ind],momk[ind])
  } else{
    dad_hap<-dat.hap[sort(c(seq(1,(4*n),4),seq(2,(4*n),4))),,drop=FALSE]
    ind_dad<-2*(1:n)+sample(c(0,-1),n,replace = T)
    ind_dad2<-(1:(2*n))[-ind_dad]
    dad1<-dad_hap[ind_dad,,drop=FALSE]
    dad2<-dad_hap[ind_dad2,,drop=FALSE]
    
    mom_hap<-dat.hap[sort(c(seq(3,(4*n),4),seq(4,(4*n),4))),,drop=FALSE]
    ind_mom<-2*(1:n)+sample(c(0,-1),n,replace = T)
    ind_mom2<-(1:(2*n))[-ind_mom]
    mom1<-mom_hap[ind_mom,,drop=FALSE]
    mom2<-mom_hap[ind_mom2,,drop=FALSE]
    
    dadk1<-create.MK.original.new(X=dad1,pos,M,corr_max=maxcor,maxBP.neighbor=maxbp)
    dadk2<-create.MK.original.new(X=dad2,pos,M,corr_max=maxcor,maxBP.neighbor=maxbp)
    momk1<-create.MK.original.new(X=mom1,pos,M,corr_max=maxcor,maxBP.neighbor=maxbp)
    momk2<-create.MK.original.new(X=mom2,pos,M,corr_max=maxcor,maxBP.neighbor=maxbp)
    dadk<-dadk1+dadk2
    momk<-momk1+momk2
    kidk<-array(dim=dim(dadk))
    phasing.dad<-phasing.dad%%2
    phasing.mom<-phasing.mom%%2
    for (i in 1:n){
      if (ind_dad[i]%%2==phasing.dad[i]) hap1<-dadk1[i,,,drop=FALSE] else hap1<-dadk2[i,,,drop=FALSE]
      if (ind_mom[i]%%2==phasing.mom[i]) hap2<-momk1[i,,,drop=FALSE] else hap2<-momk2[i,,,drop=FALSE]
      kidk[i,,]<-hap1+hap2
    }
  }
  dat1<-array(dim=c(3*n,nsnp,M)) #knockoff for trio, dim=(3*n,nsnp,M)
  for (i in 1:n){
    dat1[i*3-2,,]<-dadk[i,,]
    dat1[i*3-1,,]<-momk[i,,]
    dat1[i*3,,]<-kidk[i,,]
  }
  return(dat1)
}

#Real data version
knockofftrio_calculate_statistics<-function(dat,dat1=NA,pos,size=c(1,1000,5000,10000),p_value_only=FALSE,adjust_for_covariates=FALSE,y=NA,xchr=FALSE,sex=NA,external.window=FALSE){
  if (nrow(dat) %% 3!=0) stop("The number of rows of the original trio matrix must be a multiple of three.")
  #if (!all(dat %in% c(0,1,2))) stop("The original trio matrix can only contain 0, 1, or 2.")
  n<-nrow(dat)/3 #number of trios
  nsnp<-ncol(dat)
  #out1<-out1.pa<-out1.ma<-rep(0,nsnp)
  out_fbat<-fbat(dat,adjust_for_covariates,y,xchr=xchr,sex=sex)
  z<-out_fbat$additive #z scores for SNP-based FBAT in original trios
  p1<-2*pnorm(-abs(z)) #p-values for SNP-based FBAT in original trios
  p1_sign<-sign(z)
  # out1.pa<-out_fbat$paternal
  # out1.ma<-out_fbat$maternal
  # p1.pa<-2*pnorm(-abs(out1.pa))
  # p1.ma<-2*pnorm(-abs(out1.ma))
  # p1.pa_sign<-sign(out1.pa)
  # p1.ma_sign<-sign(out1.ma)
  
  #rownames(dat)<-1:nrow(dat)
  #p1.tdt<-unname(colTDT(dat)$pval)
  # tat<-colTAT(dat,stratified = T,bothHet = 1)
  # p.tat.pa<-unname(tat$pvalPaternal)
  # p.tat.ma<-unname(tat$pvalMaternal)
  #p.tat<-unname(tat$pval)
  
  index_dad<-seq(1,(3*n),3) #index for parent1
  index_mom<-seq(2,(3*n),3) #index for parent2
  index_off<-seq(3,(3*n),3) #index for offspring
  mac<-apply(dat[-index_off,,drop=F],2,sum) #minor allele count for parents
  if (!xchr) maf<-mac/(4*n) else maf<-(apply(dat[index_dad,,drop=F],2,sum)/n+apply(dat[index_mom,,drop=F],2,sum)/(2*n))/2
  
  #Prepare Z matrix for set-based fbat
  weight<-1/(sqrt(n*maf*(1-maf)))
  weight[which(weight==Inf)]<-0
  #original
  Z<-sweep(out_fbat$Z, 2, weight, '*')
  Zsq<-t(Z)%*%Z
  Zsum<-apply(Z,2,sum)
  # Z.pa<-sweep(out_fbat$Z.pa, 2, weight, '*')
  # Zsq.pa<-t(Z.pa)%*%Z.pa
  # Zsum.pa<-apply(Z.pa,2,sum)
  # Z.ma<-sweep(out_fbat$Z.ma, 2, weight, '*')
  # Zsq.ma<-t(Z.ma)%*%Z.ma
  # Zsum.ma<-apply(Z.ma,2,sum)
  
  #Prepare windows
  window<-data.frame()
  index<-which(maf>=0.01)
  for (i in 1:length(size)){
    if (size[i]==1 & length(index)>0) {
      window<-data.frame(start=pos[index],end=pos[index],ind3=index,ind4=index,n=rep(1,length(index)))
    } else {
      if (!external.window) window<-rbind(window,getslidingwindow(min(pos),max(pos),size[i],pos)) else window<-rbind(window,getslidingwindow_extwin(min(pos),max(pos),size[i],pos,ext.win))
    }
  }
  window<-window[which(window$n>0),]
  window<-window[!duplicated(window[,c(3,4)]),]
  tmp<-which(window$n==1 & !(window$ind3 %in% index))
  if (length(tmp)>0) window<-window[-tmp,]
  nwindow<-nrow(window)
  q1<-p.fbat<-z1<-rep(NA,nwindow) #window p-values for original trios q1.cv<-q1.rv<-q1.pa<-q1.ma<-
  #q1.burden.cv<-q1.burden.rv<-q1.burden.urv<-rep(NA,nwindow) #p.denovo
  #n.denovo<-rep(0,nwindow)
  direction<-rep(NA,nwindow) #direction.pa<-direction.ma<-
  
  if (!p_value_only){
    #The following three warnings only work when xchr==FALSE
    #if (dim(dat1)[1] %% 3!=0) stop("The number of rows of the knockoff trio matrix must be a multiple of three.")
    #if (nrow(dat)!=dim(dat1)[1]) stop("The number of rows of the two matrices must be equal.")
    #if (ncol(dat)!=dim(dat1)[2]) stop("The number of columns of the two matrices must be equal.")
    M<-dim(dat1)[3]
    out_fbat<-fbat(dat,adjust_for_covariates,y,dosage=TRUE,dat1,xchr,sex) #SNP-based FBAT statistics for knockoff trios
    z.ko<-out_fbat$additive #z-scores for SNP-based FBAT in knockoff trios
    Z2<-out_fbat$Z
    p2<-2*pnorm(-abs(z.ko)) #p-values for SNP-based FBAT in knockoff trios, dim=c(nsnp,M)
    # out2.pa<-out_fbat$paternal
    # out2.ma<-out_fbat$maternal
    # Z2.pa<-out_fbat$Z.pa
    # Z2.ma<-out_fbat$Z.ma
    # p2.pa<-2*pnorm(-abs(out2.pa))
    # p2.ma<-2*pnorm(-abs(out2.ma))
    Z2sq<-array(dim=c(nsnp,nsnp,M)) #Z2sq.pa<-Z2sq.ma<-
    Z2sum<-array(dim=c(nsnp,M)) #<-Z2sum.pa<-Z2sum.ma
    if (dim(Z2)[2]>1){
      for (m in 1:M){
        Z2[,,m]<-sweep(Z2[,,m], 2, weight, '*')
        Z2sq[,,m]<-t(Z2[,,m])%*%Z2[,,m]
        Z2sum[,m]<-apply(Z2[,,m],2,sum)
        # Z2.pa[,,m]<-sweep(Z2.pa[,,m], 2, weight, '*')
        # Z2sq.pa[,,m]<-t(Z2.pa[,,m])%*%Z2.pa[,,m]
        # Z2sum.pa[,m]<-apply(Z2.pa[,,m],2,sum)
        # Z2.ma[,,m]<-sweep(Z2.ma[,,m], 2, weight, '*')
        # Z2sq.ma[,,m]<-t(Z2.ma[,,m])%*%Z2.ma[,,m]
        # Z2sum.ma[,m]<-apply(Z2.ma[,,m],2,sum)
      }
    } else if (dim(Z2)[2]==1){
      for (m in 1:M){
        Z2[,,m]<-Z2[,,m]*weight
        Z2sq[,,m]<-t(Z2[,,m])%*%Z2[,,m]
        Z2sum[,m]<-sum(Z2[,,m])
        # Z2.pa[,,m]<-Z2.pa[,,m]*weight
        # Z2sq.pa[,,m]<-t(Z2.pa[,,m])%*%Z2.pa[,,m]
        # Z2sum.pa[,m]<-sum(Z2.pa[,,m])
        # Z2.ma[,,m]<-Z2.ma[,,m]*weight
        # Z2sq.ma[,,m]<-t(Z2.ma[,,m])%*%Z2.ma[,,m]
        # Z2sum.ma[,m]<-sum(Z2.ma[,,m])
      }
    }
    q2<-z2<-array(dim=c(M,nwindow)) #window p-values for knockoff trios q2.cv<-q2.rv<-q2.pa<-q2.ma<-
  }
  
  for (i in 1:nwindow) {
    ind<-window$ind3[i]:window$ind4[i]
    #n.denovo[i]<-sum(n_denovo[ind])
    #p.denovo[i]<-dpois(n.denovo[i], 1.18*10^-8*window$n[i]*n*2) #p=5.65e-05 for a variant to have one de novo in 2394 trios
    
    if (length(ind)==1){
      if (maf[ind]>=0.01){
        q1[i]<-p.fbat[i]<-p1[ind] #q1.cv[i]
        z1[i]<-z[ind]
        # q1.pa[i]<-p1.pa[ind]
        # q1.ma[i]<-p1.ma[ind]
        direction[i]<-p1_sign[ind]
        # direction.pa[i]<-p1.pa_sign[ind]
        # direction.ma[i]<-p1.ma_sign[ind]
      }
    }else{
      z1[i]<-fbat_set(Zsum,Zsq,ind)
      p.fbat[i]<-2*pnorm(-abs(z1[i]))
      
      #single-variant for all
      ind.single<-ind #[mac[ind]>=5]
      p.single<-p1[ind.single]
      # p.pa.single<-p1.pa[ind.single]
      # p.ma.single<-p1.ma[ind.single]
      
      #burden for cv
      # ind.cv<-ind[maf[ind]>=0.01]
      # stat1<-fbat_set(Zsum,Zsq,ind.cv)
      # p.burden.cv<-2*pnorm(-abs(stat1))
      # q1.burden.cv[i]<-p.burden.cv
      # stat<-stat1
      
      #burden for rv
      # ind.rv<-ind[maf[ind]<0.01 & mac[ind]>=5]
      # stat1<-fbat_set(Zsum,Zsq,ind.rv)
      # p.burden.rv<-2*pnorm(-abs(stat1))
      # q1.burden.rv[i]<-p.burden.rv
      # stat<-c(stat,stat1)
      
      #burden for ultra-rare variant (mac<5)
      # ind.urv<-ind[mac[ind]<5]
      # stat1<-fbat_set(Zsum,Zsq,ind.urv)
      # p.burden.urv<-2*pnorm(-abs(stat1))
      # q1.burden.urv[i]<-p.burden.urv
      # stat<-c(stat,stat1)
      
      #Paternal and maternal p-values
      # stat1<-fbat_set(Zsum.pa,Zsq.pa,ind.cv)
      # p.pa.burden.cv<-2*pnorm(-abs(stat1))
      # stat.pa<-stat1
      # stat1<-fbat_set(Zsum.pa,Zsq.pa,ind.rv)
      # p.pa.burden.rv<-2*pnorm(-abs(stat1))
      # stat.pa<-c(stat.pa,stat1)
      # stat1<-fbat_set(Zsum.pa,Zsq.pa,ind.urv)
      # p.pa.burden.urv<-2*pnorm(-abs(stat1))
      # stat.pa<-c(stat.pa,stat1)
      # 
      # stat1<-fbat_set(Zsum.ma,Zsq.ma,ind.cv)
      # p.ma.burden.cv<-2*pnorm(-abs(stat1))
      # stat.ma<-stat1
      # stat1<-fbat_set(Zsum.ma,Zsq.ma,ind.rv)
      # p.ma.burden.rv<-2*pnorm(-abs(stat1))
      # stat.ma<-c(stat.ma,stat1)
      # stat1<-fbat_set(Zsum.ma,Zsq.ma,ind.urv)
      # p.ma.burden.urv<-2*pnorm(-abs(stat1))
      # stat.ma<-c(stat.ma,stat1)
      
      #single-variant for cv
      #p.single.cv<-p1[ind.cv]
      #single-variant for rv
      #p.single.rv<-p1[ind.rv]
      
      #q1[i]<-ACAT(c(p.single,p.burden.cv,p.burden.rv,p.burden.urv))
      q1[i]<-ACAT(c(p.single,p.fbat[i]))
      # q1.cv[i]<-ACAT(c(p.single.cv,p.burden.cv))
      # q1.rv[i]<-ACAT(c(p.single.rv,p.burden.rv,p.burden.urv))
      
      # q1.pa[i]<-ACAT(c(p.pa.single,p.pa.burden.cv,p.pa.burden.rv,p.pa.burden.urv))
      # q1.ma[i]<-ACAT(c(p.ma.single,p.ma.burden.cv,p.ma.burden.rv,p.ma.burden.urv))
      
      #direction of q1[i]
      stat<-c(z1[i],z[ind.single])
      tmp<-which.max(abs(stat))
      if (length(tmp)>0) direction[i]<-sign(stat[tmp])
      
      # stat.pa<-c(stat.pa,out1.pa[ind.single])
      # tmp<-which.max(abs(stat.pa))
      # if (length(tmp)>0) direction.pa[i]<-sign(stat.pa[tmp])
      # 
      # stat.ma<-c(stat.ma,out1.ma[ind.single])
      # tmp<-which.max(abs(stat.ma))
      # if (length(tmp)>0) direction.ma[i]<-sign(stat.ma[tmp])
    }
    
    if (!p_value_only){
      if (length(ind)==1){
        if (maf[ind]>=0.01){
          q2[,i]<-p2[ind,] #q2.cv[,i]<-
          z2[,i]<-z.ko[ind,]
          # q2.pa[,i]<-p2.pa[ind,]
          # q2.ma[,i]<-p2.ma[ind,]
        }
      }else{
        for (m in 1:M){
          p.single.ko<-p2[ind.single,m]
          # p.pa.single.ko<-p2.pa[ind.single,m]
          # p.ma.single.ko<-p2.ma[ind.single,m]
          
          z2[m,i]<-fbat_set(Z2sum[,m],Z2sq[,,m],ind)
          p.burden.ko<-2*pnorm(-abs(z2[m,i]))
          
          # p.burden.cv.ko<-2*pnorm(-abs(fbat_set(Z2sum[,m],Z2sq[,,m],ind.cv)))
          # p.burden.rv.ko<-2*pnorm(-abs(fbat_set(Z2sum[,m],Z2sq[,,m],ind.rv)))
          # p.burden.urv.ko<-2*pnorm(-abs(fbat_set(Z2sum[,m],Z2sq[,,m],ind.urv)))
          # p.pa.burden.cv.ko<-2*pnorm(-abs(fbat_set(Z2sum.pa[,m],Z2sq.pa[,,m],ind.cv)))
          # p.pa.burden.rv.ko<-2*pnorm(-abs(fbat_set(Z2sum.pa[,m],Z2sq.pa[,,m],ind.rv)))
          # p.pa.burden.urv.ko<-2*pnorm(-abs(fbat_set(Z2sum.pa[,m],Z2sq.pa[,,m],ind.urv)))
          # p.ma.burden.cv.ko<-2*pnorm(-abs(fbat_set(Z2sum.ma[,m],Z2sq.ma[,,m],ind.cv)))
          # p.ma.burden.rv.ko<-2*pnorm(-abs(fbat_set(Z2sum.ma[,m],Z2sq.ma[,,m],ind.rv)))
          # p.ma.burden.urv.ko<-2*pnorm(-abs(fbat_set(Z2sum.ma[,m],Z2sq.ma[,,m],ind.urv)))
          # p.single.cv.ko<-p2[ind.cv,m]
          # p.single.rv.ko<-p2[ind.rv,m]
          
          #q2[m,i]<-ACAT(c(p.single.ko,p.burden.cv.ko,p.burden.rv.ko,p.burden.urv.ko))
          q2[m,i]<-ACAT(c(p.single.ko,p.burden.ko))
          #q2.cv[m,i]<-ACAT(c(p.single.cv.ko,p.burden.cv.ko))
          #q2.rv[m,i]<-ACAT(c(p.single.rv.ko,p.burden.rv.ko,p.burden.urv.ko))
          
          #q2.pa[m,i]<-ACAT(c(p.pa.single.ko,p.pa.burden.cv.ko,p.pa.burden.rv.ko,p.pa.burden.urv.ko))
          #q2.ma[m,i]<-ACAT(c(p.ma.single.ko,p.ma.burden.cv.ko,p.ma.burden.rv.ko,p.ma.burden.urv.ko))
        }
      }
    }
  }
  window$ind3<-pos[window$ind3]
  window$ind4<-pos[window$ind4]
  colnames(window)[3]<-"actual_start"
  colnames(window)[4]<-"actual_end"
  if (!p_value_only){
    out<-calculate_w_kappatau(q1,q2)
    w<-out$w
    #w.raw<-out$w.raw
    kappatau<-out$kappatau
    
    # out<-calculate_w_kappatau(q1.pa,q2.pa)
    # w.pa<-out$w
    # w.raw.pa<-out$w.raw
    # kappatau.pa<-out$kappatau
    # colnames(kappatau.pa)<-c("kappa.pa","tau.pa")
    
    # out<-calculate_w_kappatau(q1.ma,q2.ma)
    # w.ma<-out$w
    # w.raw.ma<-out$w.raw
    # kappatau.ma<-out$kappatau
    # colnames(kappatau.ma)<-c("kappa.ma","tau.ma")
    
    # out<-calculate_w_kappatau(q1.cv,q2.cv)
    # w.cv<-out$w
    # kappatau.cv<-out$kappatau
    # colnames(kappatau.cv)<-c("kappa.cv","tau.cv")
    
    # out<-calculate_w_kappatau(q1.rv,q2.rv)
    # w.rv<-out$w
    # kappatau.rv<-out$kappatau
    # colnames(kappatau.rv)<-c("kappa.rv","tau.rv")
    
    rownames(q2)<-paste0("p_",1:M)
    q2<-t(q2)
    rownames(z2)<-paste0("z_",1:M)
    z2<-t(z2)
    # rownames(q2.cv)<-paste0("p_",1:M,"_cv")
    # q2.cv<-t(q2.cv)
    # rownames(q2.rv)<-paste0("p_",1:M,"_rv")
    # q2.rv<-t(q2.rv)
    
    # window<-cbind(chr,window,n.denovo,p.denovo,dir=direction,dir.pa=direction.pa,dir.ma=direction.ma,w,w.raw,w.pa,w.ma,p=q1,p.pa=q1.pa,p.ma=q1.ma,p.burden=p.fbat,w.cv,p.cv=q1.cv,p.burden.cv=q1.burden.cv,w.rv,p.rv=q1.rv,p.burden.rv=q1.burden.rv,p.burden.urv=q1.burden.urv,kappatau,kappatau.pa,kappatau.ma,q2)
    # window<-cbind(chr,window,n.denovo,p.denovo,dir=direction,dir.pa=direction.pa,dir.ma=direction.ma,p=q1,p.pa=q1.pa,p.ma=q1.ma,p.burden=p.fbat,p.cv=q1.cv,p.burden.cv=q1.burden.cv,p.rv=q1.rv,p.burden.rv=q1.burden.rv,p.burden.urv=q1.burden.urv)
    #p: acat p, z: z score for single-variant or set FBAT (depends on if a window contains 1 or more variants), p.burden: p value for single-variant or set FBAT
    window<-cbind(chr,window,dir=direction,w,p=q1,z=z1,p.burden=p.fbat,kappatau,q2,z2)
  } else window<-cbind(chr,window,dir=direction,p=q1,z=z1,p.burden=p.fbat)
  
  # if (file.exists(fname)){
  #   write.table(window,fname,quote=F,row.name=F,col.names = F,append = T)
  # }else{
  #   write.table(window,fname,quote=F,row.name=F,col.names = T,append = F)
  # }
  #write.table(window,fname,quote=F,row.name=F,col.names = T,append = F)
  return(window)
}


# calculate_power_fdr<-function(kappa,tau,causal,M){
#   index_causal<-which(causal==TRUE)
#   q<-MK.q.byStat(kappa=kappa,tau=tau,M=M)
#   fdr.target<-seq(0,0.2,0.01)
#   fdr.observed<-power<-rep(0,21)
#   current<-1
#   for (i in seq(0.01,0.2,0.01)){
#     current<-current+1
#     index<-which(q<=i)
#     len<-length(index)
#     if (len>0) {
#       ndetect<-sum(index %in% index_causal)
#       fdr.observed[current]<-(len-ndetect)/len
#       if (length(index_causal)>0) power[current]<-ndetect/length(index_causal)
#     }
#   }
#   return(data.frame(fdr.target,fdr.observed,power))
# }

# calculate_power_fdr_curve<-function(window=NA){
#   result.fbat<-calculate_power_fdr(kappa = window$kappa,tau=window$tau,causal=window$causal,M=M)
#   #result.smmat<-calculate_power_fdr(kappa = window$kappa.smmat,tau=window$tau.smmat,causal=window$causal,M=M)
#   colnames(result.fbat)<-c("fdr.target","fdr.observed.fbat","power.fbat")
#   #colnames(result.smmat)<-c("fdr.target","fdr.observed.smmat","power.smmat")
#   #result<-cbind(result.fbat,result.smmat)
#   #return(result)
#   return(result.fbat)
# }

MK.q.byStat<-function (kappa,tau,M,Rej.Bound=10000){
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

#original and knockoff
check_exchangeability<-function(dat,dat1,nsnp,member=3,fname){
  library(GGally)
  cov1<-cov2<-cov3<-c()
  ind<-seq(member,nrow(dat),3) #member=1=father, member=2=mother, member=3=kid
  if (nsnp>ncol(dat)) nsnp<-ncol(dat)
  cov1<-cov2<-cov3<-rep(0,nsnp*(nsnp-1)/2)
  current<-0
  for (i in 1:(nsnp-1)){
    for (j in (i+1):nsnp){
      current<-current+1
      cov1[current]<-cov(dat[ind,i],dat[ind,j])
      cov2[current]<-cov(dat[ind,i],dat1[ind,j,1])
      cov3[current]<-cov(dat1[ind,i,1],dat1[ind,j,1])
    }
  }
  covdata<-data.frame(cov.X_X=cov1,cov.X_Xk=cov2,cov.Xk_Xk=cov3)
  #exchangeability holds if cov.X_X=cov.X_Xk=cov.Xk_Xk; cov.X_X=cov.X_Xk is more important
  pm<-ggpairs(covdata,
              columns = 1:ncol(covdata), title = "",  
              axisLabels = "show", columnLabels = colnames(covdata))
  pm2 <- pm
  for(i in 2:pm$nrow) {
    for(j in 1:(i-1)) {
      pm2[i,j] <- pm[i,j] +
        scale_x_continuous(limits = c(-0.8, 0.8)) +
        scale_y_continuous(limits = c(-0.8, 0.8))
    }
  }
  #pdf(paste0(fname,".pdf"), height = 10, width = 10)
  png(paste0(fname,nsnp,"_",member,".png"), height = 700, width = 700)
  print(pm2)
  dev.off()
  #return(covdata)
}

#hap and geno
check_exchangeability<-function(dat1,dat2,nsnp,member=3,fname){
  library(GGally)
  cov1<-cov2<-cov3<-c()
  ind<-seq(member,dim(dat1)[1],3) #member=1=father, member=2=mother, member=3=kid
  if (nsnp>dim(dat1)[2]) nsnp<-dim(dat1)[2]
  cov1<-cov2<-cov3<-rep(0,nsnp*(nsnp-1)/2)
  current<-0
  for (i in 1:(nsnp-1)){
    for (j in (i+1):nsnp){
      current<-current+1
      cov1[current]<-cov(dat1[ind,i,1],dat1[ind,j,1])
      cov2[current]<-cov(dat1[ind,i,1],dat2[ind,j,1])
      cov3[current]<-cov(dat2[ind,i,1],dat2[ind,j,1])
    }
  }
  covdata<-data.frame(cov.X_X=cov1,cov.X_Xk=cov2,cov.Xk_Xk=cov3)
  #exchangeability holds if cov.X_X=cov.X_Xk=cov.Xk_Xk; cov.X_X=cov.X_Xk is more important
  pm<-ggpairs(covdata,
              columns = 1:ncol(covdata), title = "",  
              axisLabels = "show", columnLabels = colnames(covdata))
  pm2 <- pm
  for(i in 2:pm$nrow) {
    for(j in 1:(i-1)) {
      pm2[i,j] <- pm[i,j] +
        scale_x_continuous(limits = c(-0.8, 0.8)) +
        scale_y_continuous(limits = c(-0.8, 0.8))
    }
  }
  #pdf(paste0(fname,".pdf"), height = 10, width = 10)
  png(paste0(fname,"_snp",nsnp,"_member",member,".png"), height = 700, width = 700)
  print(pm2)
  dev.off()
  #return(covdata)
}

overlap_distance<-function(window,cwindow,cwindow_actual,thr){
  proportion.ko<-distance.ko<-proportion.con<-distance.con<-NA #con=conventional
  proportion.actual.ko<-distance.actual.ko<-proportion.actual.con<-distance.actual.con<-NA
  index<-which(window$w>=thr)
  if (length(index)>0){
    n<-length(index)
    overlapped<-window$start[index]<=cwindow[1] & window$end[index]>=cwindow[1] | window$start[index]<=cwindow[2] & window$end[index]>=cwindow[2] | window$start[index]>=cwindow[1] & window$end[index]<=cwindow[2]
    proportion.ko<-sum(overlapped)/n
    distance<-abs((window$start+window$end)/2-(cwindow[1]+cwindow[2])/2)
    distance[overlapped]<-0
    distance.ko<-max(distance)
    
    overlapped<-window$actual_start[index]<=cwindow_actual[1] & window$actual_end[index]>=cwindow_actual[1] | window$actual_start[index]<=cwindow_actual[2] & window$actual_end[index]>=cwindow_actual[2] | window$actual_start[index]>=cwindow_actual[1] & window$actual_end[index]<=cwindow_actual[2]
    proportion.actual.ko<-sum(overlapped)/n
    distance<-abs((window$actual_start+window$actual_end)/2-(cwindow_actual[1]+cwindow_actual[2])/2)
    distance[overlapped]<-0
    distance.actual.ko<-max(distance)
  }
  
  index<-which(window$p<=(0.05/nrow(window)))
  if (length(index)>0){
    n<-length(index)
    overlapped<-window$start[index]<=cwindow[1] & window$end[index]>=cwindow[1] | window$start[index]<=cwindow[2] & window$end[index]>=cwindow[2] | window$start[index]>=cwindow[1] & window$end[index]<=cwindow[2]
    proportion.con<-sum(overlapped)/n
    distance<-abs((window$start+window$end)/2-(cwindow[1]+cwindow[2])/2)
    distance[overlapped]<-0
    distance.con<-max(distance)
    
    overlapped<-window$actual_start[index]<=cwindow_actual[1] & window$actual_end[index]>=cwindow_actual[1] | window$actual_start[index]<=cwindow_actual[2] & window$actual_end[index]>=cwindow_actual[2] | window$actual_start[index]>=cwindow_actual[1] & window$actual_end[index]<=cwindow_actual[2]
    proportion.actual.con<-sum(overlapped)/n
    distance<-abs((window$actual_start+window$actual_end)/2-(cwindow_actual[1]+cwindow_actual[2])/2)
    distance[overlapped]<-0
    distance.actual.con<-max(distance)  
  }
  return(data.frame(seed=para,proportion.ko,distance.ko,proportion.con,distance.con,proportion.actual.ko,distance.actual.ko,proportion.actual.con,distance.actual.con))
}

add_causal<-function(window,pos_causal){
  causal<-rep(FALSE,nrow(window))
  if (length(pos_causal)>0){
    for (i in 1:length(pos_causal)) 
      causal[which(window$actual_start<=pos_causal[i] & window$actual_end>=pos_causal[i])]<-TRUE
  }
  return(cbind(window,causal))
}

#random order
meta_analysis<-function(window1,window2,n1,n2,M){
  window1$win<-paste0(window1$start,window1$end)
  window2$win<-paste0(window2$start,window2$end)
  ind1<-match(window1$win,window2$win)
  ind2<-match(window2$win,window1$win)
  w1<-sqrt(n1)
  w2<-sqrt(n2)
  denom<-sqrt(w1^2+w2^2)
  for (i in 1:length(ind1))
    if (!is.na(ind1[i])){
      z1<-qnorm(window1$p[i]/2,lower.tail=F)*window1$dir[i]
      z2<-qnorm(window2$p[ind1[i]]/2,lower.tail=F)*window2$dir[ind1[i]]
      z<-(z1*w1+z2*w2)/denom
      window1$p[i]<-2*pnorm(-abs(z))
      #ord1<-order(window1[i,paste0("p_",1:M)])
      #ord2<-order(window2[ind1[i],paste0("p_",1:M)])
      for (m in 1:M){
        #index<-ord2[which(ord1==m)]
        coln<-paste0("p_",m)
        z1<-qnorm(window1[i,coln]/2,lower.tail=F)*window1$dir[i]
        z2<-qnorm(window2[ind1[i],coln]/2,lower.tail=F)*window2$dir[ind1[i]]
        z<-(z1*w1+z2*w2)/denom
        window1[i,coln]<-2*pnorm(-abs(z))
      }
    }
  
  out<-calculate_w_kappatau(window1$p,t(window1[,paste0("p_",1:M)]))
  window1$w<-out$w
  window1$w.raw<-out$w.raw
  window1$kappa<-out$kappatau[,1]
  window1$tau<-out$kappatau[,2]
  
  return(rbind(window1,window2[which(is.na(ind2)),]))
}

meta_analysis_zscore<-function(window,n,M){
  window1<-window[[1]]
  #z1<-qnorm(window1$p/2,lower.tail=F)*window1$dir*sqrt(n[1])
  #z1.ko<-matrix(NA,nrow=nrow(window1),ncol=M)
  #coln<-paste0("p_",1:M)
  #for (m in 1:M) z1.ko[,m]<-qnorm(window1[,coln[m]]/2,lower.tail=F)*window1$dir*sqrt(n[1])
  coln<-paste0("z_",1:M)
  window1$z<-window1$z*sqrt(n[1])
  window1[,coln]<-window1[,coln]*sqrt(n[1])
  #window1[,coln][which(is.na(window1[,coln]))]<-0
  meta<-rep(list(1),nrow(window1)) #meta[[i]] records the indexes for the windows that are meta-analyzed for the ith locus in window1
  window1$win<-paste0(window1$start,window1$end)
  uniq<-c()
  for (j in 2:length(window)){
    window2<-window[[j]]
    n2<-n[j]
    window2$win<-paste0(window2$start,window2$end)
    ind1<-match(window1$win,window2$win)
    ind2<-match(window2$win,window1$win)
    w2<-sqrt(n2)
    ind<-which(!is.na(ind1))
    for (i in ind){
      window1$z[i]<-window1$z[i]+window2$z[ind1[i]]*w2
      window1[i,coln]<-window1[i,coln]+window2[ind1[i],coln]*w2
      meta[[i]]<-c(meta[[i]],j)
    }
    uniq<-rbind(uniq,window2[which(is.na(ind2)),])
  }
  colnp<-paste0("p_",1:M)
  for (i in 1:nrow(window1)){
    denom<-sqrt(sum(n[meta[[i]]]))
    window1$z[i]<-window1$z[i]/denom
    window1[i,coln]<-window1[i,coln]/denom
    window1$p[i]<-2*pnorm(-abs(window1$z[i]))
    window1[i,colnp]<-2*pnorm(-abs(as.matrix(window1[i,coln])))
  }
  out<-calculate_w_kappatau(window1$p,t(window1[,colnp]))
  window1$w<-out$w
  #window1$w.raw<-out$w.raw
  window1$kappa<-out$kappatau[,1]
  window1$tau<-out$kappatau[,2]
  out<-rbind(window1,uniq)
  out<-out[,-ncol(out)]
  return(out)
}

#Backup functions
MK.threshold.byStat<-function (kappa,tau,M,fdr = 0.1,Rej.Bound=10000){
  b<-order(tau,decreasing=T)
  c_0<-kappa[b]==0
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    #if(i==1){temp_0=c_0[i]}
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  ok<-which(ratio<=fdr)
  if(length(ok)>0){
    #ok<-ok[which(ok-ok[1]:(ok[1]+length(ok)-1)<=0)]
    return(tau[b][ok[length(ok)]])
  }else{return(Inf)}
}

MK.threshold<-function (T_0,T_k, fdr = 0.1,method='median',Rej.Bound=10000){
  stat<-MK.statistic(T_0,T_k,method=method)
  kappa<-stat[,1];tau<-stat[,2]
  t<-MK.threshold.byStat(kappa,tau,M=ncol(T_k),fdr=fdr,Rej.Bound=Rej.Bound)
  return(t)
}

targetfdr_M<-function(t1,t2,M,jth,method='median'){
  t2_med<-apply(t2,2,median)
  t2_max<-apply(t2,2,max)
  if (method=='median') w<-(t1-t2_med)*(t1>=t2_max) else if (method=='max') w<-(t1-t2_max)*(t1>=t2_max)
  target<-seq(0,0.2,0.01)
  observed<-power<-rep(0,21)
  current<-0
  #falsepos<-c()
  for (q in seq(0,0.2,0.01)){
    current<-current+1
    t<-MK.threshold(t1,t(t2),fdr=q,method=method)
    index<-which(w>=t)
    len<-length(index)
    ndetect<-sum(index %in% jth)
    if (len>0) observed[current]<-(len-ndetect)/len
    if (length(jth)>0) power[current]<-ndetect/length(jth)
    #falsepos<-c(falsepos,index[which(!index %in% jth)]) #false positives
  }
  #falsepos<-sort(unique(falsepos))
  #return(list(data.frame(target,observed,power),falsepos))
  return(data.frame(target,observed,power))
}

targetfdr_M_q<-function(t1,t2,M,jth,method='median'){
  kappatau<-MK.statistic(t1,t(t2),method=method)
  q<-MK.q.byStat(kappatau[,1],kappatau[,2],M=M)
  target<-seq(0,0.2,0.01)
  observed<-power<-rep(0,21)
  current<-1
  for (i in seq(0.01,0.2,0.01)){
    current<-current+1
    index<-which(q<=i)
    len<-length(index)
    if (len>0) {
      ndetect<-sum(index %in% jth)
      observed[current]<-(len-ndetect)/len
      if (length(jth)>0) power[current]<-ndetect/length(jth)
    }
  }
  return(data.frame(target,observed,power))
}

create.MK.original.new <- function(X,pos,M=5,corr_max=0.75,maxN.neighbor=Inf,maxBP.neighbor=100000) {
  
  X <- as.matrix(X)
  sparse.fit <- sparse.cor(X)
  cor.X <- sparse.fit$cor
  cor.X[is.na(cor.X)] <- 0
  cor.X[is.infinite(cor.X)] <- 0
  cov.X <- sparse.fit$cov
  cov.X[is.na(cov.X)] <- 0
  cov.X[is.infinite(cov.X)] <- 0
  Sigma.distance = as.dist(1 - abs(cor.X))
  if(ncol(X)>1){
    fit = hclust(Sigma.distance, method="single")
    corr_max = corr_max
    clusters = cutree(fit, h=1-corr_max)
  }else{clusters<-1}
  
  #temp.M<-matrix(1,M+1,M+1)
  #cov.M<-kronecker(temp.M,cov.X)
  
  #X_k<-matrix(0,nrow(X),ncol(X));index.exist<-c()
  X_k<-array(0,dim=c(nrow(X),ncol(X),M));index.exist<-c()
  for (k in unique(clusters)){
    cluster.fitted<-cluster.residuals<-matrix(NA,nrow(X),sum(clusters==k))
    for(i in which(clusters==k)){
      #print(i)
      index.pos<-which(pos>=max(pos[i]-maxBP.neighbor,pos[1]) & pos<=min(pos[i]+maxBP.neighbor,pos[length(pos)]))
      temp<-abs(cor.X[i,]);temp[which(clusters==k)]<-0;temp[-index.pos]<-0
      
      index<-order(temp,decreasing=T)
      index<-setdiff(index[1:min(length(index),sum(temp>0.05),floor((nrow(X))^(1/3)),maxN.neighbor)],i)
      
      y<-X[,i]
      if(length(index)==0){fitted.values<-mean(y)}else{
        
        x<-X[,index,drop=F];temp.xy<-rbind(mean(y),crossprod(x,y)/length(y)-colMeans(x)*mean(y))
        x.exist<-c()
        for(j in 1:M){
          x.exist<-cbind(x.exist,X_k[,intersect(index,index.exist),j])
        }
        temp.xy<-rbind(temp.xy,crossprod(x.exist,y)/length(y)-colMeans(x.exist)*mean(y))
        
        temp.cov.cross<-sparse.cov.cross(x,x.exist)$cov
        temp.cov<-sparse.cor(x.exist)$cov
        temp.xx<-cov.X[index,index]
        temp.xx<-rbind(cbind(temp.xx,temp.cov.cross),cbind(t(temp.cov.cross),temp.cov))
        
        temp.xx<-cbind(0,temp.xx)
        temp.xx<-rbind(c(1,rep(0,ncol(temp.xx)-1)),temp.xx)
        
        pca.fit<-princomp(covmat=temp.xx)
        v<-pca.fit$loadings
        cump<-cumsum(pca.fit$sdev^2)/sum(pca.fit$sdev^2)
        n.pc<-which(cump>=0.999)[1]#nrow(temp.xx)#nrow(temp.xx)#
        pca.index<-intersect(1:n.pc,which(pca.fit$sdev!=0))#which(cump<=0.99)
        #calculate
        #inverse ZZ matrix
        temp.inv<-v[,pca.index,drop=F]%*%(pca.fit$sdev[pca.index]^(-2)*t(v[,pca.index,drop=F]))
        #beta coefficients
        temp.beta<-temp.inv%*%temp.xy
        
        temp.j<-1
        fitted.values<-temp.beta[1]+crossprod(t(x),temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F])-sum(colMeans(x)*temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F])
        temp.j<-temp.j+ncol(x)
        for(j in 1:M){
          temp.x<-as.matrix(X_k[,intersect(index,index.exist),j])
          if(ncol(temp.x)>=1){
            fitted.values<-fitted.values+crossprod(t(temp.x),temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),,drop=F])-sum(colMeans(temp.x)*temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),,drop=F])
          }
          temp.j<-temp.j+ncol(temp.x)
        }
      }
      residuals<-y-fitted.values
      cluster.fitted[,match(i,which(clusters==k))]<-as.vector(fitted.values)
      cluster.residuals[,match(i,which(clusters==k))]<-as.vector(residuals)
      
      index.exist<-c(index.exist,i)
    }
    #sample mutiple knockoffs
    cluster.sample.index<-sapply(1:M,function(x)sample(1:nrow(X)))
    for(j in 1:M){
      X_k[,which(clusters==k),j]<-cluster.fitted+cluster.residuals[cluster.sample.index[,j],,drop=F]
    }
  }
  return(X_k)
}

ld_diagnosis<-function(falsepos,jth,cor1,cor2){
  out<-data.frame()
  # cor1<-cor(Y2)
  # cor2<-cor(Y2k)
  for (i in 1:length(falsepos)){
    result<-c(1,cor1[falsepos[i],jth],cor2[falsepos[i],jth])
    out<-rbind(out,result)
  }
  trueneg<-which(!1:ncol(cor1) %in% c(falsepos,jth))
  for (i in 1:length(trueneg)){
    result<-c(0,cor1[trueneg[i],jth],cor2[trueneg[i],jth])
    out<-rbind(out,result)
  }
  colnames(out)<-c("falsepos",paste0("Y2_",1:length(jth)),paste0("Y2k_",1:length(jth)))
  return(out)
}

conditional_var<-function(x){
  m<-nrow(x)/3 #x is a matrix, m families
  n<-ncol(x) #n snps
  V<-array(dim=c(n,n))
  for (i in 1:n){
    for (j in i:n){
      s<-0
      for (k in 1:m){
        i1<-k*3-2
        i2<-k*3-1
        i3<-k*3
        e1<-(x[i1,i]+x[i2,i])/2
        e2<-(x[i1,j]+x[i2,j])/2
        s<-s+(x[i3,i]-e1)*(x[i3,j]-e2)
      }
      V[i,j]<-V[j,i]<-s/m
    }
  }
  return (V)
}

prune<-function(X,index,k,ld){
  #index: a random order for X's columns; k: total number of SNPs to be sampled
  index2<-index[1]
  current<-2
  len<-length(index)
  while (current<=len){
    tmp<-cor(X[,c(index2,index[current])])
    diag(tmp)<-0
    if (!(any(tmp>ld))) {
      index2<-c(index2,index[current])
      #print(length(index2))
    }
    current<-current+1
    if (length(index2)==k) break
  }
  return(index2)
}

samplesnp_ldblock<-function(nsnp,maxld,corX0){
  available<-1:ncol(corX0)
  current<-sample(available,1)
  out<-current
  available<-available[-which(available==current)]
  while (length(out)<nsnp){
    tmp<-corX0[current,available]
    available2<-available
    index<-which(tmp>=maxld)
    if (length(index)>0) {
      available2<-available2[-index]
      tmp<-tmp[-index]
    }
    tmp2<-sort(tmp,decreasing=T)[1:3]
    weight<-tmp2/sum(tmp2)
    tmp2<-order(tmp,decreasing=T)[1:3]
    current<-sample(available2[tmp2],1,prob=weight)
    out<-c(out,current)
    available<-available[-which(available==current)]
    #print(length(out))
  }
  return(out)
}

samplesnp<-function(nsnp,maxld,corX0){
  available<-1:ncol(corX0)
  current<-sample(available,1)
  out<-current
  available<-available[-which(available==current)]
  while (length(out)<nsnp & length(available)>0){
    current<-sample(available,1)
    if (all(corX0[current,out]<maxld)){
      out<-c(out,current)
    }
    available<-available[-which(available==current)]
  }
  return(sort(out))
}
generate_sim_data_linear<-function(effectsize = 0.6, CVorRV="cv",n=10000,sigwin=5000,para=1,xchr=FALSE,quan=FALSE,null=FALSE,maxld=0.7,path="/ifs/scratch/msph/eigen/yy3136/P1KO/dbgap/dbGaP-29895/phg000294.v1.NIMH_AutismGenomeProject_v2.genotype-calls-matrixfmt.c1.ARD/"){
  # library(SKAT)
  # data(SKAT.haplotypes)
  #maxld<-0.7
  if (null) maxld<-1
  # if (CVorRV=='cvrv'){
  #   X0<-SKAT.haplotypes$Haplotype
  #   maf0<-SKAT.haplotypes$SNPInfo$FREQ1
  #   pos0<-SKAT.haplotypes$SNPInfo$CHROM_POS
  # } else{
  #   if (CVorRV=="cv") index<-which(SKAT.haplotypes$SNPInfo$FREQ1>=0.01) 
  #   if (CVorRV=="rv") index<-which(SKAT.haplotypes$SNPInfo$FREQ1<0.01)
  #   X0<-SKAT.haplotypes$Haplotype[,index]
  #   maf0<-SKAT.haplotypes$SNPInfo$FREQ1[index]
  #   pos0<-SKAT.haplotypes$SNPInfo$CHROM_POS[index]
  # }
  
  posname<-"pos.str.eur.cv.chr"
  triolistname<-"trio_list.str.eur.txt"
  datname<-"AGP.phs000267.str.eur.cv.chr"
  hapname<-"AGP.phs000267.str.eur.cv.hap.chr"
  chr<-20
  posinfo<-read.table(paste0(path,posname,chr,".txt"),header = T)
  #MACROD2 in hg18: chr20:13924162-15981842
  start<-15981843
  end<-15981843+1000000-1
  flank<-0
  start.flank<-start-flank
  end.flank<-end+flank
  ind<-which(posinfo$pos>=start.flank & posinfo$pos<=end.flank)
  trio<-read.table(paste0(path,triolistname),header=T)
  trio<-trio[-seq(3,nrow(trio),3),] #remove offspring
  pos<-posinfo$pos[ind]
  sex<-trio$SEX[seq(3,nrow(trio),3)] #original coding: 1 male, 2 female
  sex[which(sex==2)]<-0 #new coding: 1 male, 0 female
  
  g0.hap<-read.table(file=paste0(path,hapname,chr,".haps"),skip=min(ind)-1,nrows=length(ind),header = FALSE, stringsAsFactors = FALSE) # rows are SNPs, 7581 columns (5 columns for SNP info + 2*3788)
  g0.hap<-as.matrix(t(g0.hap[,-(1:5)]));dimnames(g0.hap)<-NULL
  iid<-read.table(file=paste0(path,datname,chr,".fam"),stringsAsFactors = F) #can also use AGP.phs000267.str.eur.hap.chr20.sample
  colnames(iid)<-c("FID","IID","father","mother","sex","pheno")
  g0.hap.ind<-match(trio$IID,iid$IID)
  g0.hap.ind1<-g0.hap.ind*2-1
  g0.hap.ind2<-g0.hap.ind*2
  ind.hap<-as.vector(rbind(g0.hap.ind1,g0.hap.ind2))
  g0.hap<-g0.hap[ind.hap,,drop=FALSE]
  g0<-matrix(nrow=nrow(g0.hap)/2,ncol=ncol(g0.hap))
  for (i in 1:nrow(g0)) g0[i,]<-g0.hap[i*2-1,]+g0.hap[i*2,]
  qc<-quality_control(g0=g0,g0.hap=g0.hap,error_number=1,missing_rate=0.05,xchr=xchr,sex=sex,maf.thr=0.01) #error_number: allowable # of erroneous trios, missing_rate: allowable % of trios with missingness
  rm(g0.hap);rm(g0)
  X0<-qc$g0.hap
  maf0<-maf1<-apply(X0,2,sum)/(nrow(X0))
  pos0<-pos
  if (length(qc$removed)>0) pos0<-pos0[-qc$removed]
  
  corX0<-cor(X0)
  # write.table(corX0,"agpcor_cv.txt",quote=F,col.names = F,row.names = F)
  # write.table(pos0,"agppos_cv.txt",quote=F,col.names = F,row.names = F)
  #corX0<-as.matrix(read.table("agpcor_cv.txt",header = F))
  #dimnames(corX0)<-NULL
  #rm(SKAT.haplotypes)
  #corX0<-read.table(paste0("skatcor_",CVorRV,".txt"),header = F)
  totalnsnp<-ncol(X0)
  if (CVorRV=='cvrv') {
    nsnp<-nsnp0<-500
    ncausal<-5
    effect<-0.02
  }
  if (CVorRV=="cv") {
    nsnp<-nsnp0<-500
    if (null) nsnp<-nsnp0<-200
    ncausal<-3
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
  #print(para)
  set.seed(para) #seed<-2
  # #snpindex<-sort(sample(totalnsnp,nsnp))
  # snpindex<-samplesnp(nsnp0,maxld,corX0)
  # nsnp<-length(snpindex)
  tmp<-df %>% group_by(clusters) %>% sample_n(size = 1)
  nsnp<-min(max(clusters),nsnp0)
  snpindex<-sort(sample(tmp$ind,nsnp))
  
  # snpindex<-sample(totalnsnp,totalnsnp)
  # snpindex<-prune(X0,snpindex,nsnp0,ld)
  # if (length(snpindex)<=nsnp0) {
  #   snpindex<-sort(snpindex)
  #   nsnp<-length(snpindex)
  # }else {
  #   snpindex<-sort(snpindex[1:nsnp0])
  #   nsnp<-nsnp0
  # }
  # snpindex<-samplesnp(nsnp,ld,corX0)
  X<-X0[,snpindex] #summary(as.vector(cor(X))) hist(as.vector(cor(X)))
  maf<-maf0[snpindex]
  maf2<-maf1[snpindex]
  pos<-pos0[snpindex]
  rm(X0)
  
  #within 10kb signal window
  conti<-TRUE
  #len<-(max(pos)-min(pos))/5
  len<-0
  while (conti){ 
    window<-getslidingwindow(min(pos)+len,max(pos)-len,sigwin,pos)
    index<-which(window$n>=ncausal)
    if (length(index)>0){
      window<-window[index,]
      window<-window[sample(nrow(window),1),]
      #index<-which(pos>=window$ind1 & pos<=window$ind2)
      index<-window$ind3:window$ind4
      jth<-sort(sample(index,ncausal))
      cwindow<-c(window$start,window$end)
      cwindow_actual<-c(min(pos[jth]),max(pos[jth]))
      conti<-FALSE
    } else sigwin<-sigwin+1000
  }
  
  #heritability models
  # v<-2*maf2[jth]*(1-maf2[jth])
  # tmp<-2*maf[jth]*(1-maf[jth])
  # a<-sqrt(effect/sum(v/tmp))
  # beta<-a/sqrt(tmp)
  
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
  if (CVorRV=="cv") beta<-effect_cv*abs(log10(maf[jth]))
  if (CVorRV=="rv") beta<-effect_rv*abs(log10(maf[jth]))
  if (CVorRV=="cvrv") {
    beta<-effect_rv*abs(log10(maf[jth]))
    ind<-which(maf[jth]>=0.01)
    if (length(ind)>0) beta[ind]<-effect_cv*abs(log10(maf[jth][ind]))
  }
  
  
  b0<-log(p0/(1-p0))
  ncase<-1
  #Z2<-P2<-array(dim=c(n,nsnp)) #P2: average of parents; Z2: offspring minus average of parents
  dat<-array(dim=c(3*n,nsnp))
  available<-1:nrow(X)
  phasing.dad<-phasing.mom<-rep(1,n)
  if (!quan){
    if (!xchr){
      dat.hap<-array(dim=c(4*n,nsnp))
      while (ncase<=n){
        tmp<-sample(available,4)
        hap.dad<-sample(c(1,2),1)
        hap.mom<-sample(c(1,2),1)
        Y<-X[tmp[hap.dad],]+X[tmp[2+hap.mom],]
        P1<-X[tmp[1],]+X[tmp[2],]
        P2<-X[tmp[3],]+X[tmp[4],]
        x0<-0 #rnorm(1)
        alpha<-b0+2*as.numeric(beta%*%Y[jth])
        alpha_p1<-b0+2*as.numeric(beta%*%P1[jth])
        alpha_p2<-b0+2*as.numeric(beta%*%P2[jth])
        prob<-exp(alpha)/(exp(alpha)+1)
        prob_p1<-exp(alpha_p1)/(exp(alpha_p1)+1)
        prob_p2<-exp(alpha_p2)/(exp(alpha_p2)+1)
        if (null | runif(1)<prob & runif(1)>=prob_p1 & runif(1)>=prob_p2){
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
      }
    } else{
      #sex<-c(rep(1,n/2),rep(0,n/2))
      dat.hap<-array(0,dim=c(4*n,nsnp))
      sex<-rep(c(1,0),n/2)
      phasing.dad<-rep(c(2,1),n/2)
      while (ncase<=n){
        tmp<-sample(available,4)
        hap.dad<-1
        hap.mom<-sample(c(1,2),1)
        Y<-(1-sex[ncase])*X[tmp[hap.dad],]+X[tmp[2+hap.mom],]
        P1<-X[tmp[1],]#+X[tmp[2],]
        P2<-X[tmp[3],]+X[tmp[4],]
        x0<-0 #rnorm(1)
        alpha<-b0+as.numeric(beta%*%Y[jth])
        alpha_p1<-b0+as.numeric(beta%*%P1[jth])
        alpha_p2<-b0+as.numeric(beta%*%P2[jth])
        prob<-exp(alpha)/(exp(alpha)+1)
        prob_p1<-exp(alpha_p1)/(exp(alpha_p1)+1)
        prob_p2<-exp(alpha_p2)/(exp(alpha_p2)+1)
        if (null | runif(1)<prob & runif(1)>=prob_p1 & runif(1)>=prob_p2){
          dat[ncase*3-2,]<-P1
          dat[ncase*3-1,]<-P2
          dat[ncase*3,]<-Y
          dat.hap[ncase*4-3,]<-X[tmp[1],]
          dat.hap[ncase*4-1,]<-X[tmp[3],]
          dat.hap[ncase*4,]<-X[tmp[4],]
          #phasing.dad[ncase]<-hap.dad
          phasing.mom[ncase]<-hap.mom
          ncase<-ncase+1
        }
      }
    }
  } else{
    # if (quan){
    #   y.wo.rnorm<-as.matrix(dat[seq(3,nrow(dat),3),jth])%*%as.matrix(beta)
    #   y<-y.wo.rnorm+rnorm(n)
    # }
    dat.hap<-array(dim=c(4*n,nsnp))
    y<-rep(0,3*n)
    while (ncase<=n){
      tmp<-sample(available,4)
      hap.dad<-sample(c(1,2),1)
      hap.mom<-sample(c(1,2),1)
      P1<-X[tmp[1],]+X[tmp[2],]
      P2<-X[tmp[3],]+X[tmp[4],]
      Y<-X[tmp[hap.dad],]+X[tmp[2+hap.mom],]
      y[ncase*3-2]<-2*as.numeric(beta%*%P1[jth])+rnorm(1)
      y[ncase*3-1]<-2*as.numeric(beta%*%P2[jth])+rnorm(1)
      y[ncase*3]<-2*as.numeric(beta%*%Y[jth])+rnorm(1)
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
  }
  
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
      #Z2<-Z2[,-index]
      #P2<-P2[,-index]
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
  out$cwindow<-cwindow
  out$cwindow_actual<-cwindow_actual
  if (xchr) out$sex<-sex
  if (quan) {
    out$y<-as.numeric(y)
    #out$y.wo.rnorm<-as.numeric(y.wo.rnorm)
  }
  out$phasing.dad<-phasing.dad
  out$phasing.mom<-phasing.mom
  return(out)
}
generate_sim_data_nonlinear<-function(effectsize, linear = FALSE, CVorRV="cv",n = 10000, p = 1000, region_size = 10000000, sigwin=5000,para=1,xchr=FALSE,quan=FALSE,null=FALSE,maxld=0.7,path="/ifs/scratch/msph/eigen/yy3136/P1KO/dbgap/dbGaP-29895/phg000294.v1.NIMH_AutismGenomeProject_v2.genotype-calls-matrixfmt.c1.ARD/"){
  if (null) maxld<-1
  posname<-"pos.str.eur.cv.chr"
  triolistname<-"trio_list.str.eur.txt"
  datname<-"AGP.phs000267.str.eur.cv.chr"
  hapname<-"AGP.phs000267.str.eur.cv.hap.chr"
  chr<-20
  posinfo<-read.table(paste0(path,posname,chr,".txt"),header = T)
  #MACROD2 in hg18: chr20:13924162-15981842
  start<-15981843
  end<-15981843+region_size-1
  flank<-0
  start.flank<-start-flank
  end.flank<-end+flank
  ind<-which(posinfo$pos>=start.flank & posinfo$pos<=end.flank)
  trio<-read.table(paste0(path,triolistname),header=T)
  trio<-trio[-seq(3,nrow(trio),3),] #remove offspring
  pos<-posinfo$pos[ind]
  sex<-trio$SEX[seq(3,nrow(trio),3)] #original coding: 1 male, 2 female
  sex[which(sex==2)]<-0 #new coding: 1 male, 0 female
  
  g0.hap<-read.table(file=paste0(path,hapname,chr,".haps"),skip=min(ind)-1,nrows=length(ind),header = FALSE, stringsAsFactors = FALSE) # rows are SNPs, 7581 columns (5 columns for SNP info + 2*3788)
  g0.hap<-as.matrix(t(g0.hap[,-(1:5)]));dimnames(g0.hap)<-NULL
  iid<-read.table(file=paste0(path,datname,chr,".fam"),stringsAsFactors = F) #can also use AGP.phs000267.str.eur.hap.chr20.sample
  colnames(iid)<-c("FID","IID","father","mother","sex","pheno")
  g0.hap.ind<-match(trio$IID,iid$IID)
  g0.hap.ind1<-g0.hap.ind*2-1
  g0.hap.ind2<-g0.hap.ind*2
  ind.hap<-as.vector(rbind(g0.hap.ind1,g0.hap.ind2))
  g0.hap<-g0.hap[ind.hap,,drop=FALSE]
  g0<-matrix(nrow=nrow(g0.hap)/2,ncol=ncol(g0.hap))
  for (i in 1:nrow(g0)) g0[i,]<-g0.hap[i*2-1,]+g0.hap[i*2,]
  qc<-quality_control(g0=g0,g0.hap=g0.hap,error_number=1,missing_rate=0.05,xchr=xchr,sex=sex,maf.thr=0.01) #error_number: allowable # of erroneous trios, missing_rate: allowable % of trios with missingness
  rm(g0.hap);rm(g0)
  X0<-qc$g0.hap
  maf0<-maf1<-apply(X0,2,sum)/(nrow(X0))
  pos0<-pos
  if (length(qc$removed)>0) pos0<-pos0[-qc$removed]
  
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
  maf2<-maf1[snpindex]
  pos<-pos0[snpindex]
  rm(X0)
  
  #within 10kb signal window
  conti<-TRUE
  #len<-(max(pos)-min(pos))/5
  len<-0
  while (conti){ 
    window<-getslidingwindow(min(pos)+len,max(pos)-len,sigwin,pos)
    # index<-which(window$n>=ncausal)
    index<-which(window$n>=ncausal)
    if (length(index)>0){
      window<-window[index,]
      window<-window[sample(nrow(window),1),]
      #index<-which(pos>=window$ind1 & pos<=window$ind2)
      index<-window$ind3:window$ind4
      jth<-sort(sample(index,ncausal))
      cwindow<-c(window$start,window$end)
      cwindow_actual<-c(min(pos[jth]),max(pos[jth]))
      conti<-FALSE
    } else sigwin<-sigwin+1000
  }
  
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
    # beta<-effect_cv/(sqrt(2*maf[jth]*(1-maf[jth])))
  }
  if (CVorRV=="rv") beta<-effect_rv*abs(log10(maf[jth]))
  if (CVorRV=="cvrv") {
    beta<-effect_rv*abs(log10(maf[jth]))
    ind<-which(maf[jth]>=0.01)
    if (length(ind)>0) beta[ind]<-effect_cv*abs(log10(maf[jth][ind]))
  }
  
  b0<-log(p0/(1-p0))
  ncase<-1
  #Z2<-P2<-array(dim=c(n,nsnp)) #P2: average of parents; Z2: offspring minus average of parents
  dat<-array(dim=c(3*n,nsnp))
  available<-1:nrow(X)
  phasing.dad<-phasing.mom<-rep(1,n)
  ##
  if (!quan){
    if (!xchr){
      dat.hap<-array(dim=c(4*n,nsnp))
      y <- rep(0, 3 * n)
      count_0 <- 0
      count_1 <- 0
      print(maf[jth])
      if (!linear){
        while (ncase<=n){
          tmp<-sample(available,4)
          hap.dad<-sample(c(1,2),1)
          hap.mom<-sample(c(1,2),1)
          Y<-X[tmp[hap.dad],]+X[tmp[2+hap.mom],]
          P1<-X[tmp[1],]+X[tmp[2],]
          P2<-X[tmp[3],]+X[tmp[4],]
          ##nonlinear
          potential_p1 = sin(beta %*%P1[jth])*exp(beta %*%P1[jth]) + b0
          potential_p2 = sin(beta %*%P2[jth])*exp(beta %*%P2[jth]) + b0
          potential_Y = sin(beta %*%Y[jth])*exp(beta %*%Y[jth]) + b0
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
          while (ncase<=n){
            tmp<-sample(available,4)
            hap.dad<-sample(c(1,2),1)
            hap.mom<-sample(c(1,2),1)
            Y<-X[tmp[hap.dad],]+X[tmp[2+hap.mom],]
            P1<-X[tmp[1],]+X[tmp[2],]
            P2<-X[tmp[3],]+X[tmp[4],]
            x0<-0 #rnorm(1)
            alpha<-b0+as.numeric(beta%*%Y[jth])
            alpha_p1<-b0+as.numeric(beta%*%P1[jth])
            alpha_p2<-b0+as.numeric(beta%*%P2[jth])
            prob<-exp(alpha)/(exp(alpha)+1)
            prob_p1<-exp(alpha_p1)/(exp(alpha_p1)+1)
            prob_p2<-exp(alpha_p2)/(exp(alpha_p2)+1)
            if (null | runif(1)<prob & runif(1)>=prob_p1 & runif(1)>=prob_p2){
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
            }}}}
    else{
      #sex<-c(rep(1,n/2),rep(0,n/2))
      dat.hap<-array(0,dim=c(4*n,nsnp))
      sex<-rep(c(1,0),n/2)
      phasing.dad<-rep(c(2,1),n/2)
      while (ncase<=n){
        tmp<-sample(available,4)
        hap.dad<-1
        hap.mom<-sample(c(1,2),1)
        Y<-(1-sex[ncase])*X[tmp[hap.dad],]+X[tmp[2+hap.mom],]
        P1<-X[tmp[1],]#+X[tmp[2],]
        P2<-X[tmp[3],]+X[tmp[4],]
        x0<-0 #rnorm(1)
        alpha<-b0+as.numeric(beta%*%Y[jth])
        alpha_p1<-b0+as.numeric(beta%*%P1[jth])
        alpha_p2<-b0+as.numeric(beta%*%P2[jth])
        prob<-exp(alpha)/(exp(alpha)+1)
        prob_p1<-exp(alpha_p1)/(exp(alpha_p1)+1)
        prob_p2<-exp(alpha_p2)/(exp(alpha_p2)+1)
        if (null | runif(1)<prob & runif(1)>=prob_p1 & runif(1)>=prob_p2){
          dat[ncase*3-2,]<-P1
          dat[ncase*3-1,]<-P2
          dat[ncase*3,]<-Y
          dat.hap[ncase*4-3,]<-X[tmp[1],]
          dat.hap[ncase*4-1,]<-X[tmp[3],]
          dat.hap[ncase*4,]<-X[tmp[4],]
          #phasing.dad[ncase]<-hap.dad
          phasing.mom[ncase]<-hap.mom
          ncase<-ncase+1
        }
      }
    }
  } else{
    dat.hap<-array(dim=c(4*n,nsnp))
    y<-rep(0,3*n)
    while (ncase<=n){
      tmp<-sample(available,4)
      hap.dad<-sample(c(1,2),1)
      hap.mom<-sample(c(1,2),1)
      P1<-X[tmp[1],]+X[tmp[2],]
      P2<-X[tmp[3],]+X[tmp[4],]
      Y<-X[tmp[hap.dad],]+X[tmp[2+hap.mom],]
      potential_p1 = sin(beta %*%P1[jth])*exp(beta %*%P1[jth])+b0
      potential_p2 = sin(beta %*%P2[jth])*exp(beta %*%P2[jth])+b0
      potential_Y = sin(beta %*%Y[jth])*exp(beta %*%Y[jth])+b0
      y_p1 <- exp(potential_p1)/(1+exp(potential_p1))
      y_p2 <- exp(potential_p2)/(1+exp(potential_p2))
      y_Y <- exp(potential_Y)/(1+exp(potential_Y))
      if (runif(1)<y_Y & runif(1)>=y_p1 & runif(1)>=y_p2){
      y[ncase*3-2]<-potential_p1+rnorm(1,0,2)
      y[ncase*3-1]<-potential_p2+rnorm(1,0,2)
      y[ncase*3]<-potential_Y+rnorm(1,0,2)
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
    }}
  }
  
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
  out$cwindow<-cwindow
  out$cwindow_actual<-cwindow_actual
  if (xchr) out$sex<-sex
  out$y<-as.numeric(y)
  out$phasing.dad<-phasing.dad
  out$phasing.mom<-phasing.mom
  return(out)
}


generate_sim_data_nonlinear_control<-function(effectsize, linear = FALSE, CVorRV="cv",n = 10000, p = 1000, region_size = 10000000, sigwin=5000,para=1,xchr=FALSE,quan=FALSE,null=FALSE,maxld=0.7,path="/ifs/scratch/msph/eigen/yy3136/P1KO/dbgap/dbGaP-29895/phg000294.v1.NIMH_AutismGenomeProject_v2.genotype-calls-matrixfmt.c1.ARD/"){
  if (null) maxld<-1
  posname<-"pos.str.eur.cv.chr"
  triolistname<-"trio_list.str.eur.txt"
  datname<-"AGP.phs000267.str.eur.cv.chr"
  hapname<-"AGP.phs000267.str.eur.cv.hap.chr"
  chr<-20
  posinfo<-read.table(paste0(path,posname,chr,".txt"),header = T)
  #MACROD2 in hg18: chr20:13924162-15981842
  start<-15981843
  end<-15981843+region_size-1
  flank<-0
  start.flank<-start-flank
  end.flank<-end+flank
  ind<-which(posinfo$pos>=start.flank & posinfo$pos<=end.flank)
  trio<-read.table(paste0(path,triolistname),header=T)
  trio<-trio[-seq(3,nrow(trio),3),] #remove offspring
  pos<-posinfo$pos[ind]
  sex<-trio$SEX[seq(3,nrow(trio),3)] #original coding: 1 male, 2 female
  sex[which(sex==2)]<-0 #new coding: 1 male, 0 female
  
  g0.hap<-read.table(file=paste0(path,hapname,chr,".haps"),skip=min(ind)-1,nrows=length(ind),header = FALSE, stringsAsFactors = FALSE) # rows are SNPs, 7581 columns (5 columns for SNP info + 2*3788)
  g0.hap<-as.matrix(t(g0.hap[,-(1:5)]));dimnames(g0.hap)<-NULL
  iid<-read.table(file=paste0(path,datname,chr,".fam"),stringsAsFactors = F) #can also use AGP.phs000267.str.eur.hap.chr20.sample
  colnames(iid)<-c("FID","IID","father","mother","sex","pheno")
  g0.hap.ind<-match(trio$IID,iid$IID)
  g0.hap.ind1<-g0.hap.ind*2-1
  g0.hap.ind2<-g0.hap.ind*2
  ind.hap<-as.vector(rbind(g0.hap.ind1,g0.hap.ind2))
  g0.hap<-g0.hap[ind.hap,,drop=FALSE]
  g0<-matrix(nrow=nrow(g0.hap)/2,ncol=ncol(g0.hap))
  for (i in 1:nrow(g0)) g0[i,]<-g0.hap[i*2-1,]+g0.hap[i*2,]
  qc<-quality_control(g0=g0,g0.hap=g0.hap,error_number=1,missing_rate=0.05,xchr=xchr,sex=sex,maf.thr=0.01) #error_number: allowable # of erroneous trios, missing_rate: allowable % of trios with missingness
  rm(g0.hap);rm(g0)
  X0<-qc$g0.hap
  maf0<-maf1<-apply(X0,2,sum)/(nrow(X0))
  pos0<-pos
  if (length(qc$removed)>0) pos0<-pos0[-qc$removed]
  
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
  maf2<-maf1[snpindex]
  pos<-pos0[snpindex]
  rm(X0)
  
  #within 10kb signal window
  conti<-TRUE
  #len<-(max(pos)-min(pos))/5
  len<-0
  while (conti){ 
    window<-getslidingwindow(min(pos)+len,max(pos)-len,sigwin,pos)
    # index<-which(window$n>=ncausal)
    index<-which(window$n>=ncausal)
    if (length(index)>0){
      window<-window[index,]
      window<-window[sample(nrow(window),1),]
      #index<-which(pos>=window$ind1 & pos<=window$ind2)
      index<-window$ind3:window$ind4
      jth<-sort(sample(index,ncausal))
      cwindow<-c(window$start,window$end)
      cwindow_actual<-c(min(pos[jth]),max(pos[jth]))
      conti<-FALSE
    } else sigwin<-sigwin+1000
  }
  
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
    # beta<-effect_cv/(sqrt(2*maf[jth]*(1-maf[jth])))
  }
  if (CVorRV=="rv") beta<-effect_rv*abs(log10(maf[jth]))
  if (CVorRV=="cvrv") {
    beta<-effect_rv*abs(log10(maf[jth]))
    ind<-which(maf[jth]>=0.01)
    if (length(ind)>0) beta[ind]<-effect_cv*abs(log10(maf[jth][ind]))
  }
  
  b0<-log(p0/(1-p0))
  ncase<-1
  #Z2<-P2<-array(dim=c(n,nsnp)) #P2: average of parents; Z2: offspring minus average of parents
  dat<-array(dim=c(3*n,nsnp))
  available<-1:nrow(X)
  phasing.dad<-phasing.mom<-rep(1,n)
  ##
  if (!quan){
    if (!xchr){
      dat.hap<-array(dim=c(4*n,nsnp))
      y <- rep(0, 3 * n)
      count_0 <- 0
      count_1 <- 0
      print(maf[jth])
      if (!linear){
        while (ncase<=n){
          tmp<-sample(available,4)
          hap.dad<-sample(c(1,2),1)
          hap.mom<-sample(c(1,2),1)
          Y<-X[tmp[hap.dad],]+X[tmp[2+hap.mom],]
          P1<-X[tmp[1],]+X[tmp[2],]
          P2<-X[tmp[3],]+X[tmp[4],]
          ##nonlinear
          potential_p1 = sin(beta %*%P1[jth])*exp(beta %*%P1[jth]) + b0
          potential_p2 = sin(beta %*%P2[jth])*exp(beta %*%P2[jth]) + b0
          potential_Y = sin(beta %*%Y[jth])*exp(beta %*%Y[jth]) + b0
          y_p1 <- exp(potential_p1)/(1+exp(potential_p1))
          y_p2 <- exp(potential_p2)/(1+exp(potential_p2))
          y_Y <- exp(potential_Y)/(1+exp(potential_Y))
          if (runif(1)<y_Y & runif(1)>=y_p1 & runif(1)>=y_p2){
            print(dim(P1))
            print(dim(dat))
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
          while (ncase<=n){
            tmp<-sample(available,4)
            hap.dad<-sample(c(1,2),1)
            hap.mom<-sample(c(1,2),1)
            Y<-X[tmp[hap.dad],]+X[tmp[2+hap.mom],]
            P1<-X[tmp[1],]+X[tmp[2],]
            P2<-X[tmp[3],]+X[tmp[4],]
            x0<-0 #rnorm(1)
            alpha<-b0+as.numeric(beta%*%Y[jth])
            alpha_p1<-b0+as.numeric(beta%*%P1[jth])
            alpha_p2<-b0+as.numeric(beta%*%P2[jth])
            prob<-exp(alpha)/(exp(alpha)+1)
            prob_p1<-exp(alpha_p1)/(exp(alpha_p1)+1)
            prob_p2<-exp(alpha_p2)/(exp(alpha_p2)+1)
            if (null | runif(1)<prob & runif(1)>=prob_p1 & runif(1)>=prob_p2){
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
            }}}}
    else{
      #sex<-c(rep(1,n/2),rep(0,n/2))
      dat.hap<-array(0,dim=c(4*n,nsnp))
      sex<-rep(c(1,0),n/2)
      phasing.dad<-rep(c(2,1),n/2)
      while (ncase<=n){
        tmp<-sample(available,4)
        hap.dad<-1
        hap.mom<-sample(c(1,2),1)
        Y<-(1-sex[ncase])*X[tmp[hap.dad],]+X[tmp[2+hap.mom],]
        P1<-X[tmp[1],]#+X[tmp[2],]
        P2<-X[tmp[3],]+X[tmp[4],]
        x0<-0 #rnorm(1)
        alpha<-b0+as.numeric(beta%*%Y[jth])
        alpha_p1<-b0+as.numeric(beta%*%P1[jth])
        alpha_p2<-b0+as.numeric(beta%*%P2[jth])
        prob<-exp(alpha)/(exp(alpha)+1)
        prob_p1<-exp(alpha_p1)/(exp(alpha_p1)+1)
        prob_p2<-exp(alpha_p2)/(exp(alpha_p2)+1)
        if (null | runif(1)<prob & runif(1)>=prob_p1 & runif(1)>=prob_p2){
          dat[ncase*3-2,]<-P1
          dat[ncase*3-1,]<-P2
          dat[ncase*3,]<-Y
          dat.hap[ncase*4-3,]<-X[tmp[1],]
          dat.hap[ncase*4-1,]<-X[tmp[3],]
          dat.hap[ncase*4,]<-X[tmp[4],]
          #phasing.dad[ncase]<-hap.dad
          phasing.mom[ncase]<-hap.mom
          ncase<-ncase+1
        }
      }
    }
  } else{
    dat.hap<-array(dim=c(4*n,nsnp))
    y<-rep(0,3*n)
    while (ncase<=n){
      tmp<-sample(available,4)
      hap.dad<-sample(c(1,2),1)
      hap.mom<-sample(c(1,2),1)
      P1<-X[tmp[1],]+X[tmp[2],]
      P2<-X[tmp[3],]+X[tmp[4],]
      Y<-X[tmp[hap.dad],]+X[tmp[2+hap.mom],]
      potential_p1 = sin(beta %*%P1[jth])*exp(beta %*%P1[jth])+b0
      potential_p2 = sin(beta %*%P2[jth])*exp(beta %*%P2[jth])+b0
      potential_Y = sin(beta %*%Y[jth])*exp(beta %*%Y[jth])+b0
      y_p1 <- exp(potential_p1)/(1+exp(potential_p1))
      y_p2 <- exp(potential_p2)/(1+exp(potential_p2))
      y_Y <- exp(potential_Y)/(1+exp(potential_Y))
      if (runif(1)<y_Y & runif(1)>=y_p1 & runif(1)>=y_p2){
        y[ncase*3-2]<-potential_p1+rnorm(1,0,2)
        y[ncase*3-1]<-potential_p2+rnorm(1,0,2)
        y[ncase*3]<-potential_Y+rnorm(1,0,2)
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
      }}
  }
  
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
  out$cwindow<-cwindow
  out$cwindow_actual<-cwindow_actual
  if (xchr) out$sex<-sex
  out$y<-as.numeric(y)
  out$phasing.dad<-phasing.dad
  out$phasing.mom<-phasing.mom
  out_all = list()
  out_all$case = out
  
  ncase<-1
  #Z2<-P2<-array(dim=c(n,nsnp)) #P2: average of parents; Z2: offspring minus average of parents
  dat<-array(dim=c(3*n,nsnp))
  available<-1:nrow(X)
  phasing.dad<-phasing.mom<-rep(1,n)
  ##
  print('generating control trios')
  if (!quan){
    if (!xchr){
      dat.hap<-array(dim=c(4*n,nsnp))
      y <- rep(0, 3 * n)
      count_0 <- 0
      count_1 <- 0
      if (!linear){
        while (ncase<=n){
          tmp<-sample(available,4)
          hap.dad<-sample(c(1,2),1)
          hap.mom<-sample(c(1,2),1)
          Y<-X[tmp[hap.dad],]+X[tmp[2+hap.mom],]
          P1<-X[tmp[1],]+X[tmp[2],]
          P2<-X[tmp[3],]+X[tmp[4],]
          ##nonlinear
          potential_p1 = sin(beta %*%P1[jth])*exp(beta %*%P1[jth]) + b0
          potential_p2 = sin(beta %*%P2[jth])*exp(beta %*%P2[jth]) + b0
          potential_Y = sin(beta %*%Y[jth])*exp(beta %*%Y[jth]) + b0
          y_p1 <- exp(potential_p1)/(1+exp(potential_p1))
          y_p2 <- exp(potential_p2)/(1+exp(potential_p2))
          y_Y <- exp(potential_Y)/(1+exp(potential_Y))
          if (runif(1)>=y_Y & runif(1)>=y_p1 & runif(1)>=y_p2){
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
          while (ncase<=n){
            tmp<-sample(available,4)
            hap.dad<-sample(c(1,2),1)
            hap.mom<-sample(c(1,2),1)
            Y<-X[tmp[hap.dad],]+X[tmp[2+hap.mom],]
            P1<-X[tmp[1],]+X[tmp[2],]
            P2<-X[tmp[3],]+X[tmp[4],]
            x0<-0 #rnorm(1)
            alpha<-b0+as.numeric(beta%*%Y[jth])
            alpha_p1<-b0+as.numeric(beta%*%P1[jth])
            alpha_p2<-b0+as.numeric(beta%*%P2[jth])
            prob<-exp(alpha)/(exp(alpha)+1)
            prob_p1<-exp(alpha_p1)/(exp(alpha_p1)+1)
            prob_p2<-exp(alpha_p2)/(exp(alpha_p2)+1)
            if (null | runif(1)<prob & runif(1)>=prob_p1 & runif(1)>=prob_p2){
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
            }}}}
    else{
      #sex<-c(rep(1,n/2),rep(0,n/2))
      dat.hap<-array(0,dim=c(4*n,nsnp))
      sex<-rep(c(1,0),n/2)
      phasing.dad<-rep(c(2,1),n/2)
      while (ncase<=n){
        tmp<-sample(available,4)
        hap.dad<-1
        hap.mom<-sample(c(1,2),1)
        Y<-(1-sex[ncase])*X[tmp[hap.dad],]+X[tmp[2+hap.mom],]
        P1<-X[tmp[1],]#+X[tmp[2],]
        P2<-X[tmp[3],]+X[tmp[4],]
        x0<-0 #rnorm(1)
        alpha<-b0+as.numeric(beta%*%Y[jth])
        alpha_p1<-b0+as.numeric(beta%*%P1[jth])
        alpha_p2<-b0+as.numeric(beta%*%P2[jth])
        prob<-exp(alpha)/(exp(alpha)+1)
        prob_p1<-exp(alpha_p1)/(exp(alpha_p1)+1)
        prob_p2<-exp(alpha_p2)/(exp(alpha_p2)+1)
        if (null | runif(1)<prob & runif(1)>=prob_p1 & runif(1)>=prob_p2){
          dat[ncase*3-2,]<-P1
          dat[ncase*3-1,]<-P2
          dat[ncase*3,]<-Y
          dat.hap[ncase*4-3,]<-X[tmp[1],]
          dat.hap[ncase*4-1,]<-X[tmp[3],]
          dat.hap[ncase*4,]<-X[tmp[4],]
          #phasing.dad[ncase]<-hap.dad
          phasing.mom[ncase]<-hap.mom
          ncase<-ncase+1
        }
      }
    }
  } else{
    dat.hap<-array(dim=c(4*n,nsnp))
    y<-rep(0,3*n)
    while (ncase<=n){
      tmp<-sample(available,4)
      hap.dad<-sample(c(1,2),1)
      hap.mom<-sample(c(1,2),1)
      P1<-X[tmp[1],]+X[tmp[2],]
      P2<-X[tmp[3],]+X[tmp[4],]
      Y<-X[tmp[hap.dad],]+X[tmp[2+hap.mom],]
      potential_p1 = sin(beta %*%P1[jth])*exp(beta %*%P1[jth])+b0
      potential_p2 = sin(beta %*%P2[jth])*exp(beta %*%P2[jth])+b0
      potential_Y = sin(beta %*%Y[jth])*exp(beta %*%Y[jth])+b0
      y_p1 <- exp(potential_p1)/(1+exp(potential_p1))
      y_p2 <- exp(potential_p2)/(1+exp(potential_p2))
      y_Y <- exp(potential_Y)/(1+exp(potential_Y))
      if (runif(1)>=y_Y & runif(1)>=y_p1 & runif(1)>=y_p2){
        y[ncase*3-2]<-potential_p1+rnorm(1,0,2)
        y[ncase*3-1]<-potential_p2+rnorm(1,0,2)
        y[ncase*3]<-potential_Y+rnorm(1,0,2)
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
      }}
  }
  
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
  out$cwindow<-cwindow
  out$cwindow_actual<-cwindow_actual
  if (xchr) out$sex<-sex
  out$y<-as.numeric(y)
  out$phasing.dad<-phasing.dad
  out$phasing.mom<-phasing.mom
  out_all$control = out
  
  return(out_all)
}