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
      index.exist <- c(index.exist, i)
    }
  }
  return(X_k)
}

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

maxcor<-0.7
maxbp<-500000