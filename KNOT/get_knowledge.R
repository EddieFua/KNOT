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


get_knowledge <- function(dat1, dat, y, path, quan = TRUE) {
  p_matrix <- matrix(nrow = dim(dat1)[3] + 1, ncol = dim(dat)[2])
  for (i in 2:(dim(dat1)[3] + 1)) {
    if (quan) {
      res_fbat <- fbat(
        dat1[, , i - 1],
        adjust_for_covariates = TRUE,
        y = y[seq(3, (3 * nrow(dat) / 3), 3)] - mean(y),
        dosage = FALSE,
        dat1 = NA,
        xchr = FALSE,
        sex = NA
      )
    } else {
      res_fbat <- fbat(
        dat1[, , i - 1],
        adjust_for_covariates = FALSE,
        y = NA,
        dosage = FALSE,
        dat1 = NA,
        xchr = FALSE,
        sex = NA
      )
    }
    
    p <- 2 * pnorm(-abs(res_fbat$additive))
    Z <- abs(res_fbat$additive)
    p_matrix[i, ] <- -log10(p) + Z
  }
  if (quan) {
    res_fbat <- fbat(
      dat,
      adjust_for_covariates = TRUE,
      y = y[seq(3, (3 * nrow(dat) / 3), 3)] - mean(y),
      dosage = FALSE,
      dat1 = NA,
      xchr = FALSE,
      sex = NA
    )
  } else {
    res_fbat <- fbat(
      dat,
      adjust_for_covariates = FALSE,
      y = NA,
      dosage = FALSE,
      dat1 = NA,
      xchr = FALSE,
      sex = NA
    )
  }
  
  p <- 2 * pnorm(-abs(res_fbat$additive))
  Z <- abs(res_fbat$additive)
  p_matrix[1, ] <- -log10(p) + Z
  out_file <- file.path(path, "weight.csv")
  write.csv(p_matrix, out_file, row.names = FALSE)
  
  return(p_matrix)
}