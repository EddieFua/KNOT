rm(list = ls())
for (linear in c(FALSE)){
  for (quan in c(FALSE)){
for (times in 1:10){
  one_class = TRUE
  if (linear & !quan){
    effectsize<-0.25  ####0.5
  }else if(linear & quan) {
    effectsize <- 0.2 #######Original is 0.4
  }else if(!linear & !quan) {
    effectsize<- 0.6 #0.6
  }else if(!linear & quan) {
    effectsize <- 0.6  ####0.6
  }
  SSC = FALSE # If False AGP
  p = 500
  n = 3000
  con = paste0('no_control_',n,'_',p,'_',ifelse(quan,'quan_',''),times)
  if (SSC){
    source(paste0('/Users/fuyinghao/Documents/KONet-Trio/Knockoff_R/function/simulation_function','_one_class_SSC','.R'))
    folder = paste0("sim_one_class_",ifelse(linear,'linear_','nonlinear_'),con)
  }else{
    source(paste0('/Users/fuyinghao/Documents/KONet-Trio/Knockoff_R/function/simulation_function','_one_class','.R'))
    folder = paste0("sim_one_class_",ifelse(linear,'linear_','nonlinear_'),con)}
  CVorRV<-"cv"
  sigwin<-1000000
  para = times
  xchr<-FALSE
  null<-FALSE
  p0 = 0.01
  path="/Users/fuyinghao/Documents/KONet-Trio/simulation_10_replicates/"
  if (!dir.exists(paste0(path,folder))) dir.create(paste0(path,folder))
  set.seed(para)
  name<-"1030"
  if (null) name<-paste0(name,"null")
  if (CVorRV=="cv") nsnp<-p else if (CVorRV=="rv") nsnp<-400 else if (CVorRV=="cvrv") nsnp<-500
  if (null) nsnp<-200
  p_value_only<-FALSE #TRUE
  size<-c(1,1000,5000,10000,20000,50000)
  if (!null) M.vec<-c(10) else M.vec<-10 #M.vec<-c(10,8,6,4)
  chrom<-"1" ###20
  if (chrom=="X") xchr<-TRUE else xchr<-FALSE
  meta<-FALSE
  original<-FALSE #if meta==TRUE, whether print original windows for two sub-studies
  mega<-TRUE #Figure 2 in KnockoffTrio paper
  causalwindow<-TRUE
  if (null) causalwindow<-FALSE
  yout<-FALSE #if quan==TRUE, whether print y (with and without adding rnorm) for calculating WGS y
  ywgs<-FALSE #if quan==TRUE, whether use y from WGS (still testing) or from current sim
  maxcor<-0.7 #0.7
  if (linear == FALSE){                                                                                       
    sim<-generate_sim_data_nonlinear(effectsize, linear = FALSE, CVorRV=CVorRV, n = n, p = nsnp, region_size = 1000000, sigwin=sigwin,para=para,xchr=xchr,quan=quan,null=null,maxld=maxcor,path="/Users/fuyinghao/Documents/KONet-Trio/simulation_10_replicates/")}else{
      sim<-generate_sim_data_linear(effectsize, CVorRV=CVorRV, n = n, sigwin=sigwin,para=para,xchr=xchr,quan=quan,null=null,maxld=maxcor,path="/Users/fuyinghao/Documents/KONet-Trio/simulation_10_replicates/")
    }
  n_denovo<-rep(0,ncol(sim$dat))
  maxbp<-500000
  chr <- 1


  if (mega){
    if (p_value_only) dat1<-NA else {
      #dat1<-knockofftrio_create_knockoff(sim$dat,sim$pos,max(M.vec))
      dat1<-knockofftrio_create_knockoff(dat=sim$dat, pos=sim$pos, M=max(M.vec),hap=T, dat.hap = sim$dat.hap, xchr=xchr, sex=sim$sex, phasing.dad=sim$phasing.dad, phasing.mom=sim$phasing.mom)
    }
    if (quan) {
      if (ywgs){
        if (xchr) fname<-paste0("kotrio.y sim",name,".x.quan ",CVorRV," n",n," nsnp",nsnp," sigwin",sigwin," seed1_200.txt") else fname<-paste0("kotrio.y sim",name,".auto.quan ",CVorRV," n",n," nsnp",nsnp," sigwin",sigwin," seed1_200.txt")
        y<-read.table(fname);y<-y$V1
      } else y<-sim$y[seq(3,(3*n),3)]
    }
    for (M in M.vec){
      if (quan) window.mega<-knockofftrio_calculate_statistics(sim$dat,dat1[,,1:M,drop=F],sim$pos,size,p_value_only,adjust_for_covariates=TRUE,y=y-mean(sim$y),xchr=xchr,sex=sim$sex) else window.mega<-knockofftrio_calculate_statistics(sim$dat,dat1[,,1:M,drop=F],sim$pos,size,p_value_only,adjust_for_covariates=FALSE,y=NA,xchr=xchr,sex=sim$sex)
      window.mega<-add_causal(window.mega,sim$pos_causal)
      if (xchr){
        if (quan) fname<-paste0("kotrio sim",name,".x.quan ",CVorRV," M",M," n",n," nsnp",nsnp," sigwin",sigwin," seed",para,".txt") else fname<-paste0("kotrio sim",name,".x.dich ",CVorRV," M",M," n",n," nsnp",nsnp," sigwin",sigwin," seed",para,".txt")
      } else{
        if (quan) fname<-paste0("kotrio sim",name,".auto.quan ",CVorRV," M",M," n",n," nsnp",nsnp," sigwin",sigwin," seed",para,".txt") else fname<-paste0("kotrio sim",name,".auto.dich ",CVorRV," M",M," n",n," nsnp",nsnp," sigwin",sigwin," seed",para,".txt")
      }
      write.table(window.mega,paste0(path,fname),quote=F,row.name=F,col.names=T,append=F)
    }
  }
  save(dat1,file = paste0('/Users/fuyinghao/Documents/KONet-Trio/simulation_10_replicates/',folder,'/knockoffs.RData'))
  save(sim,file = paste0('/Users/fuyinghao/Documents/KONet-Trio/simulation_10_replicates/',folder,'/original.RData'))


  # window.mega<-knockofftrio_calculate_statistics(sim$dat,dat1,sim$pos,size<-c(1,1000,5000),p_value_only = FALSE,adjust_for_covariates=TRUE,y=sim$y[seq(3,(3*nrow(sim$dat)/3),3)]-mean(sim$y),xchr = FALSE,sex=sim$sex)
  # KnockoffTrio_res = calculate_power_fdr(window.mega$kappa,window.mega$tau,sim$pos%in%sim$pos_causal,10)


  parent_matrix <- Compute_expected(sim$dat)
  index_off <- seq(3, dim(sim$dat)[1], 3)
  index_dad <- seq(1, dim(sim$dat)[1]-2, 3)
  index_mom <- seq(2, dim(sim$dat)[1]-1, 3)
  child_array <- array(dim=c(dim(parent_matrix)[1],dim(sim$dat)[2],max(M.vec)+1))
  child_array[,,1] <- sim$dat[index_off,]
  for (i in 2:(max(M.vec)+1)){
    child_array[,,i] <- dat1[index_off,,i-1]
  }

  dad_array <- array(dim=c(dim(parent_matrix)[1],dim(parent_matrix)[2],max(M.vec)+1))
  dad_array[,,1] <- sim$dat[index_dad,]
  for (i in 2:(max(M.vec)+1)){
    dad_array[,,i] <- dat1[index_dad,,i-1]
  }

  mom_array <- array(dim=c(dim(parent_matrix)[1],dim(parent_matrix)[2],max(M.vec)+1))
  mom_array[,,1] <- sim$dat[index_mom,]
  for (i in 2:(max(M.vec)+1)){
    mom_array[,,i] <- dat1[index_mom,,i-1]
  }

  save(dad_array,file = paste0('/Users/fuyinghao/Documents/KONet-Trio/simulation_10_replicates/',folder,'/dad_array.RData'))
  save(mom_array,file = paste0('/Users/fuyinghao/Documents/KONet-Trio/simulation_10_replicates/',folder,'/mom_array.RData'))
  save(child_array,file = paste0('/Users/fuyinghao/Documents/KONet-Trio/simulation_10_replicates/',folder,'/child_array.RData'))


  original_matrix =  array(dim=c(dim(sim$dat)[1],dim(sim$dat)[2],1))
  original_matrix[,,1] = sim$dat
  save(original_matrix,file = paste0('/Users/fuyinghao/Documents/KONet-Trio/simulation_10_replicates/',folder,'/original_matrix.RData'))
  Y <- cbind(sim$y[index_dad], sim$y[index_mom], sim$y[index_off])
  save(Y,file = paste0('/Users/fuyinghao/Documents/KONet-Trio/simulation_10_replicates/',folder,'/y.RData'))
  print(which(sim$pos%in%sim$pos_causal))

  p_matrix = matrix(nrow = dim(dat1)[3]+1,ncol = dim(sim$dat)[2])
  for (i in 2:(dim(dat1)[3]+1)){
    if (quan){
      res_fbat = fbat(dat1[,,i-1],adjust_for_covariates=TRUE,y=sim$y[seq(3,(3*nrow(sim$dat)/3),3)]-mean(sim$y),dosage=FALSE,dat1=NA,xchr=FALSE,sex=NA)
      p_matrix[i,] = -log10(2*pnorm(-abs(res_fbat$additive))) + abs(res_fbat$additive)
    }else{
    res_fbat = fbat(dat1[,,i-1],adjust_for_covariates=FALSE,y=NA,dosage=FALSE,dat1=NA,xchr=FALSE,sex=NA)}
    p_matrix[i,] = -log10(2*pnorm(-abs(res_fbat$additive))) + abs(res_fbat$additive)
  }
  if (quan){
    res_fbat = fbat(sim$dat,adjust_for_covariates=TRUE,y=sim$y[seq(3,(3*nrow(sim$dat)/3),3)]-mean(sim$y),dosage=FALSE,dat1=NA,xchr=FALSE,sex=NA)
  }else{
    res_fbat = fbat(sim$dat,adjust_for_covariates=FALSE,y=NA,dosage=FALSE,dat1=NA,xchr=FALSE,sex=NA)}
  p_matrix[1,] = -log10(2*pnorm(-abs(res_fbat$additive))) + abs(res_fbat$additive)
  print((2*pnorm(-abs(res_fbat$additive)))[sim$pos%in%sim$pos_causal])
  write.csv(p_matrix,paste0('/Users/fuyinghao/Documents/KONet-Trio/simulation_10_replicates/',folder,'/weight.csv'), row.names = FALSE)}}}






