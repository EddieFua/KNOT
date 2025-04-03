library(dplyr)
library(ggplot2)
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

add_causal<-function(window,pos_causal){
  causal<-rep(FALSE,nrow(window))
  if (length(pos_causal)>0){
    for (i in 1:length(pos_causal)) 
      causal[which(window$actual_start<=pos_causal[i] & window$actual_end>=pos_causal[i])]<-TRUE
  }
  return(cbind(window,causal))
}

calculate_power_fdr<-function(kappa,tau,causal,M){
  index_causal<-which(causal==TRUE)
  q<-MK.q.byStat(kappa=kappa,tau=tau,M=M)
  fdr.target<-seq(0,0.2,0.01)
  fdr.observed<-power<-rep(0,21)
  current<-1
  for (i in seq(0.01,0.2,0.01)){
    current<-current+1
    index<-which(q<=i)
    len<-length(index)
    if (len>0) {
      ndetect<-sum(index %in% index_causal)
      fdr.observed[current]<-(len-ndetect)/len
      if (length(index_causal)>0) power[current]<-ndetect/length(index_causal)
    }
  }
  return(data.frame(fdr.target,fdr.observed,power))
}

get_target_fdr_idx<-function(kappa,tau,M,fdr){
  q<-MK.q.byStat(kappa=kappa,tau=tau,M=M)
  idx <- which(q<=fdr)
  return(idx)
}

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


library(lightgbm)
library(xgboost)
library(pROC)

benchmark_function <- function(U_matrix, y, M, quan) {
  # Purpose: Train LightGBM and XGBoost models on a 3D feature matrix, compute feature importance,
  #          and evaluate performance on a validation set.
  # Args:
  #   U_matrix: 3D array (n_samples x n_features x (M+1)) of input features
  #   y: Vector of target values (length n_samples)
  #   M: Integer, number of additional feature sets (total sets = M+1)
  #   quan: Character, "quan" for regression, "binary" for classification
  # Returns:
  #   List with two elements (lgb and xgb), each containing feature_importance (matrix) and metrics (list)

  # --- Hyperparameters to prevent overfitting ---
  MAX_DEPTH <- 3          # Shallower trees
  MIN_SAMPLES <- 200      # More samples per leaf
  SUBSAMPLE <- 1          # Subsample ratio
  REG_STRENGTH <- 20      # L1 and L2 regularization
  EARLY_STOP <- 50        # Early stopping rounds
  LEARNING_RATE <- 0.0001 # Slower learning rate
  NROUNDS <- 1000         # Number of boosting rounds

  # --- Input validation ---
  stopifnot(
    is.array(U_matrix),
    length(dim(U_matrix)) == 3,
    dim(U_matrix)[1] == length(y),
    dim(U_matrix)[3] == M + 1,
    quan %in% c("quan", "binary")
  )
  if (quan == "binary") stopifnot(all(y %in% c(0, 1)))

  # --- Feature matrix construction ---
  num_features <- dim(U_matrix)[2]
  num_samples <- dim(U_matrix)[1]
  # Flatten 3D array into 2D matrix: (n_samples x (n_features * (M+1)))
  X <- do.call(cbind, lapply(1:(M + 1), function(i) as.matrix(U_matrix[, , i])))
  colnames(X) <- paste0("feature_", 1:ncol(X))  # Assign unique feature names

  # --- Data splitting ---
  set.seed(123)  # For reproducibility
  shuffled_idx <- sample(num_samples)  # Shuffle data
  X <- X[shuffled_idx, ]
  y <- y[shuffled_idx]
  if (quan == "binary") {
    # Stratified sampling for binary classification
    class_prop <- table(y) / length(y)
    train_idx <- unlist(lapply(names(class_prop), function(cls) {
      idx <- which(y == cls)
      sample(idx, size = round(0.8 * length(idx)))
    }))
  } else {
    # Random sampling for regression
    train_idx <- sample(num_samples, 0.8 * num_samples)
  }
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_valid <- X[-train_idx, ]
  y_valid <- y[-train_idx]

  # --- Model training function ---
  train_model <- function(method) {
    if (method == "lgb") {
      # LightGBM configuration
      params <- list(
        boosting_type = "gbdt",
        objective = ifelse(quan == "quan", "regression", "binary"),
        metric = ifelse(quan == "quan", "rmse", "auc"),
        max_depth = MAX_DEPTH,
        min_data_in_leaf = MIN_SAMPLES,
        lambda_l1 = REG_STRENGTH,
        lambda_l2 = REG_STRENGTH,
        learning_rate = LEARNING_RATE
      )
      dtrain <- lgb.Dataset(X_train, label = y_train)
      dvalid <- lgb.Dataset(X_valid, label = y_valid, reference = dtrain)
      model <- lgb.train(
        params = params,
        data = dtrain,
        nrounds = NROUNDS,
        valids = list(valid = dvalid),
        early_stopping_rounds = EARLY_STOP,
        verbose = -1  # Silent training
      )
      # Feature importance
      imp <- lgb.importance(model)
      feature_importance <- numeric(ncol(X))
      feature_names <- colnames(X)
      for (i in seq_along(imp$Feature)) {
        idx <- which(feature_names == imp$Feature[i])
        if (length(idx) > 0) feature_importance[idx] <- imp$Gain[i]
      }
      preds <- predict(model, X_valid)
    } else if (method == "xgb") {
      # XGBoost configuration
      params <- list(
        objective = ifelse(quan == "quan", "reg:squarederror", "binary:logistic"),
        eval_metric = ifelse(quan == "quan", "rmse", "auc"),
        max_depth = MAX_DEPTH,
        subsample = SUBSAMPLE,
        colsample_bytree = SUBSAMPLE,
        lambda = REG_STRENGTH,
        alpha = REG_STRENGTH,
        eta = LEARNING_RATE
      )
      dtrain <- xgb.DMatrix(X_train, label = y_train)
      dvalid <- xgb.DMatrix(X_valid, label = y_valid)
      model <- xgb.train(
        params = params,
        data = dtrain,
        nrounds = NROUNDS,
        watchlist = list(train = dtrain, eval = dvalid),
        early_stopping_rounds = EARLY_STOP,
        verbose = 0  # Silent training
      )
      # Feature importance
      imp <- xgb.importance(model = model)
      feature_importance <- numeric(ncol(X))
      feature_names <- colnames(X)
      for (i in seq_along(imp$Feature)) {
        idx <- which(feature_names == imp$Feature[i])
        if (length(idx) > 0) feature_importance[idx] <- imp$Gain[i]
      }
      preds <- predict(model, dvalid)
    } else {
      stop("Method must be 'lgb' or 'xgb'")
    }

    # --- Compute metrics ---
    if (quan == "quan") {
      r2 <- 1 - sum((y_valid - preds)^2) / sum((y_valid - mean(y_valid))^2)
      mse <- mean((y_valid - preds)^2)
      metrics <- list(R2 = r2, MSE = mse)
    } else {
      accuracy <- mean((preds > 0.5) == y_valid)
      auc <- roc(y_valid, preds)$auc
      metrics <- list(Accuracy = accuracy, AUC = auc)
    }

    # --- Format feature importance ---
    feature_importance <- pmax(feature_importance, 0)  # Ensure non-negative
    imp_matrix <- matrix(
      feature_importance,
      nrow = M + 1,
      ncol = num_features,
      byrow = TRUE
    )
    return(list(feature_importance = imp_matrix, metrics = metrics))
  }

  # --- Train models and return results ---
  lgb_result <- train_model("lgb")
  xgb_result <- train_model("xgb")
  return(list(
    lgb = lgb_result,
    xgb = xgb_result
  ))
}

from_FI_to_fdr = function(FI,M,sim){
  k <- apply(FI, 2, function(x) {
    if (all(x == 0)) {
      return(1)  # Return NA for columns with all zeros
    } else {
      return(which.max(x) - 1)  # Subtract 1 only for valid cases
    }
  })
  tau <- apply(FI, 2, function(i) max(i)-median(i[-which.max(i)]))
  W <- apply(FI, 2, function(i) (i[1]-median(i[2:(M+1)]))*ifelse(i[1]>=max(i[2:(M+1)]),1,0))
  ours_res = calculate_power_fdr(k,tau,sim$pos%in%sim$pos_causal,M)
  return(ours_res)
}

get_target_idx = function(FI,M,target_fdr){
  k <- apply(FI, 2, function(x) {
    if (all(x == 0)) {
      return(1)  # Return NA for columns with all zeros
    } else {
      return(which.max(x) - 1)  # Subtract 1 only for valid cases
    }
  })
  tau <- apply(FI, 2, function(i) max(i)-median(i[-which.max(i)]))
  W <- apply(FI, 2, function(i) (i[1]-median(i[2:(M+1)]))*ifelse(i[1]>=max(i[2:(M+1)]),1,0))
  ours_res = get_target_fdr_idx(k,tau,M,target_fdr)
  return(ours_res)
}