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


# calculate_power_fdr<-function(W, kappa, tau,causal,M=10){
#   index_causal<-which(causal==TRUE)
#   fdr.target<-seq(0,0.2,0.01)
#   fdr.observed<-power<-rep(0,21)
#   current<-1
#   for (i in seq(0.01,0.2,0.01)){
#     thres<-MK.threshold.byStat(kappa=kappa,tau=tau,M=M,i)
#     current<-current+1
#     index<-which(W>=thres)
#     len<-length(index)
#     if (len>0) {
#       ndetect<-sum(index %in% index_causal)
#       fdr.observed[current]<-(len-ndetect)/len
#       if (length(index_causal)>0) power[current]<-ndetect/length(index_causal)
#     }
#   }
#   return(data.frame(fdr.target,fdr.observed,power))
# }

# library(ggplot2)
# library(dplyr)
# library(cowplot)
# 
# compare_methods_plot <- function(
#     results_list,
#     quan,
#     model_type = c("Linear", "Nonlinear"),    # e.g., "linear" or "non-linear"
#     data_type  = c("Binary", "Continuous"),    # e.g., "binary" or "continuous"
#     output_file = "comparison_plot.pdf"
# ) {
#   # Match the user-specified arguments or choose default
#   model_type <- match.arg(model_type)
#   data_type  <- match.arg(data_type)
# 
#   # Define a 6-color palette (one color per method), in the order of your factor below
#   color_values <- c("#83639f", "#ea7827", "#c22f2f", "#449945", "#1f70a9", "#8c564b")
# 
#   # Define shapes, also in a vector of length 6
#   shape_values <- c(16, 17, 18, 19, 15)  # points for each method
# 
#   # Combine all results into a single data frame
#   combined_res <- bind_rows(
#     lapply(results_list, function(res) {
#       data.frame(
#         fdr_target   = res$fdr.target,
#         fdr_observed = res$fdr.observed,
#         power        = res$power,
#         method       = res$method
#       )
#     })
#   )
# 
#   # Ensure factor levels match the desired plotting order
#   combined_res$method <- factor(
#     combined_res$method,
#     levels = c("Ours_shap", "Ours_gradient", "KnockoffTrio", "LightGBM", "XGBoost")
#   )
# 
#   # Create labels that indicate the model/data type
#   fdr_title   <- paste0(model_type, ",  ", data_type)
#   power_title <- paste0(model_type, ", ", data_type)
# 
#   # A helper function to create a custom theme
#   my_theme <- theme_minimal(base_size = 14) +
#     theme(
#       panel.grid          = element_blank(),  # Remove all grid lines
#       panel.border        = element_rect(color = "black", fill = NA, size = 1),
#       plot.title          = element_text(hjust = 0.5, face = "bold", size = 16),
#       legend.text         = element_text(size = 12),
#       legend.title        = element_text(size = 12, face = "bold"),
#       legend.background   = element_blank()  # No border around legend
#     )
# 
#   # Build the FDR plot (no legend)
#   p_fdr <- ggplot(combined_res, aes(
#     x       = fdr_target,
#     y       = fdr_observed,
#     color   = method,
#     shape   = method,
#     linetype= method
#   )) +
#     geom_line(size = 1.2, alpha = 0.8) +
#     geom_point(size = 3) +
#     geom_abline(intercept = 0, slope = 1, color = "grey50", linetype = "dashed") +
#     labs(x = "FDR Target", y = "FDR Observed", title = fdr_title) +
#     scale_color_manual(values = color_values) +
#     scale_shape_manual(values = shape_values) +
#     scale_linetype_manual(values = c("solid", "solid", "solid", "solid", "solid", "solid")) +
#     scale_y_continuous(limits = c(0, 1)) +
#     my_theme +
#     theme(legend.position = "none")  # Remove legend
# 
#   # Build the Power plot (legend in bottom-right corner)
#   p_power <- ggplot(
#     combined_res[!is.na(combined_res$power), ],
#     aes(
#       x       = fdr_target,
#       y       = power,
#       color   = method,
#       shape   = method,
#       linetype= method
#     )
#   ) +
#     geom_line(size = 1.2, alpha = 0.8) +
#     geom_point(size = 3) +
#     labs(x = "FDR Target", y = "Power", title = power_title) +
#     scale_color_manual(values = color_values) +
#     scale_shape_manual(values = shape_values) +
#     scale_linetype_manual(values = c("solid", "solid", "solid", "solid", "solid", "solid")) +
#     scale_y_continuous(limits = c(0, 1)) +
#     my_theme +
#     theme(
#       legend.position = c(0.9, 0.05),  # Place legend in bottom-right corner of the plot
#       legend.justification = c(1, 0) # Align the legend to the bottom-right corner
#     )
# 
#   # Combine the two plots side-by-side
#   p_combined <- cowplot::plot_grid(
#     p_fdr, p_power, ncol = 2, align = "v", rel_widths = c(1, 1)
#   )
# 
#   # Save the final figure with adjusted dimensions
#   ggsave(output_file, plot = p_combined, width = 12, height = 6)  # Adjusted size for precision
# 
#   message("Plots saved as ", output_file)
# }



#####quan
# MAX_DEPTH <- 5             # Shallower trees
# MIN_SAMPLES <- 30          # More samples per leaf
# SUBSAMPLE <- 0.8           # Subsample ratio
# REG_STRENGTH <- 0.1        # L1 and L2 regularization
# EARLY_STOP <- 20           # Early stopping rounds
# LEARNING_RATE <- 0.001     # Slower learning rate

benchmark_function_laptop <- function(U_matrix, y, M, quan) {
  library(lightgbm)
  library(xgboost)
  library(pROC)
  
  # # 参数配置 ----------------------------------------------------------------
  # MAX_DEPTH <- 3
  # MIN_SAMPLES <- 200
  # SUBSAMPLE <- 1         # 改为子采样防止过拟合
  # REG_STRENGTH <- 20
  # EARLY_STOP <- 50
  # LEARNING_RATE <- 0.0001     # 适当提高学习率
  # NROUNDS <- 1000
  if (quan == "binary") {
    NROUNDS <- 1000            # Define number of boosting rounds
    MAX_DEPTH <- 3             # Shallower trees
    MIN_SAMPLES <- 200         # More samples per leaf
    SUBSAMPLE <- 1             # Subsample ratio
    REG_STRENGTH <- 20         # L1 and L2 regularization
    EARLY_STOP <- 50           # Early stopping rounds
    LEARNING_RATE <- 0.0001    # Slower learning rate
  } else {
    NROUNDS <- 1000            # Define number of boosting rounds
    MAX_DEPTH <- 5             # Reduce tree depth
    MIN_SAMPLES <- 50          # Increase minimum samples in leaf
    SUBSAMPLE <- 0.8           # Lower subsample ratio
    REG_STRENGTH <- 0.1        # Strengthen regularization
    EARLY_STOP <- 50           # Early stopping rounds
    LEARNING_RATE <- 0.05      # Lower learning rate
  }
  # 数据校验 ----------------------------------------------------------------
  validate_data <- function(mat, y) {
    stopifnot(
      nrow(mat) == length(y),
      all(!is.na(mat)),
      all(!is.na(y)),
      if(quan == "binary") all(y %in% c(0,1)) else TRUE
    )
  }
  
  # 特征矩阵构建 -------------------------------------------------------------
  num_features <- dim(U_matrix)[2]
  n <- dim(U_matrix)[1]       # Number of rows
  p <- dim(U_matrix)[2]       # Number of columns per slice
  m_plus_1 <- dim(U_matrix)[3] # Number of matrices (should be M+1)
  
  # Check for existing column names, or create default ones
  feature_names <- dimnames(U_matrix)[[2]]
  if (is.null(feature_names)) {
    feature_names <- paste0("Feature", 1:p)
  }
  
  # Combine matrices and assign column names
  X_list <- lapply(1:m_plus_1, function(i) as.matrix(U_matrix[,,i]))
  X <- do.call(cbind, X_list)
  
  # Generate column names: "F1_Feature1", "F1_Feature2", ..., "F2_Feature1", ...
  col_names <- unlist(lapply(1:m_plus_1, function(i) paste0("F", i, "_", feature_names)))
  colnames(X) <- col_names
  validate_data(X, y)
  
  # 数据分割 ----------------------------------------------------------------
  split_data <- function(X, y) {
    set.seed(123)
    if (quan == "binary") {
      caret::createDataPartition(y, p = 0.8, list = FALSE)[,1]
    } else {
      sample(nrow(X), 0.8 * nrow(X))
    }
  }
  
  shuffled_idx <- sample(nrow(X))
  X <- X[shuffled_idx, ]
  y <- y[shuffled_idx]
  train_idx <- split_data(X, y)
  
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_valid <- X[-train_idx, ]
  y_valid <- y[-train_idx]
  
  # 确保验证集不为空
  stopifnot(nrow(X_valid) > 0) 
  
  # 模型训练 ----------------------------------------------------------------
  train_model <- function(method) {
    if (method == "lgb") {
      # LightGBM -----------------------------
      params <- list(
        boosting_type = 'gbdt',
        objective = ifelse(quan == "quan", "regression", "binary"),
        metric = ifelse(quan == "quan", "rmse", "auc"),
        max_depth = MAX_DEPTH,
        min_data_in_leaf = MIN_SAMPLES,
        lambda_l1 = REG_STRENGTH,
        lambda_l2 = REG_STRENGTH,
        learning_rate = LEARNING_RATE,
        feature_pre_filter = FALSE
      )
      
      dtrain <- lgb.Dataset(X_train, label = y_train)
      dvalid <- lgb.Dataset(X_valid, label = y_valid)
      
      model <- lgb.train(
        params = params,
        data = dtrain,
        nrounds = NROUNDS,
        valids = list(valid = dvalid),
        early_stopping_rounds = EARLY_STOP,
        verbose = -1
      )
      
      imp <- lgb.importance(model)
      feature_importance <- setNames(rep(0, ncol(X)), colnames(X))
      feature_importance[imp$Feature] <- imp$Gain
      
    } else if (method == "xgb") {
      # XGBoost ------------------------------
      params <- list(
        objective = ifelse(quan == "quan", "reg:squarederror", "binary:logistic"),
        eval_metric = ifelse(quan == "quan", "rmse", "auc"),
        max_depth = MAX_DEPTH,
        subsample = SUBSAMPLE,
        colsample_bytree = SUBSAMPLE,
        reg_lambda = REG_STRENGTH,
        reg_alpha = REG_STRENGTH,
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
        verbose = 0
      )
      
      imp <- xgb.importance(model = model)
      feature_importance <- setNames(rep(0, ncol(X)), colnames(X))
      feature_importance[imp$Feature] <- imp$Gain
    }
    
    # 结果处理 ------------------------------
    feature_importance <- pmax(feature_importance, 0)
    feature_importance <- matrix(
      feature_importance,
      nrow = M + 1,
      ncol = num_features,
      byrow = TRUE,
      dimnames = list(NULL, colnames(U_matrix[,,1])))
    
    # 验证集预测
    preds <- if (method == "lgb") {
      predict(model, X_valid)
    } else {
      predict(model, dvalid)
    }
    
    metrics <- if (quan == "quan") {
      list(
        R2 = 1 - sum((y_valid - preds)^2)/sum((y_valid - mean(y_valid))^2),
        MSE = mean((y_valid - preds)^2)
      )
    } else {
      list(
        AUC = roc(y_valid, preds)$auc,
        Accuracy = mean(round(preds) == y_valid)
      )
    }
    
    list(feature_importance = feature_importance, metrics = metrics)
  }
  
  # 执行训练 ----------------------------------------------------------------
  results <- list(
    lgb = tryCatch(train_model("lgb"), error = function(e) NULL),
    xgb = tryCatch(train_model("xgb"), error = function(e) NULL)
  )
  
  # 结果校验
  if (is.null(results$lgb)) stop("LightGBM training failed")
  if (is.null(results$xgb)) stop("XGBoost training failed")
  
  return(results)
}

benchmark_function <- function(U_matrix, y, M, quan) {
  library(lightgbm)
  library(xgboost)
  library(pROC)
  
  # Optimized parameter settings to prevent overfitting
  MAX_DEPTH <- 3             # Shallower trees
  MIN_SAMPLES <- 200         # More samples per leaf
  SUBSAMPLE <- 1             # Subsample ratio
  REG_STRENGTH <- 20         # L1 and L2 regularization
  EARLY_STOP <- 50           # Early stopping rounds
  LEARNING_RATE <- 0.0001    # Slower learning rate
  
  num_features <- dim(U_matrix)[2]
  num_samples <- dim(U_matrix)[1]
  X <- do.call(cbind, lapply(1:(M+1), function(i) as.matrix(U_matrix[,,i])))
  colnames(X) <- paste0("feature_", 1:ncol(X))  # Unique, descriptive feature names
  
  # Data splitting
  if (quan == "binary") {
    # Use stratified sampling for classification
    class_prop <- table(y) / length(y)
    train_idx <- unlist(lapply(names(class_prop), function(cls) {
      idx <- which(y == cls)
      sample(idx, size = round(0.8 * length(idx)))
    }))
  } else {
    # Simple random sampling for regression
    train_idx <- sample(nrow(X), 0.8 * nrow(X))
  }
  shuffled_idx <- sample(nrow(X))  # Full shuffle for consistency
  X <- X[shuffled_idx, ]
  y <- y[shuffled_idx]
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_valid <- X[-train_idx, ]
  y_valid <- y[-train_idx]
  
  # Model training function
  train_model <- function(method) {
    if (method == "lgb") {
      # LightGBM setup
      params <- list(
        boosting_type = 'gbdt',
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
        nrounds = ifelse(quan == "quan", 1000, 1000),
        valids = list(valid = dvalid),
        early_stopping_rounds = EARLY_STOP,
        verbose = -1
      )
      # Feature importance
      imp <- lgb.importance(model)
      feature_importance <- rep(0, ncol(X))
      feature_names <- colnames(X)
      for (i in seq_along(imp$Feature)) {
        feature_idx <- which(feature_names == imp$Feature[i])
        if (length(feature_idx) > 0) {
          feature_importance[feature_idx] <- imp$Gain[i]
        }
      }
      # Predict on validation set and compute metrics
      if (quan == "quan") {
        preds <- predict(model, X_valid)
        r2 <- 1 - sum((y_valid - preds)^2) / sum((y_valid - mean(y_valid))^2)
        mse <- mean((y_valid - preds)^2)
        metrics <- list(R2 = r2, MSE = mse)
      } else {
        preds <- predict(model, X_valid)
        accuracy <- mean((preds > 0.5) == y_valid)
        auc <- roc(y_valid, preds)$auc
        metrics <- list(Accuracy = accuracy, AUC = auc)
      }
      
    } else if (method == "xgb") {
      # XGBoost setup
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
      evals <- list(train = dtrain, eval = dvalid)
      
      model <- xgb.train(
        params = params,
        data = dtrain,
        nrounds = ifelse(quan == "quan", 1000, 1000),
        evals = evals,
        early_stopping_rounds = EARLY_STOP,
        verbose = 0
      )
      # Feature importance
      imp <- xgb.importance(model = model)
      feature_importance <- rep(0, ncol(X))
      feature_names <- colnames(X)
      for (i in seq_along(imp$Feature)) {
        feature_idx <- which(feature_names == imp$Feature[i])
        if (length(feature_idx) > 0) {
          feature_importance[feature_idx] <- imp$Gain[i]
        }
      }
      # Predict on validation set and compute metrics
      if (quan == "quan") {
        preds <- predict(model, dvalid)
        r2 <- 1 - sum((y_valid - preds)^2) / sum((y_valid - mean(y_valid))^2)
        mse <- mean((y_valid - preds)^2)
        metrics <- list(R2 = r2, MSE = mse)
      } else {
        preds <- predict(model, dvalid)
        accuracy <- mean((preds > 0.5) == y_valid)
        auc <- roc(y_valid, preds)$auc
        metrics <- list(Accuracy = accuracy, AUC = auc)
      }
    }
    
    # Ensure non-negative importance
    feature_importance <- pmax(feature_importance, 0)
    return(list(feature_importance = feature_importance, metrics = metrics))
  }
  
  # Train models and get results
  lgb_result <- train_model("lgb")
  xgb_result <- train_model("xgb")
  aggregate_importance <- function(imp) {
    imp_matrix <- matrix(imp, nrow = M + 1, ncol = num_features, byrow = TRUE)
    return(imp_matrix)
  }
  lgb_result$feature_importance = aggregate_importance(lgb_result$feature_importance)
  xgb_result$feature_importance = aggregate_importance(xgb_result$feature_importance)
  # Return results
  return(list(
    lgb = list(feature_importance = lgb_result$feature_importance, metrics = lgb_result$metrics),
    xgb = list(feature_importance = xgb_result$feature_importance, metrics = xgb_result$metrics)
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
  # ours_res = calculate_power_fdr(W, k,tau,sim$pos%in%sim$pos_causal,M = 10)
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
  # ours_res = calculate_power_fdr(W, k,tau,sim$pos%in%sim$pos_causal,M = 10)
  return(ours_res)
}

# benchmark_function <- function(U_matrix, y, M, quan) {
#   library(caret)
#   library(doParallel)
#   library(glmnet)
#   library(e1071)
# 
#   if (quan != "quan") {
#     y <- as.factor(y)
#     levels(y) <- make.names(levels(y))
#   }
# 
#   num_features <- dim(U_matrix)[2]
#   lasso <- matrix(nrow = M + 1, ncol = num_features)
#   ridge <- matrix(nrow = M + 1, ncol = num_features)
#   # svm   <- matrix(nrow = M + 1, ncol = num_features)
#   num_cores <- max(1, detectCores() - 1)
#   cl <- makeCluster(num_cores)
#   registerDoParallel(cl)
# 
#   results <- foreach(i = 1:(M + 1), .packages = c("caret", "glmnet","e1071")) %dopar% {
#     X <- as.matrix(U_matrix[, , i])
#     colnames(X) <- paste0("V", 1:ncol(X))
# 
#     if (quan == "quan") {
#       #--------------------------------
#       # 1) Lasso (回归)
#       #--------------------------------
#       cvfit_lasso <- cv.glmnet(
#         x = X, y = y,
#         alpha = 1,
#         nfolds = 5,
#         lambda = 10^seq(1, -4, length.out = 30)
#       )
#       best_lambda_lasso <- cvfit_lasso$lambda.min
#       final_lasso <- glmnet(
#         x = X, y = y,
#         alpha = 1,
#         lambda = best_lambda_lasso
#       )
#       lasso_importance <- as.vector(coef(final_lasso))[-1]
# 
#       #--------------------------------
#       # 2) Ridge (回归)
#       #--------------------------------
#       cvfit_ridge <- cv.glmnet(
#         x = X, y = y,
#         alpha = 0,
#         nfolds = 5,
#         lambda = 10^seq(1, -4, length.out = 30)
#       )
#       best_lambda_ridge <- cvfit_ridge$lambda.min
#       final_ridge <- glmnet(
#         x = X, y = y,
#         alpha = 0,
#         lambda = best_lambda_ridge
#       )
#       ridge_importance <- as.vector(coef(final_ridge))[-1]
# 
#       #--------------------------------
#       # 3) SVM (回归)
#       #--------------------------------
#       # svm_model <- svm(X, y, kernel = "linear", probability = TRUE)
#       # svm_importance <- abs(t(svm_model$coefs) %*% svm_model$SV)
# 
#     } else {
#       #==============================================================
#       # quan != "quan" => 分类场景 (family="binomial")
#       #==============================================================
# 
#       #--------------------------------
#       # 1) Lasso (分类)
#       #--------------------------------
#       cvfit_lasso <- cv.glmnet(
#         x = X, y = y,
#         alpha = 1,
#         family = "binomial",
#         nfolds = 5,
#         type.measure = "auc",  # 以 AUC 为指标
#         lambda = 10^seq(1, -4, length.out = 30)
#       )
#       best_lambda_lasso <- cvfit_lasso$lambda.min
#       final_lasso <- glmnet(
#         x = X, y = y,
#         alpha = 1,
#         family = "binomial",
#         lambda = best_lambda_lasso
#       )
#       lasso_importance <- as.vector(coef(final_lasso))[-1]
# 
#       #--------------------------------
#       # 2) Ridge (分类)
#       #--------------------------------
#       cvfit_ridge <- cv.glmnet(
#         x = X, y = y,
#         alpha = 0,
#         family = "binomial",
#         nfolds = 5,
#         lambda = 10^seq(1, -4, length.out = 30),
#         type.measure = "auc"
#       )
#       best_lambda_ridge <- cvfit_ridge$lambda.min
#       final_ridge <- glmnet(
#         x = X, y = y,
#         alpha = 0,
#         family = "binomial",
#         lambda = best_lambda_ridge
#       )
#       ridge_importance <- as.vector(coef(final_ridge))[-1]
# 
#       #--------------------------------
#       # 3) SVM (分类)
#       #--------------------------------
#       # svm_model <- svm(X, y, kernel = "linear", probability = TRUE)
#       # svm_importance <- abs(t(svm_model$coefs) %*% svm_model$SV)
#     }
# 
#     list(
#       lasso = abs(lasso_importance),
#       ridge = abs(ridge_importance)
#       # svm   = svm_importance
#     )
#   }
# 
#   stopCluster(cl)
#   for (i in seq_along(results)) {
#     lasso[i, ] <- results[[i]]$lasso
#     ridge[i, ] <- results[[i]]$ridge
#     # svm[i, ]   <- results[[i]]$svm
#   }
#   FIs <- list(lasso = lasso, ridge = ridge, svm = svm)
#   return(FIs)
# }

# benchmark_function <- function(U_matrix, y, M, quan) {
#   library(caret)
#   library(glmnet)
#   library(e1071)
# 
#   # 如果不是 "quan"，就把 y 转成因子（分类）
#   if (quan != "quan") {
#     y <- as.factor(y)
#     levels(y) <- make.names(levels(y))
#   }
# 
#   num_features <- dim(U_matrix)[2]
#   lasso <- matrix(nrow = M + 1, ncol = num_features)
#   ridge <- matrix(nrow = M + 1, ncol = num_features)
#   # svm   <- matrix(nrow = M + 1, ncol = num_features)
#   X <- as.matrix(U_matrix[, , 1])
#   for (i in 2:(M+1)){
#     X <- cbind(X,as.matrix(U_matrix[, , i]))
#   }
#     colnames(X) <- paste0("V", 1:ncol(X))
# 
#     if (quan == "quan") {
#       #--------------------------------
#       # 1) Lasso (回归)
#       #--------------------------------
#       cvfit_lasso <- cv.glmnet(
#         x = X, y = y,
#         alpha = 1,
#         nfolds = 5,
#         lambda = 10^seq(1, -4, length.out = 30)
#       )
#       best_lambda_lasso <- cvfit_lasso$lambda.min
#       final_lasso <- glmnet(
#         x = X, y = y,
#         alpha = 1,
#         lambda = best_lambda_lasso
#       )
#       lasso_importance <- as.vector(coef(final_lasso))[-1]
# 
#       #--------------------------------
#       # 2) Ridge (回归)
#       #--------------------------------
#       cvfit_ridge <- cv.glmnet(
#         x = X, y = y,
#         alpha = 0,
#         nfolds = 5,
#         lambda = 10^seq(1, -4, length.out = 30)
#       )
#       best_lambda_ridge <- cvfit_ridge$lambda.min
#       final_ridge <- glmnet(
#         x = X, y = y,
#         alpha = 0,
#         lambda = best_lambda_ridge
#       )
#       ridge_importance <- as.vector(coef(final_ridge))[-1]
# 
#       #--------------------------------
#       # 3) SVM (回归)
#       #--------------------------------
#       # ctrl_svm <- trainControl(
#       #   method = "cv",
#       #   number = 5,
#       #   savePredictions = TRUE
#       # )
#       # svm_model <- train(
#       #   x = X,
#       #   y = y,
#       #   method = "svmLinear",
#       #   trControl = ctrl_svm,
#       #   tuneLength = 5
#       # )
#       # svm_importance <- varImp(svm_model, scale = FALSE)$importance[, 1]
# 
#     } else {
#       #==============================================================
#       # quan != "quan" => 分类场景 (family="binomial")
#       #==============================================================
# 
#       #--------------------------------
#       # 1) Lasso (分类)
#       #--------------------------------
#       cvfit_lasso <- cv.glmnet(
#         x = X, y = y,
#         alpha = 1,
#         family = "binomial",
#         nfolds = 5,
#         type.measure = "auc",  # 以 AUC 为指标
#         lambda = 10^seq(-1, -4, length.out = 30)
#       )
#       best_lambda_lasso <- cvfit_lasso$lambda.min
#       final_lasso <- glmnet(
#         x = X, y = y,
#         alpha = 1,
#         family = "binomial",
#         lambda = best_lambda_lasso
#       )
#       lasso_importance <- as.vector(coef(final_lasso))[-1]
# 
#       #--------------------------------
#       # 2) Ridge (分类)
#       #--------------------------------
#       cvfit_ridge <- cv.glmnet(
#         x = X, y = y,
#         alpha = 0,
#         family = "binomial",
#         nfolds = 5,
#         lambda = 10^seq(-1, -4, length.out = 30),
#         type.measure = "auc"
#       )
#       best_lambda_ridge <- cvfit_ridge$lambda.min
#       final_ridge <- glmnet(
#         x = X, y = y,
#         alpha = 0,
#         family = "binomial",
#         lambda = best_lambda_ridge
#       )
#       ridge_importance <- as.vector(coef(final_ridge))[-1]
# 
#       #--------------------------------
#       # 3) SVM (分类)
#       #--------------------------------
#       # svm_model <- svm(X, y, kernel = "linear", probability = TRUE)
#       # svm_importance <- abs(t(svm_model$coefs) %*% svm_model$SV)
#     }
# 
#     results =list(
#       lasso = abs(lasso_importance),
#       ridge = abs(ridge_importance)
#       # svm   = svm_importance
#     )
#   for (i in 1:(M+1)) {
#     lasso[i, ] <- results$lasso[(1+(i-1)*num_features):(i*num_features)]
#     ridge[i, ] <- results$ridge[(1+(i-1)*num_features):(i*num_features)]
#   }
#   FIs <- list(lasso = lasso, ridge = ridge)
#   return(FIs)
# }


# benchmark_function <- function(U_matrix, y, M, quan) {
#   library(caret)
#   library(doParallel)
#   library(glmnet)
#   library(e1071)
# 
#   if (quan != "quan") {
#     y <- as.factor(y)
#     levels(y) <- make.names(levels(y))
#   }
# 
#   num_features <- dim(U_matrix)[2]
#   lasso <- matrix(nrow = M + 1, ncol = num_features)
#   ridge <- matrix(nrow = M + 1, ncol = num_features)
#   # svm   <- matrix(nrow = M + 1, ncol = num_features)
#   num_cores <- max(1, detectCores() - 1)
#   cl <- makeCluster(num_cores)
#   registerDoParallel(cl)
# 
#   results <- foreach(i = 1:(M + 1), .packages = c("caret", "glmnet","e1071")) %dopar% {
#     X <- as.matrix(U_matrix[, , i])
#     colnames(X) <- paste0("V", 1:ncol(X))
# 
#     if (quan == "quan") {
#       #--------------------------------
#       # 1) Lasso (回归)
#       #--------------------------------
#       cvfit_lasso <- cv.glmnet(
#         x = X, y = y,
#         alpha = 1,
#         nfolds = 5,
#         lambda = 10^seq(1, -4, length.out = 30)
#       )
#       best_lambda_lasso <- cvfit_lasso$lambda.min
#       final_lasso <- glmnet(
#         x = X, y = y,
#         alpha = 1,
#         lambda = best_lambda_lasso
#       )
#       lasso_importance <- as.vector(coef(final_lasso))[-1]
# 
#       #--------------------------------
#       # 2) Ridge (回归)
#       #--------------------------------
#       cvfit_ridge <- cv.glmnet(
#         x = X, y = y,
#         alpha = 0,
#         nfolds = 5,
#         lambda = 10^seq(1, -4, length.out = 30)
#       )
#       best_lambda_ridge <- cvfit_ridge$lambda.min
#       final_ridge <- glmnet(
#         x = X, y = y,
#         alpha = 0,
#         lambda = best_lambda_ridge
#       )
#       ridge_importance <- as.vector(coef(final_ridge))[-1]
# 
#       #--------------------------------
#       # 3) SVM (回归)
#       #--------------------------------
#       # svm_model <- svm(X, y, kernel = "linear", probability = TRUE)
#       # svm_importance <- abs(t(svm_model$coefs) %*% svm_model$SV)
# 
#     } else {
#       #==============================================================
#       # quan != "quan" => 分类场景 (family="binomial")
#       #==============================================================
# 
#       #--------------------------------
#       # 1) Lasso (分类)
#       #--------------------------------
#       cvfit_lasso <- cv.glmnet(
#         x = X, y = y,
#         alpha = 1,
#         family = "binomial",
#         nfolds = 5,
#         type.measure = "auc",  # 以 AUC 为指标
#         lambda = 10^seq(1, -4, length.out = 30)
#       )
#       best_lambda_lasso <- cvfit_lasso$lambda.min
#       final_lasso <- glmnet(
#         x = X, y = y,
#         alpha = 1,
#         family = "binomial",
#         lambda = best_lambda_lasso
#       )
#       lasso_importance <- as.vector(coef(final_lasso))[-1]
# 
#       #--------------------------------
#       # 2) Ridge (分类)
#       #--------------------------------
#       cvfit_ridge <- cv.glmnet(
#         x = X, y = y,
#         alpha = 0,
#         family = "binomial",
#         nfolds = 5,
#         lambda = 10^seq(1, -4, length.out = 30),
#         type.measure = "auc"
#       )
#       best_lambda_ridge <- cvfit_ridge$lambda.min
#       final_ridge <- glmnet(
#         x = X, y = y,
#         alpha = 0,
#         family = "binomial",
#         lambda = best_lambda_ridge
#       )
#       ridge_importance <- as.vector(coef(final_ridge))[-1]
# 
#       #--------------------------------
#       # 3) SVM (分类)
#       #--------------------------------
#       # svm_model <- svm(X, y, kernel = "linear", probability = TRUE)
#       # svm_importance <- abs(t(svm_model$coefs) %*% svm_model$SV)
#     }
# 
#     list(
#       lasso = abs(lasso_importance),
#       ridge = abs(ridge_importance)
#       # svm   = svm_importance
#     )
#   }
# 
#   stopCluster(cl)
#   for (i in seq_along(results)) {
#     lasso[i, ] <- results[[i]]$lasso
#     ridge[i, ] <- results[[i]]$ridge
#     # svm[i, ]   <- results[[i]]$svm
#   }
#   FIs <- list(lasso = lasso, ridge = ridge, svm = svm)
#   return(FIs)
# }



quantile.aggregation <- function(pvals, gamma = NULL){
  
  #Input checks
  if (!is.matrix(pvals)){
    stop("Input pvals must be a matrix")
  }
  
  p <- ncol(pvals)
  B <- nrow(pvals)
  
  #If no gamma is provided, compute sequence
  if(is.null(gamma)){ gamma <- seq(ceiling(0.05 * B) / B, 1 - 1 / B, by = 1/B)}
  
  pvals.aggr <- numeric(p)
  
  for(j in 1:p){#loop over all features
    Qj <- quantile(pvals[,j], gamma) / gamma
    penalty <- if(length(gamma) > 1) (1 - log(min(gamma))) else 1
    pvals.pre <- min(Qj) * penalty
    pvals.aggr[j] <- pmin(pvals.pre, 1)
  }
  return(pvals.aggr)
}

agg.pKO <-  function(W.list, q = 0.2, gamma = 0.3, offset = 1, method = "BH", pvals = FALSE){
  
  #Input checks
  if (!is.list(W.list)) stop('Input W.list must be a list')
  
  if(offset!=1 && offset!=0) {
    stop('Input offset must be either 0 or 1')
  }
  
  if(q < 0 | q > 1) {
    stop('q must be between 0 and 1')
  }
  
  # Input dimensions
  p <- length(W.list[[1]])
  B <- length(W.list)
  
  #Store of all B x p p-values
  int.pval <- matrix(0, nrow = B, ncol = p)
  
  #Step1: Compute Intermediate p-values
  for(b in 1:B){#Loop over all b
    
    W <- W.list[[b]]
    
    #Intermediate p-values
    for(j in 1:p){
      if(W[j] <= 0){
        int.pval[b,j] <- 1
      }
      else{
        int.pval[b,j] <- (1/p)* (offset + sum(W <= -W[j]))
      }
    }
  }
  
  #Step2: Aggregated p-values
  agg.pval <- quantile.aggregation(int.pval, gamma)
  
  #Step3: Compute BH/BY
  if(method %in% c("BH", "BY")){
    p_corr <-  p.adjust(agg.pval,method = method)
  } else  stop("Unknown method specified")
  
  
  #Step 4: Selection Set
  Shat <- as.numeric(which(p_corr <= q))
  res <- list(Shat = Shat, B = B)
  
  if(pvals == TRUE){res <- list(Shat=Shat, B = B, pvals = agg.pval)}
  
  return(res)
}



