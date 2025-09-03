compute_shap_interaction_pvalues <- function(genotype_dat, y, gene_closer, N = 100, seed = 10) {
  set.seed(seed)
  
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required. Please install it with install.packages('xgboost').")
  }
  
  n_samples  <- nrow(genotype_dat)
  n_features <- ncol(genotype_dat)
  
  is_binary <- length(unique(y)) == 2 && all(sort(unique(y)) %in% c(0, 1))
  
  # Choose objective/metric based on outcome type
  if (is_binary) {
    objective   <- "binary:logistic"
    eval_metric <- "logloss"
  } else {
    objective   <- "reg:squarederror"
    eval_metric <- "rmse"
  }
  
  params <- list(
    objective = objective,
    eta       = 0.3,
    max_depth = 3,
    eval_metric = eval_metric
  )
  
  dtrain <- xgboost::xgb.DMatrix(data = genotype_dat, label = y)
  
  cv <- xgboost::xgb.cv(
    params = params,
    data = dtrain,
    nrounds = 200,
    nfold = 5,
    early_stopping_rounds = 10,
    verbose = 0
  )
  
  nrounds <- cv$best_iteration
  if (is.null(nrounds) || nrounds == 0) nrounds <- 200
  
  model <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = nrounds,
    verbose = 0
  )
  
  compute_shap_raw_matrix <- function(model, data_matrix, snps1, snps2) {
    dmatrix <- xgboost::xgb.DMatrix(data = data_matrix)
    shap_int_3d <- predict(model, dmatrix, predinteraction = TRUE)   
    mean_shap <- apply(shap_int_3d, c(2, 3), mean)             
    
    raw_matrix <- matrix(0, nrow = length(snps1), ncol = length(snps2))
    for (i in seq_along(snps1)) {
      for (j in seq_along(snps2)) {
        raw_matrix[i, j] <- mean_shap[snps1[i], snps2[j]]
      }
    }
    raw_matrix
  }
  
  # Indices for permutation
  if (is_binary) {
    cases    <- which(y == 1)
    controls <- which(y == 0)
  } else {
    all_idx <- seq_len(n_samples)
  }
  
  # -----------------
  # Gene–gene testing
  # -----------------
  unique_genes <- unique(gene_closer)
  gene_pairs   <- combn(unique_genes, 2, simplify = FALSE)
  n_gene_pairs <- length(gene_pairs)
  
  gene_results <- data.frame(
    gene1        = sapply(gene_pairs, `[`, 1),
    gene2        = sapply(gene_pairs, `[`, 2),
    max_raw_value = NA_real_,
    p_value_raw   = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (g in seq_len(n_gene_pairs)) {
    g1 <- gene_pairs[[g]][1]
    g2 <- gene_pairs[[g]][2]
    snps1 <- which(gene_closer == g1)
    snps2 <- which(gene_closer == g2)
    
    # Observed max absolute raw interaction
    original_raw <- compute_shap_raw_matrix(model, genotype_dat, snps1, snps2)
    obs_max_raw  <- max(abs(original_raw), na.rm = TRUE)
    
    # Permutation: shuffle within case/control for binary; shuffle all for continuous
    null_max_raw <- numeric(N)
    for (k in seq_len(N)) {
      perm_geno <- genotype_dat
      
      if (length(snps1) > 0) {
        if (is_binary) {
          perm_geno[cases,   snps1] <- perm_geno[sample(cases),   snps1]
          perm_geno[controls, snps1] <- perm_geno[sample(controls), snps1]
        } else {
          perm_geno[, snps1] <- perm_geno[sample(all_idx), snps1]
        }
      }
      if (length(snps2) > 0) {
        if (is_binary) {
          perm_geno[cases,   snps2] <- perm_geno[sample(cases),   snps2]
          perm_geno[controls, snps2] <- perm_geno[sample(controls), snps2]
        } else {
          perm_geno[, snps2] <- perm_geno[sample(all_idx), snps2]
        }
      }
      
      null_raw <- compute_shap_raw_matrix(model, perm_geno, snps1, snps2)
      null_max_raw[k] <- max(abs(null_raw), na.rm = TRUE)
    }
    
    gene_results$max_raw_value[g] <- obs_max_raw
    gene_results$p_value_raw[g]   <- (sum(null_max_raw >= obs_max_raw) + 1) / (N + 1)
  }
  
  # -----------------
  # SNP–SNP testing
  # -----------------
  snp_pairs    <- utils::combn(seq_len(n_features), 2)
  n_snp_pairs  <- ncol(snp_pairs)
  
  snp_results <- data.frame(
    snp1 = if (is.null(colnames(genotype_dat)))
      paste0("SNP", snp_pairs[1, ]) else colnames(genotype_dat)[snp_pairs[1, ]],
    snp2 = if (is.null(colnames(genotype_dat)))
      paste0("SNP", snp_pairs[2, ]) else colnames(genotype_dat)[snp_pairs[2, ]],
    raw_value   = numeric(n_snp_pairs),
    p_value_raw = numeric(n_snp_pairs),
    stringsAsFactors = FALSE
  )
  
  # Observed full matrix (raw)
  original_full_raw <- {
    dmatrix <- xgboost::xgb.DMatrix(data = genotype_dat)
    shap_int_3d <- predict(model, dmatrix, predinteraction = TRUE)
    apply(shap_int_3d, c(2, 3), mean)  # [p, p]
  }
  
  # Null array for raw interactions
  null_raw_array <- array(0, dim = c(N, n_features, n_features))
  
  for (k in seq_len(N)) {
    perm_geno <- genotype_dat
    if (is_binary) {
      for (col in seq_len(n_features)) {
        perm_geno[cases,   col] <- sample(perm_geno[cases,   col])
        perm_geno[controls, col] <- sample(perm_geno[controls, col])
      }
    } else {
      perm_geno <- perm_geno[sample(seq_len(n_samples)), , drop = FALSE]
    }
    
    dmatrix_null <- xgboost::xgb.DMatrix(data = perm_geno)
    shap_int_3d_null <- predict(model, dmatrix_null, predinteraction = TRUE)
    null_raw_array[k, , ] <- apply(shap_int_3d_null, c(2, 3), mean)
  }
  
  # Per-pair statistics and p-values (raw only)
  for (p in seq_len(n_snp_pairs)) {
    i <- snp_pairs[1, p]
    j <- snp_pairs[2, p]
    
    obs_raw     <- original_full_raw[i, j]
    obs_abs_raw <- abs(obs_raw)
    null_abs_raw <- abs(null_raw_array[, i, j])
    
    snp_results$raw_value[p]   <- obs_raw
    snp_results$p_value_raw[p] <- (sum(null_abs_raw >= obs_abs_raw) + 1) / (N + 1)
  }
  
  list(
    snp_results  = snp_results,
    gene_results = gene_results,
    null_raw_array = null_raw_array
  )
}