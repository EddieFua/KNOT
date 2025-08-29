compute_shap_interaction_pvalues <- function(genotype_dat, y, gene_closer, N = 100, cores = 10, seed = 10) {
  set.seed(seed)
  n_samples <- nrow(genotype_dat)
  n_features <- ncol(genotype_dat)
  if (length(unique(y)) == 2 && all(y %in% c(0, 1))) {
    objective <- "binary:logistic"
    eval_metric <- "logloss"
  } else {
    objective <- "reg:squarederror"
    eval_metric <- "rmse"
  }
  params <- list(
    objective = objective,
    eta = 0.3,
    max_depth = 3,
    eval_metric = eval_metric
  )
  dtrain <- xgboost::xgb.DMatrix(data = genotype_dat, label = y)
  cv <- xgboost::xgb.cv(params = params, data = dtrain, nrounds = 200, nfold = 5,
                        early_stopping_rounds = 10, verbose = 0)
  nrounds <- cv$best_iteration
  if (is.null(nrounds) || nrounds == 0) nrounds <- 200
  model <- xgboost::xgb.train(params = params, data = dtrain, nrounds = nrounds)
  compute_shap_matrices <- function(model, data_matrix, snps1, snps2) {
    dmatrix <- xgboost::xgb.DMatrix(data = data_matrix)
    shap_int_3d <- predict(model, dmatrix, predinteraction = TRUE)
    mean_shap <- apply(shap_int_3d, c(2, 3), mean)
    raw_matrix <- matrix(0, nrow = length(snps1), ncol = length(snps2))
    std_matrix <- matrix(0, nrow = length(snps1), ncol = length(snps2))
    for (i in seq_along(snps1)) {
      for (j in seq_along(snps2)) {
        idx1 <- snps1[i]
        idx2 <- snps2[j]
        raw_matrix[i, j] <- mean_shap[idx1, idx2]
        denom <- sqrt(abs(mean_shap[idx1, idx1]) * abs(mean_shap[idx2, idx2]) + 1e-6)
        std_matrix[i, j] <- raw_matrix[i, j] / denom
      }
    }
    return(list(raw_matrix = raw_matrix, std_matrix = std_matrix))
  }
  
  unique_genes <- unique(gene_closer)
  gene_pairs <- combn(unique_genes, 2, simplify = FALSE)
  n_gene_pairs <- length(gene_pairs)
  gene_results <- data.frame(
    gene1 = sapply(gene_pairs, `[`, 1),
    gene2 = sapply(gene_pairs, `[`, 2),
    max_std_value = NA_real_,
    max_raw_value = NA_real_,
    p_value_std = NA_real_,
    p_value_raw = NA_real_
  )
  cases <- which(y == 1)
  controls <- which(y == 0)
  cl <- makePSOCKcluster(cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  for (g in 1:n_gene_pairs) {
    g1 <- gene_pairs[[g]][1]
    g2 <- gene_pairs[[g]][2]
    snps1 <- which(gene_closer == g1)
    snps2 <- which(gene_closer == g2)
    original <- compute_shap_matrices(model, genotype_dat, snps1, snps2)
    original_std <- original$std_matrix
    original_raw <- original$raw_matrix
    obs_max_std <- max(abs(original_std), na.rm = TRUE)
    obs_max_raw <- max(abs(original_raw), na.rm = TRUE)
    null_max <- foreach(k = 1:N, .combine = 'rbind', .packages = "xgboost") %dopar% {
      perm_geno <- genotype_dat
      if (length(snps1) > 0) {
        perm_indices_cases <- sample(length(cases))
        perm_indices_controls <- sample(length(controls))
        perm_geno[cases, snps1] <- perm_geno[cases[perm_indices_cases], snps1]
        perm_geno[controls, snps1] <- perm_geno[controls[perm_indices_controls], snps1]
      }
      if (length(snps2) > 0) {
        perm_indices_cases <- sample(length(cases))
        perm_indices_controls <- sample(length(controls))
        perm_geno[cases, snps2] <- perm_geno[cases[perm_indices_cases], snps2]
        perm_geno[controls, snps2] <- perm_geno[controls[perm_indices_controls], snps2]
      }
      null <- compute_shap_matrices(model, perm_geno, snps1, snps2)
      c(max(abs(null$std_matrix), na.rm = TRUE), max(abs(null$raw_matrix), na.rm = TRUE))
    }
    null_max_std <- null_max[, 1]
    null_max_raw <- null_max[, 2]
    
    gene_results$max_std_value[g] <- obs_max_std
    gene_results$max_raw_value[g] <- obs_max_raw
    gene_results$p_value_std[g] <- (sum(null_max_std >= obs_max_std) + 1) / (N + 1)
    gene_results$p_value_raw[g] <- (sum(null_max_raw >= obs_max_raw) + 1) / (N + 1)
  }
  snp_pairs <- combn(1:n_features, 2)
  n_snp_pairs <- ncol(snp_pairs)
  snp_results <- data.frame(
    snp1 = ifelse(is.null(colnames(genotype_dat)), paste0("SNP", snp_pairs[1, ]), colnames(genotype_dat)[snp_pairs[1, ]]),
    snp2 = ifelse(is.null(colnames(genotype_dat)), paste0("SNP", snp_pairs[2, ]), colnames(genotype_dat)[snp_pairs[2, ]]),
    std_value = numeric(n_snp_pairs),
    raw_value = numeric(n_snp_pairs),
    p_value_std = numeric(n_snp_pairs),
    p_value_raw = numeric(n_snp_pairs)
  )
  original_full <- compute_shap_matrices(model, genotype_dat, 1:n_features, 1:n_features)
  original_full_std <- original_full$std_matrix
  original_full_raw <- original_full$raw_matrix
  null_std_array <- array(0, dim = c(N, n_features, n_features))
  null_raw_array <- array(0, dim = c(N, n_features, n_features))
  for (k in 1:N) {
    perm_geno <- genotype_dat
    for (col in 1:n_features) {
      perm_geno[cases, col] <- sample(perm_geno[cases, col])
      perm_geno[controls, col] <- sample(perm_geno[controls, col])
    }
    null <- compute_shap_matrices(model, perm_geno, 1:n_features, 1:n_features)
    null_std_array[k, , ] <- null$std_matrix
    null_raw_array[k, , ] <- null$raw_matrix
  }
  for (p in 1:n_snp_pairs) {
    i <- snp_pairs[1, p]
    j <- snp_pairs[2, p]
    obs_std <- original_full_std[i, j]
    obs_abs_std <- abs(obs_std)
    null_abs_std <- abs(null_std_array[, i, j])
    snp_results$std_value[p] <- obs_std
    snp_results$p_value_std[p] <- (sum(null_abs_std >= obs_abs_std) + 1) / (N + 1)
    obs_raw <- original_full_raw[i, j]
    obs_abs_raw <- abs(obs_raw)
    null_abs_raw <- abs(null_raw_array[, i, j])
    snp_results$raw_value[p] <- obs_raw
    snp_results$p_value_raw[p] <- (sum(null_abs_raw >= obs_abs_raw) + 1) / (N + 1)
  }
  
  return(list(snp_results = snp_results, gene_results = gene_results, null_raw_array = null_raw_array, null_std_array = null_std_array))
}