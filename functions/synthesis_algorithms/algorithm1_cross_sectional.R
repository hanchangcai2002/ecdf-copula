library(reshape2)
library(copula)
library(matrixStats)
# Inverse eCDF transformation
InverseCDF.Sim <- function(Fx.matrix, x.matrix, scale.rank.matrix) {
  n_sim <- nrow(Fx.matrix)
  p <- ncol(Fx.matrix)
  
  Sim_dat <- matrix(0, n_sim, p)
  
  for (j in 1:p) {
    Fx <- Fx.matrix[, j]
    x <- x.matrix[, j]
    scale.rank.x <- scale.rank.matrix[, j]
    
    # Remove NA values
    valid_idx <- !is.na(x) & !is.na(scale.rank.x)
    x <- x[valid_idx]
    scale.rank.x <- scale.rank.x[valid_idx]
    
    if (length(x) == 0) {
      Sim_dat[, j] <- runif(n_sim, 0, 1)
      next
    }
    
    for (i in 1:n_sim) {
      if (is.na(Fx[i])) {
        Sim_dat[i, j] <- sample(x, 1)
        next
      }
      
      valid_ranks <- which(Fx[i] <= scale.rank.x)
      
      if (length(valid_ranks) == 0) {
        Sim_dat[i, j] <- max(x)
      } else {
        floorrk <- scale.rank.x[valid_ranks]
        floorx <- x[valid_ranks]
        diff_ranks <- abs(floorrk - Fx[i])
        min_idx <- which.min(diff_ranks)
        Sim_dat[i, j] <- floorx[min_idx]
      }
    }
  }
  return(Sim_dat)
}
Sim.cross.sectional.feature <- function(dat, 
                                        feature_cols, 
                                        n_sim) {
  
  p <- length(feature_cols)
  cat("Processing", p, "features\n")

  cat("Step 1: Computing Spearman correlation matrix...\n")
  
  dat_feat <- dat[, feature_cols, drop = FALSE]
  
  # Remove rows with all NA
  valid_rows <- apply(dat_feat, 1, function(x) !all(is.na(x)))
  dat_feat_clean <- dat_feat[valid_rows, , drop = FALSE]
  
  if (nrow(dat_feat_clean) < 5) {
    stop("Insufficient data: need at least 5 observations")
  }
  
  # spearman cor
  cor_matrix <- cor(dat_feat_clean, use = "pairwise.complete.obs", method = "spearman")
  cor_matrix[is.na(cor_matrix)] <- 0
  diag(cor_matrix) <- 1
  
  # Ensure positive definiteness
  eigenvals <- eigen(cor_matrix, only.values = TRUE)$values
  min_eigenval <- min(eigenvals)
  
  if (min_eigenval <= 0) {
    cat("Correlation matrix not positive definite, adjusting...\n")
    diag(cor_matrix) <- diag(cor_matrix) + abs(min_eigenval) + 0.01
  }
  
  cat("Step 2: Extracting copula parameters...\n")
  cormat_lower <- cor_matrix
  cormat_lower[upper.tri(cor_matrix)] <- NA
  diag(cormat_lower) <- NA
  
  melted <- melt(cormat_lower)
  param.vec <- melted$value[!is.na(melted$value)]
  
  cat("Parameter vector length:", length(param.vec), "\n")
  

  cat("Step 3: Generating correlated uniform samples via Gaussian copula...\n")
  tryCatch({
    mycop <- normalCopula(param = param.vec, dim = p, dispstr = "un")
    el <- list(min = 0, max = 1)
    dups <- list(el)[rep(1, p)]
    mydist <- mvdc(mycop, rep("unif", p), dups)
    Mv_unif <- rMvdc(n_sim, mydist)}, error = function(e) {
    cat("Error in copula creation:", e$message, "\n")
    cat("Using independent uniform variables instead...\n")
    Mv_unif <<- matrix(runif(n_sim * p), nrow = n_sim, ncol = p)
  })
  
  cat("Step 4: Computing empirical CDF ranks...\n")
  original_matrix <- as.matrix(dat_feat_clean)
  n_ori <- nrow(original_matrix)
  
  scale.rank.matrix <- matrix(0, nrow = n_ori, ncol = p)
  
  for (j in 1:p) {
    col_data <- original_matrix[, j]
    valid_data <- !is.na(col_data)
    
    if (sum(valid_data) > 0) {
      ranks <- rank(col_data[valid_data], ties.method = "first")
      scale.rank.matrix[valid_data, j] <- ranks / sum(valid_data)
    }
  }
  
  cat("Step 5: Performing inverse eCDF transformation...\n")
  Sim_feature <- InverseCDF.Sim(Mv_unif, original_matrix, scale.rank.matrix)
  colnames(Sim_feature) <- feature_cols
  cat("Simulation completed!\n")
  return(list(
    simulated_data = Sim_feature,
    correlation_matrix = cor_matrix,
    original_data = original_matrix,
    features = feature_cols,
    n_sim = n_sim
  ))
}


# Example Usage
# set.seed(123)
# base_path <- "/Users/hazel/Dropbox/jinyuan/Project2/sce2/impt"
# dat <- readRDS(file.path(base_path, "impt2.RDS"))
# dat <- subset(dat, visit == "baseline")
# feature_cols <- c("ma_tot", "hars_score")
# result <- Sim.cross.sectional.feature(
#   dat = dat_baseline,
#   feature_cols = feature_cols,
#   n_sim = 200)
# 
# head(result$simulated_data)
# print(cor(result$original_data, use = "complete.obs", method = "spearman"))
# print(cor(result$simulated_data, method = "spearman"))
