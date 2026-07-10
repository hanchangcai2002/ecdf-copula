get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

InverseCDF.Sim2 <- function(Fx.matrix,           ## Correlated Unif matrix: n_sim*(t*p)
                            dat.cor,             ## Original longitudinal data
                            feature_cols,        ## Feature column names
                            t_vals) {            ## Time points
  
  n_sim <- dim(Fx.matrix)[1]
  t <- length(t_vals)
  p <- length(feature_cols)
  total_dim <- t * p
  
  if (ncol(Fx.matrix) != total_dim) {
    stop("Fx.matrix dimensions don't match expected total dimensions")
  }
  
  Sim_dat <- matrix(as.numeric(0), n_sim, total_dim)
  
  # For each feature at each time point
  col_idx <- 1
  for (i in 1:t) {
    time_point <- t_vals[i]
    
    for (j in 1:p) {
      feature_name <- feature_cols[j]
      
      # Extract data for this specific feature and time point
      time_data <- subset(dat.cor, months_to_base == time_point)
      
      if (nrow(time_data) == 0) {
        # If no data for this time point, use random uniform values
        Sim_dat[, col_idx] <- runif(n_sim, 0, 1)
        col_idx <- col_idx + 1
        next
      }
      
      # Get the feature values (remove NA)
      x <- time_data[[feature_name]]
      valid_idx <- !is.na(x)
      x_clean <- x[valid_idx]
      
      if (length(x_clean) == 0) {
        # If no valid data for this feature at this time, use random uniform
        Sim_dat[, col_idx] <- runif(n_sim, 0, 1)
        col_idx <- col_idx + 1
        next
      }
      
      # Calculate empirical CDF ranks for this feature at this time point
      ranks <- rank(x_clean, ties.method = "first")
      scale_ranks <- ranks / length(x_clean)
      
      # Get the uniform values for this column
      Fx_col <- Fx.matrix[, col_idx]
      
      # Perform inverse CDF transformation
      for (sim_idx in 1:n_sim) {
        if (is.na(Fx_col[sim_idx])) {
          Sim_dat[sim_idx, col_idx] <- sample(x_clean, 1)
          next
        }
        
        # Find values where Fx_col[sim_idx] <= scale_ranks
        valid_ranks <- which(Fx_col[sim_idx] <= scale_ranks)
        
        if (length(valid_ranks) == 0) {
          # If no valid ranks, use the maximum value
          Sim_dat[sim_idx, col_idx] <- max(x_clean)
        } else {
          # Find the closest rank
          candidate_ranks <- scale_ranks[valid_ranks]
          candidate_values <- x_clean[valid_ranks]
          
          diff_ranks <- abs(candidate_ranks - Fx_col[sim_idx])
          min_idx <- which.min(diff_ranks)
          Sim_dat[sim_idx, col_idx] <- candidate_values[min_idx]
        }
      }
      
      col_idx <- col_idx + 1
    }
  }
  
  return(Sim_dat)
}

Reshape.to.visits <- function(sim_result) {
  simulated_data <- sim_result$simulated_data
  t_vals <- sim_result$time_points
  features <- sim_result$features
  
  t <- length(t_vals)
  p <- length(features)
  n_sim <- nrow(simulated_data)
  
  visit_list <- list()
  
  for (i in 1:t) {
    col_start <- (i-1) * p + 1
    col_end <- i * p
    
    visit_data <- simulated_data[, col_start:col_end, drop = FALSE]
    colnames(visit_data) <- features
    
    visit_list[[paste0("visit_", t_vals[i])]] <- visit_data
  }
  
  return(visit_list)
}

# Step 1: Get within-visit correlation matrices (simplified as per Algorithm 2)
Get.within.visit.cor <- function(dat.cor, feature_cols) {
  t_vals <- sort(unique(dat.cor$months_to_base))
  t_vals <- t_vals[t_vals != 0]  # Remove baseline if needed
  t <- length(t_vals)
  p <- length(feature_cols)
  
  matrix_list <- list()
  
  for (time_point in t_vals) {
    tbl_cor <- subset(dat.cor, months_to_base == time_point)
    t_id <- tbl_cor$record_id
    
    mat_name <- paste0("matrix_visit_", time_point)
    
    if (length(t_id) < 5) {
      # Use identity matrix with small off-diagonal correlations
      cor_mat <- diag(p) + matrix(0.1, p, p) - diag(0.1, p)
      message(paste(mat_name, ": <5 subjects — using default correlation"))
    } else {
      dat_feat <- tbl_cor[, feature_cols, drop = FALSE]
      cor_mat <- cor(dat_feat, use = "pairwise.complete.obs", method = "spearman")
      cor_mat[is.na(cor_mat)] <- 0
      diag(cor_mat) <- 1  # Ensure diagonal is 1
    }
    
    matrix_list[[mat_name]] <- cor_mat
  }
  
  return(matrix_list)
}
# Step 2: Get between-visit correlation matrices
Get.between.visit.cor <- function(dat.cor, feature_cols) {
  t_vals <- sort(unique(dat.cor$months_to_base))
  t_vals <- t_vals[t_vals != 0]  # Remove baseline if needed
  t <- length(t_vals)
  p <- length(feature_cols)
  
  combinations <- t(combn(t_vals, 2))
  comb_dat <- data.frame(t1 = combinations[,1], t2 = combinations[,2])
  
  matrix_list <- list()
  
  for (n in 1:nrow(comb_dat)) {
    time_point_1 <- comb_dat$t1[n]
    time_point_2 <- comb_dat$t2[n]
    
    tbl_cor_1 <- subset(dat.cor, months_to_base == time_point_1)
    tbl_cor_2 <- subset(dat.cor, months_to_base == time_point_2)
    
    t1_id <- tbl_cor_1$record_id
    t2_id <- tbl_cor_2$record_id
    t12_id <- intersect(t1_id, t2_id)
    
    mat_name <- paste0("matrix_", time_point_1, "_", time_point_2)
    
    if (length(t12_id) < 5) {
      # Use small positive correlations for between-visit
      cor_mat <- matrix(0.2, p, p)
      message(paste(mat_name, ": <5 common subjects — using default correlation"))
    } else {
      dat1 <- tbl_cor_1[tbl_cor_1$record_id %in% t12_id, feature_cols, drop = FALSE]
      dat2 <- tbl_cor_2[tbl_cor_2$record_id %in% t12_id, feature_cols, drop = FALSE]
      
      cor_mat <- cor(dat1, dat2, use = "pairwise.complete.obs", method = "spearman")
      cor_mat[is.na(cor_mat)] <- 0.1  # Small positive correlation for missing
    }
    
    matrix_list[[mat_name]] <- cor_mat
  }
  
  return(matrix_list)
}
# Step 3: Assemble correlation matrices into block correlation matrix
Assemble.correlation.matrix <- function(within_visit_cors, between_visit_cors, t_vals, p) {
  t <- length(t_vals)
  total_dim <- t * p
  
  # Initialize the large correlation matrix
  big_cor_matrix <- matrix(0, nrow = total_dim, ncol = total_dim)
  
  # Fill diagonal blocks (within-visit correlations)
  for (i in 1:t) {
    time_point <- t_vals[i]
    mat_name <- paste0("matrix_visit_", time_point)
    
    row_start <- (i-1) * p + 1
    row_end <- i * p
    col_start <- (i-1) * p + 1
    col_end <- i * p
    
    if (mat_name %in% names(within_visit_cors)) {
      big_cor_matrix[row_start:row_end, col_start:col_end] <- within_visit_cors[[mat_name]]
    } else {
      # Default to identity if missing
      diag(big_cor_matrix[row_start:row_end, col_start:col_end]) <- 1
    }
  }
  
  # Fill off-diagonal blocks (between-visit correlations)
  for (i in 1:(t-1)) {
    for (j in (i+1):t) {
      time_point_1 <- t_vals[i]
      time_point_2 <- t_vals[j]
      mat_name <- paste0("matrix_", time_point_1, "_", time_point_2)
      
      row_start <- (i-1) * p + 1
      row_end <- i * p
      col_start <- (j-1) * p + 1
      col_end <- j * p
      
      if (mat_name %in% names(between_visit_cors)) {
        # Upper triangle
        big_cor_matrix[row_start:row_end, col_start:col_end] <- between_visit_cors[[mat_name]]
        # Lower triangle (transpose for symmetry)
        big_cor_matrix[col_start:col_end, row_start:row_end] <- t(between_visit_cors[[mat_name]])
      } else {
        # Default small correlation
        default_cor <- matrix(0.1, p, p)
        big_cor_matrix[row_start:row_end, col_start:col_end] <- default_cor
        big_cor_matrix[col_start:col_end, row_start:row_end] <- t(default_cor)
      }
    }
  }
  
  # Ensure positive definiteness
  eigenvals <- eigen(big_cor_matrix, only.values = TRUE)$values
  min_eigenval <- min(eigenvals)
  
  if (min_eigenval <= 0) {
    cat("Matrix not positive definite, adjusting...\n")
    # Add small positive value to diagonal
    diag(big_cor_matrix) <- diag(big_cor_matrix) + abs(min_eigenval) + 0.01
  }
  
  return(big_cor_matrix)
}
# Step 4: Extract parameter vector from correlation matrix (lower triangle)
Extract.param.vector <- function(big_cor_matrix) {
  # Get lower triangle
  cormat_lower <- as.matrix(get_lower_tri(big_cor_matrix))
  
  # Convert diagonal to NA
  diag(cormat_lower) <- NA
  
  # Melt and extract non-NA values
  melted_cormat_lower <- melt(cormat_lower)
  param.vec <- melted_cormat_lower$value[!is.na(melted_cormat_lower$value)]
  
  return(param.vec)
}

# main function
Sim.longitudinal.feature <- function(dat.cor, 
                                     feature_cols, 
                                     n_sim) {
  
  # Get unique time points
  t_vals <- sort(unique(dat.cor$months_to_base))
  t <- length(t_vals)
  p <- length(feature_cols)
  total_dim <- t * p
  
  cat("Processing", t, "time points and", p, "features\n")
  cat("Total dimensions:", total_dim, "\n")
  
  # Step 1: Get within-visit and between-visit correlations
  cat("Step 1: Computing correlations...\n")
  within_visit_cors <- Get.within.visit.cor(dat.cor, feature_cols)
  between_visit_cors <- Get.between.visit.cor(dat.cor, feature_cols)
  
  # Step 2: Assemble big correlation matrix
  cat("Step 2: Assembling correlation matrix...\n")
  big_cor_matrix <- Assemble.correlation.matrix(within_visit_cors, between_visit_cors, t_vals, p)
  param.vec <- Extract.param.vector(big_cor_matrix)
  
  cat("Parameter vector length:", length(param.vec), "\n")
  cat("Expected length:", total_dim * (total_dim - 1) / 2, "\n")
  
  # Step 3: Create copula and simulate
  cat("Step 3: Creating copula and simulating uniform variables...\n")
  tryCatch({
    mycop <- normalCopula(param = param.vec, dim = total_dim, dispstr = "un")
    
    el <- list(min = 0, max = 1)
    dups <- list(el)[rep(1, total_dim)]
    
    mydist <- mvdc(mycop, c(rep("unif", total_dim)), dups)
    Mv_unif <- rMvdc(n_sim, mydist)
    
  }, error = function(e) {
    cat("Error in copula creation:", e$message, "\n")
    cat("Using independent uniform variables instead...\n")
    Mv_unif <<- matrix(runif(n_sim * total_dim), nrow = n_sim, ncol = total_dim)
  })
  
  # Step 4: Inverse CDF transformation using the updated function
  cat("Step 4: Performing inverse CDF transformation...\n")
  Sim_feature <- InverseCDF.Sim2(Mv_unif,        ## Correlated Unif matrix: n_sim*total_dim
                                 dat.cor,        ## Original longitudinal data
                                 feature_cols,   ## Feature column names  
                                 t_vals)         ## Time points
  
  # Add column names (same format as before)
  col_names <- c()
  for (i in 1:t) {
    for (j in 1:p) {
      col_names <- c(col_names, paste0(feature_cols[j], "_t", t_vals[i]))
    }
  }
  colnames(Sim_feature) <- col_names
  
  cat("Simulation completed!\n")
  
  # Return results as a list (same format as before)
  return(list(
    simulated_data = Sim_feature,
    correlation_matrix = big_cor_matrix,
    within_visit_cors = within_visit_cors,
    between_visit_cors = between_visit_cors,
    original_data = dat.cor,  # Store original long format data instead
    time_points = t_vals,
    features = feature_cols
  ))
}