library(tidyverse)
library(gbmt)
library(mclust)
paths <- list(
  real_train = "../../../../data/raw/rembrandt/0_Depression_Data_Real_Train.csv",
  real_test  = "../../../../data/raw/rembrandt/0_Depression_Data_Real_Test.csv",
  ecdf_syn   = "../../../../data/processed/rembrandt/0_Depression_Data_Synthetic_ECDF.csv",
  sdv_syn    = "../../../../data/processed/rembrandt/0_Depression_Data_Synthetic_SDV.csv"
)
prefix.date <- paste0(format(Sys.time(), "%b%d"), "_")
prep_long <- function(dat) {
  dat <- as.data.frame(dat)
  
  # Use row number as Subject
  dat$Subject <- seq_len(nrow(dat))

  # Pivot ma1~ma14, hars1~hars14 into long format
  dat_long <- dat %>%
    pivot_longer(
      cols      = matches("^(ma|hars)\\d+$"),
      names_to  = c(".value", "Time"),
      names_pattern = "^(ma|hars)(\\d+)$"
    ) %>%
    mutate(
      Time  = as.numeric(Time),
      Group = group          # keep the original group column, renamed to Group
    ) %>%
    select(Subject, Time, Group, ma, hars)
  
  return(dat_long)
}
assign_test_subjects <- function(fit, test_dat, x_vars = c("ma","hars")) {
  test_dat  <- as.data.frame(test_dat)
  subjects  <- unique(test_dat$Subject)
  n_subj    <- length(subjects)
  ng        <- length(fit$prior)
  
  # log density of N(0, Sigma) for a vector x (length p), no external package
  log_dmvnorm0 <- function(x, Sigma) {
    p <- length(x)
    # numerically stable inverse + logdet via chol
    R <- chol(Sigma)
    logdet <- 2 * sum(log(diag(R)))
    z <- backsolve(R, x, transpose = TRUE)  # solves R' z = x
    quad <- sum(z^2)
    -0.5 * (p * log(2*pi) + logdet + quad)
  }
  
  log_post <- matrix(NA_real_, nrow = n_subj, ncol = ng)
  
  for (i in seq_len(n_subj)) {
    subj_dat <- test_dat[test_dat$Subject == subjects[i], ]
    time_vec <- subj_dat$Time
    X_poly   <- cbind(1, time_vec, time_vec^2)  # 14x3
    
    # observed matrix Y: 14x2 (ma,hars) in the given order
    Y <- as.matrix(subj_dat[, x_vars, drop = FALSE])
    
    for (k in seq_len(ng)) {
      beta_k  <- fit$beta[[k]]   # 3x2
      Sigma_k <- fit$Sigma[[k]]  # 2x2
      
      # predicted mean: 14x2
      Mu <- X_poly %*% beta_k
      
      # residuals: 14x2
      E <- Y - Mu
      
      # sum loglik over time points, skipping rows with NA
      log_lik_k <- 0
      for (t in seq_len(nrow(E))) {
        e_t <- E[t, ]
        ok <- is.finite(e_t)
        if (!all(ok)) next   # skip this time point if ma/hars has any NA
        log_lik_k <- log_lik_k + log_dmvnorm0(e_t, Sigma_k)
      }
      
      log_post[i, k] <- log(pmax(fit$prior[k], 1e-12)) + log_lik_k
    }
  }
  
  # softmax -> posterior
  post_prob <- exp(log_post - apply(log_post, 1, max))
  post_prob <- post_prob / rowSums(post_prob)
  
  pred_cluster <- apply(post_prob, 1, which.max)
  
  data.frame(
    Subject = subjects,
    Predicted_Cluster = pred_cluster,
    Post_Prob_Cl1 = round(post_prob[, 1], 4),
    Post_Prob_Cl2 = round(post_prob[, 2], 4)
  )
}
best_mapping_ari <- function(pred_labels, true_labels) {
  # mapping A: pred 1 -> label 1, pred 2 -> label 2 (original)
  ari_a <- adjustedRandIndex(pred_labels, true_labels)

  # mapping B: predicted labels flipped
  flipped <- ifelse(pred_labels == 1, 2, 1)
  ari_b   <- adjustedRandIndex(flipped, true_labels)
  
  if (ari_a >= ari_b) {
    list(ari = ari_a, mapped_pred = pred_labels,
         mapping = "1=Group1, 2=Group2")
  } else {
    list(ari = ari_b, mapped_pred = flipped,
         mapping = "1=Group2, 2=Group1 (flipped)")
  }
}

# ── Core function: fit GBMT on training data, then assign test subjects via posterior probability ─────────────────
run_gbmt_trts <- function(train_dat, test_dat, dataset_name,
                          x_vars = c("ma", "hars"), seed = 123) {
  cat("\n===", dataset_name, "===\n")
  train_dat <- as.data.frame(train_dat)
  test_dat  <- as.data.frame(test_dat)
  set.seed(seed)
  
  # 1. Fit k=2 GBMT on the training data
  cat("Fitting GBMT on training data (k=2)...\n")
  fit <- gbmt(
    x.names = x_vars,
    unit     = "Subject",
    time     = "Time",
    ng       = 2,
    d        = 2,
    data     = train_dat,
    scaling  = 0,
    nstart   = 3,
    maxit    = 1000,
    quiet    = FALSE
  )
  
  # 2. Training set cluster assignments (given directly by fit$assign)
  train_subjects  <- unique(train_dat$Subject)
  train_true      <- as.numeric(as.factor(train_dat$Group[!duplicated(train_dat$Subject)]))
  train_pred      <- fit$assign
  
  train_map       <- best_mapping_ari(train_pred, train_true)
  train_result_df <- data.frame(
    Subject           = train_subjects,
    True_Group        = train_dat$Group[!duplicated(train_dat$Subject)],
    Predicted_Cluster = train_map$mapped_pred
  )
  
  cat("Train ARI:", round(train_map$ari, 4),
      " | Mapping:", train_map$mapping, "\n")
  print(table(train_result_df$True_Group, train_result_df$Predicted_Cluster))
  
  # 3. Assign test data to clusters via posterior probability, using the fixed model parameters
  cat("Assigning test subjects via posterior probabilities...\n")
  test_post <- assign_test_subjects(fit, test_dat, x_vars)
  
  test_true      <- as.numeric(as.factor(test_dat$Group[!duplicated(test_dat$Subject)]))
  test_map       <- best_mapping_ari(test_post$Predicted_Cluster, test_true)
  
  test_result_df <- data.frame(
    Subject           = unique(test_dat$Subject),
    True_Group        = test_dat$Group[!duplicated(test_dat$Subject)],
    Predicted_Cluster = test_map$mapped_pred,
    Post_Prob_G1      = test_post$Post_Prob_Cl1,
    Post_Prob_G2      = test_post$Post_Prob_Cl2
  )
  
  cat("Test  ARI:", round(test_map$ari, 4),
      " | Mapping:", test_map$mapping, "\n")
  print(table(test_result_df$True_Group, test_result_df$Predicted_Cluster))
  
  list(
    dataset       = dataset_name,
    fit           = fit,
    train_result  = train_result_df,
    test_result   = test_result_df,
    train_ARI     = train_map$ari,
    test_ARI      = test_map$ari
  )
}

# ── Data preprocessing: wide -> long, row number used as Subject ─────────────
# ── Read and preprocess all datasets ─────────────────────────────────────────
real_train <- prep_long(read_csv(paths$real_train))
real_test  <- prep_long(read_csv(paths$real_test))
ecdf_syn   <- prep_long(read_csv(paths$ecdf_syn))
sdv_syn    <- prep_long(read_csv(paths$sdv_syn))

# Quick sanity check
cat("Train columns:", paste(colnames(real_train), collapse = ", "), "\n")
cat("Train n subjects:", n_distinct(real_train$Subject), "\n")
cat("Time points:", sort(unique(real_train$Time)), "\n")
print(head(real_train))
head(real_train)
head(real_test)
head(ecdf_syn)
head(sdv_syn)
# Quick check of column names (ensure Subject / Time / Group / ma / hars are present)
cat("Train columns:", paste(colnames(real_train), collapse = ", "), "\n")
cat("Test  columns:", paste(colnames(real_test),  collapse = ", "), "\n")

# ── Run the 3 scenarios ──────────────────────────────────────────────────────
# Scenario 1: ECDF synthetic train → Real test
result_ecdf <- run_gbmt_trts(
  train_dat    = ecdf_syn,
  test_dat     = real_test,
  dataset_name = "ECDF_synthetic → Real_test",
  seed         = 42
)

# Scenario 2: SDV/SVM synthetic train → Real test
result_sdv <- run_gbmt_trts(
  train_dat    = sdv_syn,
  test_dat     = real_test,
  dataset_name = "SDV_synthetic → Real_test",
  seed         = 42
)

# Scenario 3: TRTR — Real train → Real test
result_trtr <- run_gbmt_trts(
  train_dat    = real_train,
  test_dat     = real_test,
  dataset_name = "TRTR (Real_train → Real_test)",
  seed         = 42
)

# ── Summary comparison ───────────────────────────────────────────────────────
gbmt.comparison <- data.frame(
  Scenario      = c("ECDF → Real_test", "SDV → Real_test", "TRTR"),
  Train_ARI     = c(result_ecdf$train_ARI, result_sdv$train_ARI, result_trtr$train_ARI),
  Test_ARI      = c(result_ecdf$test_ARI,  result_sdv$test_ARI,  result_trtr$test_ARI)
)

cat("\n====== ARI Summary ======\n")
print(gbmt.comparison)

# ── Save results ──────────────────────────────────────────────────────────────
# The `fit` component of each result (the full gbmt model object) is not saved:
# it is not CSV-representable and is not needed to reproduce the ARI values
# below, which are already fully determined by train_result/test_result.
write.csv(result_ecdf$train_result, paste0(prefix.date, "gbmt_ecdf_train_result.csv"), row.names = FALSE)
write.csv(result_ecdf$test_result,  paste0(prefix.date, "gbmt_ecdf_test_result.csv"),  row.names = FALSE)
write.csv(result_sdv$train_result,  paste0(prefix.date, "gbmt_sdv_train_result.csv"),  row.names = FALSE)
write.csv(result_sdv$test_result,   paste0(prefix.date, "gbmt_sdv_test_result.csv"),   row.names = FALSE)
write.csv(result_trtr$train_result, paste0(prefix.date, "gbmt_trtr_train_result.csv"), row.names = FALSE)
write.csv(result_trtr$test_result,  paste0(prefix.date, "gbmt_trtr_test_result.csv"),  row.names = FALSE)
write.csv(gbmt.comparison,          paste0(prefix.date, "gbmt_comparison.csv"),        row.names = FALSE)

cat("\nAll results saved!\n")