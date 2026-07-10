library(tidyverse)
library(gbmt)
library(mclust)

# ── Paths ──────────────────────────────────────────────────────────────────────
run_log_path <- "../../../../data/processed/chap/ensemble_50run/metadata/run_log_ecdf.csv"

# Output prefix — uses a distinct subdirectory "pooled/" to avoid clashing with
# outputs from x_gbmt_dta1_0220.R
prefix.date <- paste0(
  format(Sys.time(), "%b%d"), "/pooled_",
  format(Sys.time(), "%b%d"), "_"
)

# ── prep_long for chap data (sbp_bin / dbp_bin, 10 time points) ───────────────
prep_long_chap <- function(dat, is_synthetic = FALSE) {
  dat <- as.data.frame(dat)

  if (is_synthetic || !"Subject" %in% names(dat)) {
    dat$Subject <- seq_len(nrow(dat))
    if ("group" %in% names(dat) && !"Group" %in% names(dat))
      dat$Group <- dat$group
  } else {
    if (!"Group" %in% names(dat) && "group" %in% names(dat))
      dat$Group <- dat$group
  }

  dat_long <- dat %>%
    pivot_longer(
      cols          = matches("^(sbp_bin|dbp_bin)\\d+$"),
      names_to      = c(".value", "Time"),
      names_pattern = "^(sbp_bin|dbp_bin)(\\d+)$"
    ) %>%
    mutate(Time = as.numeric(Time)) %>%
    select(Subject, Time, Group, sbp_bin, dbp_bin)

  return(dat_long)
}

# ── Pool 5 imputed real datasets: average imputed values per Subject ───────────
# Same subjects appear in all 5 imputations; only the imputed (continuous)
# values differ. Taking the row mean per Subject is the standard way to
# collapse multiple imputations into one completed dataset when Rubin's rules
# cannot be applied directly (e.g. clustering).
pool_real_wide <- function(paths) {
  dfs <- lapply(seq_along(paths), function(i) {
    df <- read_csv(paths[i], show_col_types = FALSE)
    df$.imp <- i
    df
  })
  stacked <- bind_rows(dfs)
  num_cols <- names(stacked)[sapply(stacked, is.numeric) &
                               names(stacked) != ".imp"]

  stacked %>%
    group_by(Subject) %>%
    summarise(
      Group = first(Group),
      across(all_of(num_cols), \(x) mean(x, na.rm = TRUE)),
      .groups = "drop"
    )
}

# ── Pool 5 ecdf synthetic datasets: row-wise mean ─────────────────────────────
# Synthetic subjects have no real IDs. Each of the 5 files is an independent
# ECDF draw from the same training distribution, so averaging row-by-row gives
# one smoothed representative synthetic dataset of the same size.
pool_syn_wide <- function(paths) {
  dfs <- lapply(paths, \(p) read_csv(p, show_col_types = FALSE))

  dfs <- lapply(dfs, function(df) {
    if ("group" %in% names(df) && !"Group" %in% names(df))
      df$Group <- df$group
    df
  })

  num_cols <- names(dfs[[1]])[sapply(dfs[[1]], is.numeric)]

  num_mat <- Reduce("+", lapply(dfs, \(df) as.matrix(df[, num_cols]))) /
               length(dfs)
  result           <- dfs[[1]]
  result[, num_cols] <- as.data.frame(num_mat)
  result
}

# ── assign_test_subjects ───────────────────────────────────────────────────────
assign_test_subjects <- function(fit, test_dat,
                                  x_vars = c("sbp_bin", "dbp_bin")) {
  test_dat <- as.data.frame(test_dat)
  subjects <- unique(test_dat$Subject)
  n_subj   <- length(subjects)
  ng       <- length(fit$prior)

  log_dmvnorm0 <- function(x, Sigma) {
    p      <- length(x)
    R      <- chol(Sigma)
    logdet <- 2 * sum(log(diag(R)))
    z      <- backsolve(R, x, transpose = TRUE)
    quad   <- sum(z^2)
    -0.5 * (p * log(2 * pi) + logdet + quad)
  }

  log_post <- matrix(NA_real_, nrow = n_subj, ncol = ng)

  for (i in seq_len(n_subj)) {
    subj_dat <- test_dat[test_dat$Subject == subjects[i], ]
    time_vec <- subj_dat$Time
    X_poly   <- cbind(1, time_vec, time_vec^2)
    Y        <- as.matrix(subj_dat[, x_vars, drop = FALSE])

    for (k in seq_len(ng)) {
      beta_k    <- fit$beta[[k]]
      Sigma_k   <- fit$Sigma[[k]]
      E         <- Y - X_poly %*% beta_k
      log_lik_k <- 0
      for (t in seq_len(nrow(E))) {
        e_t <- E[t, ]
        if (!all(is.finite(e_t))) next
        log_lik_k <- log_lik_k + log_dmvnorm0(e_t, Sigma_k)
      }
      log_post[i, k] <- log(pmax(fit$prior[k], 1e-12)) + log_lik_k
    }
  }

  post_prob    <- exp(log_post - apply(log_post, 1, max))
  post_prob    <- post_prob / rowSums(post_prob)
  pred_cluster <- apply(post_prob, 1, which.max)

  data.frame(
    Subject           = subjects,
    Predicted_Cluster = pred_cluster,
    Post_Prob_Cl1     = round(post_prob[, 1], 4),
    Post_Prob_Cl2     = round(post_prob[, 2], 4)
  )
}

# ── best_mapping_ari ───────────────────────────────────────────────────────────
best_mapping_ari <- function(pred_labels, true_labels) {
  ari_a   <- adjustedRandIndex(pred_labels, true_labels)
  flipped <- ifelse(pred_labels == 1, 2, 1)
  ari_b   <- adjustedRandIndex(flipped, true_labels)

  if (ari_a >= ari_b) {
    list(ari = ari_a, mapped_pred = pred_labels,
         mapping = "1=Group0, 2=Group1")
  } else {
    list(ari = ari_b, mapped_pred = flipped,
         mapping = "1=Group1, 2=Group0 (flipped)")
  }
}

# ── run_gbmt_trts ──────────────────────────────────────────────────────────────
run_gbmt_trts <- function(train_dat, test_dat, dataset_name,
                           x_vars = c("sbp_bin", "dbp_bin"), seed = 123) {
  cat("\n===", dataset_name, "===\n")
  train_dat <- as.data.frame(train_dat)
  test_dat  <- as.data.frame(test_dat)
  set.seed(seed)

  cat("Fitting GBMT on training data (ng=2)...\n")
  fit <- gbmt(
    x.names = x_vars,
    unit    = "Subject",
    time    = "Time",
    ng      = 2,
    d       = 2,
    data    = train_dat,
    scaling = 0,
    nstart  = 3,
    maxit   = 1000,
    quiet   = FALSE
  )

  train_subjects  <- unique(train_dat$Subject)
  train_true      <- as.numeric(as.factor(
    train_dat$Group[!duplicated(train_dat$Subject)]
  ))
  train_map       <- best_mapping_ari(fit$assign, train_true)
  train_result_df <- data.frame(
    Subject           = train_subjects,
    True_Group        = train_dat$Group[!duplicated(train_dat$Subject)],
    Predicted_Cluster = train_map$mapped_pred
  )

  cat("Train ARI:", round(train_map$ari, 4),
      " | Mapping:", train_map$mapping, "\n")
  print(table(train_result_df$True_Group, train_result_df$Predicted_Cluster))

  cat("Assigning test subjects via posterior probabilities...\n")
  test_post  <- assign_test_subjects(fit, test_dat, x_vars)
  test_true  <- as.numeric(as.factor(
    test_dat$Group[!duplicated(test_dat$Subject)]
  ))
  test_map   <- best_mapping_ari(test_post$Predicted_Cluster, test_true)

  test_result_df <- data.frame(
    Subject           = unique(test_dat$Subject),
    True_Group        = test_dat$Group[!duplicated(test_dat$Subject)],
    Predicted_Cluster = test_map$mapped_pred,
    Post_Prob_G0      = test_post$Post_Prob_Cl1,
    Post_Prob_G1      = test_post$Post_Prob_Cl2
  )

  cat("Test  ARI:", round(test_map$ari, 4),
      " | Mapping:", test_map$mapping, "\n")
  print(table(test_result_df$True_Group, test_result_df$Predicted_Cluster))

  list(
    dataset      = dataset_name,
    fit          = fit,
    train_result = train_result_df,
    test_result  = test_result_df,
    train_ARI    = train_map$ari,
    test_ARI     = test_map$ari
  )
}

# ── Read run log ───────────────────────────────────────────────────────────────
run_log <- read_csv(run_log_path, show_col_types = FALSE)
cat("Run log loaded:", nrow(run_log), "rows\n")

# ── Collect paths across all 5 imputations ────────────────────────────────────
# NOTE: train/test split files do NOT contain a Subject column, so we pool
# from the full imputed files (stage == "imputed") which DO have Subject,
# then re-split here to match the original 80/20 ratio.
imp_paths <- run_log %>% filter(stage == "imputed") %>%
  arrange(imp_id) %>% pull(file_path)

# One ecdf synthetic per imputation (syn_id == 1)
ecdf_paths <- run_log %>% filter(stage == "ecdf_syn", syn_id == 1) %>%
  arrange(imp_id) %>% pull(file_path)

cat("Imputed files:", length(imp_paths),
    "| ECDF files:", length(ecdf_paths), "\n")

# ── Pool full imputed datasets → single pooled wide dataset ───────────────────
cat("\nPooling", length(imp_paths),
    "full imputed datasets (mean per Subject)...\n")
full_pooled_wide <- pool_real_wide(imp_paths)
cat("Pooled full dataset: ", nrow(full_pooled_wide), "subjects\n")

# ── 80/20 train/test split on pooled data (mirrors original split ratio) ──────
# Original: 2890 train / 722 test out of 3612 total ≈ 80 %
set.seed(42)
train_idx        <- sample(nrow(full_pooled_wide),
                           size  = round(nrow(full_pooled_wide) * 0.8))
train_pooled_wide <- full_pooled_wide[ train_idx, ]
test_pooled_wide  <- full_pooled_wide[-train_idx, ]
cat("Train:", nrow(train_pooled_wide), "| Test:", nrow(test_pooled_wide), "\n")

# ── Pool ecdf synthetic datasets (row-wise mean) ──────────────────────────────
cat("Pooling", length(ecdf_paths), "ecdf synthetic datasets (row-wise mean)...\n")
ecdf_pooled_wide <- pool_syn_wide(ecdf_paths)

# ── Convert to long format ────────────────────────────────────────────────────
train_dat <- prep_long_chap(train_pooled_wide)
test_dat  <- prep_long_chap(test_pooled_wide)
ecdf_dat  <- prep_long_chap(ecdf_pooled_wide, is_synthetic = TRUE)

cat("\nPooled train subjects:", n_distinct(train_dat$Subject),
    "| Time points:", paste(sort(unique(train_dat$Time)), collapse = " "), "\n")
cat("Pooled test  subjects:", n_distinct(test_dat$Subject), "\n")
cat("Pooled ecdf  subjects:", n_distinct(ecdf_dat$Subject), "\n")

# ── GBMT: Scenario A — TRTR (pooled real train → pooled real test) ────────────
result_real <- run_gbmt_trts(
  train_dat    = train_dat,
  test_dat     = test_dat,
  dataset_name = "TRTR (pooled real)",
  x_vars       = c("sbp_bin", "dbp_bin"),
  seed         = 42
)

# ── GBMT: Scenario B — ECDF synthetic → pooled real test ─────────────────────
result_syn <- run_gbmt_trts(
  train_dat    = ecdf_dat,
  test_dat     = test_dat,
  dataset_name = "ECDF synthetic → pooled real test",
  x_vars       = c("sbp_bin", "dbp_bin"),
  seed         = 42
)

# ── Summary comparison ─────────────────────────────────────────────────────────
gbmt_comparison <- data.frame(
  Scenario  = c("TRTR (Real)", "ECDF (Synthetic)"),
  Train_ARI = c(result_real$train_ARI, result_syn$train_ARI),
  Test_ARI  = c(result_real$test_ARI,  result_syn$test_ARI)
)

cat("\n====== ARI Comparison (pooled imputation) ======\n")
print(gbmt_comparison)

# ── Save results ───────────────────────────────────────────────────────────────
saveRDS(result_real,     paste0(prefix.date, "gbmt_pooled_real_result.RDS"))
saveRDS(result_syn,      paste0(prefix.date, "gbmt_pooled_syn_result.RDS"))
saveRDS(gbmt_comparison, paste0(prefix.date, "gbmt_pooled_comparison.RDS"))

cat("\nAll results saved!\n")
