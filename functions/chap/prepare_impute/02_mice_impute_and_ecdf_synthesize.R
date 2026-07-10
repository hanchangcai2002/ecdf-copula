# =============================================================================
# CHAP: mice × 5  →  ECDF × 10 Pipeline (the M=5 x L=10 = 50-dataset ensemble used for Sec 3.3 "Variability")
# Run AFTER: 01_bin_visits.R
#
# Output folder structure (this repo's pre-generated copy lives under
# data/processed/chap/ensemble_50run/, with the "synthetic/" folder below
# named "synthetic_ecdf/" there to sit alongside "synthetic_ctgan/"):
#   imputed/   imp_01.csv … imp_05.csv
#   train/     train_imp_01.csv … train_imp_05.csv
#   test/      test_imp_01.csv  … test_imp_05.csv
#   synthetic/ imp_01/syn_01.csv … imp_05/syn_10.csv
#   metadata/  run_log.csv
#
# Seed strategy (fully deterministic from 3 constants below):
#   mice  : SEED_MICE
#   split : SEED_SPLIT_BASE + imp_id          (e.g. 1001 … 1005)
#   ecdf  : SEED_ECDF_BASE  + (imp_id-1)*N_SYN + syn_id  (e.g. 2001 … 2050)
# =============================================================================

rm(list = ls())

# ── 0. Libraries & Source ─────────────────────────────────────────────────────
library(dplyr)
library(mice)
library(reshape2)
library(copula)
source("../../synthesis_algorithms/algorithm2_longitudinal.R")

# ── 1. Constants ──────────────────────────────────────────────────────────────
SEED_MICE       <- 42
SEED_SPLIT_BASE <- 1000   # split seed for imp i  →  1000 + i
SEED_ECDF_BASE  <- 2000   # ecdf  seed for (imp i, syn j) → 2000 + (i-1)*N_SYN + j

N_IMP        <- 5
N_SYN        <- 10
SPLIT_RATIO  <- 0.8

feature_cols <- c(paste0("sbp_bin", 1:10), paste0("dbp_bin", 1:10))

# sbp/dbp are 2 longitudinal features observed at 10 time points (bin1..bin10);
# Sim.longitudinal.feature() expects long-format input (record_id,
# months_to_base, sbp, dbp) and returns time-major columns (sbp_t1, dbp_t1,
# sbp_t2, dbp_t2, ...). These helpers convert to/from that format.
to_long_bins <- function(df) {
  record_id <- seq_len(nrow(df))
  do.call(rbind, lapply(1:10, function(t) {
    data.frame(
      record_id      = record_id,
      months_to_base = t,
      sbp            = df[[paste0("sbp_bin", t)]],
      dbp            = df[[paste0("dbp_bin", t)]]
    )
  }))
}

reorder_longitudinal_sim <- function(sim_result) {
  sim_data <- as.data.frame(sim_result$simulated_data)
  out <- sim_data[, c(paste0("sbp_t", sim_result$time_points),
                      paste0("dbp_t", sim_result$time_points))]
  colnames(out) <- c(paste0("sbp_bin", 1:10), paste0("dbp_bin", 1:10))
  out
}

# ── 2. Setup Output Directories ───────────────────────────────────────────────
base_out <- "../../../data/processed/chap/ensemble_50run"

dirs <- c(
  file.path(base_out, "imputed"),
  file.path(base_out, "train"),
  file.path(base_out, "test"),
  file.path(base_out, "metadata"),
  file.path(base_out, "synthetic", sprintf("imp_%02d", seq_len(N_IMP)))
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── 3. Load Pre-imputation Data (from 1_binned.R) ───────────────
chap_flat <- readRDS("../../../data/raw/chap/chap_flat_preimpute.RDS")
cat("Loaded chap_flat:", nrow(chap_flat), "subjects\n")
cat("Missing cells before imputation:",
    sum(is.na(chap_flat[, feature_cols])), "\n")

# ── 4. Run mice (m = N_IMP, single call → all imputations share one seed) ─────
cat("\nRunning mice (m =", N_IMP, ", seed =", SEED_MICE, ") ...\n")
mids_obj <- mice(
  chap_flat[, feature_cols],
  method    = "pmm",
  m         = N_IMP,
  seed      = SEED_MICE,
  printFlag = FALSE
)
cat("mice done.\n")

# Initialise run_log
run_log <- data.frame(
  stage     = "mice",
  imp_id    = NA_integer_,
  syn_id    = NA_integer_,
  seed      = SEED_MICE,
  n_rows    = nrow(chap_flat),
  n_group0  = NA_integer_,
  n_group1  = NA_integer_,
  file_path = NA_character_,
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  stringsAsFactors = FALSE
)

# ── 5. Loop: imputation → split → ECDF × N_SYN ────────────────────────────────
for (i in seq_len(N_IMP)) {

  imp_tag <- sprintf("imp_%02d", i)
  cat("\n══ Imputation", i, "════════════════════════════════════════\n")

  # 5a. Extract imputed dataset --------------------------------------------------
  imp_feat              <- complete(mids_obj, i)        # feature columns only
  imp_df                <- chap_flat
  imp_df[, feature_cols] <- imp_feat

  imp_df_out        <- imp_df
  imp_df_out$Group  <- ifelse(imp_df_out$Group == "0", "Group0", "Group1")

  write.csv(
    imp_df_out,
    file.path(base_out, "imputed", paste0(imp_tag, ".csv")),
    row.names = FALSE
  )
  cat("  Saved imputed dataset:", imp_tag, "\n")

  run_log <- rbind(run_log, data.frame(
    stage = "imputed", imp_id = i, syn_id = NA_integer_,
    seed = SEED_MICE, n_rows = nrow(imp_df),
    n_group0 = sum(imp_df$Group == "0"),
    n_group1 = sum(imp_df$Group == "1"),
    file_path = normalizePath(file.path(base_out, "imputed", paste0(imp_tag, ".csv")), mustWork = FALSE),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  ))

  # 5b. 80:20 train/test split ---------------------------------------------------
  split_seed <- SEED_SPLIT_BASE + i
  set.seed(split_seed)
  # Stratified split: sample 80% within each group independently
  train_ids <- unlist(lapply(split(imp_df$Subject, imp_df$Group), function(ids) {
    sample(ids, round(SPLIT_RATIO * length(ids)))
  }))
  test_ids  <- setdiff(imp_df$Subject, train_ids)

  train_df  <- imp_df[imp_df$Subject %in% train_ids, ]
  test_df   <- imp_df[imp_df$Subject %in% test_ids,  ]

  cat("  Split (seed", split_seed, "): train =", nrow(train_df),
      "| test =", nrow(test_df), "\n")

  # Save train and test (group + features; no Subject ID for downstream eval)
  fmt <- function(df) {
    df$Group <- ifelse(df$Group == "0", "Group0", "Group1")
    df[, c("Group", feature_cols)]
  }
  write.csv(fmt(train_df),
            file.path(base_out, "train", paste0("train_", imp_tag, ".csv")),
            row.names = FALSE)
  write.csv(fmt(test_df),
            file.path(base_out, "test",  paste0("test_",  imp_tag, ".csv")),
            row.names = FALSE)

  run_log <- rbind(run_log, data.frame(
    stage = "train", imp_id = i, syn_id = NA_integer_,
    seed = split_seed, n_rows = nrow(train_df),
    n_group0 = sum(train_df$Group == "0"),
    n_group1 = sum(train_df$Group == "1"),
    file_path = normalizePath(file.path(base_out, "train", paste0("train_", imp_tag, ".csv")), mustWork = FALSE),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  ))
  run_log <- rbind(run_log, data.frame(
    stage = "test", imp_id = i, syn_id = NA_integer_,
    seed = split_seed, n_rows = nrow(test_df),
    n_group0 = sum(test_df$Group == "0"),
    n_group1 = sum(test_df$Group == "1"),
    file_path = normalizePath(file.path(base_out, "test", paste0("test_", imp_tag, ".csv")), mustWork = FALSE),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  ))

  # 5c. ECDF synthesis × N_SYN --------------------------------------------------
  g0_train <- train_df[train_df$Group == "0", ]
  g1_train <- train_df[train_df$Group == "1", ]

  for (j in seq_len(N_SYN)) {

    ecdf_seed <- SEED_ECDF_BASE + (i - 1) * N_SYN + j
    set.seed(ecdf_seed)

    sim_g0 <- Sim.longitudinal.feature(
      dat.cor      = to_long_bins(g0_train),
      feature_cols = c("sbp", "dbp"),
      n_sim        = nrow(g0_train)
    )
    sim_g1 <- Sim.longitudinal.feature(
      dat.cor      = to_long_bins(g1_train),
      feature_cols = c("sbp", "dbp"),
      n_sim        = nrow(g1_train)
    )

    syn_g0 <- reorder_longitudinal_sim(sim_g0)
    syn_g1 <- reorder_longitudinal_sim(sim_g1)

    syn_df <- rbind(
      cbind(group = "Group0", syn_g0),
      cbind(group = "Group1", syn_g1)
    )

    syn_file <- file.path(base_out, "synthetic", imp_tag,
                          sprintf("syn_%02d.csv", j))
    write.csv(syn_df, syn_file, row.names = FALSE)

    run_log <- rbind(run_log, data.frame(
      stage = "ecdf_syn", imp_id = i, syn_id = j,
      seed = ecdf_seed, n_rows = nrow(syn_df),
      n_group0 = sum(syn_df$group == "Group0"),
      n_group1 = sum(syn_df$group == "Group1"),
      file_path = normalizePath(syn_file, mustWork = FALSE),
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      stringsAsFactors = FALSE
    ))

    cat("  syn", sprintf("%02d", j), "(seed", ecdf_seed, ") →", basename(syn_file), "\n")
  }
}

# ── 6. Save run_log ───────────────────────────────────────────────────────────
log_path <- file.path(base_out, "metadata", "run_log.csv")
write.csv(run_log, log_path, row.names = FALSE)

cat("\n══ Done ══════════════════════════════════════════════════════\n")
cat("Imputed datasets  :", N_IMP, "\n")
cat("Synthetic datasets:", N_IMP * N_SYN, "(", N_IMP, "imp ×", N_SYN, "syn )\n")
cat("run_log saved →", log_path, "\n")