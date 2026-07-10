# =============================================================================
# CHAP Binned Average Pipeline (no FPCA)
# Dataset: Treatment for Mild Chronic Hypertension during Pregnancy (CHAP trial)
# Steps:
#   realworld.RDS
#     → Bin weeks 1-40 into 10 × 4-week bins
#       - Average sbp/dbp within each bin if ≥1 obs exists; NA if none
#     → 80:20 train/test split (subject-level)
#     → ECDF simulation on training set (separately by Group)
# Features: sbp_bin1 … sbp_bin10, dbp_bin1 … dbp_bin10  (20 total)
# =============================================================================

rm(list = ls())

# ── 0. Libraries & Source ─────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(mice)
source("../../synthesis_algorithms/algorithm1_cross_sectional.R")

# ── 1. Load Data ──────────────────────────────────────────────────────────────
# chap.csv columns: Subject, Time (week), Group (= chtndx_any), sbp, dbp.
# sbp/dbp are already z-scored.
chap_long <- read.csv("../../../data/raw/chap/chap.csv")
chap_long$Subject <- as.character(chap_long$Subject)

# Keep only weeks 2-40 (Time=0 baseline excluded; ceiling(0/4)=0 otherwise)
chap_long <- chap_long %>% filter(Time >= 2, Time <= 40)

# ── 2. Bin into 4-Week Intervals ──────────────────────────────────────────────
# bin 1: weeks 1-4, bin 2: weeks 5-8, ..., bin 10: weeks 37-40
chap_long <- chap_long %>%
  mutate(bin = ceiling(Time / 4))   # 1 … 10

# Per subject × bin: average sbp/dbp (ignore within-bin NAs)
# If a subject has no observations in a bin, that bin is NA in the wide table
chap_binned <- chap_long %>%
  group_by(Subject, bin) %>%
  summarise(
    sbp = if (all(is.na(sbp))) NA_real_ else mean(sbp, na.rm = TRUE),
    dbp = if (all(is.na(dbp))) NA_real_ else mean(dbp, na.rm = TRUE),
    .groups = "drop"
  )

# ── 3. Pivot to Wide (one row per subject) ────────────────────────────────────
chap_wide <- chap_binned %>%
  pivot_wider(
    id_cols     = Subject,
    names_from  = bin,
    values_from = c(sbp, dbp),
    names_glue  = "{.value}_bin{bin}"
  )

# Guarantee all 20 columns exist (fill absent bins with NA)
feature_cols <- c(paste0("sbp_bin", 1:10), paste0("dbp_bin", 1:10))
for (col in feature_cols) {
  if (!col %in% names(chap_wide)) chap_wide[[col]] <- NA_real_
}
chap_wide <- chap_wide[, c("Subject", feature_cols)]


# Attach Group label
subj_group <- unique(chap_long[, c("Subject", "Group")])
chap_flat  <- merge(chap_wide, subj_group, by = "Subject", all.x = TRUE)
chap_flat  <- chap_flat[, c("Subject", "Group", feature_cols)]

# ── 3b. Missing Value Summary per Subject ─────────────────────────────────────
chap_flat$n_missing <- rowSums(is.na(chap_flat[, feature_cols]))
hist(chap_flat$n_missing)
cat("Missing bin count distribution by Group:\n")
miss_summary <- chap_flat %>%
  group_by(Group, n_missing) %>%
  summarise(n_subjects = n(), .groups = "drop") %>%
  as.data.frame()
print(miss_summary)

# Drop subjects with too many missing bins (fill in threshold below)
miss_threshold <- 10   # ← adjust: drop if n_missing >= this value
chap_flat <- chap_flat[chap_flat$n_missing < miss_threshold, ]
chap_flat$n_missing <- NULL

cat("After dropping (>= ", miss_threshold, " missing bins):",
    nrow(chap_flat), "subjects\n")
cat("Group distribution:\n"); print(table(chap_flat$Group))
cat("Missing counts per feature:\n")
print(colSums(is.na(chap_flat[, feature_cols])))

# ── Handoff: save pre-imputation flat data for 2_mice_ecdf.R ──────────────────
dir.create("../../../data/raw/chap", showWarnings = FALSE)
saveRDS(chap_flat, "../../../data/raw/chap/chap_flat_preimpute.RDS")
cat("Saved pre-imputation data → ../../../data/raw/chap/chap_flat_preimpute.RDS\n")