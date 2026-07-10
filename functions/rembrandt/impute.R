# MICE (PMM) imputation of the raw REMBRANDT longitudinal data followed by a long -> wide reshape.
# Input:  data/raw/rembrandt/rembrandt.csv          (raw, pre-imputation)
# Output: data/raw/rembrandt/real_long.csv          (long: Subject/Time/TimeLabel/Group/ma/hars)
#         data/raw/rembrandt/0_Depression_Data_Real.csv  (wide: subject/group/ma1/hars1/.../ma14/hars14)
rm(list = ls())
library(tidyverse)
library(mice)
library(purrr)

RAW_DIR <- "../../data/raw/rembrandt"

# ── Stage 1: MICE imputation ─────────────────────────────────────────────────
visit_levels <- c("screening", "baseline",
                  sprintf("month_%02d", seq(2, 24, by = 2)))
data <- read.csv(file.path(RAW_DIR, "rembrandt.csv"))
data$visit <- factor(data$visit, levels = visit_levels)

full_dat <- expand.grid(
  subj_num = unique(data$subj_num),
  visit = levels(data$visit)
) %>%
  left_join(data, by = c("subj_num", "visit")) %>%
  arrange(subj_num, visit)
# keep demographic info
demo_dat <- full_dat %>%
  group_by(subj_num) %>%
  slice(1) %>%
  select(subj_num, condition, age, sex, race, ethnicity)

impute_one_var <- function(varname, full_dat, demo_dat) {
  visit_levels <- levels(full_dat$visit)

  # Reshape to wide: one column per visit for this variable
  wide <- full_dat %>%
    select(subj_num, visit, !!sym(varname)) %>%
    pivot_wider(
      names_from = visit,
      values_from = !!sym(varname),
      names_prefix = paste0(varname, "_")
    )

  # Attach demographic covariates
  wide <- left_join(wide, demo_dat, by = "subj_num")

  # Multiple imputation via mice (single completed dataset, m = 1)
  imp <- mice(wide, method = "pmm", m = 1, maxit = 10, seed = 123)
  completed <- complete(imp, 1)

  # Extract the completed variable columns (named varname_visit) back to long
  long <- completed %>%
    select(subj_num, starts_with(paste0(varname, "_"))) %>%
    pivot_longer(
      cols = starts_with(paste0(varname, "_")),
      names_to = "visit",
      names_prefix = paste0(varname, "_"),
      values_to = varname
    )

  return(long)
}
ma_vars <- paste0("ma_", 1:10)
hars_vars <- paste0("hars_", 1:14, "_sev")
all_vars <- c(ma_vars, hars_vars)

# Impute each item variable separately, then merge results
all_imputed <- map(all_vars, ~ impute_one_var(.x, full_dat, demo_dat)) %>%
  reduce(full_join, by = c("subj_num", "visit"))

# Compute total scores (ma_tot, hars_score) from the imputed items
impt1 <- all_imputed %>%
  mutate(
    ma_tot = rowSums(across(all_of(ma_vars)), na.rm = TRUE),
    hars_score = rowSums(across(all_of(hars_vars)), na.rm = TRUE)
  ) %>%
  select(subj_num, visit, ma_tot, hars_score, all_of(ma_vars), all_of(hars_vars)) %>%
  arrange(subj_num, visit)

# Re-attach the condition (treatment group) column
condition_df <- data %>%
  select(subj_num, condition) %>%
  distinct()
impt1 <- impt1 %>%
  left_join(condition_df, by = "subj_num")

stopifnot(!any(impt1$ma_tot != rowSums(select(impt1, ma_1:ma_10))))
stopifnot(!any(impt1$hars_score != rowSums(select(impt1, hars_1_sev:hars_14_sev))))

# Reshape to real_long's schema: standardize the total scores (ma, hars),
# rename to Subject/Time/Group, and map visit labels to a numeric time index
# (screening = -1, baseline = 0, month_02 = 2, ...) while keeping the
# original visit label in TimeLabel for provenance.
visit_time_map <- setNames(c(-1, 0, seq(2, 24, by = 2)), visit_levels)

real_long <- impt1 %>%
  select(subj_num, visit, condition, ma_tot, hars_score) %>%
  mutate(
    ma        = as.vector(scale(ma_tot)),
    hars      = as.vector(scale(hars_score)),
    TimeLabel = as.character(visit),
    Time      = visit_time_map[as.character(visit)]
  ) %>%
  select(-ma_tot, -hars_score, -visit) %>%
  rename(Subject = subj_num, Group = condition) %>%
  select(Subject, Time, TimeLabel, Group, ma, hars)

dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(real_long, file.path(RAW_DIR, "real_long.csv"), row.names = FALSE)

# ── Stage 2: long -> wide reshape (14 visits per subject) ───────────────────
to_wide_14 <- function(df) {
  df <- df %>%
    mutate(
      Subject = as.character(Subject),
      Group   = as.character(Group)
    )

  # Enforce within-subject time order and create 1..14 index
  df <- df %>%
    group_by(Subject) %>%
    arrange(Time, .by_group = TRUE) %>%
    mutate(t_idx = row_number()) %>%
    ungroup()

  # Wide: group, ma1,hars1,...,ma14,hars14 (no Subject column per downstream format)
  wide <- df %>%
    select(Subject, Group, t_idx, ma, hars) %>%
    pivot_wider(
      id_cols = c(Subject, Group),
      names_from = t_idx,
      values_from = c(ma, hars),
      names_glue = "{.value}{t_idx}"
    ) %>%
    rename(subject = Subject, group = Group)

  # Ensure column order: subject, group, ma1,hars1,ma2,hars2,...
  ord <- c(
    "subject", "group",
    as.vector(rbind(paste0("ma", 1:14), paste0("hars", 1:14)))
  )
  wide %>% select(any_of(ord), everything())
}

real_wide <- to_wide_14(real_long)
write.csv(real_wide, file.path(RAW_DIR, "0_Depression_Data_Real.csv"), row.names = FALSE)
