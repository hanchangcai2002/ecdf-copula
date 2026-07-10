# Algorithm 2 (longitudinal eCDF-copula) applied to REMBRANDT training set.
# Input:  data/raw/rembrandt/0_Depression_Data_Real_Train_Long.csv (long format, 14 visits, only the training split is used per the paper's design)
# Output: data/processed/rembrandt/0_Depression_Data_Synthetic_ECDF.csv
rm(list = ls())
library(reshape2)
library(copula)
library(matrixStats)
library(tidyverse)
source("../../synthesis_algorithms/algorithm2_longitudinal.R")
impt1 <- read.csv("../../../data/raw/rembrandt/0_Depression_Data_Real_Train_Long.csv")
table(impt1$group)
head(impt1)
# Example head(impt1):
#   subject   group timestamp         ma       hars
# 1       1 Control         1 -0.8723685 -0.9433202
# 2       1 Control         2 -0.8723685 -0.9433202
# 3       1 Control         3 -0.8723685 -1.1302705
feature_cols <- c("ma", "hars")
prepare_data_for_sim <- function(data_long) {
  data_long %>%
    mutate(
      record_id    = as.numeric(factor(subject)),
      months_to_base = timestamp   # timestamp 1~14 as timepoint
    )
}

convert_sim_result_to_wide <- function(sim_result, group_name) {
  sim_data   <- sim_result$simulated_data   # n_sim × (14*2) matrix
  t_vals     <- sim_result$time_points      # 1:14
  features   <- sim_result$features        # c("ma","hars")
  
  n_sim      <- nrow(sim_data)
  n_t        <- length(t_vals)
  n_feat     <- length(features)
  
  # sort t from small to large (Sim function already sorted by t_vals)
  t_order    <- order(t_vals)
  
  result <- data.frame(group = group_name)
  
  for (ti in t_order) {
    col_start <- (ti - 1) * n_feat + 1
    col_end   <-  ti      * n_feat
    block     <- sim_data[, col_start:col_end, drop = FALSE]
    # block column names are originally ma_t1, hars_t1, etc. rename as ma1, hars1
    colnames(block) <- paste0(features, t_vals[ti])
    result <- cbind(result, block)
  }
  
  # make sure order as：group, ma1,hars1, ma2,hars2,...
  ord <- c("group",
           as.vector(rbind(paste0("ma",   1:14),
                           paste0("hars", 1:14))))
  result[, ord]
}

impt1_prepared <- prepare_data_for_sim(impt1)
control_data <- impt1_prepared %>% filter(group == "Control")
remit_data   <- impt1_prepared %>% filter(group == "Remitted")


set.seed(123)
sim1_con_result <- Sim.longitudinal.feature(
  dat.cor = control_data,
  feature_cols = feature_cols,
  n_sim = 31
)
sim_con_wide <- convert_sim_result_to_wide(sim1_con_result, "Control")
sim1_rem_result <- Sim.longitudinal.feature(
  dat.cor = remit_data,
  feature_cols = feature_cols,
  n_sim = 65
)
sim_rem_wide <- convert_sim_result_to_wide(sim1_rem_result, "Remitted")

sim_combined_wide <- rbind(sim_con_wide, sim_rem_wide)
syn_out_dir <- "../../../data/processed/rembrandt"
dir.create(syn_out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(sim_combined_wide,
          file = file.path(syn_out_dir, "0_Depression_Data_Synthetic_ECDF.csv"),
          row.names = FALSE)