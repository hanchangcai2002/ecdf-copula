library(ggplot2)
library(tidyr)
library(dplyr)

# ── Load results ───────────────────────────────────────────────────────────────
gbmt_comparison <- readRDS(
  "Mar27_gbmt_comparison.RDS"
)

# ── Pivot to long (Mean ARI only) ─────────────────────────────────────────────
df_long <- gbmt_comparison %>%
  select(Scenario, Mean_Train_ARI, Mean_Test_ARI) %>%
  pivot_longer(
    cols      = c("Mean_Train_ARI", "Mean_Test_ARI"),
    names_to  = "Split",
    values_to = "ARI"
  )

# ── Corresponding SD for annotation ───────────────────────────────────────────
df_sd <- gbmt_comparison %>%
  select(Scenario, SD_Train_ARI, SD_Test_ARI) %>%
  pivot_longer(
    cols      = c("SD_Train_ARI", "SD_Test_ARI"),
    names_to  = "Split",
    values_to = "SD"
  ) %>%
  mutate(Split = recode(Split,
    "SD_Train_ARI" = "Mean_Train_ARI",
    "SD_Test_ARI"  = "Mean_Test_ARI"
  ))

df_plot <- left_join(df_long, df_sd, by = c("Scenario", "Split"))

# ── Plot ───────────────────────────────────────────────────────────────────────
ggplot(df_plot, aes(x = Split, y = Scenario, fill = ARI)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = paste0(round(ARI, 3), "\n(±", round(SD, 3), ")")),
    color = "white", size = 3.8, lineheight = 1.2
  ) +
  scale_fill_gradient(low = "#d4e8f5", high = "#08306b") +
  scale_x_discrete(
    labels = c("Mean_Train_ARI" = "Train ARI", "Mean_Test_ARI" = "Test ARI")
  ) +
  labs(
    title = "GBMT ARI Comparison: Real vs ECDF vs CTGAN Synthetic\n(Mean ± SD across 5 imputations)",
    x = NULL, y = NULL, fill = "Mean ARI"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text        = element_text(color = "black"),
    panel.grid       = element_blank(),
    plot.title       = element_text(hjust = 0.5)
  )

ggsave(
  filename = "plot_0327_gbmt_ari.png",
  width = 7, height = 4.5, dpi = 300
)

cat("Plot saved to plot_0327_gbmt_ari.png\n")
