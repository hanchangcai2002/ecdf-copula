# Plots Figure 6 (Membership Inference Attack: precision/accuracy vs. the
# proportion of the dataset known by an attacker, by Hamming-distance
# threshold) from the results table produced by
# 1_Membership_Inference_Dataset0.ipynb.
library(tidyverse)

results <- read.csv("INFERENCE TESTS RESULTS/membership_inference_results.csv")

synth_levels <- c("GM", "ECDF", "CTGAN", "SDV", "WGANGP")
# matplotlib tab:blue/orange/green/red equivalents, keyed to threshold value
# (matches the original notebook's color order for thresholds 0.4/0.3/0.2/0.1)
threshold_colors <- c("0.4" = "#1f77b4", "0.3" = "#ff7f0e",
                      "0.2" = "#2ca02c", "0.1" = "#d62728")

plot_dat <- results %>%
  mutate(
    synthesizer = factor(synthesizer, levels = synth_levels),
    threshold   = factor(threshold, levels = c(0.4, 0.3, 0.2, 0.1))
  ) %>%
  pivot_longer(
    cols      = c(precision, accuracy),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = c("accuracy", "precision"),
                          labels = c("acc", "prec")))

p <- ggplot(plot_dat, aes(x = proportion, y = value, color = threshold)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", alpha = 0.7) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_grid(metric ~ synthesizer) +
  scale_color_manual(values = threshold_colors, name = "Threshold") +
  scale_x_continuous(breaks = sort(unique(plot_dat$proportion))) +
  scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(x = "Proportion of dataset known by an attacker", y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
    strip.background = element_rect(fill = "gray90"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("INFERENCE TESTS RESULTS/MEMBERSHIP_INFERENCE_TESTS_RESULTS.svg",
       plot = p, width = 13, height = 5.5)
ggsave("INFERENCE TESTS RESULTS/MEMBERSHIP_INFERENCE_TESTS_RESULTS.png",
       plot = p, width = 13, height = 5.5, dpi = 300)

cat("Saved plot to INFERENCE TESTS RESULTS/MEMBERSHIP_INFERENCE_TESTS_RESULTS.svg and .png\n")
