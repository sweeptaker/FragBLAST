library(ggplot2)

# Plot 1: human and environmental samples under the same 80/80 standard
counts_80 <- data.frame(
  dataset = factor(
    c("IHSMGC", "IHSMGC", "IGC", "IGC", "Tara Ocean", "Tara Ocean"),
    levels = c("IHSMGC", "IGC", "Tara Ocean")
  ),
  type = factor(
    c("Existing ARG", "New ARG", "Existing ARG", "New ARG", "Existing ARG", "New ARG"),
    levels = c("Existing ARG", "New ARG")
  ),
  count = c(5534, 857, 1176, 138, 118, 96)
)

labels_80 <- c(
  "IHSMGC\nHuman skin",
  "IGC\nHuman gut",
  "Tara Ocean\nEnvironment"
)

p1 <- ggplot(counts_80, aes(x = dataset, y = count, fill = type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_text(
    aes(label = count, group = type),
    position = position_dodge(width = 0.7),
    vjust = -0.5,
    size = 4.5,
    fontface = "bold"
  ) +
  scale_x_discrete(labels = labels_80) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12)),
    limits = c(0, 6500)
  ) +
  scale_fill_manual(
    values = c("Existing ARG" = "#4C78A8", "New ARG" = "#F58518"),
    name = NULL
  ) +
  labs(
    title = "ARG detection with identity >80% and qcov >80%",
    subtitle = "Human skin (IHSMGC), human gut (IGC) and environmental (Tara Ocean) samples",
    x = NULL,
    y = "Number of ARG sequences",
    caption = "Existing: DIAMOND hit above threshold; New: additionally confirmed by FragBLAST"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(face = "bold")
  )

# Plot 2: Tara Ocean under the two alternative criteria
counts_tara <- data.frame(
  standard = factor(
    c("Identity >70%\nqcov >40%", "Identity >70%\nqcov >40%",
      "Identity >80%\nqcov >80%", "Identity >80%\nqcov >80%"),
    levels = c("Identity >70%\nqcov >40%", "Identity >80%\nqcov >80%")
  ),
  type = factor(
    c("Existing ARG", "New ARG", "Existing ARG", "New ARG"),
    levels = c("Existing ARG", "New ARG")
  ),
  count = c(1157, 3375, 118, 96)
)

p2 <- ggplot(counts_tara, aes(x = standard, y = count, fill = type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_text(
    aes(label = count, group = type),
    position = position_dodge(width = 0.7),
    vjust = -0.5,
    size = 4.5,
    fontface = "bold"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12)),
    limits = c(0, 4000)
  ) +
  scale_fill_manual(
    values = c("Existing ARG" = "#4C78A8", "New ARG" = "#F58518"),
    name = NULL
  ) +
  labs(
    title = "Effect of detection criteria on Tara Ocean ARGs",
    x = NULL,
    y = "Number of ARG sequences",
    caption = "Human-like samples were screened with 80%/80%; environmental samples used 70%/40%"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold")
  )

dir.create("figures", showWarnings = FALSE)
ggsave(
  "figures/arg_counts_80_80.png",
  plot = p1,
  width = 9,
  height = 6.5,
  dpi = 300
)
ggsave(
  "figures/tara_standard_comparison.png",
  plot = p2,
  width = 8,
  height = 6,
  dpi = 300
)
