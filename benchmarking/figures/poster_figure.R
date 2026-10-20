library(ggplot2)
library(ggsignif)
library(readr)

results <- read_delim("../out/results.csv")

results <- results[!(results$duration > 2e-04 & results$length == 100 & results$repeat_type == "str_3"),]
results <- results[!(results$duration > 0.01 & results$length == 10000 & results$repeat_type == "str_3"),]
results <- results[!(results$duration > 0.01 & results$length == 100 & results$repeat_type == "vntr_50"),]

results$length_label <- paste("length:", results$length)

# rows = vars(as.character(length)), cols = vars(repeat_type),

plot <- ggplot(data = results, mapping = aes(x = caller, y = duration)) +
  facet_wrap(vars(repeat_type, length_label), scales = "free") +
  geom_boxplot(outlier.alpha = 0.2) +
  theme_minimal() +
  theme(text = element_text(size = 24)) +
  labs(x = "Decomposition tool", y = "Duration (s)") +
  ylim(0, NA)

ggsave("../out/poster_figure.pdf", width = 26, height = 12, dpi = 300)
