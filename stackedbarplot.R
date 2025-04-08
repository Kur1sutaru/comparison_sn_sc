setwd()

library(dplyr)
library(ggplot2)

# Use full meta.data, not a single column
cell_counts <- dmmbatch_assigned@meta.data %>%
  group_by(orig.ident, scpred_prediction) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

# Rename columns for clarity
colnames(cell_counts) <- c("Sample", "CellType", "Count", "Proportion")
library(ggplot2)

ggplot(cell_counts, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    x = "Sample", 
    y = "Proportion", 
    title = "Cell Type Proportions by Sample (scPred)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))