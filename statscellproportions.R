#statistical analysis to verify differences in cell proportions - DRG single nuclei (mouse)
library(reshape2)
library(ggplot2)
library(dplyr)

# -----------------------------
# Data: counts per sample
# -----------------------------
celltypes <- c("Satglia", "Fibroblast", "Endothelial", "Schwann_M", 
               "Schwann_N", "Pericyte", "Immune")

sham0220 <- c(10694, 3579, 2450, 2175, 907, 749, 574)
exer0220 <- c(7345, 3509, 2040, 1975, 452, 409, 644)
dmm0220  <- c(8887, 4205, 2439, 2005, 937, 609, 450)

# Combine into a matrix
counts <- matrix(c(sham0220, exer0220, dmm0220), 
                 nrow = length(celltypes), byrow = FALSE,
                 dimnames = list(celltypes, c("Sham0220", "Exer0220", "DMM0220")))

# -----------------------------
# 1. Chi-squared test
# -----------------------------
chisq_test <- chisq.test(counts)
print("Chi-squared test across all cell types:")
print(chisq_test)

# -----------------------------
# 2. Fisher's exact tests per cell type
# -----------------------------
get_star <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("")
}

comparison_stars <- data.frame(
  CellType = character(),
  Comparison = character(),
  p_value = numeric(),
  Significance = character(),
  stringsAsFactors = FALSE
)

for (cell in celltypes) {
  total_sham <- sum(counts[, "Sham0220"])
  total_exer <- sum(counts[, "Exer0220"])
  total_dmm  <- sum(counts[, "DMM0220"])
  
  mat_exer <- matrix(c(counts[cell, "Sham0220"], total_sham - counts[cell, "Sham0220"],
                       counts[cell, "Exer0220"], total_exer - counts[cell, "Exer0220"]),
                     nrow = 2, byrow = TRUE)
  mat_dmm <- matrix(c(counts[cell, "Sham0220"], total_sham - counts[cell, "Sham0220"],
                      counts[cell, "DMM0220"], total_dmm - counts[cell, "DMM0220"]),
                    nrow = 2, byrow = TRUE)
  
  p_exer <- fisher.test(mat_exer)$p.value
  p_dmm  <- fisher.test(mat_dmm)$p.value
  
  comparison_stars <- rbind(comparison_stars,
                            data.frame(CellType = cell, Comparison = "Exer0220", 
                                       p_value = p_exer, Significance = get_star(p_exer)),
                            data.frame(CellType = cell, Comparison = "DMM0220", 
                                       p_value = p_dmm, Significance = get_star(p_dmm)))
}

# -----------------------------
# 3. Plot with significance asterisks
# -----------------------------
library(reshape2)
library(ggplot2)
library(dplyr)

# Melt counts for plotting
prop_df <- as.data.frame(counts)
prop_df$CellType <- rownames(prop_df)
prop_df_melt <- melt(prop_df, id.vars = "CellType", variable.name = "Group", value.name = "Count")

# Calculate percentages
total_counts <- tapply(prop_df_melt$Count, prop_df_melt$Group, sum)
prop_df_melt$Percentage <- mapply(function(count, group) {
  count / total_counts[[group]] * 100
}, prop_df_melt$Count, prop_df_melt$Group)

# Merge significance stars
prop_df_melt <- left_join(prop_df_melt, comparison_stars, by = c("CellType", "Group" = "Comparison"))

# Plot
ggplot(prop_df_melt, aes(x = CellType, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = Significance), 
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 4, na.rm = TRUE) +
  ylab("Cell Type Percentage") +
  xlab("Cell Type") +
  ggtitle("Cell Type Proportions (Significance vs Sham0220)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")
