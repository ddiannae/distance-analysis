## #############################################################
## This file builds a normal and cancer heatmap of the p-values 
## from the wilcoxon tests that compare mi distributions of bins  
## in a single condition. 
## Its input comes from the binTest.R script
################################################################
library(vroom)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

tissue <- snakemake@params[["tissue"]]
tissue <- paste0(toupper(substring(tissue, 1, 1)), substring(tissue, 2))

getHeatmap <- function(ttest_file, title) {
  htcolors <- circlize::colorRamp2(breaks = c(0, 1),
                                   colors = c("violetred4", "lightyellow1"))
  cat("Reading file ", ttest_file, "\n")
  ttest_df <- vroom::vroom(ttest_file)
  cat("Building matrix\n")
  ttest_matrix <- ttest_df %>% 
    pivot_wider(id_cols = bin1, names_from = bin2, values_from = p.value) %>%
    as.matrix()
  cat("Setting col and rownames\n")
  rownames(ttest_matrix) <- paste0("bin", ttest_matrix[, 1])
  ttest_matrix <- ttest_matrix[, -1]
  colnames(ttest_matrix) <- c(colnames(ttest_matrix)[1], paste0("bin", colnames(ttest_matrix)[-1]))
  cat("Copying upper triangular\n")
  ttest_matrix[lower.tri(ttest_matrix, diag = F)] = t(ttest_matrix)[lower.tri(ttest_matrix, diag = F)]
  cat("Building heatmap\n")
  return(Heatmap(ttest_matrix, cluster_rows = FALSE, cluster_columns =  FALSE, show_row_names = F, 
                 show_column_names = F, column_title = title, col = htcolors, name = "p-val", 
                 heatmap_legend_param = list(legend_height = unit(6, "cm"), 
                                             title_gp = gpar(fontsize = 18), 
                                             labels_gp = gpar(fontsize = 18), 
                                             direction = "vertical"),
                 column_title_gp = grid::gpar(fontsize = 22)))
  
}

h1 <- getHeatmap(snakemake@input[["normal"]], "Normal")
h2 <- getHeatmap(snakemake@input[["cancer"]], "Cancer")
ht_opt("legend_gap" = unit(2, "cm"))
png(snakemake@output[[1]], width = 1000, height = 400)
draw(h1+h2, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legends = TRUE,
     column_title =  , column_title_gp = gpar(fontsize = 24))
dev.off()