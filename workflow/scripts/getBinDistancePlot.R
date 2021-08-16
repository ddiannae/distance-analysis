log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)

cat("Loading data\n")
color_pal <- c("#e3a098", "#a32e27")
fitted_data <- read_tsv(snakemake@input[[1]])
TISSUE <- snakemake@params[["tissue"]]
substring(TISSUE, 1, 1) <- toupper(substring(TISSUE, 1, 1))

cat("Building plot\n")
g <- ggplot(fitted_data) + 
  geom_line(aes(x = dist_mean/1e6, y = mean_fitted, color=cond)) +
  geom_ribbon(aes(x = dist_mean/1e6, ymin = mean_fitted - mi_sd, ymax = mean_fitted + mi_sd, fill = cond), 
              alpha = .2) +
  facet_wrap(~cond, nrow = 1) + 
  xlab("Distance (Mbp)") + 
  ylab("Mutual Information") + 
  theme_few(base_size = 25) +
  scale_fill_manual(values = color_pal) + 
  scale_color_manual(values = color_pal) + 
  theme(legend.position = "none", strip.text.x = element_text(size = 30)) +
  ggtitle(TISSUE)

cat("Saving plot\n")
png(snakemake@output[[1]], width = 1200, height = 600)
print(g)
dev.off()  