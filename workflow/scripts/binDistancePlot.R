## #############################################################
## This file builds a line plot for the bin mean or loess fitted
## mean MI vs distance with variance ribbons.
## Its input comes from the binFittingByChr.R and binStats.R script
################################################################

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)

cat("Loading data\n")
color_pal <- c("#e3a098", "#a32e27")
labels <- c( "Normal", "Cancer")

fitted_data <- read_tsv(snakemake@input[[1]])
TISSUE <- snakemake@params[["tissue"]]
STAT <- snakemake@params[["stat"]]
substring(TISSUE, 1, 1) <- toupper(substring(TISSUE, 1, 1))

fitted_data <- fitted_data %>% mutate(
  min_plot =  pmax(get(STAT) - mi_sd, 0), 
  max_plot = pmin(get(STAT) + mi_sd, 0.1), 
  cond = factor(cond, levels = labels, labels = labels))

cat("Building plot\n")
g <- ggplot(fitted_data) + 
  geom_line(aes(x = dist_mean/1e6, y = get(STAT), color=cond), size = 0.7) +
  geom_ribbon(aes(x = dist_mean/1e6, ymin = min_plot, 
                  ymax = max_plot, fill = cond), 
              alpha = .2) +
  facet_wrap(~cond, nrow = 1) + 
  xlab("Distance (Mbp)") + 
  ylab("Mutual Information") + 
  theme_few(base_size = 30) +
  ylim(c(0,0.1)) +
  scale_fill_manual(values = color_pal) + 
  expand_limits(x = 0, y = 0) +
  scale_color_manual(values = color_pal) + 
  theme(legend.position = "none", strip.text.x = element_text(size = 30), 
        plot.title = element_text(size = 40, face = "bold")) +
  ggtitle(TISSUE)

cat("Saving plot\n")
png(snakemake@output[[1]], width =1200, height = 600)
print(g)
dev.off()  