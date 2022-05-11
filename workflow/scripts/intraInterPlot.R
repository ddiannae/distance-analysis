## #############################################################
## This script plots the intra-chromosomal fractions for cancer
## and normal tissue at different MI values thresholds.
## Input comes from the intraInterCount.R script
################################################################

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(ggplot2)
library(ggthemes)
library(dplyr)

cat("Reading files\n")
files <- list(snakemake@input[["normal"]], snakemake@input[["cancer"]])
BINTYPE <- snakemake@params[["bintype"]]

TISSUE <- snakemake@params[["tissue"]]
substring(TISSUE, 1, 1) <- toupper(substring(TISSUE, 1, 1))

mi_data <- lapply(files, function(file) {
  read_tsv(file)
})
mi_data <- bind_rows(mi_data) 
mi_data$cond <- factor(mi_data$cond, levels = c("normal", "cancer"), labels = c("Normal", "Cancer"))

if(BINTYPE == "onek") {
  mi_data$bin <- mi_data$bin * 1000
} 

colors <- c("#e3a098", "#a32e27")

cat("Building plot\n")
p <- ggplot(mi_data, aes(y = intra_fraction, x = bin, color=cond)) +
  geom_point(size = 6) +
  scale_x_log10(name="Total Interactions", breaks = c(1e3,1e4,1e5,1e6,1e7)) +
  theme_base(base_size = 30) +
  theme(legend.position = c(0.85, 0.8), axis.title.y = element_text(size=28), 
        plot.background=element_blank()) +
  scale_y_continuous(name = "Fraction of intra-chromosomal", breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     limits = c(0, 1)) +
  scale_color_manual(name="Condition", values = colors) +
  ggtitle(TISSUE)

if(BINTYPE == "log") {
  p <- p + geom_line(size = 3)
}

cat("Saving plot\n")
png(snakemake@output[[1]], width = 1200, height = 600)
print(p)
dev.off()  
