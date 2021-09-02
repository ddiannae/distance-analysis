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
mi_data$cond <- factor(mi_data$cond, levels = c("normal", "cancer"), labels = c("Healthy", "Cancer"))

colors <- c("#e3a098", "#a32e27")

cat("Building plot\n")
p <- ggplot(mi_data, aes(y = intra_frac, x = bin, color=cond)) +
  geom_point(size = 2) +
  scale_x_log10(name="Number of Interactions", ) +
  theme_few(base_size = 20) +
  theme(legend.position = c(0.85, 0.8)) +
  ylab("Fraction of -cis interactions") +
  scale_color_manual(name="Condition", values = colors) +
  ggtitle(TISSUE)

if(BINTYPE == "log") {
  p <- p + geom_line()
}

cat("Saving plot\n")
png(snakemake@output[[1]], width = 800, height = 400)
print(p)
dev.off()  