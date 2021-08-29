log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)

NORMAL_NETWORK <- snakemake@input[["normal"]]
CANCER_NETWORK <- snakemake@input[["cancer"]]
TISSUE <- snakemake@params[["tissue"]]
substring(TISSUE, 1, 1) <- toupper(substring(TISSUE, 1, 1))

colors <- c("#e3a098", "#a32e27")
labels <- c( "Healthy", "Cancer")
names(colors) <- labels
  
healthy <- read_tsv(NORMAL_NETWORK)
cancer <- read_tsv(CANCER_NETWORK)

DT <- bind_rows(healthy, cancer)
DT$cond <- factor(DT$cond,   levels = c("normal", "cancer"), labels = labels)

p <- ggplot(DT) + 
  geom_density(aes(x = mi, y=..scaled.., fill = cond, color = cond), alpha=1/2)  + 
  xlab("Mutual Information") +
  ylab("Density") +
  scale_fill_manual(name = "Condition", values = colors) +
  scale_color_manual(name = "Condition", values = colors) +
  guides(color = FALSE) +
  theme_minimal(base_size = 30)  +
  ggtitle(TISSUE)

png(filename=snakemake@output[["density"]], width = 1000, height = 500)
print(p)
dev.off()  
  
p <- ggplot(DT, aes(x = cond, y = mi,  fill = cond)) + 
  geom_boxplot() +
  theme_minimal(base_size = 30) +
  scale_fill_manual(name = "Condition", values = colors) +
  scale_color_manual(name = "Condition", values = colors) +
  xlab("") +
  ylab("Mutual Information") +
  theme(legend.position = "none")  +
  ggtitle(TISSUE)
    
png(snakemake@output[["boxplot"]], width = 750, height = 750)
print(p)
dev.off()  
