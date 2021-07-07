library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 4 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  NORMAL_NETWORK <- args[1]
  CANCER_NETWORK <- args[2]
  FIGDIR <- args[3]
  SUFFIX <- args[4]
}

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
  theme_minimal(base_size = 20) 

png(paste0(FIGDIR, "/mi_density_network_", SUFFIX, ".png"), width = 750, height = 500)
print(p)
dev.off()  
  
p <- ggplot(DT, aes(x = cond, y = mi,  fill = cond)) + 
  geom_boxplot() +
  theme_minimal(base_size = 20) +
  scale_fill_manual(name = "Condition", values = colors) +
  scale_color_manual(name = "Condition", values = colors) +
  xlab("") +
  theme(legend.position = "none") 
    
png(paste0(FIGDIR, "/mi_boxplot_network_", SUFFIX, ".png"), width = 750, height = 750)
print(p)
dev.off()  
