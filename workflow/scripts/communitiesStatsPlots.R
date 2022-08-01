## #############################################################
## This file builds plots for size, order and density 
## distribution plots for normal and cancer communities. 
## It requires the community summary file from communities.R as input
################################################################
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(dplyr)
library(ggplot2)

TISSUE <- snakemake@params[["tissue"]]
substring(TISSUE, 1, 1) <- toupper(substring(TISSUE, 1, 1))

cat("Reading files\n")
comms_normal <- read_tsv(snakemake@input[["comm_info_normal"]], 
                         col_types = cols(chr = col_character()))
comms_normal$cond <- "normal"
comms_cancer <- read_tsv(snakemake@input[["comm_info_cancer"]], 
                         col_types = cols(chr = col_character()))
comms_cancer$cond <- "cancer"
comms <- bind_rows(comms_normal, comms_cancer)

cat("Filtering communities with more than 5 nodes\n")
comms <- comms %>% filter(order >= 5)

colors <- c("#e3a098", "#a32e27")
labels <- c( "Normal", "Cancer")
comms$cond <- factor(comms$cond,   levels = c("normal", "cancer"), labels = labels)


cat("Building order boxplot network\n")
p <- ggplot(comms) +
  geom_boxplot(aes(x = cond, y = order, fill = cond)) +
  theme_minimal(base_size = 30) +
  scale_fill_manual(name = "Condition", values = colors) +
  xlab("") +
  ylab("Community order") +
  theme(legend.position = "none")  +
  ggtitle(TISSUE)

png(snakemake@output[[1]], width = 750, height = 750)
print(p)
dev.off()  

cat("Building order histogram network\n")
p <- ggplot(comms) +
  geom_histogram(aes(x = order,  fill = cond, color = cond),
                 binwidth = 1, position = "identity") +
  theme_minimal(base_size = 30) +
  facet_wrap(~cond, nrow = 1) +
  scale_fill_manual(name = "Condition", values = colors) +
  scale_color_manual(name = "Condition", values = colors) +
  ylab("Frequency") +
  xlab("Community order") +
  theme(legend.position = "none") +
  ggtitle(TISSUE)

png(snakemake@output[[2]], width = 1000, height = 500)
print(p)
dev.off()  

cat("Building size boxplot network\n")
p <- ggplot(comms) +
  geom_boxplot(aes(x = cond, y = size, fill = cond)) +
  theme_minimal(base_size = 30) +
  scale_fill_manual(name = "Condition", values = colors) +
  xlab("") +
  ylab("Community size") +
  theme(legend.position = "none")  +
  ggtitle(TISSUE)

png(snakemake@output[[3]], width = 750, height = 750)
print(p)
dev.off()  

cat("Building size histogram network\n")
p <- ggplot(comms) +
  geom_histogram(aes(x = size,  fill = cond, color = cond),
                binwidth = 1, position = "identity") +
  theme_minimal(base_size = 30) +
  facet_wrap(~cond, nrow = 1) +
  scale_fill_manual(name = "Condition", values = colors) +
  scale_color_manual(name = "Condition", values = colors) +
  ylab("Frequency") +
  xlab("Community size") +
  theme(legend.position = "none") +
  ggtitle(TISSUE)

png(snakemake@output[[4]], width = 1000, height = 500)
print(p)
dev.off()  

cat("Building density boxplot network\n")
p <- ggplot(comms) +
  geom_boxplot(aes(x = cond, y = (2*size)/(order * (order-1)), fill = cond)) +
  theme_minimal(base_size = 30) +
  scale_fill_manual(name = "Condition", values = colors) +
  xlab("") +
  ylab("Community density") +
  theme(legend.position = "none")  +
  ggtitle(TISSUE)

png(snakemake@output[[5]], width = 750, height = 750)
print(p)
dev.off()  

cat("Building density histogram network\n")
p <- ggplot(comms) +
  geom_histogram(aes(x = (2*size)/(order * (order-1)),  fill = cond, color = cond),
                 bins = 100, position = "identity") +
  theme_minimal(base_size = 30) +
  facet_wrap(~cond, nrow = 1) +
  scale_fill_manual(name = "Condition", values = colors) +
  scale_color_manual(name = "Condition", values = colors) +
  ylab("Frequency") +
  xlab("Community density") +
  theme(legend.position = "none") +
  ggtitle(TISSUE)

png(snakemake@output[[6]], width = 1000, height = 500)
print(p)
dev.off()  