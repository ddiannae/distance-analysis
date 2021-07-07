library(readr)
library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 4 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  CI_NORMAL <- args[1]
  CI_CANCER <- args[2]
  FIGDIR <- args[3]
  CUTOFF <- args[4]
}

comms_normal <- read_tsv(CI_NORMAL)
comms_normal$cond <- "normal"
comms_cancer <- read_tsv(CI_CANCER)
comms_cancer$cond <- "cancer"
comms <- bind_rows(comms_normal, comms_cancer)
comms <- comms %>% filter(order >= 5)

colors <- c("#e3a098", "#a32e27")
labels <- c( "Healthy", "Cancer")
comms$cond <- factor(comms$cond,   levels = c("normal", "cancer"), labels = labels)

p <- ggplot(comms) +
  geom_boxplot(aes(x = cond, y = size, fill = cond)) +
  theme_minimal(base_size = 30) +
  scale_fill_manual(name = "Condition", values = colors) +
  xlab("") +
  ylab("Community size") +
  theme(legend.position = "none") 

png(paste0(FIGDIR, "/comm_size_boxplot_network_", CUTOFF, ".png"), width = 750, height = 750)
print(p)
dev.off()  

p <- ggplot(comms) +
  geom_histogram(aes(x = size,  fill = cond, color = cond),
                binwidth = 1, position = "identity") +
  theme_minimal(base_size = 30) +
  facet_wrap(~cond, nrow = 1) +
  scale_fill_manual(name = "Condition", values = colors) +
  scale_color_manual(name = "Condition", values = colors) +
  ylab("Frequency") +
  xlab("Community size") +
  theme(legend.position = "none")

png(paste0(FIGDIR, "/comm_size_histogram_network_", CUTOFF, ".png"), width = 1000, height = 500)
print(p)
dev.off()  

p <- ggplot(comms) +
  geom_boxplot(aes(x = cond, y = order, fill = cond)) +
  theme_minimal(base_size = 30) +
  scale_fill_manual(name = "Condition", values = colors) +
  xlab("") +
  ylab("Community order") +
  theme(legend.position = "none") 

png(paste0(FIGDIR, "/comm_order_boxplot_network_", CUTOFF, ".png"), width = 750, height = 750)
print(p)
dev.off()  

p <- ggplot(comms) +
  geom_histogram(aes(x = order,  fill = cond, color = cond),
                 binwidth = 1, position = "identity") +
  theme_minimal(base_size = 30) +
  facet_wrap(~cond, nrow = 1) +
  scale_fill_manual(name = "Condition", values = colors) +
  scale_color_manual(name = "Condition", values = colors) +
  ylab("Frequency") +
  xlab("Community order") +
  theme(legend.position = "none")

png(paste0(FIGDIR, "/comm_order_histogram_network_", CUTOFF, ".png"), width = 1000, height = 500)
print(p)
dev.off()  

p <- ggplot(comms) +
  geom_boxplot(aes(x = cond, y = (2*size)/(order * (order-1)), fill = cond)) +
  theme_minimal(base_size = 30) +
  scale_fill_manual(name = "Condition", values = colors) +
  xlab("") +
  ylab("Community density") +
  theme(legend.position = "none") 

png(paste0(FIGDIR, "/comm_density_boxplot_network_", CUTOFF, ".png"), width = 750, height = 750)
print(p)
dev.off()  

p <- ggplot(comms) +
  geom_histogram(aes(x = (2*size)/(order * (order-1)),  fill = cond, color = cond),
                 bins = 100, position = "identity") +
  theme_minimal(base_size = 30) +
  facet_wrap(~cond, nrow = 1) +
  scale_fill_manual(name = "Condition", values = colors) +
  scale_color_manual(name = "Condition", values = colors) +
  ylab("Frequency") +
  xlab("Community density") +
  theme(legend.position = "none")

png(paste0(FIGDIR, "/comm_density_histogram_network_", CUTOFF, ".png"), width = 1000, height = 500)
print(p)
dev.off()  