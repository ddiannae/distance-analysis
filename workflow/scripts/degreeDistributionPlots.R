################################################################################
## Script get a degree distribution and cumulative degree distribution plot of  
## cancer an normal networks.
## It requires the RData files from networStats.R as input
###############################################################################
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(igraph)
library(ggplot2)
library(ggthemes)

cat("Loading networks\n")
load(snakemake@input[["network_cancer"]])
net_cancer <- net
load(snakemake@input[["network_normal"]])

TISSUE <- snakemake@params[["tissue"]]
substring(TISSUE, 1, 1) <- toupper(substring(TISSUE, 1, 1))

cat("Getting cumulative degree distribution\n")
ddis_cancer <- igraph::degree_distribution(net_cancer, cumulative = T)
ddis_normal <- igraph::degree_distribution(net, cumulative = T)
# We don't have any 0 degree nodes
ddis_cancer <- ddis_cancer[-1]
ddis_normal <- ddis_normal[-1]

maxd_normal <- max(igraph::degree(net))
maxd_cancer <- max(igraph::degree(net_cancer))
maxd <- max(maxd_normal, maxd_cancer)

ddis_cancer <- c(ddis_cancer, rep(0, abs(maxd - maxd_cancer)))
ddis_normal <- c(ddis_normal, rep(0, abs(maxd - maxd_normal)))

cat("Building cumulative degree distribution plot\n")
p <- ggplot() +
  geom_line(aes(x=1:maxd, y=1-ddis_normal), color="#e3a098") +
  geom_line(aes(x=1:maxd, y=1-ddis_cancer), color="#a32e27") +
  theme_minimal(base_size = 30) +
  xlab("Node degree") +
  ylab("Fraction of nodes") +
  ggtitle(TISSUE)

cat("Saving plot\n")
png(filename=snakemake@output[["cdd"]], width = 1000, height = 500)
print(p)
dev.off()  

cat("Getting degree distribution\n")
ddis_cancer <- igraph::degree_distribution(net_cancer, cumulative = F)
ddis_normal <- igraph::degree_distribution(net, cumulative = F)
# We don't have any 0 degree nodes
ddis_cancer <- ddis_cancer[-1]
ddis_normal <- ddis_normal[-1]

ddis_cancer <- c(ddis_cancer, rep(0, abs(maxd - maxd_cancer)))
ddis_normal <- c(ddis_normal, rep(0, abs(maxd - maxd_normal)))

cat("Building degree distribution plot\n")
p <- ggplot() +
  geom_line(aes(x=1:maxd, y=ddis_normal), color="#e3a098") +
  geom_line(aes(x=1:maxd, y=ddis_cancer), color="#a32e27") +
  theme_minimal(base_size = 30) +
  xlab("Node degree") +
  ylab("Fraction of nodes") +
  ggtitle(TISSUE)

cat("Saving plot\n")
png(filename=snakemake@output[["dd"]], width = 1000, height = 500)
print(p)
dev.off()  
