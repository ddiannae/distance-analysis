library(EnvStats)
library(data.table)
setwd("/mnt/ddisk/transpipeline-data/breast-data/subtypes/networks/aracne-cluster/")

cond <- "basal"
all.networks <- list.files(path = paste0("samples/", cond), pattern = NULL, all.files = FALSE,
                           full.names = T, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)

network1 <- fread(all.networks[3],  col.names = c("source", "MI", "target"))
network2 <- fread(all.networks[80],  col.names = c("source", "MI", "target"))
q1 <- quantile(network1$MI, seq(1:5000)/5000)
q2 <- quantile(network2$MI, seq(1:5000)/5000)
plot(q1, q2)
abline(a = 0, b = 1, col= "red")
qqPlot(network1$MI, network2$MI, line.col = "red", add.line = T)
