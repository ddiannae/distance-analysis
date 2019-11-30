library(data.table)
setwd("/media/ddisk/transpipeline-data/kidney-data/networks")

cutoffval <-  0.105668
all_inter <- fread("kidney-healthy.sif", sep ="\t", header = F)
colnames(all_inter) <- c("source", "MI", "target")
inter_e8 <- all_inter[all_inter$MI > cutoffval, ]
fwrite(inter_e8, file = "1e8/kidney-healthy-1e8.sif", sep ="\t")
