### Manual process. 
### Complete sif contains all possible interactions. 
### The cutoff value is found by running aracne again 
### with a p value = 1e-8
library(data.table)
setwd("/media/ddisk/transpipeline-data/utero/networks")

cutoffval <- 0.244248 
all_inter <- fread("utero-healthy.sif", sep ="\t", header = F)
colnames(all_inter) <- c("source", "MI", "target")
inter_e8 <- all_inter[all_inter$MI > cutoffval, ]
inter_e8 <- inter_e8[order(-inter_e8$MI),]
fwrite(inter_e8, file = "1e8/utero-healthy-1e8.sif", sep ="\t")
