library(readr)
library(dplyr)
library(tidyr)
matrix <- "/home/diana/Workspace/regulaciontrans-data-server-take1/kidney/networks/normal_network.adj"
annot <- "/home/diana/Workspace/regulaciontrans-data-server-take1/kidney/rdata/annot.RData"
load(annot)
MImatrix <- vroom::vroom(matrix)
MImatrix <- rbind(MImatrix,  rep(NA, ncol(MImatrix)))
MImatrix$genes <- colnames(MImatrix)
mill <- MImatrix %>% pivot_longer(cols = starts_with("ENSG"))
