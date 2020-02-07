### Script to find communities (hot-spots) in the subnetworks spanned by 
### intra-chromosomal interactions. 
### Requires the type/intra/cond-all-distance-mi.txt file from the 
### 1_getIntraInteractions.R script
library(igraph)
library(data.table)

types <- c("utero", "kidney", "colon", "tiroides", "lung")

for(type in types) {
  
  setwd(paste0("/media/ddisk/transpipeline-data/", type))
  
  conds <- c("healthy", "cancer")
  
  m <- lapply(conds, function(cond) {
    mi_data <- fread(file = paste0("networks/network-tables/", type, "-", cond, "-interactions.tsv"), 
                    header = T, sep = "\t")
    mi_data <- mi_data[mi_data$interaction_type == "Intra-Cytoband" | 
                         mi_data$interaction_type == "Inter-Cytoband", ]
    colnames(mi_data)[1:2] <- c("from", "to")
    vertices <- union(mi_data$from, mi_data$to)
      net <- graph_from_data_frame(mi_data, 
                                 directed=F, vertices = vertices)
    comm <- cluster_infomap(graph = net, nb.trials = 10)
    df_comm <- data.frame(comm$names, comm$membership)
    colnames(df_comm) <- c("gene", "community")
    fwrite(df_comm, paste0("networks/network-tables/", type, "-", cond, 
                           "-communities.tsv"), sep = "\t")
    
  })
  
}
