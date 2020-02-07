### Script to find communities (hot-spots) in the subnetworks spanned by 
### intra-chromosomal interactions. 
### Requires the type/intra/cond-all-distance-mi.txt file from the 
### 1_getIntraInteractions.R script
library(igraph)
library(data.table)
library(dplyr)

types <- c("utero", "kidney", "colon", "tiroides", "lung")

for(type in types) {
  
  setwd(paste0("/media/ddisk/transpipeline-data/", type))
  
  conds <- c("healthy", "cancer")
  
  m <- lapply(conds, function(cond) {
    interactions <- fread(file = paste0("networks/network-tables/", type, "-", cond, "-interactions.tsv"), 
                          header = T, sep = "\t")
    interactions <- interactions %>% filter(interaction_type == "Intra-Cytoband" | 
                                              interaction_type == "Inter-Cytoband")
    
    vertices <- fread(file = paste0("networks/network-tables/", type, "-", cond, "-vertices.tsv"), 
                      header = T, sep = "\t")
    vertices <- vertices %>% filter(ensemblID %in% union(interactions$source, interactions$target))
    
    colnames(interactions)[1:2] <- c("from", "to")
    net <- graph_from_data_frame(interactions, 
                                 directed=F, vertices = vertices)
    comm <- cluster_infomap(graph = net, nb.trials = 10)
    names(comm$membership) <- comm$names
    
    df_comm <- data.frame(comm$names, comm$membership)
    colnames(df_comm) <- c("ensemblID", "community")
    fwrite(df_comm, paste0("networks/network-tables/", type, "-", cond, 
                           "-communities.tsv"), sep = "\t")
    
    comm_info <- getComInfo(comm$membership, net)
    fwrite(comm_info, paste0("networks/network-tables/", type, "-", cond, 
                           "-communities-info.tsv"), sep = "\t")
  })
}

getComInfo <- function(cmembership, network){
  comp.info <- lapply(unique(cmembership), function(idc){
    mem <- names(cmembership[cmembership == idc])
    com <- induced.subgraph(network, mem)
    prs <- page.rank(com)
    chrs <- table(V(com)$chr)
    return(data.frame(com_id = idc,
                      pg_gene = names(which.max(prs$vector))[1], 
                      chr = names(which.max(chrs))[1], order = length(V(com)), 
                      size = length(E(com))))
  })
  comp.info <- plyr::compact(comp.info)
  comp.info <- plyr::ldply(comp.info)
  comp.info <- comp.info %>% arrange(desc(size))
  return(comp.info)
}
