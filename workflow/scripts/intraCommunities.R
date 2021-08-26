log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(igraph)
library(readr)
library(dplyr)

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

cat("Reading files\n")
interactions <- read_tsv(snakemake@input[["interactions"]])
vertices <- read_tsv(snakemake@input[["vertices"]])

interactions <- interactions %>% filter(interaction_type == "Intra")

## Keep only vertices in intra-chromosomal interactions
vertices <- vertices %>% filter(ensembl %in% union(interactions$source_ensembl, 
                                                    interactions$target_ensembl))

colnames(interactions)[1:2] <- c("from", "to")
net <- graph_from_data_frame(interactions, 
                             directed=F, vertices = vertices)

cat("Getting intra-interactions communities\n")
comm <- cluster_louvain(graph = net)
names(comm$membership) <- comm$names

df_comm <- data.frame(comm$names, comm$membership)
colnames(df_comm) <- c("ensembl", "community")

comm_info <- getComInfo(comm$membership, net)

cat("Saving files\n")
write_tsv(df_comm, snakemake@output[["comm"]])
write_tsv(comm_info, snakemake@output[["comm_info"]])

