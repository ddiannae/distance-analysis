library(igraph)
library(readr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 4 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  INTERACTIONS <- args[1]
  VERTICES <- args[2]
  COMM <- args[3]
  COMM_INFO <- args[4]
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

interactions <- read_tsv(INTERACTIONS)
interactions <- interactions %>% filter(interaction_type == "Intra")

vertices <- read_tsv(VERTICES)

## Keep only vertices in intra-chromosomal interactions
vertices <- vertices %>% filter(ensembl %in% union(interactions$source_ensembl, 
                                                    interactions$target_ensembl))

colnames(interactions)[1:2] <- c("from", "to")
net <- graph_from_data_frame(interactions, 
                             directed=F, vertices = vertices)
comm <- cluster_louvain(graph = net)
names(comm$membership) <- comm$names

df_comm <- data.frame(comm$names, comm$membership)
colnames(df_comm) <- c("ensembl", "community")

write_tsv(df_comm, COMM)

comm_info <- getComInfo(comm$membership, net)
write_tsv(comm_info, COMM_INFO)

