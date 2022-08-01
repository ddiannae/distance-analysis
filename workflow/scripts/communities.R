################################################################################
## Script to get network communities using an specific algorithm specified in
## the params. It gets a membership file and a file with community features.
## It requires interactions and vertices files from networkTables.R as input
###############################################################################

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(igraph)
library(readr)
library(dplyr)

type <- snakemake@params[["comm_type"]]
algorithm <- snakemake@params[["commalg"]]

getComInfo <- function(cmembership, network){
  comp.info <- lapply(unique(cmembership), function(idc){
    mem <- names(cmembership[cmembership == idc])
    com <- induced.subgraph(network, mem)
    intras <- igraph::as_data_frame(com, what = "edges") %>% 
      filter(interaction_type == "Intra") %>% nrow()
    prs <- page.rank(com)
    chrs <- table(V(com)$chr)
    edges <- length(E(com))
    return(data.frame(com_id = idc,
                      pg_gene = names(which.max(prs$vector))[1], 
                      chr = names(which.max(chrs))[1], order = length(V(com)), 
                      size = edges,
                      intra_fraction = intras/edges))
  })
  comp.info <- plyr::compact(comp.info)
  comp.info <- plyr::ldply(comp.info)
  comp.info <- comp.info %>% arrange(desc(size))
  return(comp.info)
}

cat("Reading files\n")
interactions <- read_tsv(snakemake@input[["interactions"]])  %>% 
  dplyr::rename("from" = "source_ensembl", "to" = "target_ensembl",
                "weight" = "mi")

vertices <- read_tsv(snakemake@input[["vertices"]]) %>% 
  dplyr::rename("name" = "ensembl")

if(type == "intra") {
  interactions <- interactions %>% filter(interaction_type == "Intra")
  
  ## Keep only vertices in intra-chromosomal interactions
  vertices <- vertices %>% filter(name %in% union(interactions$from, 
                                                  interactions$to))
}

net <- graph_from_data_frame(interactions, 
                             directed=F, vertices = vertices)

cat("Getting communities\n")
# Weights used by default

switch(algorithm,
       infomap = {comm <- cluster_infomap(graph = net)},
       fgreedy = {comm <- cluster_fast_greedy(graph = net)},
       louvain = {comm <- comm <- cluster_louvain(graph = net)},
       leadeigen = {comm <- cluster_leading_eigen(graph = net)},
       stop("Communities identified")
)

cat(algorithm, " modularity: ", modularity(comm), "\n")
cat(algorithm, " communities: ", length(comm))

names(comm$membership) <- comm$names

df_comm <- data.frame(comm$names, comm$membership)
colnames(df_comm) <- c("ensembl", "community")

comm_info <- getComInfo(comm$membership, net)

cat("Saving files\n")
write_tsv(df_comm, snakemake@output[["comm"]])
write_tsv(comm_info, snakemake@output[["comm_info"]])

