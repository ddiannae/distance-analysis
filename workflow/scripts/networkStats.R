################################################################################
## Script to calculate network, vertices and interaction attributes, such as 
## centrality, density, transitivity, number of components, etc using igraph
## It requires interactions and vertices files from networkTables.R as input
###############################################################################
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(igraph)
library(readr)
library(tidyr)
library(dplyr)

cat("Reading files\n")
interactions <- read_tsv(snakemake@input[["interactions"]]) %>% 
  dplyr::rename("from" = "source_ensembl", "to" = "target_ensembl",
                "weight" = "mi")

vertices <- read_tsv(snakemake@input[["vertices"]]) %>% 
  dplyr::rename("name" = "ensembl")

cat("Building weighted graph\n")
net <- igraph::graph_from_data_frame(interactions, directed=FALSE, 
                                     vertices = vertices)

cat("Getting graph statistics\n")
cd <- igraph::centr_degree(net)
cb <- igraph::centr_betw(net)
ce <- igraph::centr_eigen(net)
inter <- (interactions%>% filter(interaction_type == "Inter") %>% nrow)/nrow(interactions)

centr <- tibble(name = names(V(net)), centr_degree = cd$res, 
                centr_between = cb$res, centr_eigen = ce$vector)

stats <- tibble(statistic = c("density", "transitivity", "no_components", "mean_distance", 
                "diameter", "degree_centralization", "between_centralization",
                "eigen_centralization", "vertices", "interactions", "inter_fraction"),
       value = c(igraph::edge_density(net), igraph::transitivity(net), 
                igraph::components(net)$no, igraph::mean_distance(net), 
                igraph::diameter(net), cd$centralization, 
                cb$centralization, ce$centralization, nrow(vertices), nrow(interactions), 
                inter))

stats <- stats %>% arrange(statistic)

cat("Getting node and edges statistics\n")
node_attrs <- tibble::enframe(igraph::degree(net)) %>% 
  rename("degree" = "value") %>% 
  dplyr::inner_join(tibble::enframe(igraph::strength(net)) %>% 
            dplyr::rename("strength" = "value")) %>%
  dplyr::inner_join(tibble::enframe(igraph::coreness(net, mode="all")) %>%
            dplyr::rename("coreness" = "value")) %>%
  dplyr::inner_join(tibble::enframe(igraph::hub_score(net)$vector) %>%
            dplyr::rename("hub_score" = "value")) %>%
  dplyr::inner_join(centr)

edges <- igraph::as_data_frame(net, what="edges") %>% dplyr::select(from, to)
edges$edge_between <- edge_betweenness(net, directed=T)

cat("Saving files\n")
save(net, file=snakemake@output[["network"]])
write_tsv(stats, file=snakemake@output[["network_stats"]])
write_tsv(node_attrs, file=snakemake@output[["node_attributes"]])
write_tsv(edges, file=snakemake@output[["edge_attributes"]])