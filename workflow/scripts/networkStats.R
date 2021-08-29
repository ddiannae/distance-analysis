log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(igraph)
library(tidyr)
library(readr)

cat("Reading files\n")
interactions <- read_tsv(snakemake@input[["interactions"]]) %>% 
  dplyr::rename("from" = "source_ensembl", "to" = "target_ensembl",
                "weight" = "mi")

vertices <- read_tsv(snakemake@input[["vertices"]]) %>% 
  dplyr::rename("name" = "ensembl")

cat("Building weighted graph")
net <- igraph::graph_from_data_frame(interactions, directed=FALSE, 
                                     vertices = vertices)

cat("Getting graph statistics\m")
cd <- igraph::centr_degree(net)
cb <- igraph::centr_betw(net)
ce <- igraph::centr_eigen(net)

centr <- tibble(name = names(V(net)), centr_degree = cd$res, 
                centr_between = cb$res, centr_eigen = ce$vector)

stats <- tibble(statistic = c("density", "transitivity", "no_components", "mean_distance", 
                "diameter", "degree_centralization", "between_centralization",
                "eigen_centralization"),
       value = c(igraph::edge_density(net), igraph::transitivity(net), 
                igraph::components(net)$no, igraph::mean_distance(net), 
                igraph::diameter(net), cd$centralization, 
                cb$centralization, ce$centralization))

stats <- stats %>% arrange(statistic)

cat("Getting node and edges statistics\m")
node_attrs <- tibble::enframe(igraph::degree(net)) %>% 
  rename("degree" = "value") %>% 
  inner_join(tibble::enframe(igraph::strength(net)) %>% 
               rename("strength" = "value")) %>%
  inner_join(tibble::enframe(igraph::coreness(net, mode="all")) %>%
               rename("coreness" = "value")) %>%
  inner_join(tibble::enframe(igraph::hub_score(net)$vector) %>%
               rename("hub_score" = "value")) %>%
  inner_join(centr)

edges <- as_data_frame(net, what="edges") %>% select(from, to)
edges$edge_between <- edge_betweenness(net, directed=T)

cat("Saving files\n")
save(net, file=snakemake@output[["network"]])
write_tsv(stats, file=snakemake@output[["network_stats"]])
write_tsv(node_attrs, file=snakemake@output[["node_attributes"]])
write_tsv(edges, file=snakemake@output[["edge_attributes"]])

