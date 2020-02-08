## Script to get histogram plots of size and distance 
## of intra-communities
library(data.table)
library(dplyr)
library(ggplot2)
types <- c("utero", "kidney", "colon", "tiroides", "lung")

for(type in types) {
  
  setwd(paste0("/media/ddisk/transpipeline-data/", type))
  
  conds <- c("healthy", "cancer")
  comms <- lapply(conds, function(cond) {
    gen_comm <- fread(file = paste0("networks/network-tables/", type, "-", cond, "-communities.tsv"), 
                     header = T, sep = "\t")
    interactions <- fread(file = paste0("networks/network-tables/", type, "-", cond, "-interactions.tsv"), 
                      header = T, sep = "\t")
    interactions <- interactions %>% filter(interaction_type == "Intra-Cytoband" | 
                                              interaction_type == "Inter-Cytoband")
    
    vertices <- fread(file = paste0("networks/network-tables/", type, "-", cond, "-vertices.tsv"), 
                          header = T, sep = "\t")
    vertices <- vertices %>% filter(ensemblID %in% union(interactions$source, interactions$target))
    vertices <- vertices %>% inner_join(gen_comm) 
    
    comm_dist <- lapply(unique(gen_comm$community), function(idc){
      v_comm <- vertices %>% filter(community == idc)
      if(nrow(v_comm) > 5) {
        e_comm <- distinct(bind_rows(interactions %>% semi_join(v_comm, by = c("source" = "ensemblID")),
                                     interactions %>% semi_join(v_comm, by = c("target" = "ensemblID"))))
        mean_dist <- e_comm %>% summarise(mean(distance))  %>% unlist() %>% unname()
        diameter <- max(v_comm$start) - min(v_comm$start)
        return(list(diameter = diameter, mean_dist = mean_dist))
      }
      return(NULL)
    })
    comm_dist <- bind_rows(comm_dist)
    comm_dist$cond <- cond
    return(comm_dist)
  })
  
  comms <- bind_rows(comms)
  colors <- c("#e3a098", "#a32e27")
  labels <- c( "Healthy", "Cancer")
  comms$cond <- factor(comms$cond,   levels = conds, labels = labels)
  
  p <- ggplot(comms) +
    geom_boxplot(aes(x = cond, y = diameter, fill = cond)) +
    theme_minimal(base_size = 30) +
    scale_fill_manual(name = "Condition", values = colors) +
    xlab("") +
    ylab("Community diameter") +
    theme(legend.position = "none") 
  
  png("figures/comm-diameter-boxplot-min-network.png", width = 750, height = 750)
  print(p)
  dev.off()  
  
  p <- ggplot(comms) +
    geom_histogram(aes(x = diameter,  fill = cond, color = cond),
                   bins = 100, position = "identity") +
    theme_minimal(base_size = 30) +
    facet_wrap(~cond, nrow = 1) +
    scale_fill_manual(name = "Condition", values = colors) +
    scale_color_manual(name = "Condition", values = colors) +
    ylab("Frequency") +
    xlab("Community diameter") +
    theme(legend.position = "none")
  
  png("figures/comm-diameter-histogram-min-network.png", width = 1000, height = 500)
  print(p)
  dev.off()  
  
  p <- ggplot(comms) +
    geom_boxplot(aes(x = cond, y = mean_dist, fill = cond)) +
    theme_minimal(base_size = 30) +
    scale_fill_manual(name = "Condition", values = colors) +
    xlab("") +
    ylab("Community mean distance") +
    theme(legend.position = "none") 
  
  png("figures/comm-meandist-boxplot-min-network.png", width = 750, height = 750)
  print(p)
  dev.off()  
  
  p <- ggplot(comms) +
    geom_histogram(aes(x = mean_dist,  fill = cond, color = cond),
                   bins = 100, position = "identity") +
    theme_minimal(base_size = 30) +
    facet_wrap(~cond, nrow = 1) +
    scale_fill_manual(name = "Condition", values = colors) +
    scale_color_manual(name = "Condition", values = colors) +
    ylab("Frequency") +
    xlab("Community mean distance") +
    theme(legend.position = "none")
  
  png("figures/comm-meandist-histogram-min-network.png", width = 1000, height = 500)
  print(p)
  dev.off()  
}

