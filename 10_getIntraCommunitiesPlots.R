## Script to get histogram of communities sizes, 
## orders and densities
library(data.table)
library(dplyr)
library(ggplot2)

types <- c("utero", "kidney", "colon", "tiroides", "lung")

for(type in types) {
  
  setwd(paste0("/media/ddisk/transpipeline-data/", type))
  
  conds <- c("healthy", "cancer")
  
  comms <- lapply(conds, function(cond) {
    comm_info <- fread(paste0("networks/network-tables/", type, "-", cond, 
                             "-communities-info.tsv"), sep = "\t")
    comm_info <- comm_info %>% filter(order > 5)
    comm_info$cond <- cond
    return(comm_info)
  })
  comms <- bind_rows(comms)
  
  colors <- c("#e3a098", "#a32e27")
  labels <- c( "Healthy", "Cancer")
  comms$cond <- factor(comms$cond,   levels = conds, labels = labels)
  
  p <- ggplot(comms) +
    geom_boxplot(aes(x = cond, y = size, fill = cond)) +
    theme_minimal(base_size = 30) +
    scale_fill_manual(name = "Condition", values = colors) +
    xlab("") +
    ylab("Community size") +
    theme(legend.position = "none") 
  
  png("figures/comm-size-boxplot-min-network.png", width = 750, height = 750)
  print(p)
  dev.off()  

  p <- ggplot(comms) +
    geom_histogram(aes(x = size,  fill = cond, color = cond),
                  binwidth = 1, position = "identity") +
    theme_minimal(base_size = 30) +
    facet_wrap(~cond, nrow = 1) +
    scale_fill_manual(name = "Condition", values = colors) +
    scale_color_manual(name = "Condition", values = colors) +
    ylab("Frequency") +
    xlab("Community size") +
    theme(legend.position = "none")
  
  png("figures/comm-size-histogram-min-network.png", width = 1000, height = 500)
  print(p)
  dev.off()  
  
  p <- ggplot(comms) +
    geom_boxplot(aes(x = cond, y = order, fill = cond)) +
    theme_minimal(base_size = 30) +
    scale_fill_manual(name = "Condition", values = colors) +
    xlab("") +
    ylab("Community order") +
    theme(legend.position = "none") 
  
  png("figures/comm-order-boxplot-min-network.png", width = 750, height = 750)
  print(p)
  dev.off()  
  
  p <- ggplot(comms) +
    geom_histogram(aes(x = order,  fill = cond, color = cond),
                   binwidth = 1, position = "identity") +
    theme_minimal(base_size = 30) +
    facet_wrap(~cond, nrow = 1) +
    scale_fill_manual(name = "Condition", values = colors) +
    scale_color_manual(name = "Condition", values = colors) +
    ylab("Frequency") +
    xlab("Community order") +
    theme(legend.position = "none")
  
  png("figures/comm-order-histogram-min-network.png", width = 1000, height = 500)
  print(p)
  dev.off()  
  
  p <- ggplot(comms) +
    geom_boxplot(aes(x = cond, y = (2*size)/(order * (order-1)), fill = cond)) +
    theme_minimal(base_size = 30) +
    scale_fill_manual(name = "Condition", values = colors) +
    xlab("") +
    ylab("Community density") +
    theme(legend.position = "none") 
  
  png("figures/comm-density-boxplot-min-network.png", width = 750, height = 750)
  print(p)
  dev.off()  
  
  p <- ggplot(comms) +
    geom_histogram(aes(x = (2*size)/(order * (order-1)),  fill = cond, color = cond),
                   bins = 100, position = "identity") +
    theme_minimal(base_size = 30) +
    facet_wrap(~cond, nrow = 1) +
    scale_fill_manual(name = "Condition", values = colors) +
    scale_color_manual(name = "Condition", values = colors) +
    ylab("Frequency") +
    xlab("Community density") +
    theme(legend.position = "none")
  
  png("figures/comm-density-histogram-min-network.png", width = 1000, height = 500)
  print(p)
  dev.off()  
}