## Script to get histogram plots of size and distance 
## of intra-communities
library(data.table)
library(dplyr)
library(ggplot2)
types <- c("utero", "kidney", "colon", "tiroides", "lung")

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

dropLeadingZero <- function(l){
  #cat(l)
  lnew <- c()
  for(i in l){
    if(!is.na(i)) {
      if(i==0){ #zeros stay zero
        lnew <- c(lnew,"0")
      } else if (i>1){ #above one stays the same
        lnew <- c(lnew, as.character(i))
      } else
        lnew <- c(lnew, gsub("(?<![0-9])0+", "", i, perl = TRUE))
        lnew <- sub("\\+0?", "", lnew)
        lnew <- c(sub("e\\+", "e", lnew))
    }else{
      lnew <- c(lnew, "")
    }
  }
  #cat(lnew)
  as.character(lnew)
}

for(type in types) {
  
  conds <- c("healthy", "cancer")
  comms <- lapply(conds, function(cond) {
    gen_comm <- fread(file = paste0("../regulaciontrans-data/", type, "/", type, "-", cond, "-louvain-communities.tsv"), 
                     header = T, sep = "\t")
    interactions <- fread(file = paste0("../regulaciontrans-data/", type, "/", type, "-", cond, "-interactions.tsv"), 
                          header = T, sep = "\t")
    interactions <- interactions %>% filter(interaction_type == "Intra-Cytoband" | 
                                              interaction_type == "Inter-Cytoband")
    
    vertices <- fread(file = paste0("../regulaciontrans-data/", type, "/", type, "-", cond, "-vertices.tsv"), 
                      header = T, sep = "\t")
    vertices <- vertices %>% filter(ensemblID %in% union(interactions$source, interactions$target))
    vertices <- vertices %>% inner_join(gen_comm) 
    
    comm_dist <- lapply(unique(gen_comm$community), function(idc){
      v_comm <- vertices %>% filter(community == idc)
      if(nrow(v_comm) > 2) {
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
  
  png(paste0("../regulaciontrans-data/", type, "/","figures/comm-diameter-boxplot-network-2.png"), width = 750, height = 750)
  print(p)
  dev.off()  
  
  p <- ggplot(comms) +
    geom_histogram(aes(x = diameter,  fill = cond, color = cond),
                   bins = 100, position = "identity") +
    theme_minimal(base_size = 30) +
    facet_wrap(~cond, nrow = 1) +
    scale_fill_manual(name = "Condition", values = colors) +
    scale_color_manual(name = "Condition", values = colors) +
    scale_x_continuous(labels = dropLeadingZero) +
    ylab("Frequency") +
    xlab("Community diameter") +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 20))
  
  png(paste0("../regulaciontrans-data/", type, "/","figures/comm-diameter-histogram-network-2.png"), width = 1000, height = 500)
  print(p)
  dev.off()  
  
  p <- ggplot(comms) +
    geom_boxplot(aes(x = cond, y = mean_dist, fill = cond)) +
    theme_minimal(base_size = 30) +
    scale_fill_manual(name = "Condition", values = colors) +
    xlab("") +
    ylab("Community mean distance") +
    theme(legend.position = "none") 
  
  png(paste0("../regulaciontrans-data/", type, "/","figures/comm-meandist-boxplot-network-2.png"), width = 750, height = 750)
  print(p)
  dev.off()  
  
  p <- ggplot(comms) +
    geom_histogram(aes(x = mean_dist,  fill = cond, color = cond),
                   bins = 100, position = "identity") +
    theme_minimal(base_size = 30) +
    facet_wrap(~cond, nrow = 1) +
    scale_fill_manual(name = "Condition", values = colors) +
    scale_color_manual(name = "Condition", values = colors) +
    scale_x_continuous(labels = dropLeadingZero) +
    ylab("Frequency") +
    xlab("Community mean distance") +
    theme(legend.position = "none", axis.text.x = element_text(size = 20))
  
  png(paste0("../regulaciontrans-data/", type, "/","figures/comm-meandist-histogram-network-2.png"), width = 1000, height = 500)
  print(p)
  dev.off()  
}

