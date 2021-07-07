library(readr)
library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 4 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  NORMAL_INTER <- args[1]
  NORMAL_VER <- args[2]
  NORMAL_COMMS <- args[3]
  CANCER_INTER <- args[4]
  CANCER_VER <- args[5]
  CANCER_COMMS <- args[6]
  FIGDIR <- args[7]
  CUTOFF <- args[8]
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
  as.character(lnew)
}

getDistances <- function(interactions, vertices, communities, cond) {
  
  interactions <- interactions %>% filter(interaction_type == "Intra")
  vertices <- vertices %>% filter(ensembl %in% union(interactions$source_ensembl, 
                                                     interactions$target_ensembl))
  vertices <- vertices %>% inner_join(communities, by ="ensembl") 
  
  comm_dist <- lapply(unique(communities$community), function(idc){
    v_comm <- vertices %>% filter(community == idc)
    if(nrow(v_comm) > 2) {
      e_comm <- distinct(bind_rows(interactions %>% semi_join(v_comm, by = c("source_ensembl" = "ensembl")),
                                   interactions %>% semi_join(v_comm, by = c("target_ensembl" = "ensembl"))))
      mean_dist <- e_comm %>% summarise(mean(distance))  %>% unlist() %>% unname()
      diameter <- max(v_comm$start) - min(v_comm$start)
      return(list(diameter = diameter, mean_dist = mean_dist))
    }
    return(NULL)
  })
  comm_dist <- bind_rows(comm_dist)
  comm_dist$cond <- cond
  return(comm_dist) 
}

normal_inter <- read_tsv(NORMAL_INTER)
normal_ver <- read_tsv(NORMAL_VER)
normal_comms <- read_tsv(NORMAL_COMMS)

cancer_inter <- read_tsv(CANCER_INTER)
cancer_ver <- read_tsv(CANCER_VER)
cancer_comms <- read_tsv(CANCER_COMMS)

normal_distances <- getDistances(normal_inter, normal_ver, normal_comms, "normal")
cancer_distances <- getDistances(cancer_inter, cancer_ver, cancer_comms, "cancer")

comms <- bind_rows(normal_distances, cancer_distances)
colors <- c("#e3a098", "#a32e27")
labels <- c( "Healthy", "Cancer")
comms$cond <- factor(comms$cond,   levels = c("normal", "cancer"), labels = labels)

p <- ggplot(comms) +
  geom_boxplot(aes(x = cond, y = diameter, fill = cond)) +
  theme_minimal(base_size = 30) +
  scale_fill_manual(name = "Condition", values = colors) +
  xlab("") +
  ylab("Community diameter") +
  theme(legend.position = "none") 

png(paste0(FIGDIR, "/comm_diameter_boxplot_network_", CUTOFF, ".png"), width = 750, height = 750)
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

png(paste0(FIGDIR, "/comm_diameter_histogram_network_", CUTOFF, ".png"), width = 1000, height = 500)
print(p)
dev.off()  

p <- ggplot(comms) +
  geom_boxplot(aes(x = cond, y = mean_dist, fill = cond)) +
  theme_minimal(base_size = 30) +
  scale_fill_manual(name = "Condition", values = colors) +
  xlab("") +
  ylab("Community mean distance") +
  theme(legend.position = "none") 

png(paste0(FIGDIR, "/comm_meandistance_boxplot_network_", CUTOFF, ".png"), width = 750, height = 750)
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

png(paste0(FIGDIR, "/comm_meandistance_histogram_network_", CUTOFF, ".png"), width = 1000, height = 500)
print(p)
dev.off()  


