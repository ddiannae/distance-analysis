## Script to get MI density plot and MI boxplot for networks.
## It currently works with the min networks in *.sif files
## and it requires the files to have headers: source-MI-target
library(data.table)
library(ggplot2)
library(ggthemes)

setwd("/media/ddisk/transpipeline-data/")

types <- c("utero", "kidney", "tiroides", "lung", "colon")

getMIDistribution <- function(types) {
  conds <- c("healthy", "cancer")
  colors <- c("#e3a098", "#a32e27")
  labels <- c( "Healthy", "Cancer")
  names(colors) <- labels
  
  for(type in types) {
    mitables <- lapply(X = conds, FUN = function(cond) {
      cat("Reading network ", cond, "\n")
      MIvals <- fread(paste0(type, "/networks/1e8/", type, "-", cond, "-min.sif"), 
                      header = T, sep="\t", nThread = 3)
      MIvals$cond <- cond
      MIvals
    })
    
    DT <- rbindlist(mitables)
    DT$cond <- factor(DT$cond,   levels = conds, labels = labels)

    p <- ggplot(DT) + 
      geom_density(aes(x = MI, y=..scaled.., fill = cond, color = cond), alpha=1/2)  + 
      xlab("Mutual Information") +
      ylab("Density") +
      scale_fill_manual(name = "Condition", values = colors) +
      scale_color_manual(name = "Condition", values = colors) +
      guides(color = FALSE) +
      theme_minimal(base_size = 20) 
    
    png(paste0(type, "/figures/mi-density-min-network.png"), width = 750, height = 500)
    print(p)
    dev.off()  
  }
}

getMIDistribution(types)

getMIBox <- function(types) {
  conds <- c("healthy", "cancer")
  colors <- c("#e3a098", "#a32e27")
  labels <- c( "Healthy", "Cancer")
  names(colors) <- labels
  for(type in types) {
    mitables <- lapply(X = conds, FUN = function(cond) {
      cat("Reading network ", cond, "\n")
      MIvals <- fread(paste0(type, "/networks/1e8/", type, "-", cond, "-min.sif"), 
                      header = T, sep="\t", nThread = 3)
      MIvals$cond <- cond
      MIvals
    })
    
    DT <- rbindlist(mitables)
    DT$cond <- factor(DT$cond,   levels = conds, labels = labels)

    p <- ggplot(DT, aes(x = cond, y = MI,  fill = cond)) + 
      geom_boxplot() +
      theme_minimal(base_size = 20) +
      scale_fill_manual(name = "Condition", values = colors) +
      scale_color_manual(name = "Condition", values = colors) +
      xlab("") +
      theme(legend.position = "none") +
    png(paste0(type, "/figures/mi-boxplot-min-network.png"), width = 750, height = 750)
    print(p)
    dev.off()  
  }
}

getMIBox(types)
