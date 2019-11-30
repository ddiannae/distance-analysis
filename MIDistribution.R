library(data.table)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(tidyr)

setwd("/mnt/ddisk/transpipeline-data/breast-data/subtypes/")

conds <- c("healthy", "luma", "lumb", "basal", "her2")

#### Save annotated networks

mitables <- lapply(X = conds, FUN = function(network) {
  cat("Reading network ", network, "\n")
  MIvals <- fread(paste0("networks/aracne-cluster/1e8/13319/", network, ".sif"), 
                  header = F, sep="\t", nThread = 1, col.names =  c("Source", "Target", "MI"))
  
  colnames(annot) <-  c("Source", "Source.Chr", "Source.Start", "Source.End")
  MIvals <- merge(annot, MIvals, by.x = "Source", by.y = "Source")
  colnames(annot) <-  c("Target", "Target.Chr", "Target.Start", "Target.End")
  MIvals <- merge(annot, MIvals, by.x = "Target", by.y = "Target")
  MIvals$Inter <- MIvals$Target.Chr != MIvals$Source.Chr
  MIvals$Cond <- network
  write.table(MIvals, file = paste0("networks/aracne-cluster/1e8/13319/", network, "-13319-annotated.tsv"), 
              sep ="\t",  row.names = F, quote = F)
  return(MIvals)
})


#### Networks are in networks/aracne-cluster/1e8/13319/ o networks/aracne-cluster/1e8/13319/
getMIDistribution <- function(conds, path) {
  mitables <- lapply(X = conds, FUN = function(network) {
    cat("Reading network ", network, "\n")
    MIvals <- fread(paste0(path, "/", network, ".sif"), header = F, sep="\t", nThread = 1, col.names =  c("Source", "Target", "MI"))
    MIvals$Cond <- network
    MIvals
  })
  
  DT <- rbindlist(mitables)
  
  p <- ggplot(DT, aes(x = MI, color = Cond)) + geom_density()
  png("mi-distribution.png", width = 1500, height = 1500)
  print(p)
  dev.off()
}

load("rdata/annot.RData")
annot <- annot[, c("EnsemblID", "Chr", "Start", "End")]

getMIBox <- function(conds, path) {
  mitables <- lapply(X = conds, FUN = function(network) {
    cat("Reading network ", network, "\n")
    MIvals <- fread(paste0("networks/aracne-cluster/1e8/13319/", network, ".sif"), 
                    header = F, sep="\t", nThread = 1, col.names =  c("Source", "Target", "MI"))
    
    colnames(annot) <-  c("Source", "Source.Chr", "Source.Start", "Source.End")
    MIvals <- merge(annot, MIvals, by.x = "Source", by.y = "Source")
    colnames(annot) <-  c("Target", "Target.Chr", "Target.Start", "Target.End")
    MIvals <- merge(annot, MIvals, by.x = "Target", by.y = "Target")
    MIvals$Inter <- MIvals$Target.Chr != MIvals$Source.Chr
    MIvals <- MIvals[, c("Source", "Target", "MI", "Inter")]
    MIvals <- MIvals[!MIvals$Inter,  ]
    MIvals$Cond <- network
    maxMI <- max(MIvals$MI)
    MIvals$NormMI <- MIvals$MI/maxMI
    MIvals
  })
  
  DT <- rbindlist(mitables)
  p <- ggplot(DT, aes(x = Cond, y = NormMI,  fill = Cond)) + 
        geom_boxplot() +
        theme_minimal() +
        scale_fill_tableau()
  
  png("mi-distribution.png", width = 1500, height = 1500)
  print(p)
  dev.off()
}

