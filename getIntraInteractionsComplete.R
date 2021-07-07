library(data.table)
library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = T)

if (length(args) < 2 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  ANNOT_RDATA <- args[1]
  OUTDIR 
  MATRIX <- args[2]
  MCCORES <- args[3]
}

load(ANNOT_RDATA)

chrs <- c(as.character(1:22), "X", "Y")
MImatrix <- fread(file = MATRIX, header = T,  sep = "\t", nThread = MCCORES, na.strings = "NA")
MImatrix <- data.matrix(MImatrix)

rownames(MImatrix) = colnames(MImatrix)[1: nrow(MImatrix)]

annot <- annot %>% filter(gene_id %in% colnames(MImatrix))
rownames(annot) <- annot$gene_id


## Intra-chromosomal interaction pairs
## calculates distance between each pair of genes in the form:
## gene1[start] - gene2[start]
condvals <- lapply(X = chrs, FUN = function(c) {
  cat("Working with chromosome", c, "\n")
  genes_annot <- annot[annot$chr == c, ]
  genes <- genes_annot$gene_id
  ngenes <- length(genes)
  chrvals <- parallel::mclapply(X = 1:(ngenes-1), mc.cores = MCCORES,  mc.cleanup = T, FUN = function(i) {
    gene1 <- genes[i]
    other.genes <- genes[(i+1):ngenes]
    mivals <- lapply(X = other.genes, 
                     FUN = function(gene2) {
                       data.frame(source = gene1, target = gene2, 
                                  distance = max(genes_annot[gene2, "start"], genes_annot[gene1, "start"]) - 
                                    min(genes_annot[gene1, "start"], genes_annot[gene2, "start"]), 
                                  mi = max(MImatrix[gene1, gene2], MImatrix[gene2, gene1],  na.rm = T))
                     })
    return(bind_rows(mivals))
  })
  chrdf <- bind_rows(chrvals)
  chrdf$chr <- c
  return(chrdf)
})

conddf <- bind_rows(condvals)
conddf$cond <- cond
cat("Saving file.\n")
fwrite(conddf, file =  paste0("../intra/", cond,"-all-distance-mi.txt"), 
       row.names = F, col.names = T, sep = "\t", nThread = 5)  

## Function that gets the mean MI value for a set of gene-pairs
## size of the set in binsize variable. 
## Bins are sorted according to the distance between genes
## It requires the files from getMIDistanceListByCondAndChrs
## File structure required:
## path
##  /type
##    /intra
##      cond-all-distance-mi.txt
##      /intra-fixed-bin-size (results will be here)
getAllMIDistanceMeans <- function(path, binsize, types) {
  chrs <- c(as.character(1:22), "X", "Y")
  
  for (type in types){
    cat("Working with type ", type, "\n")
    setwd(paste0(path, type, "/intra"))
    conds <- c("cancer", "healthy")
    conditiondfs <- lapply(X = conds, FUN = function(cond){
      cat("Working with condition ", cond, "\n")
      dist.mi <- fread(file=paste0(cond, "-all-distance-mi.txt"), header = T,
                          nThread = 5, sep = "\t")
      
      meansbych <- parallel::mclapply(X = chrs, mc.cores = 7,  mc.cleanup = FALSE, FUN = function(c){
        dist.mi.df <- dist.mi[dist.mi$chr == c, ]
        dist.mi.df <- dist.mi.df[order(dist.mi.df$distance), ]
        rownames(dist.mi.df) <- NULL
        dist.mi.df$bin <- ((as.numeric(rownames(dist.mi.df)) - 1)%/%binsize) + 1
        dfmeans <- aggregate(cbind(distance, mi)~bin, data=dist.mi.df, FUN=mean, na.rm=TRUE)
        dfmeans$cond <- cond
        dfmeans$chr <- c
        return(dfmeans)
      })
      df <- plyr::ldply(meansbych)
      write.table(df, file = paste0("intra-fixed-bin-size/", cond, "-", binsize, ".txt"), sep="\t",
                  col.names = T, row.names = F, quote = F)
      return(df)
    })
    
    cdf <- plyr::ldply(conditiondfs)
    write.table(cdf, file = paste0("intra-fixed-bin-size/", binsize, "-all.txt"), sep="\t",
                col.names = T, row.names = F, quote = F)
    }
  }
 
#getAllMIDistanceMeans("/media/ddisk/transpipeline-data/", 100, c("utero"))
