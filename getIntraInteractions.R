library(data.table)
library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = T)

if (length(args) < 4 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  ANNOT_RDATA <- args[1]
  MATRIX <- args[2]
  OUTFILE <- args[3]
  MCCORES <- as.integer(args[4])
}

load(ANNOT_RDATA)

chrs <- c(as.character(1:22), "X", "Y")
MImatrix <- fread(file = MATRIX, header = T,  sep = ",", nThread = MCCORES)
MImatrix <- data.matrix(MImatrix)

MImatrix <- rbind(MImatrix,  rep(NA, ncol(MImatrix)))

rownames(MImatrix) = colnames(MImatrix)

annot <- annot %>% filter(gene_id %in% colnames(MImatrix))
rownames(annot) <- annot$gene_id
cond <- strsplit(strsplit(MATRIX, split = "/")[[1]][8], split="_")[[1]][1]

## Intra-chromosomal interaction pairs
## calculates distance between each pair of genes in the form:
## gene1[start] - gene2[start]
condvals <- lapply(X = chrs, FUN = function(ch) {
  cat("Working with chromosome", ch, "\n")
  genes_annot <- annot %>% filter(chr == ch)
  genes <- genes_annot$gene_id
  ngenes <- length(genes)
  
  if(ngenes > 0) {
    chrvals <- parallel::mclapply(X = 1:(ngenes-1), mc.cores = MCCORES,  mc.cleanup = T, FUN = function(i) {
      gene1 <- genes[i]
      other.genes <- genes[(i+1):ngenes]
      mivals <- lapply(X = other.genes, 
                       FUN = function(gene2) {
                         list(source = gene1, target = gene2, 
                              distance = max(genes_annot[gene2, "start"], genes_annot[gene1, "start"]) - 
                                min(genes_annot[gene1, "start"], genes_annot[gene2, "start"]), 
                              mi = max(MImatrix[gene1, gene2], MImatrix[gene2, gene1],  na.rm = T))
                       })
      return(bind_rows(mivals))
    })
    chrdf <- bind_rows(chrvals)
    chrdf$chr <- ch
    return(chrdf)
  }
})

conddf <- bind_rows(condvals)
conddf$cond <- cond
cat("Saving file.\n")
fwrite(conddf, file = OUTFILE, row.names = F, col.names = T, sep = "\t", nThread = MCCORES)