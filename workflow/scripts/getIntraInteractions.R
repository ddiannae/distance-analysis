log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(data.table)
library(readr)
library(dplyr)

cat("Reading files \n")
ANNOT_RDATA <- snakemake@params[["annot"]]
MCCORES <- as.integer(snakemake@threads[[1]])
COND <- snakemake@params[["cond"]]

load(ANNOT_RDATA)

chrs <- c(as.character(1:22), "X", "Y")
MImatrix <- fread(file = snakemake@input[["network"]], header = T,  sep = ",", nThread = MCCORES)
MImatrix <- data.matrix(MImatrix)
MImatrix <- rbind(MImatrix,  rep(NA, ncol(MImatrix)))
rownames(MImatrix) = colnames(MImatrix)
cat("MI matrix with ", nrow(MImatrix), " rows and ", ncol(MImatrix), " columns loaded \n")
annot <- annot %>% filter(gene_id %in% colnames(MImatrix))

## Intra-chromosomal interaction pairs
## calculates distance between each pair of genes in the form:
## gene1[start] - gene2[start]
distvals <- lapply(X = chrs, FUN = function(ch) {
  cat("Working with chromosome", ch, "\n")
  genes_annot <- annot %>% filter(chr == ch)
  genes <- genes_annot$gene_id
  ngenes <- length(genes)

  if(ngenes > 0) {
    chrvals <- parallel::mclapply(X = 1:(ngenes-1), mc.cores = MCCORES,  mc.cleanup = T, FUN = function(i) {
      gene1 <- genes[i]
      s1 <- genes_annot %>% filter(gene_id == gene1) %>% pull(start)
      other.genes <- genes[(i+1):ngenes]
      mivals <- lapply(X = other.genes,
                       FUN = function(gene2) {
                         s2 <- genes_annot %>% filter(gene_id == gene2) %>% pull(start)
                         list(source = gene1, target = gene2,
                              distance = max(s2, s1) - min(s2, s1),
                              mi = max(MImatrix[gene1, gene2], MImatrix[gene2, gene1],  na.rm = T))
                       })
      return(bind_rows(mivals))
    })
    chrdf <- bind_rows(chrvals)
    chrdf$chr <- ch
    return(chrdf)
  }
})

distvals <- bind_rows(distvals)
distvals$cond <- COND
cat("Saving file.\n")
fwrite(distvals, file = snakemake@output[[1]], row.names = F, col.names = T, sep = "\t", nThread = MCCORES)
