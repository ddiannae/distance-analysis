library(data.table)

# Intra-chromosomal interactions


## Function that gets intra-chromosomal interaction pairs
## calculates distance between each pair of genes in the form:
## gene1[start] - gene2[start]
## It requires an adjacency matrix with gene EnsemblIds as columns
## A diagonal matrix would work
## File structure required:
## path
##  /type
##    /networks
##      type-cond-matrix.adj
##    /intra (results will be here)
getMIDistanceListByCondAndChrs <- function(path, types, threads) {
  annot <- read.delim(paste0(path, "biomarts/Biomart_EnsemblG94_GRCh38_p12.txt"), stringsAsFactors = F, 
                      col.names = c("ensemblID", "chr", "start", "end", "GC",  "type", "symbol") )
  annot <- annot[!is.na(annot$ensemblID), ]
  annot <- annot[!duplicated(annot$ensemblID), ]
  rownames(annot) <- annot$ensemblID
  chrs <- c(as.character(1:22), "X")
  
  conditions <- c("healthy", "cancer")
  for(type in types) {
    cat("Working with type", type, "\n")
    setwd(paste0(path, type, "/networks"))
    
    for(cond in conditions) {
      cat("Working with condition", cond, "\n")
      MImatrix <- fread(file = paste0(type, "-", cond, "-matrix.adj"), header = T, 
                        sep = "\t", nThread = threads)
      MImatrix <- data.matrix(MImatrix)
      rownames(MImatrix) = colnames(MImatrix)
      annot.ok <- annot[colnames(MImatrix),]
      condvals <- lapply(X = chrs, FUN = function(c) {
        cat("Working with chromosome", c, "\n")
        genes.annot <- annot.ok[annot.ok$chr == c, ]
        genes <- genes.annot$ensemblID
        ngenes <- length(genes)
        chrvals <- parallel::mclapply(X = 1:(ngenes-1), mc.cores = threads,  mc.cleanup = T, FUN = function(i) {
          gene1 <- genes[i]
          other.genes <- genes[(i+1):ngenes]
          mivals <- lapply(X = other.genes, 
                           FUN = function(gene2) {
                             data.frame(source = gene1, target = gene2, 
                                        distance = max(genes.annot[gene2, "start"], genes.annot[gene1, "start"]) - min(genes.annot[gene1, "start"], genes.annot[gene2, "start"]), 
                                        mi = max(MImatrix[gene1, gene2], MImatrix[gene2, gene1],  na.rm = T))
                           })
          plyr::ldply(mivals)
        })
        chrdf <- plyr::ldply(chrvals)
        chrdf$chr <- c
        return(chrdf)
      })
      rm(annot.ok)
      conddf <- plyr::ldply(condvals)
      conddf$cond <- cond
      cat("Saving file.\n")
      fwrite(conddf, file =  paste0("../intra/", cond,"-all-distance-mi.txt"), 
             row.names = F, col.names = T, sep = "\t", nThread = 5)  
    }
  }
}

#getMIDistanceListByCondAndChrs("/media/ddisk/transpipeline-data/", c("utero"), 6)


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
  chrs <- c(as.character(1:22), "X")
  
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
