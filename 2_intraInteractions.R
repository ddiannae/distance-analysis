library(data.table)

getMIDistanceListByCondAndChrs <- function(types) {
  annot <- read.delim("/labs/csbig/regulaciontrans/distance-analysis/Biomart_EnsemblG94_GRCh38_p12.txt", stringsAsFactors = F, 
                      col.names = c("ensemblID", "chr", "start", "end", "GC",  "type", "symbol") )
  annot <- annot[!is.na(annot$ensemblID), ]
  annot <- annot[!duplicated(annot$ensemblID), ]
  rownames(annot) <- annot$ensemblID
  chrs <- c(as.character(1:22), "X")
  
  conditions <- c("healthy", "cancer")
  for(type in types) {
    cat("Working with type", type, "\n")
    setwd(paste0("/labs/csbig/regulaciontrans/", type, "/networks"))
    
    for(cond in conditions) {
      cat("Working with condition", cond, "\n")
      MImatrix <- fread(file = paste0(type, "-", cond, "-matrix.adj"), header = T, 
                        sep = "\t", nThread = 5)
      MImatrix <- data.matrix(MImatrix)
      rownames(MImatrix) = colnames(MImatrix)
      annot.ok <- annot[colnames(MImatrix),]
      condvals <- lapply(X = chrs, FUN = function(chr) {
        cat("Working with chromosome", chr, "\n")
        genes.annot <- annot.ok[annot.ok$chr == chr, ]
        genes <- genes.annot$ensemblID
        ngenes <- length(genes)
        chrvals <- parallel::mclapply(X = 1:(ngenes-1), mc.cores = 6,  mc.cleanup = T, FUN = function(i) {
          #chrvals <- lapply(X = 1:(ngenes-1),  FUN = function(i) {
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
        chrdf$chr <- chr
        write.table(chrdf, file =  paste0("../intra/", cond, "/chr-", chr, "-distance-mi.txt"), 
                    row.names = F, col.names = T, sep = "\t", quote = F)  
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

getMIDistanceListByCondAndChrs(c("tiroides"))

getAllMIDistanceMeans <- function(binsize, types) {
  chrs <- c(as.character(1:22), "X")
  for (type in types){
    setwd(paste0("/labs/csbig/regulaciontrans/", type, "/intra"))
    conds <- c("cancer")
    conditiondfs <- lapply(X = conds, FUN = function(cond){
      cat("Working with condition ", cond, "\n")
      
      meansbych <- parallel::mclapply(X = chrs, mc.cores = 7,  mc.cleanup = FALSE, FUN = function(chr){
        dist.mi.df <- read.delim(file=paste0(cond, "/chr_", chr, "_distance_mi.txt"), header = T)
        dist.mi.df <- dist.mi.df[order(dist.mi.df$distance), ]
        rownames(dist.mi.df) <- NULL
        dist.mi.df$bin <- ((as.numeric(rownames(dist.mi.df)) - 1)%/%binsize) + 1
        dfmeans <- aggregate(cbind(distance, mi)~bin, data=dist.mi.df, FUN=mean, na.rm=TRUE)
        dfmeans$cond <- cond
        dfmeans$ch <- chr
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
 
#  
# getAllMIDistanceMeans(100)