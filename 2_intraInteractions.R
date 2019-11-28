library(data.table)

annot <- read.delim("/media/ddisk/transpipeline-data/biomarts/Biomart_EnsemblG94_GRCh38_p12.txt", stringsAsFactors = F, 
                    col.names = c("EnsemblID", "Chr", "Start", "End", "GC",  "Type", "Symbol") )
annot <- annot[!is.na(annot$EnsemblID), ]
annot <- annot[!duplicated(annot$EnsemblID), ]
rownames(annot) <- annot$EnsemblID
chrs <- c(as.character(1:22), "X")

type <- "kidney"
conditions <- c("healthy", "tumor")
setwd(paste0("/media/ddisk/transpipeline-data/", type, "-data/networks"))

cond <- conditions[1]

network <- fread(file = paste0(type, "_", cond, ".sif"), header = F, sep = "\t", nThread = 5)
colnames(network) <- c("source", "target", "MI")
colnames(annot) <-  c("source", "source_chr", "source_start", "source_end", "source_band")
network <- merge(annot, network, by = "source")
colnames(annot) <-  c("Target", "Target.Chr", "Target.Start", "Target.End",  "Target.Band")
net <- merge(annot, net)


MImatrix <- data.matrix(MImatrix)
rownames(MImatrix) = colnames(MImatrix)
annot.ok <- annot[rownames(MImatrix),]
annot.ok <- annot.ok[!is.na(annot.ok$EnsemblID), ]
condvals <- lapply(X = chrs, FUN = function(chr) {
  cat("Working with chromosome", chr, "\n")
  genes.annot <- annot.ok[annot.ok$Chr == chr, ]
  genes <- genes.annot$EnsemblID
  ngenes <- length(genes)
  chrvals <- parallel::mclapply(X = 1:(ngenes-1), mc.cores = 5,  mc.cleanup = T, FUN = function(i) {
  #chrvals <- lapply(X = 1:(ngenes-1),  FUN = function(i) {
    gene1 <- genes[i]
    other.genes <- genes[(i+1):ngenes]
    mivals <- lapply(X = other.genes, 
                     FUN = function(gene2) {
                       data.frame(source = gene1, target = gene2, 
                                  distance = max(genes.annot[gene2, "Start"], genes.annot[gene1, "Start"]) - min(genes.annot[gene1, "Start"], genes.annot[gene2, "Start"]), 
                                  mi = max(MImatrix[gene1, gene2], MImatrix[gene2, gene1], na.rm = T))
                     })
    plyr::ldply(mivals)
  })
  chrdf <- plyr::ldply(chrvals)
  chrdf$chr <- chr
  write.table(chrdf, file =  paste("intra/", cond, "/chr_", chr, "_distance_norm_mi.txt", sep = ""), 
              row.names = F, col.names = T, sep = "\t", quote = F)  
  return(chrdf)
})
rm(MImatrix)
rm(annot.ok)
conddf <- plyr::ldply(condvals)
conddf$cond <- cond
cat("Saving file.\n")
fwrite(conddf, file =  paste("intra/", cond,"_all_distance_norm_mi.txt", sep = ""), 
       row.names = F, col.names = T, sep = "\t", nThread = 5)  

getMIDistanceListByCondAndChrs("healthy")



getAllMIDistanceMeans <- function(binsize) {
  conditiondfs <- lapply(X = conds, FUN = function(cond){
    cat("Working with condition ", cond, "\n")

    meansbych <- parallel::mclapply(X = chrs, mc.cores = 7,  mc.cleanup = FALSE, FUN = function(chr){
      dist.mi.df <- read.delim(file=paste("intra/", cond,"/", "chr_", chr, "_distance_norm_mi.txt",
                                          sep=""), header = T)
      dist.mi.df <- dist.mi.df[order(dist.mi.df$distance), ]
      rownames(dist.mi.df) <- NULL
      dist.mi.df$bin <- ((as.numeric(rownames(dist.mi.df)) - 1)%/%binsize) + 1
      dfmeans <- aggregate(cbind(distance, mi)~bin, data=dist.mi.df, FUN=mean, na.rm=TRUE)
      dfmeans$cond <- cond
      dfmeans$ch <- chr
      return(dfmeans)
    })
    df <- plyr::ldply(meansbych)
    write.table(df, file = paste0("intra/intra-fixed-bin-size/",  binsize, "-norm/", cond, "_all.txt"), sep="\t",
                col.names = T, row.names = F, quote = F)
    return(df)
  })

  cdf <- plyr::ldply(conditiondfs)
  write.table(cdf, file = paste("intra/intra-fixed-bin-size/", binsize, "-norm/all.txt", sep=""), sep="\t",
              col.names = T, row.names = F, quote = F)
}
 
getAllMIDistanceMeans(100)