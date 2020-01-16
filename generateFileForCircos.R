library(data.table)

setwd("/media/ddisk/transpipeline-data/")


getMILinks <- function(types) {
  for(type in types) {
    conds <- c("healthy", "cancer")
    for(cond in conds) {
      MIvals <- fread(paste0(type, "/networks/network-tables/", type, "-", cond, "-interactions.tsv"), 
                      header = T, sep="\t", nThread = 3)
      MIvals <- MIvals[order(-MI),]
      MIvals <- MIvals[1:25000,]
      vertices <- fread(paste0(type, "/networks/network-tables/", type, "-", cond, "-vertices.tsv"), 
                        header = T, sep="\t", nThread = 3)
      vertices <- vertices[, c("ensemblID", "chr", "start", "end")]
      colnames(vertices) <-  c("source", "source_chr", "source_start", "source_end")
      MIvals <- merge(vertices, MIvals, by.x = "source", by.y = "source")
      colnames(vertices) <-  c("target", "target_chr", "target_start", "target_end")
      MIvals <- merge(vertices, MIvals, by.x = "target", by.y = "target")
      MIvals$source_chr <- paste0("hs", MIvals$source_chr)
      MIvals$target_chr <- paste0("hs", MIvals$target_chr)
      MIvals <- MIvals[, c("source_chr", "source_start", "source_end",
                           "target_chr", "target_start", "target_end")]
      fwrite(MIvals, paste0(type, "/circos/", cond, "-circos-data.tsv"),
             sep = "\t", col.names = F,row.names = F)
    }
  }
}

types <- c("utero")
getMILinks(types)

