setwd("/media/ddisk/transpipeline-data/")
library(dplyr)

annot <- read.delim("/media/ddisk/transpipeline-data/biomarts/Biomart_EnsemblG94_GRCh38_p12_karyotype.txt",
                    col.names = c("ensemblID", "chr", "start", "end", "band"), stringsAsFactors = F)
annot <- annot[!duplicated(annot$ensemblID), ]

types <- c("utero")
for (type in types) {
  cat("Working with type", type, "\n")
  
  conds <- c("cancer", "healthy")
  net <- lapply(conds, function(cond){
    cat("Reading network ", cond, "\n")
    network <- read.delim(paste0(type, "/networks/1e8/", type, "-", cond, "-min.sif"), 
                          stringsAsFactors = F)
    cat("Is not unsorted: ", !is.unsorted(rev(network$MI)), "\n")
    network <- network %>% arrange(desc(MI))
    network$cond <- cond
    
    network$row_num <- 1:nrow(network)
    return(network)
  })
  net <- bind_rows(net)
  colnames(annot) <-  c("source", "source_chr", "source_start", "source_end", "source_band")
  net <- merge(annot, net)
  colnames(annot) <-  c("target", "target_chr", "target_start", "target_end",  "target_band")
  net <- merge(annot, net)
  net$inter <- ifelse(net$source_chr == net$target_chr, F, T)
  net$interaction_type <- ifelse(net$inter, "Trans", "Inter-Cytoband")
  net$interaction_type <- ifelse(net$source_band == net$target_band & !net$inter, "Intra-Cytoband", 
                                 net$interaction_type)
  
  annot.symbol <- read.delim("/media/ddisk/transpipeline-data/biomarts/Biomart_EnsemblG94_GRCh38_p12.txt",
                             colClasses = c("character",  "NULL",  "NULL",  "NULL", "NULL", "NULL",  "character" ), 
                             stringsAsFactors = F, 
                             col.names = c("ensemblID",  "NULL",  "NULL",  "NULL", "NULL", "NULL",  "symbol"))
  
  annot.symbol <- annot.symbol[!duplicated(annot.symbol$ensemblID), ]
  targets <- net %>% select(target, target_chr, target_start, target_end, target_band)
  sources <- net %>% select(source, source_chr, source_start, source_end, source_band)
  colnames(targets) <- c("ensemblID", "chr", "start", "end", "band")
  colnames(sources)  <- c("ensemblID", "chr", "start", "end", "band")
  vertices <- rbind(targets, sources)
  vertices <- vertices[!duplicated(vertices$ensemblID), ]
  vertices <- vertices %>% inner_join(annot.symbol, by = "ensemblID")
  
  
  housekeeping <- read.delim("/media/ddisk/transpipeline-data/annotations/housekeeping.txt", header = F, 
                             colClasses = c("character", "NULL"))
  housekeeping <- housekeeping[, 1, drop= T]
  housekeeping <- stringr::str_trim(housekeeping)
  housekeeping <- tolower(housekeeping)
  vertices$isHK <- ifelse(tolower(vertices$symbol) %in% housekeeping, T, F)
  vertices %>% filter(isHK == T)
  
  tfs <- read.delim("/media/ddisk/transpipeline-data/annotations/tfs_2018.txt", header = F,
                    col.names = c("ensemblID"))
  vertices$isTF <- ifelse(vertices$ensemblID %in% tfs$ensemblID, T, F)
  
  kariotype <- read.delim("/media/ddisk/transpipeline-data/annotations/chromosome.band.hg38.txt",
                          header = T, stringsAsFactors = F)
  names(kariotype)[1] <- "chr"
  tbands <- kariotype %>% mutate(arm = unlist(lapply(as.vector(strsplit(name, split = "")), "[[", 1)))
  tbands <- tbands %>% group_by(chr, arm) %>% summarise(t = max(name))
  tbands <- tbands %>% group_by(chr) %>% summarise(t1 = min(t), t2 = max(t))
  tbands$chr <- stringr::str_replace(tbands$chr, "chr", "")
  vertices <- vertices %>% inner_join(tbands, by = "chr")
  vertices$isExtreme <- vertices$band == vertices$t1 | vertices$band == vertices$t2
  vertices$isExtreme <- as.integer(vertices$isExtreme)
  vertices$isHK <- as.integer(vertices$isHK)
  vertices$isTF <- as.integer(vertices$isTF)
  
  for(i in 1:length(conds)){
    net.cond <- net%>% filter(cond == conds[i]) %>% select(source, target, MI, row_num, interaction_type) %>%
      arrange(desc(MI)) 
  
    vertices.cond <- vertices %>% 
      filter(vertices$ensemblID %in% net.cond$source | vertices$ensemblID %in% net.cond$target) 
    vertices.cond <- vertices.cond %>% select(-t1, -t2)
    write.table(vertices.cond, file = paste0(type, "/networks/network-tables/",
                                             type,"-", conds[i], "-vertices.tsv") , 
                quote = F, row.names = F, col.names = T, sep = "\t")
    write.table(net.cond, file = paste0(type, "/networks/network-tables/", 
                                        type ,"-", conds[i], "-interactions.tsv"),
                quote = F, row.names = F, col.names = T, sep = "\t")
  }
}
