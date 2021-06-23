library(readr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 4 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  ANNOT <- args[1]
  NETWORK <- args[2]
  INTERACTIONS <- args[3]
  VERTICES <- args[4]
}

cat("Loading files\n")
load(ANNOT)
net <- read_tsv(NETWORK)

cat("Merging annotations\n")
cond <- strsplit(strsplit(NETWORK, split = "/")[[1]][8], split="_")[[1]][1]
net$cond <- cond
net$row_num <- 1:nrow(net)

annot <- annot %>% select(gene_id, chr, start, end, ensembl_id, gene_name)
colnames(annot) <-  c("source", "source_chr", "source_start", "source_end", "source_ensembl", "source_name")
net <- merge(annot, net)
colnames(annot) <-  c("target", "target_chr", "target_start", "target_end",  "target_ensembl", "target_name")
net <- merge(annot, net)
net <- net %>% mutate(inter = if_else(source_chr == target_chr,  F, T), 
                      interaction_type = if_else(inter == T, "Inter", "Intra"),
                      distance = if_else(inter == F, pmax(source_start, target_start) - 
                                               pmin(source_start, target_start), as.integer(NaN)))
targets <- net %>% select(target_ensembl, target_chr, target_start, target_end, target_name)
sources <- net %>% select(source_ensembl, source_chr, source_start, source_end, source_name)
colnames(targets) <- c("ensembl", "chr", "start", "end",  "name")
colnames(sources)  <-  c("ensembl", "chr", "start", "end", "name")
vertices <- bind_rows(targets, sources)
vertices <- vertices[!duplicated(vertices$id), ]

net <- net %>% 
        select(source_ensembl, target_ensembl, mi, distance, row_num, interaction_type) %>%
        arrange(row_num) 

cat("Saving files\n")
write_tsv(vertices, file = VERTICES)
write_tsv(net, file = INTERACTIONS)

