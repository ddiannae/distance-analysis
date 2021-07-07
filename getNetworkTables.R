library(data.table)
library(readr)
library(parallel)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 5 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  ANNOT <- args[1]
  MATRIX <- args[2]
  INTERCUT <- as.numeric(args[3])
  OUTDIR <- args[4]
  MCCORES <- as.numeric(args[5])
}

cat("Loading files\n")
load(ANNOT)
MImatrix <- fread(file = MATRIX, header = T,  sep = ",", nThread = MCCORES)
MImatrix <- as.matrix(MImatrix)
#max <- nrow(MImatrix)
max <-  nrow(MImatrix)
MIvals <- mclapply(X = 1:max, mc.cores = MCCORES, FUN = function(i) {
  return(unlist(MImatrix[i, (i+1):ncol(MImatrix)], use.names = F))
})

MIvals <- unlist(MIvals)
MIvals <- sort(MIvals, decreasing = T)
cutoff <- MIvals[INTERCUT]
genes <- colnames(MImatrix)

MIvals <- which(MImatrix>=cutoff, arr.ind = T)
mis <- apply(X = MIvals, MARGIN = 1, FUN = function(r){
  return(MImatrix[r["row"], r["col"]])
})

MIvals <- as_tibble(MIvals) %>% mutate(source = genes[row], target = genes[col])
MIvals$mi <- mis
MIvals <- MIvals %>% arrange(desc(mi)) %>% mutate(row_num = 1:nrow(MIvals))

cat("Merging annotations\n")
cond <- strsplit(strsplit(MATRIX, split = "/")[[1]][8], split="_")[[1]][1]
MIvals$cond <- cond

annot <- annot %>% select(gene_id, chr, start, end, ensembl_id, gene_name)
colnames(annot) <-  c("source", "source_chr", "source_start", "source_end", "source_ensembl", "source_name")
MIvals <- merge(annot, MIvals)
colnames(annot) <-  c("target", "target_chr", "target_start", "target_end",  "target_ensembl", "target_name")
MIvals <- merge(annot, MIvals)
MIvals <- MIvals %>% mutate(inter = if_else(source_chr == target_chr,  F, T), 
                      interaction_type = if_else(inter == T, "Inter", "Intra"),
                      distance = if_else(inter == F, pmax(source_start, target_start) - 
                                           pmin(source_start, target_start), as.integer(NaN)))
targets <- MIvals %>% select(target_ensembl, target_chr, target_start, target_end, target_name)
sources <- MIvals %>% select(source_ensembl, source_chr, source_start, source_end, source_name)
colnames(targets) <- c("ensembl", "chr", "start", "end",  "symbol")
colnames(sources)  <-  c("ensembl", "chr", "start", "end", "symbol")
vertices <- bind_rows(targets, sources)
vertices <- vertices[!duplicated(vertices$ensembl), ]

MIvals <- MIvals %>% 
  select(source_ensembl, target_ensembl, mi, distance, row_num, interaction_type, cond) %>%
  arrange(row_num) 

cat("Saving files\n")
write_tsv(vertices, file = paste0(OUTDIR, "/", cond, "_vertices_", INTERCUT, ".tsv"))
write_tsv(MIvals, file = paste0(OUTDIR, "/", cond, "_interactions_",  INTERCUT, ".tsv"))

