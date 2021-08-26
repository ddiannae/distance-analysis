log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(vroom)
library(tidyr)
library(dplyr)

CUTOFF <- as.numeric(snakemake@params[["cutoff"]])
COND <- snakemake@params[["cond"]]

cat("Loading files\n")
load(snakemake@params[["annot"]])
head(annot)
MImatrix <- vroom::vroom(snakemake@input[["mi_matrix"]])
genes <- colnames(MImatrix)

cat("MI matrix with ", nrow(MImatrix), " rows and ", ncol(MImatrix), " columns loaded \n")
MImatrix$source <- genes

cat("Annotating interactions\n")
MIvals <- MImatrix %>% pivot_longer(cols = starts_with("ENSG"), 
                                     names_to = "target",
                                     values_to = "mi",
                                     values_drop_na = TRUE) %>% 
  arrange(desc(mi)) %>% mutate(row_num = row_number()) %>% filter(row_num <= CUTOFF)

cat("Merging annotations\n")
MIvals$cond <- COND

annot <- annot %>% select(gene_id, chr, start, end, ensembl_id, gene_name)
colnames(annot) <-  c("source", "source_chr", "source_start", "source_end", "source_ensembl", "source_name")
MIvals <- merge(annot, MIvals)
colnames(annot) <-  c("target", "target_chr", "target_start", "target_end",  "target_ensembl", "target_name")
MIvals <- merge(annot, MIvals)
MIvals <- MIvals %>% mutate(inter = if_else(source_chr == target_chr,  F, T), 
                      interaction_type = if_else(inter == T, "Inter", "Intra"),
                      distance = if_else(inter == F, as.integer(pmax(source_start, target_start) - 
                                           pmin(source_start, target_start)), as.integer(-1)))

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
vroom_write(vertices, file = snakemake@output[["vertices"]])
vroom_write(MIvals, file = snakemake@output[["interactions"]])

