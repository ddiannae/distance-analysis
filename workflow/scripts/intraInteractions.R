## #############################################################
## This file filters the intra-chromosomal interactions and
## gets the distance in base pairs between each pair of genes
## Its input comes from the MI matrix.
################################################################

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(vroom)
library(dplyr)
library(tidyr)

cat("Reading files \n")

ANNOT_RDATA <- snakemake@params[["annot"]]
load(ANNOT_RDATA)

chrs <- c(as.character(1:22), "X", "Y")
MImatrix <- vroom::vroom(snakemake@input[["mi_matrix"]])
genes <- colnames(MImatrix)

cat("MI matrix with ", nrow(MImatrix), " rows and ", ncol(MImatrix), " columns loaded \n")
MImatrix$source <- genes

cat("Getting intra-chromosomal interactions\n")
mi_vals <- MImatrix %>% pivot_longer(cols = starts_with("ENSG"), 
                                     names_to = "target",
                                     values_to = "mi",
                                     values_drop_na = TRUE) %>% 
  arrange(desc(mi)) %>% 
  left_join(annot %>% dplyr::select(gene_id, chr, start), 
            by = c("source" = "gene_id")) %>% 
  rename("source_chr" = "chr", "source_start" = "start") %>%
  left_join(annot %>% dplyr::select(gene_id, chr, start), 
            by = c("target" = "gene_id")) %>% 
  rename("target_chr" = "chr", "target_start" = "start") %>%
  mutate(interaction = ifelse(source_chr == target_chr, "intra", "inter"),
         nrow = row_number())
rm(MImatrix)

annot <- annot %>% filter(gene_id %in% genes)

cat("Calculating distance between genes\n")
mi_vals <- mi_vals %>% filter(interaction == "intra") %>% 
  mutate(distance = pmax(source_start, target_start) - 
           pmin(source_start, target_start)) %>% 
  rename("chr" = "source_chr") %>%
  select(source, target, distance, mi, chr)

cat("Saving file.\n")
vroom::vroom_write(mi_vals, file = snakemake@output[[1]], delim = "\t")
