## #############################################################
## This file gets the fraction of intra-chromosomal (and inter-
## chromosomal) interactions for different thresholds of MI 
## values, in both accumulated and chunk manners.
## Input comes directly from the mi matrix.
################################################################

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(vroom)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(furrr)

MCCORES <- snakemake@threads[[1]]
future::plan(future::multisession, workers = MCCORES)

cat("Loading annotation\n")
load(snakemake@params[["annot"]])

cat("Reading matrix\n")
MImatrix <- vroom::vroom(snakemake@input[["mi_matrix"]])
COND <- snakemake@params[["cond"]]

genes <- colnames(MImatrix)
MImatrix$source <- genes

cat("Annotating interactions\n")
mi_vals <- MImatrix %>% pivot_longer(cols = starts_with("ENSG"), 
                                  names_to = "target",
                                  values_to = "mi",
                                  values_drop_na = TRUE) %>% 
  arrange(desc(mi)) %>%  
  left_join(annot %>% dplyr::select(gene_id, chr), 
                                 by = c("source" = "gene_id")) %>% 
  rename("source_chr" = "chr") %>%
  left_join(annot %>% dplyr::select(gene_id, chr), 
                                         by = c("target" = "gene_id")) %>% 
  rename("target_chr" = "chr") %>%
  mutate(interaction = ifelse(source_chr == target_chr, "intra", "inter"),
         nrow = row_number())

rm(MImatrix)
cat("Saving onek bins\n")
### bins by a thousand interactions

onek_chunks <- mi_vals %>%  
  mutate(bin = floor((nrow-1)/1000)) %>%
  group_by(bin) %>%
  count(interaction) %>% 
  pivot_wider(id_cols = bin, names_from = interaction, values_from = n, 
              values_fill = 0) %>%
  mutate(inter_fraction = round(inter/(intra+inter), 4),
         intra_fraction = round(intra/(intra+inter), 4),
         cond = COND) %>%
  ungroup() %T>% 
  vroom_write(file = snakemake@output[["onek_chunks"]])

bins <- unique(onek_chunks$bin)
names(bins) <- unique(onek_chunks$bin)

furrr::future_map_dfr(bins, ~ onek_chunks %>% 
                        filter(bin <= .x) %>% 
                        select(intra, inter) %>% colSums(),
                      .id = "bin") %>%
  mutate(bin = as.numeric(bin),
         inter_fraction = round(inter/(intra+inter), 4),
         intra_fraction = round(intra/(intra+inter), 4),
          cond = COND) %>%
  vroom_write(file = snakemake@output[["onek_bins"]])


cat("Saving log bins\n")
log_chunks <- mi_vals %>% mutate(bin_size = 10^floor(log10(nrow))) %>%
  mutate(bin_size = ifelse(bin_size < 1000, 1000, bin_size),
         bin = (floor((nrow-1)/bin_size)+1)*bin_size) %>% 
  group_by(bin) %>%
  count(interaction) %>% 
  pivot_wider(id_cols = bin, names_from = interaction, values_from = n, 
              values_fill = 0) %>%
  mutate(inter_fraction = round(inter/(intra+inter), 4),
         intra_fraction = round(intra/(intra+inter), 4),
         cond = COND) %>%
  ungroup() %T>%
  vroom_write(file = snakemake@output[["log_chunks"]])

bins <- unique(log_chunks$bin)
names(bins) <- unique(log_chunks$bin)

furrr::future_map_dfr(bins, ~ log_chunks %>% 
                         filter(bin <= .x) %>% 
                         select(intra, inter) %>% colSums(),
                       .id = "bin") %>%
  mutate(bin = as.numeric(bin),
         inter_fraction = round(inter/(intra+inter), 4),
         intra_fraction = round(intra/(intra+inter), 4),
         cond = COND) %>%
  vroom_write(file = snakemake@output[["log_bins"]])