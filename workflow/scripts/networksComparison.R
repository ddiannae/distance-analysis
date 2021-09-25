log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(dplyr)
library(tidyr)

cat("Reading files \n")
cancer_interactions <- read_tsv(file=snakemake@input[["inter_cancer"]])
normal_interactions <- read_tsv(file=snakemake@input[["inter_normal"]])
cancer_vertices <-  read_tsv(file=snakemake@input[["ver_cancer"]])
normal_vertices <- read_tsv(file=snakemake@input[["ver_normal"]])

cat("Getting intersections \n")
cancer_interactions <- cancer_interactions %>% 
  mutate(inter = paste0(pmin(source_ensembl,target_ensembl), "-", 
                        pmax(source_ensembl,target_ensembl))) 

normal_interactions <- normal_interactions %>% 
  mutate(inter = paste0(pmin(source_ensembl,target_ensembl), "-", 
                        pmax(source_ensembl,target_ensembl)))

shared_interactions <- cancer_interactions %>% 
  select(inter, distance, interaction_type) %>% 
  inner_join(normal_interactions %>% 
               select(inter, distance, interaction_type)) 

cancer_only <- cancer_interactions %>% 
  anti_join(shared_interactions, by="inter") %>%
  select(-inter) 

normal_only <- normal_interactions %>% 
  anti_join(shared_interactions, by="inter") %>%
  select(-inter) 

normal_vertices %>% 
  filter(ensembl %in% union(normal_only$source_ensembl, 
                                normal_only$target_ensembl)) %>%
  write_tsv(file = snakemake@output[["ver_normal_only"]])

cancer_vertices %>% 
  filter(ensembl %in% union(cancer_only$source_ensembl, 
                            cancer_only$target_ensembl)) %>%
  write_tsv(file = snakemake@output[["ver_cancer_only"]])

shared_interactions <- shared_interactions %>%
  separate(inter, into = c("source_ensembl", "target_ensembl")) 

normal_vertices %>%
  filter(ensembl %in% union(shared_interactions$source_ensembl, 
                            shared_interactions$target_ensembl)) %>%
  write_tsv(file = snakemake@output[["ver_shared"]])


cat("Saving files \n")
shared_interactions %>%
  write_tsv(file = snakemake@output[["inter_shared"]])

normal_only %>% 
  write_tsv(file=snakemake@output[["inter_normal_only"]])

cancer_only %>%
  write_tsv(file=snakemake@output[["inter_cancer_only"]])
