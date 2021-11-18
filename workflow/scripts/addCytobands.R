log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(dplyr)

cat("Reading file \n")

vertices_cancer <- read_tsv(snakemake@input[["vertices_cancer"]])
vertices_normal <- read_tsv(snakemake@input[["vertices_normal"]])

cytobands <- read_tsv(snakemake@params[["annot_cytobands"]])
colnames(cytobands) <- c("ensembl", "cytoband")

vertices_cancer <- vertices_cancer %>% left_join(cytobands, by = "ensembl")
vertices_normal <- vertices_normal %>% left_join(cytobands, by = "ensembl")

cat("Writing files \n")

write_tsv(vertices_cancer, file=snakemake@output[["vertices_cancer"]])
write_tsv(vertices_normal, file=snakemake@output[["vertices_normal"]])

doneFile <- file(snakemake@output[["done"]])
writeLines(c("Done!"), doneFile)
close(doneFile)
