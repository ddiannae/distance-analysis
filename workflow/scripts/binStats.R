log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(dplyr)

cat("Loading data \n")
COND <- snakemake@params[["cond"]]
BINSIZE <- snakemake@params[["binsize"]]
BINTYPE <- snakemake@params[["bintype"]]
MCCORES <- as.numeric(snakemake@threads[[1]])

cat("Reading data \n")
dist_df <- read_tsv(snakemake@input[[1]], col_types = cols("chr" = col_character()))
chrs <- c(as.character(1:22), "X", "Y")

## Gets the mean MI value for a set of gene-pairs
## size of the set in binsize variable.
## Bins are sorted according to the distance between genes

if(BINTYPE == "size") {
  dist_df <- dist_df %>% group_by(chr) %>% arrange(distance) %>%
    mutate(id = row_number(), bin = ((id-1)%/%as.integer(BINSIZE)) + 1)
} else if(BINTYPE == "distance") {
  dist_df <- dist_df %>% group_by(chr) %>% arrange(distance) %>%
    mutate(bin = distance%/%as.integer(BINSIZE))
}

dfmeans <- dist_df %>% group_by(chr, bin) %>% summarise(dist_mean = mean(distance), dist_sd = sd(distance),
                                                    dist_median = median(distance), dist_max = max(distance),
                                                    dist_min = min(distance), mi_mean = mean(mi), mi_sd = sd(mi),
                                                    mi_median = median(mi), mi_max = max(mi),
                                                    mi_min = min(mi))
dfmeans$cond <- COND

cat("Saving bychr file \n")
write_tsv(dfmeans, snakemake@output[["by_chr"]])

dist_df <- dist_df  %>% ungroup() %>% arrange(distance)

cat("Getting all bins stats \n")
if(BINTYPE == "size") {
  dist_df <- dist_df %>% 
    mutate(id = row_number(), bin = ((id-1)%/%as.integer(BINSIZE)) + 1)
} else if(BINTYPE == "distance") {
 dist_df <- dist_df %>%
   mutate(bin = (distance%/%as.integer(BINSIZE))+1)
}

allmeans <- dist_df %>% group_by(bin) %>% summarise(dist_mean = mean(distance), dist_sd = sd(distance),
                                                    dist_median = median(distance), dist_max = max(distance),
                                                    dist_min = min(distance), mi_mean = mean(mi), mi_sd = sd(mi),
                                                    mi_median = median(mi), mi_max = max(mi),
                                                    mi_min = min(mi))

allmeans$cond <- COND
cat("Saving all bins file \n")
write_tsv(allmeans, snakemake@output[["all"]])
