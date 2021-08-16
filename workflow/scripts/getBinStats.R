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

dist_df <- read_tsv(snakemake@input[[1]], col_types = cols("chr" = col_character()))
chrs <- c(as.character(1:22), "X", "Y")
cat(COND)
cat(BINSIZE)
cat(BINTYPE)
cat(MCCORES)
## Function that gets the mean MI value for a set of gene-pairs
## size of the set in binsize variable.
## Bins are sorted according to the distance between genes
meansbych <- parallel::mclapply(X = chrs, mc.cores = MCCORES,  mc.cleanup = FALSE, FUN = function(ch){
  cat("Working with chromosome", ch, "\n")
  if(BINTYPE == "size") {
    dist_chr <- dist_df %>% filter(chr == ch) %>% arrange(distance) %>%
      mutate(id = row_number(), bin = ((id-1)%/%as.integer(BINSIZE)) + 1)
  } else if(BINTYPE == "distance") {
    dist_chr <- dist_df %>% filter(chr == ch) %>% arrange(distance) %>%
      mutate(bin = distance%/%as.integer(BINSIZE))
  }

  cat("Getting bin stats for chromosome", ch, "\n")
  dfmeans <- dist_chr %>% group_by(bin) %>% summarise(dist_mean = mean(distance), dist_sd = sd(distance),
                                                      dist_median = median(distance), dist_max = max(distance),
                                                      dist_min = min(distance), mi_mean = mean(mi), mi_sd = sd(mi),
                                                      mi_median = median(mi), mi_max = max(mi),
                                                      mi_min = min(mi))
  if(nrow(dfmeans) > 0) {
    dfmeans$cond <- COND
    dfmeans$chr <- ch
    return(dfmeans)
  }
  return(NULL)
})

cat("Parallel calculations done \n")
meansbych <- bind_rows(meansbych)
cat("Saving bychr file \n")
write_tsv(meansbych, snakemake@output[["by_chr"]])

dist_df <- dist_df %>% arrange(distance)

cat("Getting all bins stats \n")
if(BINTYPE == "size") {
  bins <- seq(1, ceiling(nrow(dist_df)/as.integer(BINSIZE)))
  dist_df <- dist_df %>% arrange(distance) %>% mutate(id = row_number(), bin = ((id-1)%/%as.integer(BINSIZE)) + 1)
} else if(BINTYPE == "distance") {
  bins <- seq(1, ceiling( max(dist_df$distance)/as.integer(BINSIZE)))
  dist_df <- dist_df %>% arrange(distance) %>% mutate(bin = (distance%/%as.integer(BINSIZE))+1)
}

allmeans <- dist_df %>% group_by(bin) %>% summarise(dist_mean = mean(distance), dist_sd = sd(distance),
                                                    dist_median = median(distance), dist_max = max(distance),
                                                    dist_min = min(distance), mi_mean = mean(mi), mi_sd = sd(mi),
                                                    mi_median = median(mi), mi_max = max(mi),
                                                    mi_min = min(mi))

allmeans$cond <- COND
cat("Saving all bins file \n")
write_tsv(allmeans, snakemake@output[["all"]])
