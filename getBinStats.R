library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = T)

if (length(args) < 4 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  DIST_DATA <- args[1]
  OUTDIR <- args[2]
  COND <- args[3]
  BINSIZE <- args[4]
  BINTYPE <- args[5]
  MCCORES <- as.integer(args[6])
}

dist_df <- read_tsv(DIST_DATA, col_types = cols("chr" = col_character()))
chrs <- c(as.character(1:22), "X", "Y")

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

meansbych <- bind_rows(meansbych)
write_tsv(meansbych, paste0(OUTDIR, COND, "-fixed-", BINTYPE, "-bychr-", BINSIZE, ".tsv" ))

dist_df <- dist_df %>% arrange(distance) 

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
write_tsv(allmeans, paste0(OUTDIR, COND, "-fixed-", BINTYPE, "-all-", BINSIZE, ".tsv" ))