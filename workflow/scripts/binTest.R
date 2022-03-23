log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(vroom)
library(dplyr)
library(rstatix)

cat("Loading data \n")
COND <- snakemake@params[["cond"]]
BINSIZE <- snakemake@params[["binsize"]]
BINTYPE <- snakemake@params[["bintype"]]
MCCORES <- as.numeric(snakemake@threads[[1]])

cat("Reading data \n")
dist_df <- vroom::vroom(snakemake@input[[1]], 
                    col_types = cols("chr" = col_character()))
chrs <- c(as.character(1:22), "X", "Y")

dist_df <- dist_df %>% arrange(distance)

cat("Getting all bins stats \n")
if(BINTYPE == "size") {
  dist_df <- dist_df %>% 
    mutate(id = row_number(), bin = ((id-1)%/%as.integer(BINSIZE)) + 1)
} else if(BINTYPE == "distance") {
   dist_df <- dist_df %>%
     mutate(bin = (distance%/%as.integer(BINSIZE))+1)
}

all_tests <- parallel::mclapply(X=1:(max(dist_df$bin)-2),
                                mc.cores = MCCORES,
                                FUN = function(j) {
# all_tests <-lapply(1:(max(dist_df$bin)-2), 
#                         function(j) {
  cat("Testing bin ", j, "\n")
  wtest <- lapply((j+1):(max(dist_df$bin)-1), function(i) {
    dist_df %>% 
      filter(bin %in% c(j, i)) %>% 
      wilcox.test(mi ~ bin, data =.,) %>% 
      tidy() %>%
      mutate(bin1 = j, bin2 = i)
  })
  wtest <- bind_rows(wtest)
  return(wtest)
})
all_tests <- bind_rows(all_tests)
all_tests$cond <- COND

cat("Saving all bins file \n")
vroom::vroom_write(all_tests, file = snakemake@output[[1]])
