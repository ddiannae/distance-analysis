## #############################################################
## This script performs Kolmogorov-Smirnov tests to compare the
## distribution of intra-chromosomal values in cancer and normal
## at different starting thresholds of top MI values.
## Input comes from the intraInterCount.R script
################################################################

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(broom)
library(vroom)
library(dplyr)
library(parallel)

MCCORES <- snakemake@threads[[1]]

cat("Reading files \n")
cancer <- vroom::vroom(snakemake@input[["cancer"]])
normal <- vroom::vroom(snakemake@input[["normal"]])

nn <- nrow(cancer)

cat("Performing ", nn, " KS tests\n")
ks_results <- parallel::mclapply(X = 1:nn, FUN = function(i) {
  cat("Test for bin #", i, "\n")
  broom::tidy(ks.test(cancer$intra_fraction[i:nn], normal$intra_fraction[i:nn])) %>%
    mutate(from_bin = cancer$bin[i])
}, mc.cores = MCCORES)

cat("Writing results\n")
bind_rows(ks_results) %>%
  vroom::vroom_write(snakemake@output[[1]])
