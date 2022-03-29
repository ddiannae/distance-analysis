log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(vroom)
library(dplyr)
library(tidyr)
library(janitor)
library(rstatix)
library(purrr)

load(snakemake@params[["annot"]])
MCCORES <- snakemake@threads[[1]]

cat("Reading files\n")

bin_mus <- vroom::vroom(snakemake@input[["counts"]])
mi_vals <- vroom::vroom(snakemake@input[["mi_matrix"]])
cat("MI matrix with ", nrow(mi_vals), " rows and ", ncol(mi_vals), " columns loaded \n")

genes <- colnames(mi_vals)
mi_vals$source <- genes

mi_vals <- mi_vals %>% tidyr::pivot_longer(cols = starts_with("ENSG"), 
                                     names_to = "target",
                                     values_to = "mi",
                                     values_drop_na = TRUE) %>% 
  dplyr::arrange(desc(mi)) %>% 
  dplyr::left_join(annot %>% dplyr::select(gene_id, source_chr = chr), 
            by = c("source" = "gene_id")) %>% 
  dplyr::left_join(annot %>% dplyr::select(gene_id, target_chr = chr), 
            by = c("target" = "gene_id")) %>% 
  dplyr::mutate(interaction = ifelse(source_chr == target_chr, "intra", "inter"),
         nrow = row_number()) %>%
  dplyr::mutate(bin_size = 10^floor(log10(nrow)),
          bin_size = ifelse(bin_size < 1000, 1000, bin_size),
          bin = (floor((nrow-1)/bin_size)+1)*bin_size)

cat("Got a total of", nrow(mi_vals), " interactions \n")

rm(annot)

n_tests <- 1000

bin_mus <- bin_mus %>%  
  dplyr::filter(bin <= 5e7) %>% 
  dplyr::select(bin, intra_fraction)

cat("Getting gene lists \n")
t1 <- Sys.time()
genes_bin <-  mi_vals %>% 
  dplyr::filter(bin <= 5e7) %>%
  dplyr::group_split(bin) %>%
  purrr::map(~union(.x$source, .x$target))
t2 <- Sys.time()
cat("Time to get gene list: ", t2-t1, "\n")

## No podemos asegurar que todos los genes en un chunk incluyan a los del anterior
for(i in 2:length(genes_bin)) {
  genes_bin[[i]] = union(genes_bin[[i]], genes_bin[[i-1]])
}

cat("Getting interactions lists \n")
t1 <- Sys.time()
interactions <- purrr::map(genes_bin, function(gbin) {
  if(length(gbin) == length(genes)) {
    return("all")
  }
  mi_vals %>%
    filter(source %in% gbin & target %in% gbin) %>%
    pull(interaction)
})
t2 <- Sys.time()
cat("Time to get interactions list: ", t2-t1, "\n")

cat("Getting all samples \n")
t1 <- Sys.time()
all_samples <- mapply(FUN = function(inter_set, bin) {
  cat(length(inter_set), " ", bin, "\n")
  bin_samples <- parallel::mclapply(X = 1:n_tests, function(i) {
    if(length(inter_set) == 1) {
      sum(sample(mi_vals$interaction, bin, replace = FALSE)  == "intra")
    } else {
      sum(sample(inter_set, bin, replace = FALSE)  == "intra")
    }
  },  mc.cores = MCCORES) 
}, interactions, bin_mus$bin, SIMPLIFY = FALSE) 
t2 <- Sys.time()
cat("Time to get samples list: ", t2-t1, "\n")

all_samples <- purrr::map2(all_samples, bin_mus$bin,
                           ~dplyr::tibble(n = unlist(.x)/.y))

names(all_samples) <-  bin_mus$bin

cat("Saving ttest statistics \n")
purrr::map2_df(all_samples, bin_mus$intra_fraction,
       ~ .x %>% rstatix::t_test(n~1, mu = .y, detailed = TRUE) %>% 
         janitor::clean_names() %>%
         dplyr::select(estimate, conf_low, conf_high, statistic, p) %>%
         dplyr::mutate(mu = .y), .id = "bin") %>% 
  vroom::vroom_write(snakemake@output[[1]])
