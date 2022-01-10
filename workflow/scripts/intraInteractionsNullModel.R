log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(vroom)
library(dplyr)
library(tidyr)
library(rstatix)
library(parallel)

load(snakemake@params[["annot"]])
MCCORES <- snakemake@threads[[1]]

cat("Reading files\n")

regions <- vroom::vroom(snakemake@params[["annot_cytoband"]])
colnames(regions) <- c("ensembl_id", "cytoband")
intra_count <- vroom::vroom(snakemake@input[["intra_count"]])

annot <- annot %>% inner_join(regions, by="ensembl_id")
rm(regions)
chrs <- c(as.character(1:22), "X", "Y")
MImatrix <- vroom::vroom(snakemake@input[["mi_matrix"]])
cat("MI matrix with ", nrow(MImatrix), " rows and ", ncol(MImatrix), " columns loaded \n")
genes <- colnames(MImatrix)
MImatrix$source <- genes

k <- 1000
tops <- c(10*k, 100*k)

cat("Annotating interactions\n")
mi_vals <- MImatrix %>% pivot_longer(cols = starts_with("ENSG"), 
                                     names_to = "target",
                                     values_to = "mi",
                                     values_drop_na = TRUE) %>% 
  arrange(desc(mi)) %>% 
  left_join(annot %>% dplyr::select(gene_id, chr, cytoband), 
            by = c("source" = "gene_id")) %>% 
  rename("source_chr" = "chr", "source_cytoband" = "cytoband") %>%
  left_join(annot %>% dplyr::select(gene_id, chr, cytoband), 
            by = c("target" = "gene_id")) %>% 
  rename("target_chr" = "chr", "target_cytoband" = "cytoband") %>%
  mutate(interaction = ifelse(source_chr == target_chr, "intra", "inter"),
         nrow = row_number()) %>%
  filter(nrow <= tops[length(tops)] & interaction == "intra")

rm(MImatrix)
annot <- annot %>% filter(gene_id %in% genes)

n_tests <- 1000

all_fractions <- parallel::mclapply(X = tops, mc.cores = min(MCCORES, length(tops)), 
                                    FUN = function(top) {
  all_chrs <- lapply(chrs, function(chr) {
    cat("\tWorking with chromosome", chr, "\n")
    
    ## sacamos los genes del cromosoma
    annot_chr <- annot %>% filter(chr == !!chr)
    
    inter_in_cytobands <- annot_chr %>% group_by(cytoband) %>% 
      tally(name = "n_genes") %>%
      mutate(total_inter = (n_genes*(n_genes-1))/2) %>%
      filter(total_inter > 0) %>% pull(total_inter, name = cytoband)
    
    ## Sacamos las interacciones cis en el cromosoma
    n_inter <- mi_vals %>% filter(source_chr == chr & target_chr == chr) %>% nrow()
    
    all_chr_pairs <- expand_grid(source = annot_chr$gene_id, 
                                 target = annot_chr$gene_id) %>% 
      filter(source != target)
    
    intra_count_chr <- intra_count %>% 
      filter(chr == !!chr & top == !!top) %>% pull(fraction, name = cytoband)
  
    ### se repite n.test veces
    all_samples <- lapply(1:n_tests, function(i) {
      sample_intrak <- slice_sample(all_chr_pairs, n = n_inter) 
      sample_intrak$i <- i
      return(sample_intrak)
    })
  
    all_samples <- bind_rows(all_samples) %>% 
        left_join(annot_chr %>% select(gene_id, cytoband), by = c("source" = "gene_id")) %>%
        left_join(annot_chr %>% select(gene_id, cytoband), by = c("target" = "gene_id"),
                  suffix = c("_source", "_target")) %>%
        filter(cytoband_source == cytoband_target) %>%
        group_by(cytoband_source, i) %>% tally() %>%
        rename("cytoband" = "cytoband_source", "sample_n" = "n") %>%
        pivot_wider(id_cols = i, names_from = cytoband, values_from = sample_n, values_fill = 0)
    
    t_test_results <- lapply(names(intra_count_chr), function(cyt) {
      all_samples %>% mutate(across(all_of(cyt))/inter_in_cytobands[cyt]) %>%
        t_test(as.formula(paste(cyt, "~1" )), 
               mu = intra_count_chr[cyt])
      
    })
    t_test_results <- bind_rows(t_test_results)
    t_test_results$chr <- chr
    return(t_test_results)
  })
  
  all_chrs <- bind_rows(all_chrs) %>% 
    select(chr, .y., p, statistic, n, df) %>%
    rename("cytoband" = ".y.")
  all_chrs$top <- top
  return(all_chrs)
})

bind_rows(all_fractions) %>% 
  vroom::vroom_write(snakemake@output[[1]])
