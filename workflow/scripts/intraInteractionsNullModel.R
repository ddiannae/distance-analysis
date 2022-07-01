## #############################################################
## This file executes t-tests to evaluate the number of 
## intra-interactions in a MI threshold versus the number of 
## possible interactions in a  cytoband. 
## Its input comes from the MI matrix and the 
## intraInteractionsCount.R file
################################################################

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(vroom)
library(dplyr)
library(tidyr)
library(janitor)
library(rstatix)
library(parallel)

load(snakemake@params[["annot"]])
MCCORES <- snakemake@threads[[1]]

cat("Reading files\n")

regions <- vroom::vroom(snakemake@params[["annot_cytoband"]])
colnames(regions) <- c("ensembl_id", "cytoband")
annot <- annot %>% inner_join(regions, by="ensembl_id")
rm(regions)

intra_count <- vroom::vroom(snakemake@input[["intra_count"]])
top <- as.numeric(snakemake@params[["ninter"]])

MImatrix <- vroom::vroom(snakemake@input[["mi_matrix"]])
cat("MI matrix with ", nrow(MImatrix), " rows and ", ncol(MImatrix), " columns loaded \n")

genes <- colnames(MImatrix)
MImatrix$source <- genes

## se obtienen solo las intra
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
  filter(nrow <= top & interaction == "intra")

rm(MImatrix)

## filtramos las anotaciones para quedarnos solo con las de los genes en la matriz de 
## expresión
annot <- annot %>% filter(gene_id %in% genes)

## Obtenemos el count real y total de las interacciones. 
## Solamente nos interesan las fracciones maypres a 0.5.
## sumamos las interacciones vistas y totales
intra_count <- intra_count %>%
  filter(top <= !!top) %>%
  group_by(chr, cytoband) %>%
  summarise(fraction = sum(fraction), total_inter = sum(total_inter)) %>% 
  arrange(chr, cytoband) %>%
  filter(fraction >= 0.5)

n_tests <- 10000

all_chrs <- parallel::mclapply(X = unique(intra_count$chr), FUN = function(chr) {
 
  ## interacciones reales por citobanda. Solo nos interesan las mayores de 0.5
  intra_cyto_count <- intra_count %>% 
    filter(chr == !!chr) %>% 
    pull(fraction, name = cytoband)

  if(length(intra_cyto_count) > 0) {
    cat("\tWorking with chromosome", chr, "\n")
    
    ## sacamos los genes del cromosoma
    annot_chr <- annot %>% filter(chr == !!chr)
    
    ## total de interacciones reales en el cromosoma
    chr_intra <- mi_vals %>% filter(source_chr == chr & target_chr == chr) %>% nrow()

    ## combinaciones posibles entre genes del cromosoma
    chr_gene_pairs <- expand_grid(source = annot_chr$gene_id, 
                                  target = annot_chr$gene_id) %>% 
      filter(source != target)

    ## número posible de interacciones por citobanda 
    ## dados los genes en la matriz de expresión
    total_cyto <- intra_count %>% 
      filter(chr == !!chr) %>% 
      pull(total_inter, name = cytoband)

    cat("Getting samples \n")
    ### se repite n.test veces
    all_samples <- lapply(1:n_tests, function(i) {
      sample_intrak <- slice_sample(chr_gene_pairs, n = chr_intra) 
      sample_intrak$i <- i
      return(sample_intrak)
    })
    
    cat("Grupping samples \n")
    ## distribución del número de interacciones intra-citobanda obtenidas
    ## en los tests, por cada citobanda
    all_samples <- bind_rows(all_samples) %>% 
      left_join(annot_chr %>% select(gene_id, cytoband), by = c("source" = "gene_id")) %>%
      left_join(annot_chr %>% select(gene_id, cytoband), by = c("target" = "gene_id"),
                suffix = c("_source", "_target")) %>%
      filter(cytoband_source == cytoband_target) %>%
      group_by(cytoband_source, i) %>% tally() %>%
      rename("cytoband" = "cytoband_source", "sample_n" = "n") %>%
      pivot_wider(id_cols = i, names_from = cytoband, values_from = sample_n, values_fill = 0)
    
    cat("Running T tests\n")
    t_test_results <- lapply(names(intra_cyto_count), function(cyt) {
      all_samples %>% mutate(across(all_of(cyt))/total_cyto[cyt]) %>%
        rstatix::t_test(as.formula(paste(cyt, "~1" )), mu = intra_cyto_count[cyt],   detailed = TRUE)
    })
    
    t_test_results <- bind_rows(t_test_results) %>% 
      janitor::clean_names() %>%
      rename("cytoband" = "y") %>%
      mutate(chr = chr, mu = intra_cyto_count[cytoband]) %>%
      select(chr, cytoband, mu, statistic, estimate, conf_low, conf_high, p)
    
    return(t_test_results)
  }
  NULL
}, mc.cores = 2)

bind_rows(all_chrs)  %>% 
  vroom::vroom_write(snakemake@output[[1]])
