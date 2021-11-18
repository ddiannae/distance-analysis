log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(dplyr)
library(tidyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

cat("Reading files\n")

membership <- read_tsv(file=snakemake@input[["membership"]])
gene_universe <- read_tsv(file=snakemake@input[["universe"]], col_names = c("ensembl_id")) 

gene_universe <- gene_universe %>% 
  separate(col = ensembl_id, into = c("ensembl", "version"), sep = "\\.") %>%
  pull(ensembl)

all_enrichments <- lapply(X = unique(membership$community),
                          FUN = function(com){
			
		cat("Working with community: ", com, "\n")

    gene_list <- membership %>% filter(community == com) %>%
                        pull(ensembl)
    
    if(length(gene_list) >= 5) {
      
      terms <- enrichGO(gene          = gene_list,
                      universe      = gene_universe,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENSEMBL",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.005,
                      qvalueCutoff  = 0.01,
                      minGSSize     = 10,
                      readable      = FALSE)
      
      tryCatch({
        bp <- pairwise_termsim(terms)
        bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
        simple_results <- bp2@result
        if(nrow(simple_results) > 0) {
          simple_results$commun <- com 
        }
        return(simple_results)
      }, error = function(cond) {
        return(NULL)
      })
     
    }
    return(NULL)
})

bind_rows(all_enrichments) %>% 
  janitor::clean_names() %>%
  write_tsv(snakemake@output[[1]])


