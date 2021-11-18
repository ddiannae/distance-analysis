log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(dplyr)
library(tidyr)
library(janitor)
library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)
library(DOSE)

cat("Reading files\n")

etype <-  snakemake@params[["enrich_type"]]
membership <- read_tsv(file = snakemake@input[["membership"]])
gene_universe <- read_tsv(file = snakemake@input[["universe"]], 
                          col_names = c("ensembl_id")) 

gene_universe <- gene_universe %>% 
  separate(col = ensembl_id, into = c("ensembl", "version"), sep = "\\.") %>%
  pull(ensembl)

membership_entrez <- bitr(membership$ensembl, fromType="ENSEMBL", 
                          toType="ENTREZID", OrgDb="org.Hs.eg.db") 
colnames(membership_entrez) <- c("ensembl", "entrez")

membership <- membership %>% full_join(membership_entrez, by = "ensembl")
  
gene_universe_entrez <- bitr(gene_universe, fromType="ENSEMBL", 
                          toType="ENTREZID", OrgDb="org.Hs.eg.db") %>%
  pull(ENTREZID)

if(etype == "kegg") {
  enrich_function <- function(glist, guniverse) {
    return(enrichKEGG(gene = glist,
                      universe      = guniverse,
                      keyType       = "kegg",
                      organism      = "hsa",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.1,
                      minGSSize     = 10))
  } 
} else if(etype == "onco") {
  ONCO <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
    dplyr::select(gs_name, entrez_gene)
  
  enrich_function <- function(glist, guniverse) {
    return(enricher(gene = glist,
                      universe = guniverse,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.1,
                      minGSSize     = 10,
                      TERM2GENE = ONCO))
  } 
} else {
  enrich_function <- function(glist, guniverse) {
    return(enrichNCG(gene = glist,
                    universe = guniverse,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.1,
                    minGSSize     = 10))
  }
}

all_enrichments <- lapply(X = unique(membership$community),
                          FUN = function(com){
                            
        cat("Working with community: ", com, "\n")

        gene_list <- membership %>% filter(community == com) %>%
                             pull(entrez)
        
        if(length(gene_list) >= 5) {
          
          terms <- enrich_function(gene_list, gene_universe_entrez)
          
          if(nrow(as.data.frame(terms)) > 0) {
            terms_df <- as.data.frame(terms)
            terms_df$commun <- com 
            
            return(terms_df)  
          }
        }
          
    return(NULL)
})

bind_rows(all_enrichments) %>% 
  janitor::clean_names() %>%
  write_tsv(snakemake@output[[1]])
