log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(scales)

cat("Reading files\n")
TISSUE <- snakemake@params[["tissue"]]
COND <- snakemake@params[["cond"]]
TERMS <- as.integer(snakemake@params[["terms"]])

substring(TISSUE, 1, 1) <- toupper(substring(TISSUE, 1, 1))
                                   
chr_assortativity <- read_tsv(snakemake@input[["chr_assortativity"]],
                     col_types = cols_only(community_id = col_double(),
                                           diffraction = col_double()))
if(COND == "cancer") {
  expr_assortativity <- read_tsv(snakemake@input[["expr_assortativity"]]) %>% 
    select(community_id, diffraction, mean_log_fc)
  
  cat("Joining assortativity files\n")
  assort <- chr_assortativity %>% 
    inner_join(expr_assortativity, by = "community_id", suffix = c("_chr", "_exp"))
} else {
  assort <- chr_assortativity %>% rename("diffraction_chr" = "diffraction")
}

cat("Getting communities info\n")
comm_info <- read_tsv(snakemake@input[["comm_info"]], 
                      col_types = cols_only(com_id = col_double(), order = col_double(), 
                                            pg_gene = col_character())) %>%
    left_join(read_tsv(snakemake@input[["vertices"]],
                     col_types = cols_only(ensembl = col_character(), symbol = col_character())),
            by = c("pg_gene" = "ensembl"))

cat("Getting enrichments\n")
comm_enrich <- read_tsv(snakemake@input[["enrich"]], col_types = cols( p_adjust = col_double())) %>%
  filter(p_adjust < 1e-10) %>%
  select(id, commun) %>% group_by(commun) %>% 
  tally(name = "nterms") 

cat("Joining plot data\n")
plot <- assort %>% 
  inner_join(comm_info, by = c("community_id" = "com_id")) %>%
  left_join(comm_enrich, by = c("community_id" = "commun")) 

plot[is.na(plot$nterms), "nterms"] <- 0

cat("Writing plot data\n")

plot %>% filter(nterms > TERMS) %>% 
  dplyr::select(community_id) %>% write_tsv(snakemake@output[["comm_assortativity_nterms"]])

if(COND == "cancer") {
  plot %>% select(community_id, symbol, order, diffraction_chr, diffraction_exp, mean_log_fc, nterms) %>%
    rename(size = order, name = symbol, chr_assortativity = diffraction_chr, 
           exp_assprtativity = diffraction_exp, enriched_terms = nterms) %>% 
    write_tsv(snakemake@output[["comm_assortativity_tsv"]])
} else {
  plot %>% select(community_id, symbol, order, diffraction_chr, nterms) %>%
    rename(size = order, name = symbol, chr_assortativity = diffraction_chr,
           enriched_terms = nterms) %>% 
    write_tsv(snakemake@output[["comm_assortativity_tsv"]])
}

if(COND == "cancer") {
  cat("Building plot\n")
  p <- ggplot(plot, aes(x = diffraction_chr, y = diffraction_exp)) + 
    geom_point(aes(size = order, color = nterms)) +
    geom_text(aes(label = ifelse(nterms > TERMS, as.character(symbol), NA)),
              colour = "black", size = 1.3, check_overlap = FALSE, fontface = "bold") +
    geom_hline(yintercept = 0.0, linetype="dashed", color = "gray") +
    geom_vline(xintercept = 0.0, linetype="dashed", color = "gray") +
    theme_base() +
    labs(size = "Nodes in \ncommunity",
         color = "Enriched \nterms") +
    xlab("Chromosomal assortativity") +
    ylab("Expression assortativity") +
    scale_color_gradient(low="lightblue1", high="steelblue4", breaks = c(0, 10, 20, 30),
                         limits = c(0, 30), oob = squish) +
    scale_size_continuous(range = c(1, 7), breaks = c(1, 5, 20, 50, 100, 150, 200)) +
    theme(text = element_text(size = 6), 
          axis.title = element_text(size = 6),
          legend.text = element_text(size = 4),
          plot.background=element_blank(), legend.key.size = unit(0.6, "line"),
          legend.spacing.y = unit(0.05, 'cm')) +
    ggtitle(TISSUE)
} else {
  
  p <- ggplot(plot, aes(x = community_id, y = diffraction_chr)) + 
    geom_point(aes(color = nterms, size = order)) +
    ylab("Chromosomal assortativity") +
    xlab("") +
    geom_hline(yintercept = 0.0, linetype="dashed", color = "gray") +
    labs(size = "Nodes in \ncommunity",
         color = "Enriched \nterms") +
    geom_text(aes(label = ifelse(nterms > TERMS, as.character(symbol), NA)),
              colour = "black", size = 1.3, check_overlap = FALSE, fontface = "bold") +
    scale_color_gradient(low="lightblue1", high="steelblue4", breaks = c(0, 10, 20, 30),
                         limits = c(0, 30), oob = squish) +
    scale_size_continuous(range = c(1, 7), breaks = c(1, 5, 20, 50, 100, 150, 200)) +
   theme_base() +
    theme(text = element_text(size = 6), 
          axis.title = element_text(size = 6),
          legend.text = element_text(size = 4),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.background = element_blank(), legend.key.size = unit(0.6, "line"),
          legend.spacing.y = unit(0.05, 'cm')) +
    ggtitle(TISSUE)
}
cat("Saving plot\n")
png(snakemake@output[["comm_assortativity_png"]], 
    units="in", width=4, height=2.5, res=300)
print(p)
dev.off()