log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(scales)

cat("Reading files\n")
TISSUE <- snakemake@params[["tissue"]]
COND <- snakemake@params[["cond"]]

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
                              pg_gene = col_character(), intra_fraction = col_double())) %>%
    left_join(read_tsv(snakemake@input[["vertices"]],
                     col_types = cols_only(ensembl = col_character(), symbol = col_character())),
            by = c("pg_gene" = "ensembl"))

cat("Getting enrichments\n")
comm_enrich <- read_tsv(snakemake@input[["enrich"]], col_types = cols(p_adjust = col_double())) %>%
  filter(p_adjust < 1e-10) %>%
  select(id, commun) %>% group_by(commun) %>% 
  tally(name = "nterms") 

cat("Joining plot data\n")
plot <- assort %>% 
  inner_join(comm_info, by = c("community_id" = "com_id")) %>%
  left_join(comm_enrich, by = c("community_id" = "commun")) 

plot[is.na(plot)] <- 0
cat("Writing plot data\n")


if(COND == "cancer") {
  plot %>% select(community_id, name = symbol, size = order, intra_fraction, 
                  chr_assortativity = diffraction_chr, exp_assprtativity = diffraction_exp, 
                  mean_log_fc, enriched_terms = nterms) %>%
    write_tsv(snakemake@output[["comm_summary"]])
} else {
  plot %>% select(community_id, name = symbol, size = order, intra_fraction,
                  chr_assortativity = diffraction_chr, enriched_terms = nterms) %>%
    write_tsv(snakemake@output[["comm_summary"]])
}

plot <- plot %>% filter(nterms > 0)
if(COND == "cancer") {
  cat("Building plot\n")
  p <- ggplot(plot, aes(x = diffraction_chr, y = diffraction_exp, label = symbol)) + 
    geom_point(aes(size = pmax(pmin(order, 5),200) , fill = mean_log_fc), 
               color = "darkgrey", shape = 21) +
    geom_text(color = "black", size = 1.3, check_overlap = FALSE, fontface = "bold") +
    geom_hline(yintercept = 0.0, linetype="dashed", color = "gray") +
    geom_vline(xintercept = 0.0, linetype="dashed", color = "gray") +
    theme_base() +
    xlim(c(-1, 1)) +
    ylim(c(-1, 1)) +
    labs(size = "Nodes in \ncommunity", fill = "Mean \nLogFC", 
         x = "Chromosomal assortativity", 
         y = "Expression assortativity") +
    scale_fill_gradient2(low = "#16296B", high = "#660000", limits = c(-2, 2),
                          oob = scales::squish) +
    scale_size(range = c(1, 6), breaks = c(5, 50, 100, 200),
                          limits = c(5, 200)) +
    theme(text = element_text(size = 6), 
          axis.title = element_text(size = 6),
          legend.text = element_text(size = 4),
          plot.background=element_blank(), legend.key.size = unit(0.6, "line"),
          legend.spacing.y = unit(0.1, 'cm')) +
    ggtitle(TISSUE)
} else {
  p <- ggplot(plot, aes(x = community_id, y = diffraction_chr, label = symbol)) + 
    geom_point(aes(color = nterms, size = pmax(pmin(order, 5),200))) +
    ylab("Chromosomal assortativity") +
    xlab("") +
    ylim(c(-1, 1)) +
    geom_hline(yintercept = 0.0, linetype="dashed", color = "gray") +
    labs(size = "Nodes in \ncommunity",
         color = "Enriched \nterms") +
    geom_text(colour = "black", size = 1.3, check_overlap = FALSE, fontface = "bold") +
    scale_color_gradient(low="lightblue1", high="steelblue4", breaks = c(0, 10, 20, 30),
                         limits = c(0, 30), oob = squish) +
    scale_size(range = c(1, 6), breaks = c(5, 50, 100, 200), limits = c(5, 200)) +
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
