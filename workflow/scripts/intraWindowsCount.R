log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(dplyr)
library(tidyr)
library(parallel)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)

load(snakemake@params[["annot"]])
MCCORES <- snakemake@threads[[1]]

TISSUE <- snakemake@params[["tissue"]]
substring(TISSUE, 1, 1) <- toupper(substring(TISSUE, 1, 1))

cat("Reading files\n")
chrs <- c(as.character(1:22), "X", "Y")
MImatrix <- vroom::vroom(snakemake@input[["mi_matrix"]])
cat("MI matrix with ", nrow(MImatrix), " rows and ", ncol(MImatrix), " columns loaded \n")
genes <- colnames(MImatrix)
MImatrix$source <- genes

cat("Annotating interactions\n")
mi_vals <- MImatrix %>% pivot_longer(cols = starts_with("ENSG"), 
                                     names_to = "target",
                                     values_to = "mi",
                                     values_drop_na = TRUE) %>% 
  arrange(desc(mi)) %>% 
  left_join(annot %>% dplyr::select(gene_id, chr), 
            by = c("source" = "gene_id")) %>% 
  rename("source_chr" = "chr") %>%
  left_join(annot %>% dplyr::select(gene_id, chr), 
            by = c("target" = "gene_id")) %>% 
  rename("target_chr" = "chr") %>%
  mutate(interaction = ifelse(source_chr == target_chr, "intra", "inter"),
         nrow = row_number())

intra_mi  <- mi_vals %>% filter(interaction == "intra")
rm(MImatrix)
annot <- annot %>% filter(gene_id %in% union(intra_mi$source, intra_mi$target))

k <- 1000
tops <- c(10*k, 100*k, 500*k, 1000*k, 5000*k, 10000*k, nrow(mi_vals))

# if(nrow(mi_vals) > 5000*k) {
#   tops <- c(tops, 5000*k, nrow(mi_vals))
# } else {
#   tops <- c(tops, nrow(mi_vals))
# }

cat("Getting intra interaction counts\n")
all_fractions <- lapply(1:length(tops), function(i) {
  cat("Working with ", tops[i], "top interactions\n")
  mi_top <- mi_vals %>% filter(nrow > max(tops[i-1],0) & nrow <= tops[i])
  chr_win_fractions <- parallel::mclapply(X = chrs, mc.cores = MCCORES, FUN = function(ch) {
  #chr_win_fractions <- lapply(chrs, FUN = function(ch) {
    cat("Working with chromosome", ch, "\n")
    mi_top_chr <- mi_top %>% filter(source_chr == ch & target_chr == ch)
    annot_chr <- annot %>% filter(chr == ch)
    if(nrow(annot_chr) > 0) {
      wins <- seq(from = 0, to = max(annot_chr$start), by=800000)
      win_fractions <- lapply(wins, function(win) {
        annot_chr_win <- annot_chr %>% filter(start >= win & start <= (win+1e6))
        intra_win <- mi_top_chr %>% filter(source %in% annot_chr_win$gene_id &
                                                   target %in% annot_chr_win$gene_id)
        
        nwin <- nrow(annot_chr_win)
        if(nwin > 1) {
          ntotal <- (nwin*(nwin-1))/2
          return(list(start = win, end = win+1e6, inter_total= ntotal, inter = nrow(intra_win),
                      fraction = nrow(intra_win)/ntotal))
        }
      })
      win_fractions <- bind_rows(win_fractions)
      win_fractions$chr <- ch
      return(win_fractions)
    }
  })
  chr_win_fractions <- bind_rows(chr_win_fractions)
  chr_win_fractions$top <- tops[i]
  return(chr_win_fractions)
})
cat("Done counting interactions\n")
all_fractions <- bind_rows(all_fractions)
all_fractions$top <- factor(all_fractions$top, levels = sort(tops, decreasing = T))

cat("Writing file\n")
all_fractions %>% write_tsv(file=snakemake@output[["tsv"]])

all_fractions <- all_fractions %>% mutate(chr = factor(chr, levels=chrs), 
                                          top = factor(top, levels=rev(tops)))

cat("Building plot\n")
g <- ggplot(all_fractions, aes(x=start, y=fraction, fill=top)) +
  geom_bar(position="stack", stat="identity") + 
  theme_base(base_size = 20) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 5e7, 1e8, 1.5e8, 2e8, 2.5e8)) +
  facet_grid(rows = vars(chr)) +
  theme(panel.spacing.y = unit(1, "lines")) +
  ylab("Fraction of total interacions in region") +
  xlab("") +
  scale_fill_brewer(palette = "Pastel1", name = "Top Interactions", direction=-1) +
  ggtitle(TISSUE) 

cat("Saving plot\n")
png(snakemake@output[["plot"]], width = 1500, height = 1500)
plot(g)
dev.off()

