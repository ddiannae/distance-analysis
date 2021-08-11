library(tidyverse)
library(broom)
library(ggthemes)

files <- list(snakemake@input[["normal"]], snakemake@input[["cancer"]])
fitted_data <- read_tsv(snakemake@input[["fitted"]], col_types = cols(chr = col_character()))
TISSUE <- snakemake@params[["tissue"]]

mi_data <- lapply(files, function(file) {
  read_tsv(file, col_types = cols(chr = col_character()))
})

mi_data <- bind_rows(mi_data) 
mi_data$cond <- factor(mi_data$cond, levels = c("normal", "cancer"), labels = c("Healthy", "Cancer"))
mi_data$chr <- factor(mi_data$chr, levels = as.character(c(1:22, "X", "Y")))


chrs <- as.character(c(1:22, "X", "Y"))
chr_pal <- c("#D909D1", "#0492EE", "#D2656C", "#106F35", "#5BD2AE", "#199F41", 
             "#FE0F43", "#00FFCC", "#F495C5", "#E1BF5D", "#5F166F", "#088ACA",
             "#41CFE0", "#0F0A71", "#FFFF99", "#B06645", "#651532", "#B925AE",
             "#86C91F", "#CB97E8", "#130B9E", "#EF782B", "#79A5B9", "#F7AA7C")

names(chr_pal) <- c("22","11","12","13","14","15","16","17","18","19","1" ,"2" ,"3" ,"4" ,"5" ,
                    "6" ,"7" ,"X" ,"8" ,"9" ,"20","10","21", "Y")
chr_pal <- chr_pal[chrs]

## Todos los cromosomas
g <- ggplot() +
  geom_point(data = mi_data,
             mapping = aes(x = dist_mean/1e6, y = mi_mean), color = "gray65", size = 0.1) +
  geom_line(data = fitted_data,
            mapping = aes(x = dist_mean/1e6, y = mean_fitted, color=chr)) +
  facet_grid(chr ~ cond, scales = "free_y") +
  theme_few(base_size = 20) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 30),
    strip.text.y  = element_text(size = 24),
    axis.title = element_text(size = 30),
    title = element_text(size = 35),
  ) + xlab("Distance (Mbp)") + ylab("Mutual Information") + 
  scale_color_manual(values = chr_pal, name = "Chromosome") +
  ggtitle(TISSUE)

png(snakemake@output[[1]], width = 1000, height = 4500)
plot(p)
dev.off()
