log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)

cat("Loading data\n")
color_pal <- c("#e3a098", "#a32e27")
fitted_data <- read_tsv(snakemake@input[[1]])
TISSUE <- snakemake@params[["tissue"]]
substring(TISSUE, 1, 1) <- toupper(substring(TISSUE, 1, 1))

fitted_data <- fitted_data %>% mutate(
  min_plot =  pmax(mean_fitted - mi_sd, 0), 
  max_plot = pmin(mean_fitted + mi_sd, 0.1))

cat("Building plot\n")
g <- ggplot(fitted_data) + 
  geom_line(aes(x = dist_mean/1e6, y = mean_fitted, color=cond)) +
  geom_ribbon(aes(x = dist_mean/1e6, ymin = min_plot, 
                  ymax = max_plot, fill = cond), 
              alpha = .2) +
  facet_wrap(~cond, nrow = 1) + 
  xlab("Distance (Mbp)") + 
  ylab("Mutual Information") + 
  theme_few(base_size = 25) +
  ylim(c(0,0.1)) +
  scale_fill_manual(values = color_pal) + 
  expand_limits(x = 0, y = 0) +
  scale_color_manual(values = color_pal) + 
  theme(legend.position = "none", strip.text.x = element_text(size = 30)) +
  ggtitle(TISSUE)

cat("Saving plot\n")
png(snakemake@output[[1]], width =800, height = 450)
print(g)
dev.off()  