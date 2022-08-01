## #############################################################
## This file fits a loess curve to the mean, min, and max values
## per bin per chromosome to improve visualization 
## Its input comes from the binStats.R script
################################################################

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(broom)

cat("Reading files\n")
files <- list(snakemake@input[["normal"]], snakemake@input[["cancer"]])

mi_data <- lapply(files, function(file) {
  read_tsv(file, col_types = cols(chr = col_character()))
})
mi_data <- bind_rows(mi_data)
mi_data$cond <- factor(mi_data$cond, levels = c("normal", "cancer"), labels = c("Normal", "Cancer"))
mi_data$chr <- factor(mi_data$chr, levels = as.character(c(1:22, "X", "Y")))

cat("Fitting data\n")
mi_data <- mi_data %>% filter(chr != "Y")
fitted_data <- mi_data %>% group_by(cond, chr) %>% nest() %>%
  mutate(fit_mean = map(data, ~ loess(mi_mean ~ bin, ., span = 0.2)),
         fit_max = map(data, ~ loess(mi_max ~ bin, ., span = 0.2)),
         fit_min = map(data, ~ loess(mi_min ~ bin, ., span = 0.2)),
         meanf = map(fit_mean, augment),
         maxf = map(fit_max, augment),
         minf = map(fit_min, augment)) %>%
  select(-data, -fit_mean, -fit_max,-fit_min ) %>%
  unnest(c(meanf, maxf, minf), names_sep = "_") %>% select(-minf_bin, -maxf_bin) %>%
  rename(bin = meanf_bin,  mean = meanf_mi_mean, mean_fitted = meanf_.fitted, mean_resid = meanf_.resid,
         min = minf_mi_min, min_fitted = minf_.fitted, min_resid = minf_.resid,
         max = maxf_mi_max, max_fitted = maxf_.fitted, max_resid = maxf_.resid) %>%
  ungroup() %>%
  inner_join(mi_data %>% select(bin, cond, chr, mi_sd, dist_mean))

cat("Writing data\n")
write_tsv(fitted_data, snakemake@output[[1]])
