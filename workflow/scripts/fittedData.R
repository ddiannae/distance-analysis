library(tidyverse)
library(broom)
library(ggthemes)

files <- list(snakemake@input[["normal"]], snakemake@input[["cancer"]])

mi_data <- lapply(files, function(file) {
  read_tsv(file)
})
mi_data <- bind_rows(mi_data) 
mi_data$cond <- factor(mi_data$cond, levels = c("normal", "cancer"), labels = c("Normal", "Cancer"))

fitted_data <- mi_data %>% group_by(cond) %>% nest() %>%
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
  inner_join(mi_data %>% select(bin, cond, mi_sd, dist_mean))

write_tsv(fitted_data, snakemake@output[[1]])