library(tidyverse)
library(broom)
library(ggthemes)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  DIST_DIR <- args[1]
  FIG_DIR <- args[2]
  DIST_BIN <- as.integer(args[3])
  SIZE_BIN <- as.integer(args[4])
  TISSUE <- args[5]
}

getFittedData <- function(filesuff) {
  conds <- c("normal", "cancer")
  mi_data <- lapply(conds, function(cond) {
    read_tsv(paste0(DIST_DIR, cond, "-", filesuff))
  })
  mi_data <- bind_rows(mi_data) 
  mi_data$cond <- factor(mi_data$cond, levels = c("normal", "cancer"))
  levels(mi_data$cond) <- c("Healthy", "Cancer")
  
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
           max = maxf_mi_max, max_fitted = maxf_.fitted, max_resid = maxf_.resid) %>% ungroup()
  sds <- mi_data %>% select(bin, cond, mi_sd) 
  fitted_data <- fitted_data %>% inner_join(sds)
  write_tsv(fitted_data, paste0(DIST_DIR, "fitted-", filesuff))
  return(fitted_data)
  
}

size_data <- getFittedData(paste0("fixed-size-all-", SIZE_BIN, ".tsv"))
distance_data <- getFittedData(paste0("fixed-distance-all-", DIST_BIN, ".tsv"))
color_pal <- c("#e3a098", "#a32e27")
TISSUE <- paste0(toupper(substring(TISSUE, 1, 1)), substring(TISSUE, 2))

g <- ggplot(size_data) + 
  geom_line(aes(x = bin, y = mean_fitted, color=cond)) +
  geom_ribbon(aes(x = bin, ymin = mean_fitted - mi_sd, ymax = mean_fitted + mi_sd, fill = cond), alpha = .2) +
  facet_wrap(~cond, nrow = 1) + 
  xlab("Bin number") + 
  ylab("Mutual Information") + 
  theme_few(base_size = 25) +
  scale_fill_manual(values = color_pal) + 
  scale_color_manual(values = color_pal) + 
  theme(legend.position = "none", strip.text.x = element_text(size = 30)) +
  ggtitle(TISSUE)
  

png(paste0(FIG_DIR, "all-bins-fixed-size.png"), width = 1200, height = 600)
print(g)
dev.off()  


g <- ggplot(distance_data) + 
  geom_line(aes(x = bin, y = mean_fitted, color=cond)) +
  geom_ribbon(aes(x = bin, ymin = mean_fitted - mi_sd, ymax = mean_fitted + mi_sd, fill = cond), alpha = .2) +
  facet_wrap(~cond, nrow = 1) + 
  xlab("Distance (Mbp)") + 
  ylab("Mutual Information") + 
  theme_few(base_size = 25) +
  scale_fill_manual(values = color_pal) + 
  scale_color_manual(values = color_pal) + 
  theme(legend.position = "none", strip.text.x = element_text(size = 30)) +
  ggtitle(TISSUE)


png(paste0(FIG_DIR, "all-bins-fixed-distance.png"), width = 1200, height = 600)
print(g)
dev.off()  
