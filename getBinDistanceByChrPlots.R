library(tidyverse)
library(broom)
library(ggthemes)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  DIST_DIR <- args[1]
  FIG_DIR <- args[2]
  DIST_BIN <- as.integer(args[3])
  SIZE_BIN <- as.integer(args[4])
  TISSUE <- args[5]
}

getPointData <- function(filesuff) {
  conds <- c("normal", "cancer")
  mi_data <- lapply(conds, function(cond) {
    read_tsv(paste0(DIST_DIR, cond, "-", filesuff), col_types = cols(chr = col_character()))
  })
  mi_data <- bind_rows(mi_data) 
  mi_data$cond <- factor(mi_data$cond, levels = c("normal", "cancer"), labels = c("Healthy", "Cancer"))
  mi_data$chr <- factor(mi_data$chr, levels = as.character(c(1:22, "X", "Y")))
  
  return(mi_data)
}

getFittedData <- function(mi_data) {
 
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
  
  return(fitted_data)
}

getChrMIDistancePlot <- function(point_data, fitted_data) {
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
    geom_point(data = point_data,
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
  
  return(g)
  
}

size_data <- getPointData(paste0("fixed-size-bychr-", SIZE_BIN, ".tsv"))
distance_data <- getPointData(paste0("fixed-distance-bychr-", DIST_BIN, ".tsv"))
fitted_size_data <- getFittedData(size_data)
fitted_distance_data <- getFittedData(distance_data)
write_tsv(fitted_size_data, paste0(DIST_DIR, "fitted-size-bychr-", SIZE_BIN, ".tsv"))
write_tsv(fitted_distance_data, paste0(DIST_DIR, "fitted-distance-bychr-", DIST_BIN, ".tsv"))


p <- getChrMIDistancePlot(size_data, fitted_size_data)
png(paste0(FIG_DIR,  "bychr-bins-fixed-size-", SIZE_BIN, ".png"), 
    width = 1000, height = 4500)
plot(p)
dev.off()

p <- getChrMIDistancePlot(distance_data, fitted_distance_data)
png(paste0(FIG_DIR,  "bychr-bins-fixed-distance-", DIST_BIN, ".png"), 
    width = 1000, height = 4500)
plot(p)
dev.off()
    
# FROM BREAST SUBTYPE
### Chr8 y Chr17
# datos_chr8 <- datos[datos$Chr == "8", ]
# datos_chr8$Type <- "Chr 8"
# 
# datos_chr17 <- datos[datos$Chr == "17", ]
# datos_chr17$Type <- "Chr 17"
# datos_all <- rbind(datos_chr8, datos_chr17)
# datos_all$Type <- factor(datos_all$Type, levels = c("Chr 8", "Chr 17"))

# #### Solo Chromosoma 17 o solo chromosoma 8
# datos_chr <- datos_chr8
# p <- ggplot(
#   data = datos_chr,
#   aes(
#     x = Distance/1e6,
#     y = MI
#   )
# ) +
#   geom_point(
#     color = "gray38"
#   ) +
#   facet_grid(
#     Type ~ Subtype
#   ) +   
#   geom_line(
#     data = datos_chr,
#     aes(
#       x = datos_chr$Distance/1e6,
#       y = datos_chr$fitted
#     ), size = 1.2, colour = "#B1B719"
#   ) +
#   theme_few(base_size = 16) +
#   theme(
#     legend.position = "none",
#     strip.text.x = element_text(size = 16),
#     strip.text.y = element_text(size = 16),
#     axis.title=element_text(size=16)
#   ) + xlab("Distance (Million bp)") +
#   ylab("Mutual Information") +
#   scale_x_continuous(breaks = seq(0, 300, by = 100))
# 
# g <- ggplot_gtable(ggplot_build(p))
# strips <- which(grepl('strip-', g$layout$name))
# for (i in 1:5) {
#   k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
#     g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- subtypes.pal[i]
# }
# plot(g)
# 
# 
# datos.bh <- datos[datos$Subtype == "Basal" | datos$Subtype == "Healthy", ]
# 
# p2 <- ggplot(
#   data = datos.bh,
#   aes(
#     x = Distance/1e6,
#     y = MI, 
#     color=Chr
#   ) 
# ) + geom_point(size = 1,  alpha = 0.5, stroke = 0) + 
#   scale_color_manual(values = chromosomes.pal, name = "Chromosome") + 
#   theme_few(base_size = 30) +
#   facet_wrap(
#     datos.bh$Subtype, nrow = 1
#   ) +
#   theme(
#     legend.position = "none", 
#     legend.spacing.y = unit(1.0, "cm"),
#     legend.spacing.x = unit(1.0, 'cm'),
#     strip.text.x = element_text(size = 30),
#     legend.title = element_text(size =30),
#     legend.text = element_text(size = 30),
#     axis.title=element_text(size=30)
#   ) +
#   guides(colour = guide_legend(override.aes = list(size=5), keyheight=1.0, 
#          default.unit = "cm")) +
#   xlab("Distance (Mbp)") +
#   ylab("Mutual Information") 
# p2
# pzoom <- p2
# pzoom <- pzoom +
#   theme_few(base_size = 30) +
#   theme(legend.position = "none", 
#         strip.text.x = element_blank(),
#         axis.title.x= element_blank(),
#         axis.title.y= element_blank()) +
#   scale_y_continuous(breaks = seq(0.02, 0.07, by = .025), limits = c(0.020, 0.07))+
#   scale_x_continuous(breaks = seq(0, 20, by = 10), limits = c(0, 20))
# vp <- viewport(
#   width = 0.48,
#   height = 0.38,
#   x = 0.48,
#   y = 0.58,
#   just = c("left", "bottom")
# )
# p2; print(pzoom, vp = vp)

