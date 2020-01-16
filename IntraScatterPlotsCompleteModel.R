library(data.table)
library(ggplot2)
library(grid)
library(scales)
library(plyr)
library(ggthemes)
library(RColorBrewer)

setwd("/media/ddisk/transpipeline-data/")

types <- c("utero", "lung")


getScatterByChr <- function(types) {
  chrs <- as.character(c(1:22, "X"))
  chromosomes.pal <- c("#D909D1", "#0492EE", "#5DA0CB", "#106F35", "#5BD2AE", "#199F41", 
                       "#FE0F43", "#00FFCC", "#F495C5", "#E1BF5D", "#5F166F", "#088ACA",
                       "#41CFE0", "#0F0A71", "#FFFF99", "#B06645", "#800092", "#B925AE",
                       "#B1B719", "#CB97E8", "#130B9E", "#E12B29", "#79A5B9")
  
  names(chromosomes.pal) <- c("22","11","12","13","14","15","16","17","18","19","1" ,"2" ,"3" ,"4" ,"5" ,
                              "6" ,"7" ,"X" ,"8" ,"9" ,"20","10","21")
  chromosomes.pal <- chromosomes.pal[chrs]
  conds <- c("healthy", "cancer")
  colors <- c("#C7C7C7", "#B37700")
  labels <- c("Healthy", "Cancer")
  
  for(type in types) {
    datos <- fread(file = paste0(type, "/intra/intra-fixed-bin-size/100-all.txt"), 
                  sep = "\t", nThread = 3,  header = T)
    names(datos) <- c("bin", "distance", "MI", "cond", "chr")
    datos$cond <- factor(datos$cond, levels = conds, labels = labels)
    datos$chr <- factor(datos$chr, levels = chrs)
    
    ## Todos los cromosomas
    p <- ggplot(
      data = datos
    ) +
      geom_point( aes(
        x = distance/1e6,
        y = MI
      ), color = "gray38", size = 0.1) +
      facet_grid(chr ~ cond, scales = "free_y") +
      theme_few(base_size = 20) +
      theme(
        legend.position = "none",
        strip.text.x = element_text(size = 30),
        strip.text.y  = element_text(size = 24),
        axis.title=element_text(size=30,face="bold")
      ) + xlab("Distance (Mbp)") + ylab("Mutual Information") 
    
    g <- ggplot_gtable(ggplot_build(p))
    strips <- which(grepl('strip-', g$layout$name))
    for (i in 1:2) {
      k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
      l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
      g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- colors[i]
    }
    png(paste0(type, "/figures/mi-distance-all-chromosomes.png"), width = 1000, height = 4500)
    plot(g)
    dev.off()
    
  }
}
getScatterByChr(types)
# 
# 
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
# ### Figura 3
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

