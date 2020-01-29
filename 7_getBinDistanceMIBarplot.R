## Script to get Distance vs MI barplot
## 5Mb bins are used for the distance and the mean MI 
## value for each bin is plotted. 
## It requires the cond-all-distance-mi.txt files from 
## the 1_getIntraInteractions script.
## It is usally run in the server
library(data.table)
library(ggplot2)
library(ggthemes)

types <- c("utero")

for(type in types) {
  
  setwd(paste0("/labs/csbig/regulaciontrans/", type))
  
  conds <- c("healthy", "cancer")
  
  all.mis <- lapply(conds, function(cond) {
    MI.data <-fread(file = paste0("intra/", cond, "-all-distance-mi.txt"), 
                    header = T, sep = "\t")
    MI.data$distMb <- MI.data$distance / 1000000
    
    mis <- lapply(seq(0, 245, 5), function(mbs){
      mean(MI.data[MI.data$distMb >= mbs & MI.data$distMb < mbs + 5, mi])
    })
    
    mis <- data.frame(distance = seq(0, 245, 5), mi = unlist(mis))
    mis$cond <- cond
    return(mis)
  })
  
  all.mis <- rbindlist(all.mis)
  all.mis$cond <- factor(all.mis$cond, levels = c("healthy", "cancer"))
  levels(all.mis$cond) <- c("Healthy", "Cancer")
  subtypes.pal <- c("#e3a098", "#a32e27")
  
  fwrite(all.mis, file =  paste("intra/all_mi-by_distanceMb.txt", sep = ""), 
         row.names = F, col.names = T, sep = "\t")  
  
  
  g <- ggplot(all.mis, aes(x = distance, y = mi, fill = cond)) + 
    geom_bar(stat = "identity") + 
    facet_wrap(~cond, nrow = 1) + 
    xlab("Distance (Mbp)") + 
    ylab("Mutual Information") + 
    theme_few(base_size = 25) +
    scale_fill_manual(values = subtypes.pal) + 
    scale_x_continuous(breaks = seq(0,250,by=100)) +
    theme(legend.position = "none", strip.text.x = element_text(size = 30))
  
  png("figures/mi-barplot-distance.png", width = 1000, height = 500)
  plot(g)
  dev.off()  
}
