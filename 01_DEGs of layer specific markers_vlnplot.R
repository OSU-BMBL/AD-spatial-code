library(Seurat)
library(patchwork)
library(ggplot2)

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

load("~/guoqi/output/sample6_clusters.RData")
table(sample6.combined_f$patientID)
table(sample6.combined_f$Layer)

sample6.combined_f$Layer_factor <- sample6.combined_f$Layer
sample6.combined_f$Layer_factor <- factor(sample6.combined_f$Layer_factor)
Idents(sample6.combined_f) <- sample6.combined_f$Layer_factor
DefaultAssay(sample6.combined_f) <- "SCT"

conservMarker_df <- read.csv("~/guoqi/output/allconserved_markers.csv",header = T)
conservMarker_df <- cbind.data.frame(conservMarker_df$X,conservMarker_df$minimump_p_val, conservMarker_df$layer)
conservMarker_df <- conservMarker_df[order(conservMarker_df$`conservMarker_df$layer`),]
# plot stack violin
p <- StackedVlnPlot(obj = sample6.combined_f, features = conservMarker_df$`conservMarker_df$X`)
ggsave(plot = p , filename = "~/guoqi/output/stackedvlnplot_for_markers.tiff",  
       device = "tiff",dpi = 150, width = 6, height = 40, units = "in") 


