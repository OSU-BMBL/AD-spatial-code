library(Seurat)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(dplyr)


###################################

layer_1 <- F
median.mean <- T

###################################

setwd("/users/PAS1475/liuzeyi/guoqi/output")
load("sample6_clusters.RData")
my.merge.object.rm.noise <- sample6.combined_f
rm(sample6.combined_f)
DEGs.2456.list <- read.csv("/users/PAS1475/liuzeyi/guoqi/output/conservedmarker/specific_marker_final/specific_markers_2456w.csv", 
                           row.names = 1, header = T)

if (layer_1) {
  DEGs.1.list <- read.csv("/users/PAS1475/liuzeyi/guoqi/output/conservedmarker/specific_marker_final/specific_markers_1.csv", 
                          row.names = 1, header = T)
}

length(heatmap.all.genes)

if (layer_1) {
  gene_anno <- rbind(DEGs.1.list, DEGs.2456.list) %>% arrange(Layer) %>% select(Layer)
} else {
  gene_anno <-DEGs.2456.list %>% arrange(Layer) %>% select(Layer)
}

heatmap.all.genes <- rownames(gene_anno)
dim(gene_anno)
head(gene_anno)
tail(gene_anno)


# create median
mat_in <- my.merge.object.rm.noise@assays$SCT@data[heatmap.all.genes,]
dim(mat_in)
mat_in[1:5, 1:5]
meta_sub <- my.merge.object.rm.noise@meta.data
dim(meta_sub)
meta_sub <- meta_sub[order(meta_sub$Layer, meta_sub$patient, meta_sub$patientID),]
head(meta_sub)
tail(meta_sub)
mat_in_sub <- as.matrix(mat_in)
mat_in_sub <- mat_in_sub[, rownames(meta_sub)]
dim(mat_in_sub)
mat_in_sub[1:5, 1:5]

mat_in_merge <- as.data.frame(t(mat_in_sub))
head(mat_in_merge)
identical(rownames(meta_sub), rownames(mat_in_merge))
layer_vec <- meta_sub$Layer
# layer_vec <- layer_vec[layer_vec %in% unique(gene_anno$Layer)]
layer_vec <- gsub(" ", "_", layer_vec)
unique(layer_vec)

Pseudo_bulk <- c()
for(i in 1:length(unique(layer_vec))) {
  print(unique(layer_vec)[i])
  tmp.mat <- mat_in_merge[layer_vec == unique(layer_vec)[i], ]
  
  if (median.mean) {
    tmp.vec <- apply(tmp.mat, 2, median)
  } else {
    tmp.vec <- apply(tmp.mat, 2, mean)
  }
  
  Pseudo_bulk <- rbind(Pseudo_bulk, tmp.vec)
  rownames(Pseudo_bulk)[i] <- unique(layer_vec)[i]
}
dim(Pseudo_bulk)
head(Pseudo_bulk)
sum(Pseudo_bulk[1:5, 1:5])

my_colour = list(
  gene_cat = c( Layer_1 = brewer.pal(n = 7, name = "Dark2")[1],
                Layer_2 = brewer.pal(n = 7, name = "Dark2")[2],
                Layer_3 = brewer.pal(n = 7, name = "Dark2")[3],
                Layer_4 = brewer.pal(n = 7, name = "Dark2")[4],
                Layer_5 = brewer.pal(n = 7, name = "Dark2")[5],
                Layer_6 = brewer.pal(n = 7, name = "Dark2")[6],
                White_Matter = brewer.pal(n = 7, name = "Dark2")[7]))

gene_cat_val <- as.numeric(table(gene_anno))

check.names <- unique(gene_anno$Layer) %>% gsub(" ", "_", .)
Pseudo_bulk <- Pseudo_bulk[check.names,]
p <- pheatmap(t(Pseudo_bulk),
              gaps_row = cumsum(gene_cat_val),
              cluster_rows = F,
              cluster_cols = F,
              scale = "row",
              border_color = NA,
              annotation_colors = my_colour,
              fontsize = 20, 
              annotation_row = gene_anno,
              labels_col = rownames(Pseudo_bulk),
              color = colorRampPalette(c("blue","white","red"))(100))

if (median.mean) {
  if (layer_1) {
    ggsave(plot = p, filename = "/users/PAS1475/liuzeyi/guoqi/output/picture/heatmap_median_layers_all.tiff", 
           device = "tiff",dpi = 300, width = 15, height = 30, units = "in")
    # ggsave(plot = p, filename = "/users/PAS1475/liuzeyi/guoqi/output/picture/heatmap_median_layers_all.png", 
    #        device = "png",dpi = 300, width = 15, height = 30, units = "in")
  } else {
    ggsave(plot = p, filename = "/users/PAS1475/liuzeyi/guoqi/output/picture/heatmap_median_layers_no_1_3.tiff", 
           device = "tiff",dpi = 300, width = 15, height = 30, units = "in")
    # ggsave(plot = p, filename = "/users/PAS1475/liuzeyi/guoqi/output/picture/heatmap_median_layers_no_1_3.png", 
    #        device = "png",dpi = 300, width = 15, height = 30, units = "in")
    }
} else {
  if (layer_1) {
    ggsave(plot = p, filename = "/users/PAS1475/liuzeyi/guoqi/output/picture/heatmap_mean_layers_all.tiff", 
           device = "tiff",dpi = 300, width = 15, height = 30, units = "in")
    # ggsave(plot = p, filename = "/users/PAS1475/liuzeyi/guoqi/output/picture/heatmap_mean_layers_all.png", 
    #        device = "png",dpi = 300, width = 15, height = 30, units = "in")
  } else {
    ggsave(plot = p, filename = "/users/PAS1475/liuzeyi/guoqi/output/picture/heatmap_mean_layers_no_1_3.tiff", 
           device = "tiff",dpi = 300, width = 15, height = 30, units = "in")
    # ggsave(plot = p, filename = "/users/PAS1475/liuzeyi/guoqi/output/picture/heatmap_mean_layers_no_1_3.png", 
    #        device = "png",dpi = 300, width = 15, height = 30, units = "in")
  }
}