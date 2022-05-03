library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(devtools)
devtools::install_github("jaredhuling/jcolors")
library(jcolors)
jcolors('default')
dyn.load("/apps/hdf5-serial/gnu/9.1/1.12.0/lib/libhdf5_hl.so.200")
library(hdf5r)
options(future.globals.maxSize = 8000 * 1024 ^ 2)
#-------Module for data processing,including data loading,normalization,noise annotation(keeping noise
#for spatial plot), integrating.
#---data loading(h5 and annotated layer information) and checking
setwd("/fs/ess/PCON0022/liyang/Shane/raw_data/")
#2-5
wd <- "./2-5/outs/"
object.2_5 <-
  Load10X_Spatial(data.dir = wd, filename = "filtered_feature_bc_matrix.h5")
my.2_5.label <-
  read.csv(
    "Layer_annotation/2-5 layer manuallabel new.csv",
    header = T,
    row.names = 1
  )
identical(colnames(object.2_5), rownames(my.2_5.label))
#18-64
wd <- "./18-64/outs/"
object.18_64 <-
  Load10X_Spatial(data.dir = wd, filename = "filtered_feature_bc_matrix.h5")
my.18_64.label <-
  read.csv(
    "Layer_annotation/18-64 manuallabel new.csv",
    header = T,
    row.names = 1
  )
identical(colnames(object.18_64), rownames(my.18_64.label))
#1-1
wd <- "./1-1/outs/"
object.1_1 <-
  Load10X_Spatial(data.dir = wd, filename = "filtered_feature_bc_matrix.h5")
my.1_1.label <-
  read.csv(
    "Layer_annotation/1-1 Layer manuallabel new.csv",
    header = T,
    row.names = 1
  )
identical(colnames(object.1_1), rownames(my.1_1.label))
#2-8
wd <- "./2-8/outs/"
object.2_8 <-
  Load10X_Spatial(data.dir = wd, filename = "filtered_feature_bc_matrix.h5")
my.2_8.label <-
  read.csv(
    "Layer_annotation/2-8 layer manuallabel new.csv",
    header = T,
    row.names = 1
  )
identical(colnames(object.2_8), rownames(my.2_8.label))
#T4857
wd <- "./T4857/outs/"
object_T4857 <-
  Load10X_Spatial(data.dir = wd, filename = "filtered_feature_bc_matrix.h5")
my.T4857.label <-
  read.csv(
    "Layer_annotation/T4857 layer manuallabel new.csv",
    header = T,
    row.names = 1
  )
identical(colnames(object_T4857), rownames(my.T4857.label))
#2-3
wd <- "./2-3/outs/"
object.2_3 <-
  Load10X_Spatial(data.dir = wd, filename = "filtered_feature_bc_matrix.h5")
my.2_3.label <-
  read.csv(
    "Layer_annotation/2-3 layer manuallabel new.csv",
    header = T,
    row.names = 1
  )
identical(colnames(object.2_3), rownames(my.2_3.label))

#---normalization
object.2_5 <-
  SCTransform(object.2_5, assay = "Spatial", verbose = FALSE)
object.18_64 <-
  SCTransform(object.18_64, assay = "Spatial", verbose = FALSE)
object.1_1 <-
  SCTransform(object.1_1, assay = "Spatial", verbose = FALSE)
object.2_8 <-
  SCTransform(object.2_8, assay = "Spatial", verbose = FALSE)
object_T4857 <-
  SCTransform(object_T4857, assay = "Spatial", verbose = FALSE)
object.2_3 <-
  SCTransform(object.2_3, assay = "Spatial", verbose = FALSE)

#---annotation
#category(control,AD)
object.2_5$category <- "control"
object.18_64$category <- "control"
object.1_1$category <- "control"
object.2_8$category <- "AD"
object_T4857$category <- "AD"
object.2_3$category <- "AD"
#stage(control,middle)
object.2_5$stage <- "control"
object.18_64$stage <- "control"
object.1_1$stage <- "control"
object.2_8$stage <- "Mid-AD"
object_T4857$stage <- "Mid-AD"
object.2_3$stage <- "Mid-AD"
#id
object.2_5$patientID <- "2_5"
object.18_64$patientID <- "18_64"
object.1_1$patientID <- "1_1"
object.2_8$patientID <- "2_8"
object_T4857$patientID <- "T4857"
object.2_3$patientID <- "2_3"
#noise
colnames(my.2_5.label) <- c("Layer")
colnames(my.18_64.label) <- c("Layer")
colnames(my.1_1.label) <- c("Layer")
colnames(my.2_8.label) <- c("Layer")
colnames(my.T4857.label) <- c("Layer")
colnames(my.2_3.label) <- c("Layer")
my.2_5.label$Layer[my.2_5.label$Layer == ""]  <- "Noise"
my.18_64.label$Layer[my.18_64.label$Layer == ""]  <- "Noise"
my.1_1.label$Layer[my.1_1.label$Layer == ""]  <- "Noise"
my.2_8.label$Layer[my.2_8.label$Layer == ""]  <- "Noise"
my.T4857.label$Layer[my.T4857.label$Layer == ""]  <- "Noise"
my.2_3.label$Layer[my.2_3.label$Layer == ""]  <- "Noise"
#adding layer information to h5
object.2_5 <- AddMetaData(object.2_5, metadata = my.2_5.label)
object.18_64 <- AddMetaData(object.18_64, metadata = my.18_64.label)
object.1_1 <- AddMetaData(object.1_1, metadata = my.1_1.label)
object.2_8 <- AddMetaData(object.2_8, metadata = my.2_8.label)
object_T4857 <- AddMetaData(object_T4857, metadata = my.T4857.label)
object.2_3 <- AddMetaData(object.2_3, metadata = my.2_3.label)
#adding layer information to Idents of h5
Idents(object.2_5) <- object.2_5$Layer
Idents(object.18_64) <- object.18_64$Layer
Idents(object.1_1) <- object.1_1$Layer
Idents(object.2_8) <- object.2_8$Layer
Idents(object_T4857) <- object_T4857$Layer
Idents(object.2_3) <- object.2_3$Layer
#save object including noise
save(object.2_8, file = "/users/PAS1475/liuzeyi/guoqi/output/object.2_8_n.RData")
save(object_T4857, file = "/users/PAS1475/liuzeyi/guoqi/output/object_T4857_n.RData")
save(object.2_3, file = "/users/PAS1475/liuzeyi/guoqi/output/object.2_3_n.RData")
save(object.2_5, file = "/users/PAS1475/liuzeyi/guoqi/output/object.2_5_n.RData")
save(object.18_64, file = "/users/PAS1475/liuzeyi/guoqi/output/object.18_64_n.RData")
save(object.1_1, file = "/users/PAS1475/liuzeyi/guoqi/output/object.1_1_n.RData")
setwd("/users/PAS1475/liuzeyi/guoqi/output/")
#umap
load("/users/PAS1475/liuzeyi/guoqi/output/object.2_8_n.RData")
load("/users/PAS1475/liuzeyi/guoqi/output/object_T4857_n.RData")
load("/users/PAS1475/liuzeyi/guoqi/output/object.2_3_n.RData")
load("/users/PAS1475/liuzeyi/guoqi/output/object.2_5_n.RData")
load("/users/PAS1475/liuzeyi/guoqi/output/object.18_64_n.RData")
load("/users/PAS1475/liuzeyi/guoqi/output/object.1_1_n.RData")
#------Integrating 6 samples
my.all.object.list <- list(
  object.2_5 = object.2_5,
  object.18_64 = object.18_64,
  object.2_8 = object.2_8,
  object_T4857 = object_T4857,
  object.1_1 = object.1_1,
  object.2_3 = object.2_3
)
#anchor feature
features <-
  SelectIntegrationFeatures(object.list = my.all.object.list, nfeatures = 3000)
#integrated data
brain.list <-
  PrepSCTIntegration(
    object.list = my.all.object.list,
    anchor.features = features,
    verbose = TRUE,
  )
anchors <-
  FindIntegrationAnchors(
    object.list = brain.list,
    anchor.features = features,
    normalization.method = "SCT",
    verbose = FALSE
  )
sample6.combined <-
  IntegrateData(
    anchorset = anchors,
    normalization.method = "SCT",
    verbose = FALSE
  )
#integrated analysis
DefaultAssay(sample6.combined) <- "integrated"
sample6.combined <- ScaleData(sample6.combined, verbose = FALSE)
sample6.combined <-
  RunPCA(sample6.combined, npcs = 30, verbose = FALSE)
sample6.combined <-
  RunUMAP(sample6.combined, reduction = "pca", dims = 1:30)
sample6.combined <-
  FindNeighbors(sample6.combined, reduction = "pca", dims = 1:30)
sample6.combined <-
  FindClusters(sample6.combined, resolution = 0.17)#8
meta <- sample6.combined@meta.data
meta <- meta[, -c(10)]
sample6.combined@meta.data <- meta
save(sample6.combined, file = "./sample6_clusters_n.RData")


#umap
load("./sample6_clusters_n.RData")
#define order of legend
sample6.combined$Layer <-
  factor(
    sample6.combined$Layer,
    levels = c(
      "Layer 1",
      "Layer 2",
      "Layer 3",
      "Layer 4",
      "Layer 5",
      "Layer 6",
      "White Matter",
      "Noise"
    )
  )
# my.layer_table <- table(brain.integrated.rm.noise$seurat_clusters)
# my.layer_label <- paste0(names(my.layer_table),": (",my.layer_table,")")
my_color <-
  c(as.character(jcolors::jcolors(palette = "pal7"))[-8], "#808080")
p_u <-
  DimPlot(
    sample6.combined,
    group.by = "seurat_clusters",
    pt.size = 1,
    cols = my_color
  ) +
  #size of spot
  #scale_color_jcolors(palette = "pal7") + #color
  guides(color = guide_legend(override.aes = list(size = 5))) + #legend
  ggtitle("Unsupervised clustering")
p_u
setwd("/users/PAS1475/liuzeyi/guoqi/output/picture")
ggsave(
  p_u,
  filename = "./integrating landscape/unsurpervised_umap.tiff",
  device = "tiff",
  dpi = 300 ,
  width = 10,
  height = 10
)

# my.layer_table <- table(brain.integrated.rm.noise$Layer)
# my.layer_label <- paste0(names(my.layer_table),": (",my.layer_table,")")
p_a <-
  DimPlot(
    sample6.combined,
    group.by = "Layer",
    pt.size = 1,
    cols = my_color
  ) +
  #scale_color_jcolors(palette = "pal7") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  ggtitle("Manual annotation")
p_a
ggsave(
  p_a,
  filename = "./integrating landscape/annotation_umap.tiff",
  device = "tiff",
  dpi = 300 ,
  width = 10,
  height = 10
)


#spatial map-4.14
#unsupervised spatial
plot_spatial_map <-
  function(object = Seurat.object,
           group.by = "seurat_clusters",
           pt.size = 4.8) {
    my.cat <- unique(object$patientID)
    tmp.slice.name <- switch (
      my.cat,
      '2_5' = "slice1",
      '18_64' = "slice1.1",
      '2_8' = "slice1.2",
      'T4857' = "slice1.3",
      '1_1' = "slice1.4",
      '2_3' = "slice1.5"
    )
    xy.cor <-  object@images[[tmp.slice.name]]
    xy.cor <- xy.cor@coordinates
    my.meta <- object@meta.data[, group.by]
    plot.df <- as.data.frame(cbind(xy.cor, group.by = my.meta))
    p <-
      ggplot(plot.df, aes(
        x = imagecol,
        y = -imagerow,
        color = group.by
      )) +
      geom_point(size = pt.size)  + theme_void()
    if (any(grepl("noise", unique(my.meta), ignore.case = T))) {
      p <-
        p + scale_color_manual(
          values = c(as.character(jcolors::jcolors(palette = "pal7"))[-8], "#808080"),
          breaks = c(
            "Layer 1",
            "Layer 2",
            "Layer 3",
            "Layer 4",
            "Layer 5",
            "Layer 6",
            "White Matter",
            "Noise"
          )
        ) +
        labs(col = group.by, title = group.by)
    } else{
      p <-
        p + scale_color_manual(values = as.character(jcolors::jcolors(palette = "pal7")),
                               breaks = sort(unique(my.meta))) +
        labs(col = group.by, title = group.by)
    }
    return(p)
  }
# prepare object
table(sample6.combined$patientID)
obj.list <- SplitObject(sample6.combined, split.by = "patientID")
table(obj.list[[1]]$seurat_clusters)

p.2_5_un <-
  plot_spatial_map(object = obj.list[[1]], group.by = "seurat_clusters")
p.18_64_un <-
  plot_spatial_map(object = obj.list[[2]], group.by = "seurat_clusters")
p.2_8_un <-
  plot_spatial_map(object = obj.list[[3]], group.by = "seurat_clusters")
p.T4857_un <-
  plot_spatial_map(object = obj.list[[4]], group.by = "seurat_clusters")
p.1_1_un <-
  plot_spatial_map(object = obj.list[[5]], group.by = "seurat_clusters")
p.2_3_un <-
  plot_spatial_map(object = obj.list[[6]], group.by = "seurat_clusters")
# obj.list[[1]]@images[setdiff(name.list, i)]
# name.list <- seq(1:8)
# p.2_3 <- plot_spatial_map(object = obj.list[[6]], group.by = "Layer")
plot.un <-
  list(p.2_5_un, p.18_64_un, p.2_8_un, p.T4857_un, p.1_1_un, p.2_3_un)
names(plot.un) <- names(obj.list)

for (i in 1:6) {
  ggsave(
    plot = plot.un[[i]],
    filename = paste0(
      "/users/PAS1475/liuzeyi/guoqi/output/picture/integrating landscape/",
      names(plot.un)[i],
      "_un.tiff"
    ),
    device = "tiff",
    dpi = 150,
    width = 10,
    height = 10,
    units = "in"
  )
}



#annotation of spatial,4.14
p.2_5 <-
  plot_spatial_map(object = obj.list[[1]], group.by = "Layer")
p.18_64 <-
  plot_spatial_map(object = obj.list[[2]], group.by = "Layer")
p.2_8 <-
  plot_spatial_map(object = obj.list[[3]], group.by = "Layer")
p.T4857 <-
  plot_spatial_map(object = obj.list[[4]], group.by = "Layer")
p.1_1 <-
  plot_spatial_map(object = obj.list[[5]], group.by = "Layer")
p.2_3 <-
  plot_spatial_map(object = obj.list[[6]], group.by = "Layer")
plot.list <- list(p.2_5, p.18_64, p.2_8, p.T4857, p.1_1, p.2_3)
names(plot.list) <- names(obj.list)
for (i in 1:6) {
  ggsave(
    plot = plot.list[[i]],
    filename = paste0(
      "/users/PAS1475/liuzeyi/guoqi/output/picture/integrating landscape/",
      names(plot.list)[i],
      "_an.tiff"
    ),
    device = "tiff",
    dpi = 150,
    width = 10,
    height = 10,
    units = "in"
  )
}

#2-8 plot too small
#un-4.55
ggsave(
  plot = p.2_8_un,
  filename = paste0(
    "/users/PAS1475/liuzeyi/guoqi/output/picture/integrating landscape/",
    "p.2_8_un.tiff"
  ),
  device = "tiff",
  dpi = 150,
  width = 10,
  height = 10,
  units = "in"
)
#an-4.8
ggsave(
  plot = p.2_8,
  filename = paste0(
    "/users/PAS1475/liuzeyi/guoqi/output/picture/integrating landscape/",
    "p.2_8_an.tiff"
  ),
  device = "tiff",
  dpi = 150,
  width = 10,
  height = 10,
  units = "in"
)
