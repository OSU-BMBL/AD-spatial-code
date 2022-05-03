#------Identification of conserved markers across 7 layers in 6 samples and a/b test to determine the numbers
#of sample
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
dyn.load("/apps/hdf5-serial/gnu/9.1/1.12.0/lib/libhdf5_hl.so.200")
library(hdf5r)
options(future.globals.maxSize = 8000 * 1024 ^ 2)
#-------Module for data processing,including data loading,normalization,noise removing,annotation.
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
    "Layer_annotation/1-1 Layer annotation.csv",
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
#---noise removing
object.2_5 <- subset(object.2_5, subset = Layer != "Noise")
object.18_64 <- subset(object.18_64, subset = Layer != "Noise")
object.1_1 <- subset(object.1_1, subset = Layer != "Noise")
object.2_8 <- subset(object.2_8, subset = Layer != "Noise")
object_T4857 <- subset(object_T4857, subset = Layer != "Noise")
object.2_3 <- subset(object.2_3, subset = Layer != "Noise")
save(object.2_8, file = "/users/PAS1475/liuzeyi/guoqi/output/object.2_8.RData")
save(object_T4857, file = "/users/PAS1475/liuzeyi/guoqi/output/object_T4857.RData")
save(object.2_3, file = "/users/PAS1475/liuzeyi/guoqi/output/object.2_3.RData")
save(object.2_5, file = "/users/PAS1475/liuzeyi/guoqi/output/object.2_5.RData")
save(object.18_64, file = "/users/PAS1475/liuzeyi/guoqi/output/object.18_64.RData")
save(object.1_1, file = "/users/PAS1475/liuzeyi/guoqi/output/object.1_1.RData")
#------Integrating 6 samples
options(future.globals.maxSize = 8000 * 1024 ^ 2)
load("/users/PAS1475/liuzeyi/guoqi/output/my.all.object.list.RData")
#anchor feature
my.all.object.list <- list(
  object.2_5 = object.2_5,
  object.18_64 = object.18_64,
  object.2_8 = object.2_8,
  object_T4857 = object_T4857,
  object.1_1 = object.1_1,
  object.2_3 = object.2_3
)
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
  FindClusters(sample6.combined, resolution = 0.15)#8
meta <- sample6.combined@meta.data
meta <- meta[, -c(12)]
meta <- meta[, -10]
sample6.combined@meta.data <- meta
save(sample6.combined, file = "../../sample6_clusters.RData")


#------Module for identification of conserved markers in 6 samples
#Integration of 6 h5
sample_6 <-
  merge(object.2_5,
        c(object.18_64, object.1_1, object.2_8, object_T4857, object.2_3))
save(sample_6, file = "./specific_marker_final/sample6_merge.RData")
#Define type before DEG analysis
DefaultAssay(sample_6) <- "SCT"
allconserved_markers <- data.frame()
layer <- unique(sample_6$Layer)
for (i in layer) {
  markers <-
    FindConservedMarkers(
      sample_6,
      ident.1 = i ,
      grouping.var = "patientID",
      verbose = TRUE,
      assay = "SCT",
      only.pos = TRUE
    )
  markers$layer <- as.character(rep(i, nrow(markers)))
  if (ncol(markers) != (5 * 6 + 3) |
      nrow(markers) == 0) {
    #layer1 and layer3
    markers <- data.frame()
  }
  else{
    markers <- cbind(rownames(markers), markers)
    rownames(markers) <- c(1:nrow(markers))
  }
  #binding all conservedmarker in all layers
  allconserved_markers <- rbind(allconserved_markers, markers)#47
}

setwd("/users/PAS1475/liuzeyi/guoqi/output/")
write.csv(x = allconserved_markers, file = "allconserved_markers.csv")
allconserved_marker <-
  read.csv("allconserved_markers.csv",
           header = T,
           row.names = 1)
specific_markers_2456w <- allconserved_markers[, -c(1:31)]
colnames(specific_markers_2456w)[1] <- c("Combined_p")
colnames(specific_markers_2456w)[2] <- c("Layer")
write.csv(
  specific_markers_2456w,
  "conservedmarker/specific_marker_final/specific_markers_2456w.csv"
)

#------Module for A/B test,including filtering potential marker,similarity and variance test with different p value, test in each layer
my.all.object.list <- list(
  "2-5" = object.2_5,
  "18-64" = object.18_64,
  "2-8" = object.2_8,
  "T4857" = object_T4857,
  "1-1" = object.1_1,
  "2-3" = object.2_3
)
save(my.all.object.list, file = "/users/PAS1475/liuzeyi/guoqi/output/my.all.object.list.RData")
load("/users/PAS1475/liuzeyi/guoqi/output/my.all.object.list.RData")
#---A/B test-finding conserved marker in each layer with 1 sample,2 samples,3 samples,4 samples,5 samples
setwd(
  "/users/PAS1475/liuzeyi/guoqi/output/conservedmarker/specific-marker-correctrbind/"
)
#1 sample-1409
c_s_1 <- combn(1:6, 1)
for (i in 1:ncol(c_s_1)) {
  sample <- my.all.object.list[[as.numeric(c_s_1[i])]]
  DefaultAssay(sample) <- "SCT"
  #2-8 only has 6 layers
  layer_f <- intersect(as.character(unique(Idents(sample))), layer)
  allconserved_markers_temp <- data.frame()
  #findconservedmarker
  for (j in layer_f) {
    markers <-
      FindConservedMarkers(
        sample,
        ident.1 = j ,
        grouping.var = "patientID",
        verbose = TRUE,
        assay = "SCT",
        only.pos = TRUE
      )
    markers$layer <- as.character(rep(j, nrow(markers)))
    markers <- cbind(rownames(markers), markers)
    rownames(markers) <- c(1:nrow(markers))
    allconserved_markers_temp <-
      rbind(allconserved_markers_temp, markers)
  }#2365
  sample_name <- names(my.all.object.list)[as.numeric(c_s_1[i])]
  write.csv(
    allconserved_markers_temp,
    paste("level1", sample_name, "ConservedMarkers.csv", sep = "_")
  )
}
last <- rownames(allconserved_markers_temp)
new <- allconserved_markers_temp$`rownames(markers)`#700个不同的
#2 samples-628
c_s_2 <- combn(1:6, 2)
for (i in 1:ncol(c_s_2)) {
  sample1 <- my.all.object.list[[as.numeric(c_s_2[1, i])]]
  sample2 <- my.all.object.list[[as.numeric(c_s_2[2, i])]]
  sample <- merge(sample1, sample2)
  DefaultAssay(sample) <- "SCT"
  allconserved_markers_temp <- data.frame()
  for (j in layer) {
    markers <-
      FindConservedMarkers(
        sample,
        ident.1 = j ,
        grouping.var = "patientID",
        verbose = TRUE,
        assay = "SCT",
        only.pos = TRUE
      )
    markers$layer <- as.character(rep(j, nrow(markers)))
    if (ncol(markers) != (5 * nrow(c_s_2) + 3) | nrow(markers) == 0) {
      markers <- data.frame()
    }
    else{
      markers <- cbind(rownames(markers), markers)
      rownames(markers) <- c(1:nrow(markers))
    }
    allconserved_markers_temp <-
      rbind(allconserved_markers_temp, markers)
  }
  sample1_name <- names(my.all.object.list)[as.numeric(c_s_2[1, i])]
  sample2_name <- names(my.all.object.list)[as.numeric(c_s_2[2, i])]
  write.csv(
    allconserved_markers_temp,
    paste(
      "level2",
      sample1_name,
      sample2_name,
      "ConservedMarkers.csv",
      sep = "_"
    )
  )
}
#3 samples
c_s_3 <- combn(1:6, 3)
for (i in 1:ncol(c_s_3)) {
  sample1 <- my.all.object.list[[as.numeric(c_s_3[1, i])]]
  sample2 <- my.all.object.list[[as.numeric(c_s_3[2, i])]]
  sample3 <- my.all.object.list[[as.numeric(c_s_3[3, i])]]
  sample <- merge(sample1, c(sample2, sample3))
  DefaultAssay(sample) <- "SCT"
  allconserved_markers_temp <- data.frame()
  for (j in layer) {
    markers <-
      FindConservedMarkers(
        sample,
        ident.1 = j ,
        grouping.var = "patientID",
        verbose = TRUE,
        assay = "SCT",
        only.pos = TRUE
      )
    markers$layer <- as.character(rep(j, nrow(markers)))
    if (ncol(markers) != (5 * nrow(c_s_3) + 3) | nrow(markers) == 0) {
      markers <- data.frame()
    }
    else{
      markers <- cbind(rownames(markers), markers)
      rownames(markers) <- c(1:nrow(markers))
    }
    allconserved_markers_temp <-
      rbind(allconserved_markers_temp, markers)
  }
  sample1_name <- names(my.all.object.list)[as.numeric(c_s_3[1, i])]
  sample2_name <- names(my.all.object.list)[as.numeric(c_s_3[2, i])]
  sample3_name <- names(my.all.object.list)[as.numeric(c_s_3[3, i])]
  write.csv(
    allconserved_markers_temp,
    paste(
      "level3",
      sample1_name,
      sample2_name,
      sample3_name,
      "ConservedMarkers.csv",
      sep = "_"
    )
  )
}
#4 samples
c_s_4 <- combn(1:6, 4)
for (i in 1:ncol(c_s_4)) {
  sample1 <- my.all.object.list[[as.numeric(c_s_4[1, i])]]
  sample2 <- my.all.object.list[[as.numeric(c_s_4[2, i])]]
  sample3 <- my.all.object.list[[as.numeric(c_s_4[3, i])]]
  sample4 <- my.all.object.list[[as.numeric(c_s_4[4, i])]]
  sample <- merge(sample1, c(sample2, sample3, sample4))
  DefaultAssay(sample) <- "SCT"
  allconserved_markers_temp <- data.frame()
  for (j in layer) {
    markers <-
      FindConservedMarkers(
        sample,
        ident.1 = j ,
        grouping.var = "patientID",
        verbose = TRUE,
        assay = "SCT",
        only.pos = TRUE
      )
    markers$layer <- as.character(rep(j, nrow(markers)))
    if (ncol(markers) != (5 * nrow(c_s_4) + 3) | nrow(markers) == 0) {
      markers <- data.frame()
    }
    else{
      markers <- cbind(rownames(markers), markers)
      rownames(markers) <- c(1:nrow(markers))
    }
    allconserved_markers_temp <-
      rbind(allconserved_markers_temp, markers)
  }
  sample1_name <- names(my.all.object.list)[as.numeric(c_s_4[1, i])]
  sample2_name <- names(my.all.object.list)[as.numeric(c_s_4[2, i])]
  sample3_name <- names(my.all.object.list)[as.numeric(c_s_4[3, i])]
  sample4_name <- names(my.all.object.list)[as.numeric(c_s_4[4, i])]
  write.csv(
    allconserved_markers_temp,
    paste(
      "level4",
      sample1_name,
      sample2_name,
      sample3_name,
      sample4_name,
      "ConservedMarkers.csv",
      sep = "_"
    )
  )
}

#5 samples
c_s_5 <- combn(1:6, 5)
for (i in 1:ncol(c_s_5)) {
  sample1 <- my.all.object.list[[as.numeric(c_s_5[1, i])]]
  sample2 <- my.all.object.list[[as.numeric(c_s_5[2, i])]]
  sample3 <- my.all.object.list[[as.numeric(c_s_5[3, i])]]
  sample4 <- my.all.object.list[[as.numeric(c_s_5[4, i])]]
  sample5 <- my.all.object.list[[as.numeric(c_s_5[5, i])]]
  sample <- merge(sample1, c(sample2, sample3, sample4, sample5))
  DefaultAssay(sample) <- "SCT"
  allconserved_markers_temp <- data.frame()
  for (j in layer) {
    markers <-
      FindConservedMarkers(
        sample,
        ident.1 = j ,
        grouping.var = "patientID",
        verbose = TRUE,
        assay = "SCT",
        only.pos = TRUE
      )
    markers$layer <- as.character(rep(j, nrow(markers)))
    if (ncol(markers) != (5 * nrow(c_s_5) + 3) | nrow(markers) == 0) {
      markers <- data.frame()
    }
    else{
      markers <- cbind(rownames(markers), markers)
      rownames(markers) <- c(1:nrow(markers))
    }
    allconserved_markers_temp <-
      rbind(allconserved_markers_temp, markers)
  }
  sample1_name <- names(my.all.object.list)[as.numeric(c_s_5[1, i])]
  sample2_name <- names(my.all.object.list)[as.numeric(c_s_5[2, i])]
  sample3_name <- names(my.all.object.list)[as.numeric(c_s_5[3, i])]
  sample4_name <- names(my.all.object.list)[as.numeric(c_s_5[4, i])]
  sample5_name <- names(my.all.object.list)[as.numeric(c_s_5[5, i])]
  write.csv(
    allconserved_markers_temp,
    paste(
      "level5",
      sample1_name,
      sample2_name,
      sample3_name,
      sample4_name,
      sample5_name,
      "ConservedMarkers.csv",
      sep = "_"
    )
  )
}
#6
write.csv(
  allconserved_markers,
  "./conservedmarker/specific-marker-correctrbind/level6_conservedmarkers.csv"
)

#---filtering potential marker(removing duplicate genes,p<0.05/p<0.005/p<0.001)
#--removing duplicate genes
# setwd("/users/PAS1475/liuzeyi/guoqi/output/conservedmarker/")
# all_file<-list.files(".")
# all_file<-all_file[grep("level",all_file)]
# for (file in all_file){
#   marker_data<-read.csv(file,header = T,row.names = 1)
#   #duplicated
#   if(length(grep("\\.",rownames(marker_data)))!=0){
#     #rep_t<-grep("\\.",rownames(marker_data),value=T)
#     rep_list<-strsplit(grep("\\.",rownames(marker_data),value = T),split="\\.")
#     rep_tt<-c()
#     for(i in rep_list){
#       rep_tt<-c(rep_tt,i[1])
#     }
#     if(length(na.omit(match(rep_tt,rownames(marker_data))))!=0){
#       print(file)
#       print(na.omit(match(rep_tt,rownames(marker_data))))
#     }
#   }
# }
#--p value
#p<0.05
setwd("../specific-marker-correctrbind/")
all_file <- list.files(".")
for (file in all_file) {
  marker_data <- read.csv(file, header = T, row.names = 1)
  if (grepl("^level1", file)) {
    pre_marker_df <- marker_data[marker_data[, 6] < 0.05, ]
  }
  else{
    pre_marker_df <- marker_data[marker_data$minimump_p_val < 0.05, ]
  }
  write.csv(pre_marker_df,
            paste("../allconservedmarker_0.05_new_csv/", file, sep = ""))
}

#-10
for (file in all_file) {
  marker_data <- read.csv(file, header = T, row.names = 1)
  if (grepl("^level1", file)) {
    pre_marker_df <- marker_data[marker_data[, 6] < 1e-10, ]
  }
  else{
    pre_marker_df <- marker_data[marker_data$minimump_p_val < 1e-10, ]
  }
  write.csv(pre_marker_df,
            paste("../allconservedmarker_-10_new_csv/",
                  file, sep = ""))
}

#---similar coefficient
#jaccard
jaccard <- function(marker) {
  intersect_length <- length(intersect(marker, markers_6))
  union_length <- length(union(marker, markers_6))
  coefficient <- intersect_length / union_length
  return(coefficient)
}
level <-
  c(
    rep("6 samples", 1),
    rep("5 samples", 6),
    rep("4 samples", 15),
    rep("3 samples", 20),
    rep("2 samples", 15),
    rep("1 samples", 6)
  )
level <- factor(level, levels = unique(level))
#0.05
setwd("../allconservedmarker_0.05_new_csv/")
sample6 <-
  read.csv("level6_conservedmarkers.csv",
           header = T,
           row.names = 1)
markers_6 <- sample6$rownames.markers.
allfile_0.05 <- list.files(".")
t_allfile_0.05 <- rev(allfile_0.05)
jac <- c()
for (file in t_allfile_0.05) {
  marker_data <- read.csv(file, header = T, row.names = 1)
  coefficient <- jaccard(marker_data$rownames.markers.)
  jac <- c(jac, coefficient)
}
box_data_0.05 <-
  data.frame(sample = t_allfile_0.05,
             Level = level,
             Jaccard = jac)
#var
Level <- unique(level)
var_all <- c()
for (l in Level) {
  index <- which(box_data_0.05$Level == l)
  jac <- box_data_0.05$Jaccard[index]
  var <- var(jac)
  var_all <- c(var_all, var)
}
var_all[1] <- 0
var_data <- data.frame(Level = Level, Var = var_all)
#adding labels of x
ann_x <- c()
for (i in 1:6) {
  a <-
    paste(Level[i], "(", as.character(round(var_data$Var[i], 3)), ")", sep =
            "")
  ann_x <- c(ann_x, a)
}
box_data_0.05$x <- c(1, rep(2, 6), rep(3, 15), rep(4, 20), rep(0, 21))
median_box <- c()
for (i in unique(level)) {
  m_temp <-
    as.numeric(median(box_data_0.05[which(box_data_0.05$Level == i), 3]))
  m_temp <- rep(m_temp, length(which(box_data_0.05$Level == i)))
  median_box <- c(median_box, m_temp)
}
box_data_0.05$median <- median_box
box_0.05_all <- box_data_0.05
box_0.05_all_new <- box_0.05_all
g <-
  ggplot(box_0.05_all) + aes(x = Level, y = Jaccard) + geom_boxplot(outlier.size = -1) +
  geom_line(
    data = box_0.05_all[1:42,],
    aes(group = 1, y = median),
    lwd = 1.5,
    show.legend = F,
    colour = "blue"
  ) +
  geom_jitter() +
  labs(title = "A/B test_0.05_all_layers") +
  theme(
    axis.text = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(
      size = 15,
      face = "bold",
      angle = 40,
      hjust = 1
    ),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold")
  )
g_0.05_all <- g + scale_x_discrete(breaks = Level,
                                   labels = ann_x)
#each layer
#layer2
jac_layer2 <- c()
markers_6 <- sample6$rownames.markers.[sample6$layer == "Layer 2"]
for (file in t_allfile_0.05) {
  all_marker <- read.csv(file, header = T, row.names = 1)
  marker <- all_marker$rownames.markers.[all_marker$layer == "Layer 2"]
  coefficient <- jaccard(marker)
  jac_layer2 <- c(jac_layer2, coefficient)
}
#layer4
markers_6 <- sample6$rownames.markers.[sample6$layer == "Layer 4"]
jac_layer4 <- c()
for (file in t_allfile_0.05) {
  all_marker <- read.csv(file, header = T, row.names = 1)
  marker <- all_marker$rownames.markers.[all_marker$layer == "Layer 4"]
  coefficient <- jaccard(marker)
  jac_layer4 <- c(jac_layer4, coefficient)
}
#layer5
markers_6 <- sample6$rownames.markers.[sample6$layer == "Layer 5"]
jac_layer5 <- c()
for (file in t_allfile_0.05) {
  all_marker <- read.csv(file, header = T, row.names = 1)
  marker <- all_marker$rownames.markers.[all_marker$layer == "Layer 5"]
  coefficient <- jaccard(marker)
  jac_layer5 <- c(jac_layer5, coefficient)
}
#layer6
markers_6 <- sample6$rownames.markers.[sample6$layer == "Layer 6"]
jac_layer6 <- c()
for (file in t_allfile_0.05) {
  all_marker <- read.csv(file, header = T, row.names = 1)
  marker <- all_marker$rownames.markers.[all_marker$layer == "Layer 6"]
  coefficient <- jaccard(marker)
  jac_layer6 <- c(jac_layer6, coefficient)
}
#white matter
markers_6 <- sample6$rownames.markers.[sample6$layer == "White Matter"]
jac_layer_wm <- c()
for (file in t_allfile_0.05) {
  all_marker <- read.csv(file, header = T, row.names = 1)
  marker <-
    all_marker$rownames.markers.[all_marker$layer == "White Matter"]
  coefficient <- jaccard(marker)
  jac_layer_wm <- c(jac_layer_wm, coefficient)
}
#plot
jac_layer <-
  list(
    Layer2 = jac_layer2,
    Layer4 = jac_layer4,
    Layer5 = jac_layer5,
    Layer6 = jac_layer6,
    Layer_wm = jac_layer_wm
  )
var_all <- c()
b <- list()
for (i in 1:5) {
  box_data <-
    data.frame(sample = t_allfile_0.05,
               Level = level,
               Jaccard = jac_layer[[i]])
  box_data$x <- c(1, rep(2, 6), rep(3, 15), rep(4, 20), rep(0, 21))
  median_box <- c()
  for (j in unique(level)) {
    m_temp <- as.numeric(median(box_data[which(box_data$Level == j), 3]))
    m_temp <- rep(m_temp, length(which(box_data$Level == j)))
    median_box <- c(median_box, m_temp)
  }
  box_data$median <- median_box
  g <-
    ggplot(box_data) + aes(x = Level, y = Jaccard) + geom_boxplot(outlier.size = -1) +
    geom_line(
      data = box_data[1:42,],
      aes(group = 1, y = median),
      lwd = 1.5,
      show.legend = F,
      colour = "blue"
    ) +
    geom_jitter() +
    labs(title = paste("A/B test_0.05", names(jac_layer)[i], sep = "_")) +
    theme(
      axis.text = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(
        size = 15,
        face = "bold",
        angle = 40,
        hjust = 1
      ),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(size = 20, face = "bold")
    )
  #--variation test in each layer
  var_all <- c()
  for (l in Level) {
    index <- which(box_data$Level == l)
    jac <- box_data$Jaccard[index]
    var <- var(jac)
    var_all <- c(var_all, var)
  }
  var_all[1] <- 0
  var_data <- data.frame(Level = Level, Var = var_all)
  ann_x <- c()
  for (h in 1:6) {
    a <-
      paste(Level[h], "(", as.character(round(var_data$Var[h], 3)), ")", sep =
              "")
    ann_x <- c(ann_x, a)
  }
  g_f <- g + scale_x_discrete(breaks = Level,
                              labels = ann_x)
  b[[i]] <- g_f
}
# mutiple plot
g_0.05 <- g_0.05_all + b[1] + b[2] + b[3] + b[4] + b[5]
ggsave(
  plot = g_0.05,
  filename = paste("../../picture/",
                   "abtest_0.05_all_new.png", sep = ""),
  device = "png",
  dpi = 150,
  width = 30,
  height = 12,
  units = "in"
)
#1e-10

setwd("../allconservedmarker_-10_new_csv/")
sample6 <-
  read.csv("level6_conservedmarkers.csv",
           header = T,
           row.names = 1)
markers_6 <- sample6$rownames.markers.
allfile_0.05 <- list.files(".")
t_allfile_0.05 <- rev(allfile_0.05)
jac <- c()
for (file in t_allfile_0.05) {
  marker_data <- read.csv(file, header = T, row.names = 1)
  coefficient <- jaccard(marker_data$rownames.markers.)
  jac <- c(jac, coefficient)
}
box_data_0.05 <-
  data.frame(sample = t_allfile_0.05,
             Level = level,
             Jaccard = jac)
#var
Level <- unique(level)
var_all <- c()
for (l in Level) {
  index <- which(box_data_0.05$Level == l)
  jac <- box_data_0.05$Jaccard[index]
  var <- var(jac)
  var_all <- c(var_all, var)
}
var_all[1] <- 0
var_data <- data.frame(Level = Level, Var = var_all)
#adding labels of x
ann_x <- c()
for (i in 1:6) {
  a <-
    paste(Level[i], "(", as.character(round(var_data$Var[i], 3)), ")", sep =
            "")
  ann_x <- c(ann_x, a)
}
box_data_0.05$x <- c(1, rep(2, 6), rep(3, 15), rep(4, 20), rep(0, 21))
median_box <- c()
for (i in unique(level)) {
  m_temp <-
    as.numeric(median(box_data_0.05[which(box_data_0.05$Level == i), 3]))
  m_temp <- rep(m_temp, length(which(box_data_0.05$Level == i)))
  median_box <- c(median_box, m_temp)
}
box_data_0.05$median <- median_box
box_0.05_all <- box_data_0.05
box_0.05_all_new <- box_0.05_all
g <-
  ggplot(box_0.05_all) + aes(x = Level, y = Jaccard) + geom_boxplot(outlier.size = -1) +
  geom_line(
    data = box_0.05_all[1:42,],
    aes(group = 1, y = median),
    lwd = 1.5,
    show.legend = F,
    colour = "blue"
  ) +
  geom_jitter() +
  labs(title = "A/B test_1e-10_all_layers") +
  theme(
    axis.text = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(
      size = 15,
      face = "bold",
      angle = 40,
      hjust = 1
    ),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold")
  )
g_0.05_all <- g + scale_x_discrete(breaks = Level,
                                   labels = ann_x)
#each layer
#layer2
jac_layer2 <- c()
markers_6 <- sample6$rownames.markers.[sample6$layer == "Layer 2"]
for (file in t_allfile_0.05) {
  all_marker <- read.csv(file, header = T, row.names = 1)
  marker <- all_marker$rownames.markers.[all_marker$layer == "Layer 2"]
  coefficient <- jaccard(marker)
  jac_layer2 <- c(jac_layer2, coefficient)
}
#layer4
markers_6 <- sample6$rownames.markers.[sample6$layer == "Layer 4"]
jac_layer4 <- c()
for (file in t_allfile_0.05) {
  all_marker <- read.csv(file, header = T, row.names = 1)
  marker <- all_marker$rownames.markers.[all_marker$layer == "Layer 4"]
  coefficient <- jaccard(marker)
  jac_layer4 <- c(jac_layer4, coefficient)
}
#layer5
markers_6 <- sample6$rownames.markers.[sample6$layer == "Layer 5"]
jac_layer5 <- c()
for (file in t_allfile_0.05) {
  all_marker <- read.csv(file, header = T, row.names = 1)
  marker <- all_marker$rownames.markers.[all_marker$layer == "Layer 5"]
  coefficient <- jaccard(marker)
  jac_layer5 <- c(jac_layer5, coefficient)
}
#layer6
markers_6 <- sample6$rownames.markers.[sample6$layer == "Layer 6"]
jac_layer6 <- c()
for (file in t_allfile_0.05) {
  all_marker <- read.csv(file, header = T, row.names = 1)
  marker <- all_marker$rownames.markers.[all_marker$layer == "Layer 6"]
  coefficient <- jaccard(marker)
  jac_layer6 <- c(jac_layer6, coefficient)
}
#white matter
markers_6 <- sample6$rownames.markers.[sample6$layer == "White Matter"]
jac_layer_wm <- c()
for (file in t_allfile_0.05) {
  all_marker <- read.csv(file, header = T, row.names = 1)
  marker <-
    all_marker$rownames.markers.[all_marker$layer == "White Matter"]
  coefficient <- jaccard(marker)
  jac_layer_wm <- c(jac_layer_wm, coefficient)
}
#plot
jac_layer <-
  list(
    Layer2 = jac_layer2,
    Layer4 = jac_layer4,
    Layer5 = jac_layer5,
    Layer6 = jac_layer6,
    Layer_wm = jac_layer_wm
  )
var_all <- c()
b <- list()
for (i in 1:5) {
  box_data <-
    data.frame(sample = t_allfile_0.05,
               Level = level,
               Jaccard = jac_layer[[i]])
  box_data$x <- c(1, rep(2, 6), rep(3, 15), rep(4, 20), rep(0, 21))
  median_box <- c()
  for (j in unique(level)) {
    m_temp <- as.numeric(median(box_data[which(box_data$Level == j), 3]))
    m_temp <- rep(m_temp, length(which(box_data$Level == j)))
    median_box <- c(median_box, m_temp)
  }
  box_data$median <- median_box
  g <-
    ggplot(box_data) + aes(x = Level, y = Jaccard) + geom_boxplot(outlier.size = -1) +
    geom_line(
      data = box_data[1:42,],
      aes(group = 1, y = median),
      lwd = 1.5,
      show.legend = F,
      colour = "blue"
    ) +
    geom_jitter() +
    labs(title = paste("A/B test_1e-10", names(jac_layer)[i], sep = "_")) +
    theme(
      axis.text = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(
        size = 15,
        face = "bold",
        angle = 40,
        hjust = 1
      ),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(size = 20, face = "bold")
    )
  #--variation test in each layer
  var_all <- c()
  for (l in Level) {
    index <- which(box_data$Level == l)
    jac <- box_data$Jaccard[index]
    var <- var(jac)
    var_all <- c(var_all, var)
  }
  var_all[1] <- 0
  var_data <- data.frame(Level = Level, Var = var_all)
  ann_x <- c()
  for (h in 1:6) {
    a <-
      paste(Level[h], "(", as.character(round(var_data$Var[h], 3)), ")", sep =
              "")
    ann_x <- c(ann_x, a)
  }
  g_f <- g + scale_x_discrete(breaks = Level,
                              labels = ann_x)
  b[[i]] <- g_f
}
# mutiple plot
g_0.05 <- g_0.05_all + b[1] + b[2] + b[3] + b[4] + b[5]
ggsave(
  plot = g_0.05,
  filename = paste("../../picture/",
                   "abtest_1e-10_all_new.png", sep = ""),
  device = "png",
  dpi = 150,
  width = 30,
  height = 12,
  units = "in"
)
