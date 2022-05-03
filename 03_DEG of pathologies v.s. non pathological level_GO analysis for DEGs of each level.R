library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
dyn.load("/apps/hdf5-serial/gnu/9.1/1.12.0/lib/libhdf5_hl.so.200")
library(hdf5r)
options(future.globals.maxSize = 8000 * 1024 ^ 2)
setwd("/fs/ess/PCON0022/liyang/Shane/pathology")
ann_2_3 <-
  read.csv(
    "./confirmed_Abeta_annotation/2-3 Aβ annotation.csv",
    header = T,
    row.names = 1
  )
ann_2_8 <-
  read.csv(
    "confirmed_Abeta_annotation/2-8 Aβ annotation.csv",
    header = T,
    row.names = 1
  )
ann_T4857 <-
  read.csv(
    "confirmed_Abeta_annotation/T4857 Aβ annotation.csv",
    header = T,
    row.names = 1
  )

setwd("/users/PAS1475/liuzeyi/guoqi/output")
load("object.2_3.RData")
load("object.2_8.RData")
load("object_T4857.RData")

##Aβ
#--------2-3
identical(colnames(object.2_3), rownames(ann_2_3))#false
dim(object.2_3)
dim(ann_2_3)
length(which(ann_2_3$Aβ == ""))
intersect <-
  intersect(colnames(object.2_3), rownames(ann_2_3))#spot from annotation inclue spats from object
#check if this ann is identical with original layer annotation(noise cause fause)
label <-
  read.csv(
    "Layer_annotation/2-3 layer manuallabel new.csv",
    header = T,
    row.names = 1
  )
identical(rownames(label), rownames(ann_2_3))#T
#removing noise spot in ann of pathology
ann_2_3_f <-
  as.data.frame(ann_2_3[match(colnames(object.2_3), rownames(ann_2_3)), ])
rownames(ann_2_3_f) <-
  rownames(ann_2_3)[match(colnames(object.2_3), rownames(ann_2_3))]
colnames(ann_2_3_f) <- "condition"
setwd("/users/PAS1475/liuzeyi/guoqi/output/pathology")
write.csv(ann_2_3_f, "2-3 Aβ annotation.csv")
Ab2_3 <- read.csv("2-3 Aβ annotation.csv",
                  header = T,
                  row.names = 1)
#2-8
identical(colnames(object.2_8), rownames(ann_2_8))#T
dim(object.2_8)
dim(ann_2_8)
write.csv(ann_2_8, "2-8 Aβ annotation.csv")
Ab2_8 <- read.csv("2-8 Aβ annotation.csv",
                  header = T,
                  row.names = 1)
colnames(Ab2_8) <- colnames(T4857)
#T
identical(colnames(object_T4857), rownames(ann_T4857))#T
dim(object_T4857)
dim(ann_T4857)
intersect <- intersect(colnames(object_T4857), rownames(ann_T4857))
dif_T4857 <-
  as.data.frame(ann_T4857[-match(colnames(object_T4857), rownames(ann_T4857)), ])
#same as 2_3
label <-
  read.csv(
    "Layer_annotation/T4857 layer manuallabel new.csv",
    header = T,
    row.names = 1
  )
identical(rownames(label), rownames(ann_T4857))#T
#
ann_T4857_f <-
  as.data.frame(ann_T4857[match(colnames(object_T4857), rownames(ann_T4857)), ])
rownames(ann_T4857_f) <-
  rownames(ann_T4857)[match(colnames(object_T4857), rownames(ann_T4857))]
colnames(ann_T4857_f) <- "condition"
write.csv(ann_T4857_f, "T4857 Aβ annotation.csv")
T4857 <- read.csv("T4857 Aβ annotation.csv",
                  header = T,
                  row.names = 1)

# #addmeta infor
# object.2_3 <- AddMetaData(object.2_3, metadata = Ab2_3)
# object_T4857 <- AddMetaData(object_T4857, metadata = T4857)
# object.2_8 <- AddMetaData(object.2_8, metadata = Ab2_8)
# object_ab <- merge(object.2_3, c(object_T4857, object.2_8))
# DefaultAssay(object_ab) <- "SCT"
# ab_cat <- unique(object_ab$condition)
# ab_cat <- ab_cat[-5]
# Idents(object_ab) <- object_ab$condition
# #conserved_markers
# allconserved_markers <- data.frame()
# for (i in ab_cat[-1]) {
#   markers <-
#     FindConservedMarkers(
#       object_ab,
#       ident.1 = ab_cat[1],
#       ident.2 = i,
#       grouping.var = "patientID",
#       verbose = TRUE,
#       assay = "SCT",
#       only.pos = F
#     )
#   markers$level <- as.character(rep(i, nrow(markers)))
#   allconserved_markers <- rbind(allconserved_markers, markers)
# }
# 
# ab_down <-
#   allconserved_markers[allconserved_markers$`2_8_avg_log2FC` < 0, ]
# ab_up <-
#   allconserved_markers[allconserved_markers$`2_8_avg_log2FC` > 0, ]
# save(ab_down, file = "./result/Aβ_down.RData")
# save(ab_up, file = "./result/Aβ_up.RData")

##Abeta-specific
setwd("/fs/ess/PCON0022/liyang/Shane/pathology")
T4857_sp <-
  read.csv(
    "confirmed_Abeta_annotation_specific/T4957 Aβ specific.csv",
    header = T,
    row.names = 1
  )
ab_2_3_sp <-
  read.csv(
    "confirmed_Abeta_annotation_specific/2-3 Aβ specific.csv",
    header = T,
    row.names = 1
  )
ab_2_8_sp <-
  read.csv(
    "confirmed_Abeta_annotation_specific/2-8 Aβ specific.csv",
    header = T,
    row.names = 1
  )
#--remove noise
identical(colnames(object.2_3), rownames(ab_2_3_sp))
identical(colnames(object.2_8), rownames(ab_2_8_sp))#T
identical(colnames(object_T4857), rownames(T4857_sp))

list_object <- list(object.2_3, object_T4857)
list_data <- list(ab_2_3_sp, T4857_sp)
names(list_data) <-
  c("2-3 Aβ annotation_specific", "T4857 Aβ annotation_specific")
for (i in 1:2) {
  original <- list_data[[i]]
  f <-
    as.data.frame(original[match(colnames(list_object[[i]]), rownames(original)), ])
  rownames(f) <-
    rownames(original)[match(colnames(list_object[[i]]), rownames(original))]
  print(identical(colnames(list_object[[i]]), rownames(f)))#T
  colnames(f) <- "condition"
  write.csv(f,
            paste0(
              "/users/PAS1475/liuzeyi/guoqi/output/pathology/",
              names(list_data)[i],
              ".csv"
            ))
}
colnames(ab_2_8_sp) <- "condition"
write.csv(
  ab_2_8_sp,
  "/users/PAS1475/liuzeyi/guoqi/output/pathology/2-8 Aβ annotation_specific.csv"
)

##AT8
setwd("/fs/ess/PCON0022/liyang/Shane/pathology")
AT8_T4857 <-
  read.csv(
    "confirmed_AT8_annotation/T4857 AT8 annotation.csv",
    header = T,
    row.names = 1
  )
AT8_2_3 <-
  read.csv(
    "confirmed_AT8_annotation/2-3 AT8 annotation.csv",
    header = T,
    row.names = 1
  )
AT8_2_8 <-
  read.csv(
    "confirmed_AT8_annotation/2-8 AT8 annotation.csv",
    header = T,
    row.names = 1
  )
#--remove noise
identical(colnames(object.2_3), rownames(AT8_2_3))
identical(colnames(object.2_8), rownames(AT8_2_8))#T
identical(colnames(object_T4857), rownames(AT8_T4857))
list_data <- list(AT8_2_3, AT8_T4857)
names(list_data) <- c("2-3 AT8 annotation", "T4857 AT8 annotation")
for (i in 1:2) {
  original <- list_data[[i]]
  f <-
    as.data.frame(original[match(colnames(list_object[[i]]), rownames(original)), ])
  rownames(f) <-
    rownames(original)[match(colnames(list_object[[i]]), rownames(original))]
  print(identical(colnames(list_object[[i]]), rownames(f)))#T
  colnames(f) <- "condition"
  write.csv(f,
            paste0(
              "/users/PAS1475/liuzeyi/guoqi/output/pathology/",
              names(list_data)[i],
              ".csv"
            ))
}
colnames(AT8_2_8) <- "condition"
write.csv(AT8_2_8,
          "/users/PAS1475/liuzeyi/guoqi/output/pathology/2-8 AT8 annotation.csv")

##AT8_specific
setwd("/fs/ess/PCON0022/liyang/Shane/pathology")
AT8_T4857_sp <-
  read.csv(
    "confirmed_AT8_annotation_specific/T4957 AT8 specific.csv",
    header = T,
    row.names = 1
  )
AT8_2_3_sp <-
  read.csv(
    "confirmed_AT8_annotation_specific/2-3 AT8 specific.csv",
    header = T,
    row.names = 1
  )
AT8_2_8_sp <-
  read.csv(
    "confirmed_AT8_annotation_specific/2-8 AT8 specific.csv",
    header = T,
    row.names = 1
  )
#--remove noise
identical(colnames(object.2_3), rownames(AT8_2_3_sp))
identical(colnames(object.2_8), rownames(AT8_2_8_sp))#T
identical(colnames(object_T4857), rownames(AT8_T4857_sp))
list_data <- list(AT8_2_3_sp, AT8_T4857_sp)
names(list_data) <-
  c("2-3 AT8 annotation_specific", "T4857 AT8 annotation_specific")
for (i in 1:2) {
  original <- list_data[[i]]
  f <-
    as.data.frame(original[match(colnames(list_object[[i]]), rownames(original)), ])
  rownames(f) <-
    rownames(original)[match(colnames(list_object[[i]]), rownames(original))]
  print(identical(colnames(list_object[[i]]), rownames(f)))#T
  colnames(f) <- "condition"
  write.csv(f,
            paste0(
              "/users/PAS1475/liuzeyi/guoqi/output/pathology/",
              names(list_data)[i],
              ".csv"
            ))
}
colnames(AT8_2_8_sp) <- "condition"
write.csv(
  AT8_2_8_sp,
  "/users/PAS1475/liuzeyi/guoqi/output/pathology/2-8 AT8 annotation_specific.csv"
)

##
setwd("/fs/ess/PCON0022/liyang/Shane/pathology")
Double_T4857_sp <-
  read.csv(
    "confirmed_double_positive_annotation_specific/T4957 AT8,Aβ double positive specific.csv",
    header = T,
    row.names = 1
  )
Double_2_3_sp <-
  read.csv(
    "confirmed_double_positive_annotation_specific/2-3 AT8,Aβ double positive specific.csv",
    header = T,
    row.names = 1
  )
Double_2_8_sp <-
  read.csv(
    "confirmed_double_positive_annotation_specific/2-8 AT8,Aβ double positive specific.csv",
    header = T,
    row.names = 1
  )
Double_T4857 <-
  read.csv(
    "confirmed_double_positive_annotation/T4857 double positive.csv",
    header = T,
    row.names = 1
  )
Double_2_3 <-
  read.csv(
    "confirmed_double_positive_annotation/2-3 double positive.csv",
    header = T,
    row.names = 1
  )
Double_2_8 <-
  read.csv(
    "confirmed_double_positive_annotation/2-8 double positive.csv",
    header = T,
    row.names = 1
  )
#--remove noise
identical(colnames(object.2_3), rownames(Double_2_3_sp))
identical(colnames(object.2_3), rownames(Double_2_3))
identical(colnames(object.2_8), rownames(Double_2_8))#T
identical(colnames(object.2_8), rownames(Double_2_8_sp))#T
identical(colnames(object_T4857), rownames(Double_T4857_sp))
identical(colnames(object_T4857), rownames(Double_T4857))
list_data <-
  list(Double_2_3_sp, Double_T4857_sp, Double_2_3, Double_T4857)
list_object <- list(object.2_3, object_T4857, object.2_3, object_T4857)
names(list_data) <-
  c(
    "2-3 Double annotation_specific",
    "T4857 Double annotation_specific",
    "2-3 Double annotation",
    "T4857 Double annotation"
  )
for (i in 1:4) {
  original <- list_data[[i]]
  f <-
    as.data.frame(original[match(colnames(list_object[[i]]), rownames(original)), ])
  rownames(f) <-
    rownames(original)[match(colnames(list_object[[i]]), rownames(original))]
  print(identical(colnames(list_object[[i]]), rownames(f)))#T
  colnames(f) <- "condition"
  write.csv(f,
            paste0(
              "/users/PAS1475/liuzeyi/guoqi/output/pathology/",
              names(list_data)[i],
              ".csv"
            ))
}
colnames(Double_2_8) <- "condition"
write.csv(
  Double_2_8,
  "/users/PAS1475/liuzeyi/guoqi/output/pathology/2-8 Double annotation.csv"
)
colnames(Double_2_8_sp) <- "condition"
write.csv(
  Double_2_8_sp,
  "/users/PAS1475/liuzeyi/guoqi/output/pathology/2-8 Double annotation_specific.csv"
)

#Double WFS1_AT8
##
setwd("/fs/ess/PCON0022/liyang/Shane/pathology")
Double_W_AT8_T4857_sp <-
  read.csv(
    "confirmed_WFS1_AT8_double_positive_annotation_specific/T4857 WFS1 AT8 double positive specific.csv",
    header = T,
    row.names = 1
  )
Double_W_AT8_2_3_sp <-
  read.csv(
    "confirmed_WFS1_AT8_double_positive_annotation_specific/2-3 WFS1 AT8 double positive specific.csv",
    header = T,
    row.names = 1
  )
Double_W_AT8_2_8_sp <-
  read.csv(
    "confirmed_WFS1_AT8_double_positive_annotation_specific/2-8 WFS1 AT8 double positive specific.csv",
    header = T,
    row.names = 1
  )
Double_W_AT8_T4857 <-
  read.csv(
    "confirmed_WFS1_AT8_double_positive_annotation/T4857 WFS1 AT8 double positive nonspecific.csv",
    header = T,
    row.names = 1
  )
Double_W_AT8_2_3 <-
  read.csv(
    "confirmed_WFS1_AT8_double_positive_annotation/2-3 WFS1 AT8 double positive nonspecific.csv",
    header = T,
    row.names = 1
  )
Double_W_AT8_2_8 <-
  read.csv(
    "confirmed_WFS1_AT8_double_positive_annotation/2-8 WFS1 AT8 double positive nonspecific.csv",
    header = T,
    row.names = 1
  )
#--remove noise
identical(colnames(object.2_3), rownames(Double_W_AT8_2_3_sp))
identical(colnames(object.2_3), rownames(Double_W_AT8_2_3))
identical(colnames(object.2_8), rownames(Double_W_AT8_2_8_sp))#T
identical(colnames(object.2_8), rownames(Double_W_AT8_2_8))#T
identical(colnames(object_T4857), rownames(Double_W_AT8_T4857_sp))
identical(colnames(object_T4857), rownames(Double_W_AT8_T4857))
list_object <- list(object.2_3, object_T4857, object.2_3, object_T4857)
list_data <-
  list(Double_W_AT8_2_3_sp,
       Double_W_AT8_T4857_sp,
       Double_W_AT8_2_3,
       Double_W_AT8_T4857)
names(list_data) <-
  c(
    "2-3 Double_W_AT8 annotation_specific",
    "T4857 Double_W_AT8 annotation_specific",
    "2-3 Double_W_AT8 annotation",
    "T4857 Double_W_AT8 annotation"
  )
for (i in 1:4) {
  original <- list_data[[i]]
  f <-
    as.data.frame(original[match(colnames(list_object[[i]]), rownames(original)), ])
  rownames(f) <-
    rownames(original)[match(colnames(list_object[[i]]), rownames(original))]
  print(identical(colnames(list_object[[i]]), rownames(f)))#T
  colnames(f) <- "condition"
  write.csv(f,
            paste0(
              "/users/PAS1475/liuzeyi/guoqi/output/pathology/",
              names(list_data)[i],
              ".csv"
            ))
}
colnames(Double_W_AT8_2_8) <- "condition"
write.csv(
  Double_W_AT8_2_8,
  "/users/PAS1475/liuzeyi/guoqi/output/pathology/2-8 Double_W_AT8 annotation.csv"
)
colnames(Double_W_AT8_2_8_sp) <- "condition"
write.csv(
  Double_W_AT8_2_8_sp,
  "/users/PAS1475/liuzeyi/guoqi/output/pathology/2-8 Double_W_AT8 annotation_specific.csv"
)
#T4857 error: non_level_
temp_a <-
  read.csv(
    "T4857 Double_W_AT8 annotation_specific.csv",
    header = TRUE,
    row.names = 1
  )
temp_b <-
  read.csv("T4857 Double_W_AT8 annotation.csv",
           header = T,
           row.names = 1)
temp_a[which(temp_a$condition == "non_level_"), 1] <- c("non_level_1")
temp_b[which(temp_b$condition == "non_level_"), 1] <- c("non_level_1")
write.csv(
  temp_a,
  "/users/PAS1475/liuzeyi/guoqi/output/pathology/T4857 Double_W_AT8 annotation_specific.csv"
)
write.csv(
  temp_b,
  "/users/PAS1475/liuzeyi/guoqi/output/pathology/T4857 Double_W_AT8 annotation.csv"
)



# ##---------findconservedmarker
# #ab_specific
# setwd("/users/PAS1475/liuzeyi/guoqi/output/pathology/")
# Ab2_3_sp<-read.csv("2-3 Aβ annotation_specific.csv",header=T,row.names = 1)
# Ab2_8_sp<-read.csv("2-8 Aβ annotation_specific.csv",header=T,row.names = 1)
# Ab4857_sp<-read.csv("T4857 Aβ annotation_specific.csv",header=T,row.names = 1)
# #addmeta infor
# object.2_3 <- AddMetaData(object.2_3,metadata = Ab2_3_sp)
# object_T4857 <- AddMetaData(object_T4857,metadata = Ab4857_sp)
# object.2_8 <- AddMetaData(object.2_8,metadata = Ab2_8_sp)
# object_ab<-merge(object.2_3,c(object_T4857,object.2_8))
# DefaultAssay(object_ab) <- "SCT"
# ab_cat<-sort(unique(object_ab$condition))
# ab_cat<-ab_cat[-1]
# Idents(object_ab)<-object_ab$condition
# allconserved_markers<-data.frame()
# for(i in ab_cat[-1]){
#   markers <- FindConservedMarkers(object_ab, ident.1 =ab_cat[1],ident.2 =i, grouping.var = "patientID",
#                                   verbose = TRUE,
#                                   assay = "SCT",
#                                   only.pos=F)
#   markers$level<-as.character(rep(i,nrow(markers)))
#   allconserved_markers<-rbind(allconserved_markers,markers)
# }
# ab_down<-allconserved_markers[allconserved_markers$`2_8_avg_log2FC`<0,]
# ab_up<-allconserved_markers[allconserved_markers$`2_8_avg_log2FC`>0,]
# save(ab_down,file="./result/Aβ_specific_down.RData")
# save(ab_up,file="./result/Aβ_specific_up.RData")
#
# #AT8
# AT8_2_3<-read.csv("2-3 AT8 annotation.csv",header=T,row.names = 1)
# AT8_2_8<-read.csv("2-8 AT8 annotation.csv",header=T,row.names = 1)
# AT8_T4857<-read.csv("T4857 AT8 annotation.csv",header=T,row.names = 1)
# #addmeta infor
# object.2_3 <- AddMetaData(object.2_3,metadata = AT8_2_3)
# object_T4857 <- AddMetaData(object_T4857,metadata = AT8_T4857)
# object.2_8 <- AddMetaData(object.2_8,metadata = AT8_2_8)
# object_ab<-merge(object.2_3,c(object_T4857,object.2_8))
# DefaultAssay(object_ab) <- "SCT"
# ab_cat<-sort(unique(object_ab$condition))
# ab_cat<-ab_cat[-1]
# Idents(object_ab)<-object_ab$condition
# allconserved_markers<-data.frame()
# for(i in ab_cat[-1]){
#   markers <- FindConservedMarkers(object_ab, ident.1 =ab_cat[1],ident.2 =i, grouping.var = "patientID",
#                                   verbose = TRUE,
#                                   assay = "SCT",
#                                   only.pos=F)
#   markers$level<-as.character(rep(i,nrow(markers)))
#   allconserved_markers<-rbind(allconserved_markers,markers)
# }
# ab_down<-allconserved_markers[allconserved_markers$`2_8_avg_log2FC`<0,]
# ab_up<-allconserved_markers[allconserved_markers$`2_8_avg_log2FC`>0,]
# save(ab_down,file="./result/AT8_down.RData")
# save(ab_up,file="./result/AT8_up.RData")
#
# #AT8_specific
# AT8_2_3_sp<-read.csv("2-3 AT8 annotation_specific.csv",header=T,row.names = 1)
# AT8_2_8_sp<-read.csv("2-8 AT8 annotation_specific.csv",header=T,row.names = 1)
# AT8_T4857_sp<-read.csv("T4857 AT8 annotation_specific.csv",header=T,row.names = 1)
# #addmeta infor
# object.2_3 <- AddMetaData(object.2_3,metadata = AT8_2_3_sp)
# object_T4857 <- AddMetaData(object_T4857,metadata = AT8_T4857_sp)
# object.2_8 <- AddMetaData(object.2_8,metadata = AT8_2_8_sp)
# object_ab<-merge(object.2_3,c(object_T4857,object.2_8))
# DefaultAssay(object_ab) <- "SCT"
# ab_cat<-sort(unique(object_ab$condition))
# ab_cat<-ab_cat[-1]
# Idents(object_ab)<-object_ab$condition
# allconserved_markers<-data.frame()
# for(i in ab_cat[-1]){
#   markers <- FindConservedMarkers(object_ab, ident.1 =ab_cat[1],ident.2 =i, grouping.var = "patientID",
#                                   verbose = TRUE,
#                                   assay = "SCT",
#                                   only.pos=F)
#   markers$level<-as.character(rep(i,nrow(markers)))
#   allconserved_markers<-rbind(allconserved_markers,markers)
# }
# allconserved_markers_f<-allconserved_markers[allconserved_markers$minimump_p_val<0.05,]
# ab_down_f<-allconserved_markers[allconserved_markers$`2_8_avg_log2FC`<0,]
# ab_up_f<-allconserved_markers[allconserved_markers$`2_8_avg_log2FC`>0,]
# save(ab_down_f,file="./result/AT8_specific_down.RData")
# save(ab_up_f,file="./result/AT8_specific_up.RData")
#
# #Double_specific
# D_2_3_sp<-read.csv("2-3 Double annotation_specific.csv",header=T,row.names = 1)
# D_2_8_sp<-read.csv("2-8 Double annotation_specific.csv",header=T,row.names = 1)
# D_T4857_sp<-read.csv("T4857 Double annotation_specific.csv",header=T,row.names = 1)
# #addmeta infor
# object.2_3 <- AddMetaData(object.2_3,metadata = D_2_3_sp)
# object_T4857 <- AddMetaData(object_T4857,metadata = D_T4857_sp)
# object.2_8 <- AddMetaData(object.2_8,metadata = D_2_8_sp)
# object_ab<-merge(object.2_3,c(object_T4857,object.2_8))
# DefaultAssay(object_ab) <- "SCT"
# ab_cat<-sort(unique(object_ab$condition))
# ab_cat<-ab_cat[-1]
# Idents(object_ab)<-object_ab$condition
# allconserved_markers<-data.frame()
# for(i in ab_cat[-1]){
#   markers <- FindConservedMarkers(object_ab, ident.1 =ab_cat[1],ident.2 =i, grouping.var = "patientID",
#                                   verbose = TRUE,
#                                   assay = "SCT",
#                                   only.pos=F)
#   markers$level<-as.character(rep(i,nrow(markers)))
#   allconserved_markers<-rbind(allconserved_markers,markers)
# }
# allconserved_markers_f<-allconserved_markers[allconserved_markers$minimump_p_val<0.05,]
# ab_down<-allconserved_markers[allconserved_markers$`2_8_avg_log2FC`<0,]#0
# ab_up<-allconserved_markers[allconserved_markers$`2_8_avg_log2FC`>0,]
# #save(ab_down_f,file="./result/AT8_specific_down.RData")
# save(ab_up_f,file="./result/Double_specific_up.RData")
#
# #Double
# D_2_3<-read.csv("2-3 Double annotation.csv",header=T,row.names = 1)
# D_2_8<-read.csv("2-8 Double annotation.csv",header=T,row.names = 1)
# D_T4857<-read.csv("T4857 Double annotation.csv",header=T,row.names = 1)
# #addmeta infor
# object.2_3 <- AddMetaData(object.2_3,metadata = D_2_3)
# object_T4857 <- AddMetaData(object_T4857,metadata = D_T4857)
# object.2_8 <- AddMetaData(object.2_8,metadata = D_2_8)
# object_ab<-merge(object.2_3,c(object_T4857,object.2_8))
# DefaultAssay(object_ab) <- "SCT"
# ab_cat<-sort(unique(object_ab$condition))
# ab_cat<-ab_cat[-1]
# Idents(object_ab)<-object_ab$condition
# allconserved_markers<-data.frame()
# for(i in ab_cat[-1]){
#   markers <- FindConservedMarkers(object_ab, ident.1 =ab_cat[1],ident.2 =i, grouping.var = "patientID",
#                                   verbose = TRUE,
#                                   assay = "SCT",
#                                   only.pos=F)
#   markers$level<-as.character(rep(i,nrow(markers)))
#   allconserved_markers<-rbind(allconserved_markers,markers)
# }
# #0




###########findmarkers
#ab
setwd("/users/PAS1475/liuzeyi/guoqi/output/pathology")
Ab2_3 <- read.csv("2-3 Aβ annotation.csv",
                  header = T,
                  row.names = 1)
Ab2_8 <- read.csv("2-8 Aβ annotation.csv",
                  header = T,
                  row.names = 1)
T4857 <- read.csv("T4857 Aβ annotation.csv",
                  header = T,
                  row.names = 1)
colnames(Ab2_8) <- colnames(T4857)
object.2_3 <- AddMetaData(object.2_3, metadata = Ab2_3)
object_T4857 <- AddMetaData(object_T4857, metadata = T4857)
object.2_8 <- AddMetaData(object.2_8, metadata = Ab2_8)
object_ab <- merge(object.2_3, c(object_T4857, object.2_8))
DefaultAssay(object_ab) <- "SCT"
ab_cat <- sort(unique(object_ab$condition))
ab_cat <- ab_cat[-1]
Idents(object_ab) <- object_ab$condition
#findmarkers
allconserved_markers <- data.frame()
for (i in ab_cat[-1]) {
  markers <-
    FindMarkers(
      object_ab,
      ident.1 = ab_cat[1],
      ident.2 = i,
      group.by = "condition",
      verbose = TRUE,
      assay = "SCT",
      only.pos = F
    )
  markers <- cbind(rownames(markers), markers)
  markers$level <- as.character(rep(i, nrow(markers)))
  allconserved_markers <- rbind(allconserved_markers, markers)
}
allconserved_markers_f <-
  allconserved_markers[allconserved_markers$p_val_adj < 0.05, ]
rownames(allconserved_markers_f) <- NULL
colnames(allconserved_markers_f)[1] <- c("Marker")
ab_down <-
  allconserved_markers_f[allconserved_markers_f$avg_log2FC > 0, ]
ab_up <- allconserved_markers_f[allconserved_markers_f$avg_log2FC < 0, ]
save(ab_down, file = "./result/Aβ_down_findmarker.RData")
save(ab_up, file = "./result/Aβ_up_findmarker.RData")
##ab_specific
setwd("/users/PAS1475/liuzeyi/guoqi/output/pathology/")
Ab2_3_sp <-
  read.csv("2-3 Aβ annotation_specific.csv",
           header = T,
           row.names = 1)
Ab2_8_sp <-
  read.csv("2-8 Aβ annotation_specific.csv",
           header = T,
           row.names = 1)
Ab4857_sp <-
  read.csv("T4857 Aβ annotation_specific.csv",
           header = T,
           row.names = 1)
#addmeta infor
object.2_3 <- AddMetaData(object.2_3, metadata = Ab2_3_sp)
object_T4857 <- AddMetaData(object_T4857, metadata = Ab4857_sp)
object.2_8 <- AddMetaData(object.2_8, metadata = Ab2_8_sp)
object_ab <- merge(object.2_3, c(object_T4857, object.2_8))
DefaultAssay(object_ab) <- "SCT"
Idents(object_ab) <- object_ab$condition
ab_cat <- sort(unique(object_ab$condition))
ab_cat <- ab_cat[-1]
allconserved_markers <- data.frame()
for (i in ab_cat[-1]) {
  markers <-
    FindMarkers(
      object_ab,
      ident.1 = ab_cat[1],
      ident.2 = i,
      group.by = "condition",
      verbose = TRUE,
      assay = "SCT",
      only.pos = F
    )
  markers <- cbind(rownames(markers), markers)
  markers$level <- as.character(rep(i, nrow(markers)))
  allconserved_markers <- rbind(allconserved_markers, markers)
}
allconserved_markers_f <-
  allconserved_markers[allconserved_markers$p_val_adj < 0.05, ]
rownames(allconserved_markers_f) <- NULL
colnames(allconserved_markers_f)[1] <- c("Marker")
ab_down <-
  allconserved_markers_f[allconserved_markers_f$avg_log2FC > 0, ]
ab_up <- allconserved_markers_f[allconserved_markers_f$avg_log2FC < 0, ]
save(ab_down, file = "./result/Aβ_down_findmarker_specific.RData")
save(ab_up, file = "./result/Aβ_up_findmarker_specific.RData")

#AT8
AT8_2_3 <- read.csv("2-3 AT8 annotation.csv",
                    header = T,
                    row.names = 1)
AT8_2_8 <- read.csv("2-8 AT8 annotation.csv",
                    header = T,
                    row.names = 1)
AT8_T4857 <-
  read.csv("T4857 AT8 annotation.csv",
           header = T,
           row.names = 1)
#addmeta infor
object.2_3 <- AddMetaData(object.2_3, metadata = AT8_2_3)
object_T4857 <- AddMetaData(object_T4857, metadata = AT8_T4857)
object.2_8 <- AddMetaData(object.2_8, metadata = AT8_2_8)
object_ab <- merge(object.2_3, c(object_T4857, object.2_8))
DefaultAssay(object_ab) <- "SCT"
Idents(object_ab) <- object_ab$condition
ab_cat <- sort(unique(object_ab$condition))
ab_cat <- ab_cat[-1]
allconserved_markers <- data.frame()
for (i in ab_cat[-1]) {
  markers <-
    FindMarkers(
      object_ab,
      ident.1 = ab_cat[1],
      ident.2 = i,
      group.by = "condition",
      verbose = TRUE,
      assay = "SCT",
      only.pos = F
    )
  markers <- cbind(rownames(markers), markers)
  markers$level <- as.character(rep(i, nrow(markers)))
  allconserved_markers <- rbind(allconserved_markers, markers)
}
allconserved_markers_f <-
  allconserved_markers[allconserved_markers$p_val_adj < 0.05, ]
rownames(allconserved_markers_f) <- NULL
colnames(allconserved_markers_f)[1] <- c("Marker")
ab_down <-
  allconserved_markers_f[allconserved_markers_f$avg_log2FC > 0, ]
ab_up <- allconserved_markers_f[allconserved_markers_f$avg_log2FC < 0, ]
save(ab_down, file = "./result/AT8_down_findmarker.RData")
save(ab_up, file = "./result/AT8_up_findmarker.RData")

#AT8_specific
AT8_2_3_sp <-
  read.csv("2-3 AT8 annotation_specific.csv",
           header = T,
           row.names = 1)
AT8_2_8_sp <-
  read.csv("2-8 AT8 annotation_specific.csv",
           header = T,
           row.names = 1)
AT8_T4857_sp <-
  read.csv("T4857 AT8 annotation_specific.csv",
           header = T,
           row.names = 1)
#addmeta infor
object.2_3 <- AddMetaData(object.2_3, metadata = AT8_2_3_sp)
object_T4857 <- AddMetaData(object_T4857, metadata = AT8_T4857_sp)
object.2_8 <- AddMetaData(object.2_8, metadata = AT8_2_8_sp)
object_ab <- merge(object.2_3, c(object_T4857, object.2_8))
DefaultAssay(object_ab) <- "SCT"
ab_cat <- sort(unique(object_ab$condition))
ab_cat <- ab_cat[-1]
Idents(object_ab) <- object_ab$condition
allconserved_markers <- data.frame()
for (i in ab_cat[-1]) {
  markers <-
    FindMarkers(
      object_ab,
      ident.1 = ab_cat[1],
      ident.2 = i,
      group.by = "condition",
      verbose = TRUE,
      assay = "SCT",
      only.pos = F
    )
  markers <- cbind(rownames(markers), markers)
  markers$level <- as.character(rep(i, nrow(markers)))
  allconserved_markers <- rbind(allconserved_markers, markers)
}
allconserved_markers_f <-
  allconserved_markers[allconserved_markers$p_val_adj < 0.05, ]
rownames(allconserved_markers_f) <- NULL
colnames(allconserved_markers_f)[1] <- c("Marker")
ab_down <-
  allconserved_markers_f[allconserved_markers_f$avg_log2FC > 0, ]
ab_up <- allconserved_markers_f[allconserved_markers_f$avg_log2FC < 0, ]
save(ab_down, file = "./result/AT8_down_findmarker_specific.RData")
save(ab_up, file = "./result/AT8_up_findmarker_specific.RData")

#Double_specific
D_2_3_sp <-
  read.csv("2-3 Double annotation_specific.csv",
           header = T,
           row.names = 1)
D_2_8_sp <-
  read.csv("2-8 Double annotation_specific.csv",
           header = T,
           row.names = 1)
D_T4857_sp <-
  read.csv("T4857 Double annotation_specific.csv",
           header = T,
           row.names = 1)
#addmeta infor
object.2_3 <- AddMetaData(object.2_3, metadata = D_2_3_sp)
object_T4857 <- AddMetaData(object_T4857, metadata = D_T4857_sp)
object.2_8 <- AddMetaData(object.2_8, metadata = D_2_8_sp)
object_ab <- merge(object.2_3, c(object_T4857, object.2_8))

DefaultAssay(object_ab) <- "SCT"
ab_cat <- sort(unique(object_ab$condition))
ab_cat <- ab_cat[-1]
Idents(object_ab) <- object_ab$condition
allconserved_markers <- data.frame()
for (i in ab_cat[-1]) {
  markers <-
    FindMarkers(
      object_ab,
      ident.1 = ab_cat[1],
      ident.2 = i,
      group.by = "condition",
      verbose = TRUE,
      assay = "SCT",
      only.pos = F
    )
  markers <- cbind(rownames(markers), markers)
  markers$level <- as.character(rep(i, nrow(markers)))
  allconserved_markers <- rbind(allconserved_markers, markers)
}
allconserved_markers_f <-
  allconserved_markers[allconserved_markers$p_val_adj < 0.05, ]
rownames(allconserved_markers_f) <- NULL
colnames(allconserved_markers_f)[1] <- c("Marker")
ab_down <-
  allconserved_markers_f[allconserved_markers_f$avg_log2FC > 0, ]
ab_up <- allconserved_markers_f[allconserved_markers_f$avg_log2FC < 0, ]
save(ab_down, file = "./result/Double_down_findmarker_specific.RData")
save(ab_up, file = "./result/Double_up_findmarker_specific.RData")

#Double
D_2_3 <- read.csv("2-3 Double annotation.csv",
                  header = T,
                  row.names = 1)
D_2_8 <- read.csv("2-8 Double annotation.csv",
                  header = T,
                  row.names = 1)
D_T4857 <-
  read.csv("T4857 Double annotation.csv",
           header = T,
           row.names = 1)
#addmeta infor
object.2_3 <- AddMetaData(object.2_3, metadata = D_2_3)
object_T4857 <- AddMetaData(object_T4857, metadata = D_T4857)
object.2_8 <- AddMetaData(object.2_8, metadata = D_2_8)
object_ab <- merge(object.2_3, c(object_T4857, object.2_8))
DefaultAssay(object_ab) <- "SCT"
ab_cat <- sort(unique(object_ab$condition))
ab_cat <- ab_cat[-1]
Idents(object_ab) <- object_ab$condition
allconserved_markers <- data.frame()
for (i in ab_cat[-1]) {
  markers <-
    FindMarkers(
      object_ab,
      ident.1 = ab_cat[1],
      ident.2 = i,
      group.by = "condition",
      verbose = TRUE,
      assay = "SCT",
      only.pos = F
    )
  markers <- cbind(rownames(markers), markers)
  markers$level <- as.character(rep(i, nrow(markers)))
  allconserved_markers <- rbind(allconserved_markers, markers)
}
allconserved_markers_f <-
  allconserved_markers[allconserved_markers$p_val_adj < 0.05, ]
rownames(allconserved_markers_f) <- NULL
colnames(allconserved_markers_f)[1] <- c("Marker")
ab_down <-
  allconserved_markers_f[allconserved_markers_f$avg_log2FC > 0, ]
ab_up <- allconserved_markers_f[allconserved_markers_f$avg_log2FC < 0, ]
save(ab_down, file = "./result/Double_down_findmarker.RData")
save(ab_up, file = "./result/Double_up_findmarker.RData")


#-----------WFS1,AT8
#Double_specific
setwd("/users/PAS1475/liuzeyi/guoqi/output/pathology")
D_2_3_W_AT8_sp <-
  read.csv(
    "2-3 Double_W_AT8 annotation_specific.csv",
    header = T,
    row.names = 1
  )
D_2_8_W_AT8_sp <-
  read.csv(
    "2-8 Double_W_AT8 annotation_specific.csv",
    header = T,
    row.names = 1
  )
D_T4857_W_AT8_sp <-
  read.csv(
    "T4857 Double_W_AT8 annotation_specific.csv",
    header = T,
    row.names = 1
  )
#addmeta infor
object.2_3 <- AddMetaData(object.2_3, metadata = D_2_3_W_AT8_sp)
object_T4857 <-
  AddMetaData(object_T4857, metadata = D_T4857_W_AT8_sp)
object.2_8 <- AddMetaData(object.2_8, metadata = D_2_8_W_AT8_sp)
object_W_AT8_sp <- merge(object.2_3, c(object_T4857, object.2_8))

DefaultAssay(object_W_AT8_sp) <- "SCT"
ab_cat <- sort(unique(object_W_AT8_sp$condition))
ab_cat <- ab_cat[-1]
Idents(object_W_AT8_sp) <- object_W_AT8_sp$condition
allconserved_markers <- data.frame()
for (i in ab_cat[-4]) {
  markers <-
    FindMarkers(
      object_W_AT8_sp,
      ident.1 = ab_cat[4],
      ident.2 = i,
      group.by = "condition",
      verbose = TRUE,
      assay = "SCT",
      only.pos = F
    )
  markers <- cbind(rownames(markers), markers)
  markers$level <- as.character(rep(i, nrow(markers)))
  allconserved_markers <- rbind(allconserved_markers, markers)
}
allconserved_markers_f <-
  allconserved_markers[allconserved_markers$p_val_adj < 0.05, ]
rownames(allconserved_markers_f) <- NULL
colnames(allconserved_markers_f)[1] <- c("Marker")
ab_down <-
  allconserved_markers_f[allconserved_markers_f$avg_log2FC > 0, ]
ab_up <- allconserved_markers_f[allconserved_markers_f$avg_log2FC < 0, ]
save(ab_down, file = "./result/Double_down_findmarker_W_AT8_specific.RData")
save(ab_up, file = "./result/Double_up_findmarker_W_AT8_specific.RData")

#Double
D_2_3 <-
  read.csv("2-3 Double_W_AT8 annotation.csv",
           header = T,
           row.names = 1)
D_2_8 <-
  read.csv("2-8 Double_W_AT8 annotation.csv",
           header = T,
           row.names = 1)
D_T4857 <-
  read.csv("T4857 Double_W_AT8 annotation.csv",
           header = T,
           row.names = 1)
#addmeta infor
object.2_3 <- AddMetaData(object.2_3, metadata = D_2_3)
object_T4857 <- AddMetaData(object_T4857, metadata = D_T4857)
object.2_8 <- AddMetaData(object.2_8, metadata = D_2_8)
object_ab <- merge(object.2_3, c(object_T4857, object.2_8))
DefaultAssay(object_ab) <- "SCT"
ab_cat <- sort(unique(object_ab$condition))
ab_cat <- ab_cat[-1]
Idents(object_ab) <- object_ab$condition
allconserved_markers <- data.frame()
for (i in ab_cat[-4]) {
  markers <-
    FindMarkers(
      object_ab,
      ident.1 = ab_cat[4],
      ident.2 = i,
      group.by = "condition",
      verbose = TRUE,
      assay = "SCT",
      only.pos = F
    )
  markers <- cbind(rownames(markers), markers)
  markers$level <- as.character(rep(i, nrow(markers)))
  allconserved_markers <- rbind(allconserved_markers, markers)
}
allconserved_markers_f <-
  allconserved_markers[allconserved_markers$p_val_adj < 0.05, ]
rownames(allconserved_markers_f) <- NULL
colnames(allconserved_markers_f)[1] <- c("Marker")
ab_down <-
  allconserved_markers_f[allconserved_markers_f$avg_log2FC > 0, ]
ab_up <- allconserved_markers_f[allconserved_markers_f$avg_log2FC < 0, ]
save(ab_down, file = "./result/Double_down_findmarker_W_AT8.RData")
save(ab_up, file = "./result/Double_up_findmarker_W_AT8.RData")

#GO
RunPathway_human <- function(Seurat.DEGs = NULL) {
  library(org.Hs.eg.db)
  library(clusterProfiler)
  genes.use <- Seurat.DEGs
  #Obtain a GO object
  GO_pathway <-
    enrichGO(
      gene = genes.use,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      keyType = "SYMBOL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.4
    )
  
  if (is.null(GO_pathway)) {
    GO_simplied_res <-
      as.data.frame(cbind(NO_results = "NO_result_for_GO"))
  } else{
    #GO_res <- simplify(GO_pathway)
    GO_simplied_res <-
      GO_pathway@result[GO_pathway@result$pvalue < 0.05, ]
    dim(GO_simplied_res)
    if (dim(GO_simplied_res)[1] == 0) {
      GO_simplied_res <-
        as.data.frame(cbind(NO_results = "NO_result_for_GO"))
    }
  }
  return(GO_simplied_res)
}
##construct go file
setwd("/users/PAS1475/liuzeyi/guoqi/output/pathology/result/")
all_file <- list.files(pattern = "*find")
list <- list()
for (i in 1:length(all_file)) {
  a <- load(all_file[i])
  data <- eval(parse(text = a))
  list[[i]] <- data
}
name <- strsplit(all_file, split = "\\.")
names <- c()
for (i in 1:length(name)) {
  a <- name[[i]][1]
  names <- c(names, a)
}
names(list) <- names
#classfication based on level
list1 <- list()
list2 <- list()
list3 <- list()
for (i in 1:length(list)) {
  temp <- list[[i]]
  if (any(grepl(pattern = "1", unique(temp$level)))) {
    l1 <-
      temp$Marker[which(temp$level == unique(temp$level)[grep("1", unique(temp$level))])]
    list1[[i]] <- RunPathway_human(l1)
  }
  if (any(grepl(pattern = "2", unique(temp$level)))) {
    l2 <-
      temp$Marker[which(temp$level == unique(temp$level)[grep("2", unique(temp$level))])]
    list2[[i]] <- RunPathway_human(l2)
  }
  if (any(grepl(pattern = "3", unique(temp$level)))) {
    l3 <-
      temp$Marker[which(temp$level == unique(temp$level)[grep("3", unique(temp$level))])]
    list3[[i]] <- RunPathway_human(l3)
  }
  print(i)
}
names(list1) <- names[1:length(list1)]
names(list2) <- names[1:length(list2)]
names(list3) <- names[1:length(list3)]
save(list1, file = "level1.rdata")
save(list2, file = "level2.rdata")
save(list3, file = "level3.rdata")


# GO plot for level3
#function
GOenrichment_human <- function(Seurat.DEGs = NULL) {
  library(org.Hs.eg.db)
  library(clusterProfiler)
  genes.use <- Seurat.DEGs
  #Obtain a GO object
  GO_pathway <-
    enrichGO(
      gene = genes.use,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      keyType = "SYMBOL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.4
    )
  return(GO_pathway)
}

#load data
setwd("/users/PAS1475/liuzeyi/guoqi/output/pathology/result/")
all_file <- list.files(pattern = "*specific.RData")
list <- list()
for (i in 1:length(all_file)) {
  a <- load(all_file[i])
  data <- eval(parse(text = a))
  list[[i]] <- data
}
name <- strsplit(all_file, split = "\\.")
names <- c()
for (i in 1:length(name)) {
  a <- name[[i]][1]
  names <- c(names, a)
}
names(list) <- names
list3 <- list()
for (i in 1:length(list)){   
  temp <- list[[i]]
  if (any(grepl(pattern = "3", unique(temp$level)))) {
    l3 <-
      temp$Marker[which(temp$level == unique(temp$level)[grep("3", unique(temp$level))])]
    list3[[i]] <- GOenrichment_human(l3)
  }
  print(i)
}
names(list3) <- names[1:length(list3)]

#target gene
setwd("/users/PAS1475/liuzeyi/guoqi/output/pathology/GO for level3")
level3<-list()
for(i in 1:6){
  level3[[i]] <- read_excel("No.11 Go pathway level3 pathology.xlsx", sheet = i)
}
names(level3)<-excel_sheets(path = "No.11 Go pathway level3 pathology.xlsx")
level3_updata<-list(list3[[3]],list3[[4]],list3[[1]],list3[[2]],list3[[5]],list3[[7]])
show<-c(level3[[1]]$Description)
dotplot(level3_updata[[1]],showCategory=show)

for(i in 1:6){
  show<-c(level3[[i]]$Description)
  g<-dotplot(level3_updata[[i]],showCategory=show,font.size=20)
  ggsave(
    plot = g,
    filename = paste0(names(level3)[i],"_specific_enrichment_NO11.tiff"),
    device = "tiff",
    dpi = 150,
    width = 11,
    height = 12,
    units = "in"
  )
}
