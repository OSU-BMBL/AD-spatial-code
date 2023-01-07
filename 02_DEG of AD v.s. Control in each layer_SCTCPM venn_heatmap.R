library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
#BiocManager::install("MAST")
library("MAST")
dyn.load(
  "/usr/lib/jvm/java-1.6.0-openjdk-1.6.0.41.x86_64/jre/lib/amd64/server/libjvm.so"
)
library(xlsx, lib.loc = "/fs/scratch/PCON0022/liyang/lib")
options(future.globals.maxSize = 4000 * 1024 ^ 2)
setwd("/users/PAS1475/liuzeyi/guoqi/output/")
load("object.18_64.RData")
load("object.2_5.RData")
load("object.2_8.RData")
load("object_T4857.RData")
load("object.1_1.RData")
load("object.2_3.RData")
#merge
object_merge <-
  merge(object.18_64,
        y = list(object.2_5, object.2_8, object_T4857,
                 object.1_1, object.2_3))
# #normalize-original spatial data
# DefaultAssay(object_merge) = "Spatial"
# object_merge <- NormalizeData(object_merge)
# #variablefeatures
# my.marker.AD <-
#   FindMarkers(
#     object_merge,
#     ident.1 = "AD",
#     ident.2 = "control",
#     group.by = "category",
#     only.pos = F
#   )
# #Spatial Control-158, AD-78
# my.marker.AD_spatial<-my.marker.AD[my.marker.AD$p_val_adj < 0.05,]
# my.marker.control_spatial<-my.marker.AD_spatial[my.marker.AD_spatial$avg_log2FC<0,]
# my.marker.AD_spatial<-my.marker.AD_spatial[my.marker.AD_spatial$avg_log2FC>0,]
# write.csv(my.marker.AD_spatial,"marker.AD_spatial.csv")
# write.csv(my.marker.control_spatial,"marker.control_spatial.csv")
# my.marker.AD_spatial <- read.csv("./AD_Control_DEG/marker.AD_spatial.csv", header = T, row.names = 1)
# my.marker.control_spatial <- read.csv("./AD_Control_DEG/marker.control_spatial.csv", header = T, row.names = 1)
#SCT
DefaultAssay(object_merge) = "SCT"
my.marker.control <-
  FindMarkers(
    object_merge,
    ident.1 = "AD",
    ident.2 = "control",
    group.by = "category",
    only.pos = F,
  )
#
my.marker.control_p <-
  my.marker.control[my.marker.control$p_val_adj < 0.05, ]
my.marker.AD_sct <-
  my.marker.control_p[my.marker.control_p$avg_log2FC > 0, ]
my.marker.control_sct <-
  my.marker.control_p[my.marker.control_p$avg_log2FC < 0, ]
write.csv(my.marker.control_sct,
          "./AD_Control_DEG/marker.control_sct.csv")
write.csv(my.marker.AD_sct, "./AD_Control_DEG/marker.AD_sct.csv")


#plot
#venn
sct_ <- c(rownames(my.marker.control), rownames(my.marker.AD_sct))
spatial_gene <-
  c(rownames(my.marker.AD_spatial),
    rownames(my.marker.control_spatial))
venn <-
  list(
    sct_up = rownames(my.marker.AD_sct),
    sct_down = rownames(my.marker.control),
    cmp_up = rownames(my.marker.AD_spatial),
    cmp_down = rownames(my.marker.control_spatial)
  )

ve_up <-
  list(sct_up = rownames(my.marker.AD_sct),
       cmp_up = rownames(my.marker.AD_spatial))
ve_down <-
  list(
    sct_up = rownames(my.marker.control_sct),
    cmp_up = rownames(my.marker.control_spatial)
  )
#install.packages("VennDiagram")
library(VennDiagram)
dif_vennplot <-
  venn.diagram(
    venn,
    filename = "picture/cad_control.vennplot.tiff",
    col = "transparent",
    fill = c("orange", "blue", "grey", "yellow"),
    cex = 1,
    cat.cex = 1,
    rotation.degree = 0,
    main = "sct_cpm dif gene",
    main.cex = 2,
    sub.cex = 1,
    alpha = 0.50
  )
#up
dif_vennplot_up <-
  venn.diagram(
    ve_up,
    filename = "picture/up.vennplot.tiff",
    col = "transparent",
    fill = c("blue", "yellow"),
    cex = 1,
    cat.cex = 1,
    rotation.degree = 0,
    main = "sct_cpm up_dif gene",
    main.cex = 2,
    sub.cex = 1,
    alpha = 0.50
  )

#down
dif_vennplot_down <-
  venn.diagram(
    ve_down,
    filename = "picture/down.vennplot.tiff",
    col = "transparent",
    fill = c("orange", "blue"),
    cex = 1,
    cat.cex = 1,
    rotation.degree = 0,
    main = "sct_cpm down_dif gene",
    main.cex = 2,
    sub.cex = 1,
    alpha = 0.50
  )


category <-
  c(rep("AD", nrow(my.marker.AD.clean)), rep("control", nrow(my.marker.control.clean)))
categroy <- as.factor(category)
p <-
  -log(c(my.marker.AD.clean[, 5], my.marker.control.clean[, 5]))#?-log-inf
box_ad_control <- data.frame(categroy = categroy, p = p)
g <-
  ggplot(box_ad_control) + aes(x = categroy, y = p) + geom_boxplot() +
  geom_point(aes(x = category, y = p)) +
  labs(title = "Up_gene") + ylab("-log(p_adj)") + ylim(c(0, 300)) +
  theme(
    title = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )
ggsave(
  plot = g,
  filename = "/users/PAS1475/liuzeyi/guoqi/output/picture/AD_Control_box.tiff",
  device = "tiff",
  dpi = 150,
  width = 12,
  height = 9,
  units = "in"
)
a <- read.csv("allconserved_markers.csv",
              header = T,
              row.names = 1)

##DEG of control-AD in each layer
DefaultAssay(object_merge) = "SCT"
level <- unique(object_merge$Layer)
all_marker <- list()
for (i in 1:7) {
  object_layer <- subset(object_merge, subset = Layer == level[i])
  my.marker <-
    FindMarkers(
      object_layer,
      ident.1 = "AD",
      ident.2 = "control",
      group.by = "category",
      only.pos = F,
    )
  all_marker[[i]] <- my.marker
}

all_marker_up <- list()
for (i in 1:7) {
  marker_up <- all_marker[[i]][all_marker[[i]]$avg_log2FC > 0, ]
  marker_up <- marker_up[marker_up$p_val_adj < 0.05, ]
  all_marker_up[[i]] <- marker_up
  names(all_marker_up)[i] <- paste(level[i], "up", sep = "_")
}
all_marker_down <- list()
for (i in 1:7) {
  marker_down <- all_marker[[i]][all_marker[[i]]$avg_log2FC < 0, ]
  marker_down <- marker_down[marker_down$p_val_adj < 0.05, ]
  all_marker_down[[i]] <- marker_down
  names(all_marker_down)[i] <- paste(level[i], "down", sep = "_")
}

all_marker_output <- list()
for (i in 1:14) {
  if (i %% 2 == 0) {
    marker <- all_marker_down[[i / 2]]
    name <- names(all_marker_down)[[i / 2]]
  }
  else{
    marker <- all_marker_up[[ceiling(i / 2)]]
    name <- names(all_marker_up)[[ceiling(i / 2)]]
  }
  all_marker_output[[i]] <- marker
  names(all_marker_output)[i] <- name
}
#write xlsx file
xlsx::write.xlsx(
  all_marker_output[[1]],
  file = "/users/PAS1475/liuzeyi/guoqi/output/AD_Control_DEG/Layer_ad_control_markers.xlsx",
  sheetName = names(all_marker_output)[1],
  row.names = TRUE
)
for (sheet.idx in 2:14) {
  print(sheet.idx)
  xlsx::write.xlsx(
    all_marker_output[[sheet.idx]],
    file = "/users/PAS1475/liuzeyi/guoqi/output/AD_Control_DEG/Layer_ad_control_markers.xlsx",
    sheetName = names(all_marker_output)[sheet.idx],
    append = T,
    row.names = TRUE
  )
}
#save
save(all_marker_output, file = "/users/PAS1475/liuzeyi/guoqi/output/AD_Control_DEG/markers in each layer.RData")
load("/users/PAS1475/liuzeyi/guoqi/output/AD_Control_DEG/markers in each layer.RData")











#heatmap
load("/users/PAS1475/liuzeyi/guoqi/output/sample6_clusters.RData")
#express data
setwd("/users/PAS1475/liuzeyi/guoqi/output/AD_Control_DEG")
deg <- read.csv("Layer marker for heatmap.csv")
expr_data <- as.data.frame(sample6.combined@assays$SCT@data[deg$Marker,])
expr_data<-rbind(expr_data,layer=sample6.combined$Layer)
expr_data<-as.data.frame(t(expr_data))
temp_data<-data.frame(rep(0,7))
for(j in 1:21){
  gene_data<-c()
  for(i in sort(unique(expr_data$layer))){
    temp_gene<-mean(as.numeric(expr_data[which(expr_data$layer==i),j]))
    gene_data<-c(gene_data,temp_gene)
  }
  temp_data<-cbind(temp_data,gene_data)
}
temp_data<-temp_data[,-1]
rownames(temp_data)<-sort(unique(expr_data$layer))
colnames(temp_data)<-deg$Marker
heatmap_data<-t(temp_data)
#plot

my_colour = list(
  gene_cat = c( Layer_1 = brewer.pal(n = 7, name = "Dark2")[1],
                Layer_2 = brewer.pal(n = 7, name = "Dark2")[2],
                Layer_3 = brewer.pal(n = 7, name = "Dark2")[3],
                Layer_4 = brewer.pal(n = 7, name = "Dark2")[4],
                Layer_5 = brewer.pal(n = 7, name = "Dark2")[5],
                Layer_6 = brewer.pal(n = 7, name = "Dark2")[6],
                White_Matter = brewer.pal(n = 7, name = "Dark2")[7]))
gene_cat_val <- as.numeric(table(deg$Layer))
gene_anno <- data.frame(gene_cat = c(rep("Layer_1",as.numeric(table(deg$Layer)[1])),
                                     rep("Layer_2",as.numeric(table(deg$Layer)[2])),
                                     rep("Layer_4",as.numeric(table(deg$Layer)[3])),
                                     rep("Layer_5",as.numeric(table(deg$Layer)[4])),
                                     rep("Layer_6",as.numeric(table(deg$Layer)[5])),
                                     rep("White_Matter",as.numeric(table(deg$Layer)[6]))
) )
library("pheatmap")
p <- pheatmap(heatmap_data,
              gaps_row = cumsum(gene_cat_val),
              cluster_rows = F,
              cluster_cols = F,
              scale = "row",
              border_color = NA,
              annotation_colors = my_colour, 
              #annotation_row = gene_anno,
              labels_col =rownames(temp_data),
              color = colorRampPalette(c("blue","white","red"))(100))

ggsave(
  plot = p,
  filename = "/users/PAS1475/liuzeyi/guoqi/output/picture/integrating landscape/heatmap_revised.tiff",
  device = "tiff",
  dpi = 150,
  width = 7,
  height = 11,
  units = "in"
)

#Enrichment for DEGs between AD and control
RunPathway_human <- function(Seurat.DEGs = NULL) {
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ReactomePA)
  genes.use <- rownames(Seurat.DEGs)
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
  #convert gene type for pathway annotation method-enrichPathway based on reactome and kegg
  gene.convert <-
    bitr(
      genes.use,
      fromType = "SYMBOL",
      toType = c("ENSEMBL", "ENTREZID"),
      OrgDb = org.Hs.eg.db
    )
  reactome_analysis <-
    enrichPathway(
      gene = gene.convert$ENTREZID,
      pvalueCutoff = 0.05,
      readable = TRUE
    )
  if (is.null(reactome_analysis)) {
    reactome_res <-
      as.data.frame(cbind(NO_results = "NO_result_for_reactome"))
  } else{
    reactome_res <-
      reactome_analysis@result[reactome_analysis@result$pvalue < 0.05, ]
    dim(reactome_res)
    if (dim(reactome_res)[1] == 0) {
      reactome_res <-
        as.data.frame(cbind(NO_results = "NO_result_for_reactome"))
    }
  }
  
  KEGG_pathway <-
    enrichKEGG(
      gene = gene.convert$ENTREZID,
      organism = "hsa",
      keyType = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  #KEGG_pathway@result[1:5,]
  if (is.null(KEGG_pathway)) {
    KEGG_res <- as.data.frame(cbind(NO_results = "NO_result_for_KEGG"))
  } else{
    KEGG_res <- KEGG_pathway@result[KEGG_pathway@result$pvalue < 0.05, ]
    dim(KEGG_res)
    if (dim(KEGG_res)[1] == 0) {
      KEGG_res <- as.data.frame(cbind(NO_results = "NO_result_for_KEGG"))
    } else{
      for (j in 1:nrow(KEGG_res)) {
        tmp_name <- KEGG_res$geneID[j]
        tmp_name <- unlist(strsplit(tmp_name, "/"))
        index_entizID <- gene.convert$ENTREZID %in% tmp_name
        converted_geneSymbol <- gene.convert$SYMBOL[index_entizID]
        final_name <- paste0(converted_geneSymbol, collapse = "/")
        KEGG_res$geneID[j] <- final_name
      }
    }
  }
  Go_KEGG.bind <-
    gdata::cbindX(GO_simplied_res, KEGG_res, reactome_res)
  final.data <- gdata::cbindX(Seurat.DEGs, Go_KEGG.bind)
  return(final.data)
}

setwd("/users/PAS1475/liuzeyi/guoqi/output/")
control<-read.csv("./AD_Control_DEG/marker.control_sct.csv",header=T,row.names = 1)
ad<-read.csv("./AD_Control_DEG/marker.AD_sct.csv",header=T,row.names = 1)
ad_enrichment<-RunPathway_human(ad)
control_enrichment<-RunPathway_human(control)
control_enrichment<-control_enrichment[,-c(1:5)]
rownames(control_enrichment)<-NULL
ad_enrichment<-ad_enrichment[,-c(1:5)]
rownames(ad_enrichment)<-NULL
rownames(control_enrichment)<-NULL
write.csv(ad_enrichment,"./AD_Control_DEG/marker.ad_enrichment.csv")
write.csv(control_enrichment,"./AD_Control_DEG/marker.control_enrichment.csv")


#GO plot for given gene list
#function
GOenrichment_human_object <- function(Seurat.DEGs = NULL) {
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
dotplot_WGCNA<-function(module_gene,target_data){
  go_object<-GOenrichment_human_object(module_gene)
  target_description<-target_data$Description
  plot<-dotplot(go_object,showCategory=target_description,font.size=20)
}
#load target description
library(xlsx)
setwd("/users/PAS1475/liuzeyi/guoqi/output/AD_Control_DEG/GO plot")
up_deg <- read_excel("No.4, AD vs nonAD GO for dot plot.xlsx", sheet = 1)
down_deg<-read_excel("No.4, AD vs nonAD GO for dot plot.xlsx", sheet = 2)
up_deg<-up_deg[,-c(1,11:ncol(up_deg))]
down_deg<-down_deg[,-c(1,11:ncol(down_deg))]
data_target<-list(up_deg,down_deg)
#load original data
setwd("/users/PAS1475/liuzeyi/guoqi/output/")
control<-read.csv("./AD_Control_DEG/marker.control_sct.csv",header=T,row.names = 1)
ad<-read.csv("./AD_Control_DEG/marker.AD_sct.csv",header=T,row.names = 1)
data<-list(ad,control)
plot_list_ad_control<-list()
for (i in 1:2) {
  temp <- data[[i]]
  module_gene <- rownames(temp)
  target_data<-data_target[[i]]
  plot_list_ad_control[[i]] <- dotplot_WGCNA(module_gene,target_data)
}
up<-plot_list_ad_control[[1]]+theme(legend.key.size = unit(1.5, 'cm'),
                                                          legend.title = element_text(size=18), #change legend title font size
                                                          legend.text = element_text(size=18))
down<-plot_list_ad_control[[2]]+theme(legend.key.size = unit(1.5, 'cm'),
                                      legend.title = element_text(size=18), #change legend title font size
                                      legend.text = element_text(size=18))
ggsave(plot = up,
       filename = "AD_Control_up_enrichment_NO4.tiff",
       device = "tiff",
       dpi = 150,
       width = 13,
       height = 12,
       units = "in")
ggsave(
  plot = down,
  filename = "AD_Control_down_enrichment_NO4.tiff",
  device = "tiff",
  dpi = 150,
  width = 12,
  height = 10,
  units = "in"
)


#GO for DEGs in each layer
setwd("/users/PAS1475/liuzeyi/guoqi/output/AD_Control_DEG/GO for DEGs in each layer")
GOenrichment_human <- function(Seurat.DEGs = NULL) {
  library(org.Hs.eg.db)
  library(clusterProfiler)
  genes.use <- Seurat.DEGs$'V1'
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
#load data
list_go<-list()
library(xlsx)
for(i in 1:2){
  list_go[[i]] <- read_excel("No.15, Layer specific AD vs nonAD GO genes.xlsx", sheet = i)
}
names(list_go)<-excel_sheets(path = "No.15, Layer specific AD vs nonAD GO genes.xlsx")
#up_enrichment
list_up<-list()
up<-list_go[[1]]
for(j in 1:ncol(up)){
  temp_data<-up[,j]
  colnames(temp_data)<-'V1'
  data_go<-GOenrichment_human(temp_data)
  list_up[[j]]<-data_go
}
names(list_up)<-colnames(up)
#down_enrichment
list_down<-list()
down<-list_go[[2]]
for(j in 1:ncol(down)){
  temp_data<-down[,j]
  colnames(temp_data)<-'V1'
  data_go<-GOenrichment_human(temp_data)
  list_down[[j]]<-data_go
}
names(list_down)<-colnames(down)
#save file
dyn.load(
  "/usr/lib/jvm/java-1.6.0-openjdk-1.6.0.41.x86_64/jre/lib/amd64/server/libjvm.so"
)
qsave(list_up,"list_up.qs")
qsave(list_down,"list_down.qs")
write.xlsx2
library(xlsx, lib.loc = "/fs/scratch/PCON0022/liyang/lib")
write.xlsx2(list_up[[1]], "x.xlsx", sheetName = names(list_up)[1],
              col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx2(list_down[[1]], paste0(names(list_down)[i],"-GO.xlsx"), sheetName = names(list_down)[i],
              col.names = TRUE, row.names = T, append = FALSE)
for(j in 2:7){
  write.xlsx2(list_up[[j]], "layer specific AD vs control genes-GO.xlsx", sheetName = names(list_up)[j],
              col.names = TRUE, row.names = T, append = TRUE)
 }
library(qs)
