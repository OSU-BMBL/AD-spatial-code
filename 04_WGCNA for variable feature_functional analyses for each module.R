install.packages("BiocManager")
BiocManager::install("WGCNA")
BiocManager::install("flashClust")
BiocManager::install("methylumi")#fail
BiocManager::install("ReactomePA")#fail
library("WGCNA")
library(Seurat)
library(flashClust)
library(dplyr)
library(ggplot2)
library(methylumi)
library("genefilter")
setwd("/users/PAS1475/liuzeyi/guoqi/output/")
load("./conservedmarker/specific_marker_final/sample6_merge.RData")
load("/users/PAS1475/liuzeyi/guoqi/output/object_T4857.RData")
load("/users/PAS1475/liuzeyi/guoqi/output/object.2_8.RData")
load("/users/PAS1475/liuzeyi/guoqi/output/object.2_5.RData")
load("/users/PAS1475/liuzeyi/guoqi/output/object.1_1.RData")
load("/users/PAS1475/liuzeyi/guoqi/output/object.2_3.RData")
load("/users/PAS1475/liuzeyi/guoqi/output/object.18_64.RData")
options(stringsAsFactors = FALSE)
#obtain number of half of all genes in object as highly valuable gene for wgcna(not house keeping gene)
# nrow(sample_6)
# num <- nrow(sample_6) / 2
# VariableFeatures(sample_6) <- rownames(sample_6)
#findvariablefeatures:findvariablefeatures can not apply to the merged object after merge
HVG_18_64 <- FindVariableFeatures(
  object.18_64,
  selection.method = "vst",
  nfeatures = 10000,
  assay = "SCT"
)

HVG_2_3 <-
  Seurat::FindVariableFeatures(
    object.2_3,
    selection.method = "vst",
    nfeatures = 10000,
    assay = "SCT"
  )
HVG_2_5 <-
  Seurat::FindVariableFeatures(
    object.2_5,
    selection.method = "vst",
    nfeatures = 10000,
    assay = "SCT"
  )
HVG_1_1 <-
  Seurat::FindVariableFeatures(
    object.1_1,
    selection.method = "vst",
    nfeatures = 10000,
    assay = "SCT"
  )
HVG_T4857 <-
  Seurat::FindVariableFeatures(
    object_T4857,
    selection.method = "vst",
    nfeatures = 10000,
    assay = "SCT"
  )
HVG_2_8 <-
  Seurat::FindVariableFeatures(
    object.2_8,
    selection.method = "vst",
    nfeatures = 10000,
    assay = "SCT"
  )
H_2_8 <- HVG_2_8@assays$SCT@var.features
H_2_3 <- HVG_2_3@assays$SCT@var.features
H_1_1 <- HVG_1_1@assays$SCT@var.features
H_T4857 <- HVG_T4857@assays$SCT@var.features
H_18_64 <- HVG_18_64@assays$SCT@var.features
H_2_5 <- HVG_2_5@assays$SCT@var.features
All_HVG <- unique(c(H_2_5, H_18_64, H_T4857, H_1_1, H_2_3, H_2_8))
ad <- read.csv("./AD_Control_DEG/marker.AD_sct.csv")
control <- read.csv("./AD_Control_DEG/marker.control_sct.csv")
ad_gene <- ad$X
control_gene <- control$X
length(intersect(All_HVG, ad_gene))
length(intersect(All_HVG, control_gene))
#intersect with markers with no p filtering
counts<-as.data.frame(sample_6@assays$SCT@counts)
Idents(sample_6)<-sample_6$category
ad<-counts[,which(Idents(sample_6)=="AD")]
control<-counts[,which(Idents(sample_6)=="control")]

# library(pbmcapply)
# for(i in 1:nrow(ad)){
#   ad_exp_mean<-as.data.frame(pbmclapply(1:nrow(ad), fun(i){
#     mean(ad[i,])
#   }, mc.cores = detectCores()))
# }
ad_exp_mean<-apply(ad, 1, mean)
control_exp_mean<-apply(control, 1, mean)
ratio<-ad_exp_mean/control_exp_mean
up_gene<-names(ratio)[which(ratio>1)]
down_gene<-names(ratio)[which(ratio<1)]
length(intersect(All_HVG, up_gene))
length(intersect(All_HVG, down_gene))
#cannot fit to integrated sct
# DefaultAssay(object_T4857) <- "SCT"
# a <-
#   Seurat::FindVariableFeatures(sample6.combined, nfeatures = nrow(sample6.combined) /
#                                  2)
# a_hvg <- VariableFeatures(a)
# a_2_8 <- VariableFeatures(HVG_2_8)
# most variable features to use for integration-7260,not enough gene exist
# HVG <-
#   SelectIntegrationFeatures(object.list = my.all.object.list, nfeatures = 20055)
# length(intersect(HVG, ad_gene))
# length(intersect(HVG, control_gene))
#only used for microarray-MethyLumiM
#HVG_var <- varFilter(sample_6, var.cutoff = 0.5)


#obtain data matrix(sct-continuous data)
all_matrix_hvg <- sample_6@assays$SCT@data
all_matrix_hvg <-
  sample_6@assays$SCT@data[match(All_HVG, rownames(all_matrix_hvg)), ]
#row:spot,col:gene.
all_matrix_hvg <- t(as.matrix(all_matrix_hvg))
dim(all_matrix_hvg)
gene.names = colnames(all_matrix_hvg)
# caculate soft threshold-β，
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))  # in practice this should include powers up to 20.
sft = pickSoftThreshold(
  all_matrix_hvg,
  dataIsExpr = TRUE,
  powerVector = powers,
  corFnc = cor,
  corOptions = list(use = 'p'),
  networkType = "unsigned"
)
sft$powerEstimate
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 = 0.9
# SFT index as a function of different powers
plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n",
  main = paste("Scale independence")
)

text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "red"
)

# this line corresponds to using an R^2 cut-off of h
abline(h = 0.9, col = "red")
# Mean connectivity as a function of different powers
plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = paste("Mean connectivity")
)
text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1,
  col = "red"
)
softPower = 2
# #automatic module detection
# mergingThresh = 0
bwnet = blockwiseModules(
  all_matrix_hvg,
  corType = "spearman",
  maxBlockSize = 5000,
  mergeCutHeight = 0.25,
  networkType = "unsigned",
  power = softPower,
  minModuleSize = 30,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMFileBase = "maleMouseTOM-blockwise",
  verbose = 3
)
# table(bwnet$colors)
# plot(bwnet$dendrograms)
#
# MEs = bwnet$MEs
# MEs_col = MEs
# colnames(MEs_col) = paste0("ME", labels2colors(
#   as.numeric(stringr::str_replace_all(colnames(MEs),"ME",""))))
# MEs_col = orderMEs(MEs_col)
# plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
#                       marDendro = c(3,3,2,4),
#                       marHeatmap = c(3,4,2,2), plotDendrograms = T,
#                       xLabelsAngle = 90)

TOM = TOMsimilarityFromExpr(
  all_matrix_hvg,
  networkType = "unsigned",
  TOMType = "unsigned",
  power = softPower
)
colnames(TOM) = rownames(TOM) = gene.names

dissTOM = 1 - TOM
geneTree = flashClust(as.dist(dissTOM), method = "average")
plot(geneTree,
     xlab = "",
     sub = "",
     cex = 0.3)


# Set the minimum module size
minModuleSize = 20
# Module identification using dynamic tree cut

dynamicMods = cutreeDynamic(
  dendro = geneTree,
  method = "tree",
  minClusterSize = minModuleSize,
  deepSplit = 3
)
table(dynamicMods)

#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(
  geneTree,
  dynamicColors,
  "Dynamic Tree Cut",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)

#discard the unassigned genes, and focus on the rest
restGenes = (dynamicColors != "grey")
diss1 = 1 - TOMsimilarityFromExpr(all_matrix_hvg[, restGenes], power = softPower)

colnames(diss1) = rownames(diss1) = gene.names[restGenes]
hier1 = flashClust(as.dist(diss1), method = "average")
plotDendroAndColors(
  hier1,
  dynamicColors[restGenes],
  "Dynamic Tree Cut",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
#set the diagonal of the dissimilarity to NA
diag(diss1) = NA
#Visualize the Tom plot. Raise the dissimilarity matrix to the power of 4 to bring out the module structure
sizeGrWindow(7, 7)
TOMplot(diss1, hier1, as.character(dynamicColors[restGenes]))
# extract module
module_colors = setdiff(unique(dynamicColors), "grey")
for (color in module_colors) {
  module = gene.names[which(dynamicColors == color)]
  write.table(
    module,
    paste("WGCNA/", "module_", color, ".txt", sep = ""),
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}
#
data <- data.frame(1:8)
for (i in 1:length(module_colors)) {
  data[i, 1] <-
    length(gene.names[which(dynamicColors == module_colors[i])])
}
rownames(data) <- module_colors
colnames(data) <- c("gene_num")

####pathway enrichment
black <- read.table("./WGCNA/module_black.txt")
rownames(black) <- black$V1
a <- RunPathway_human(black)
Seurat.DEGs <- black
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
#Result--need to revised this part of code
setwd("//users//PAS1475//liuzeyi//guoqi//output//WGCNA")
all_file <- list.files("WGCNA/")
files<-all_file[grep(".txt",all_file)]
for (i in 1:length(module_colors)) {
  temp <- read.table(files[5])
  rownames(temp) <- temp$V1
  reactome <- RunPathway_human(temp)
  write.csv(a, paste0("./pathway_result/", strsplit(all_file[i], split =
                                                      "\\.")[[1]][1]))
}
## plot WGCNA
#read data
setwd("/users/PAS1475/liuzeyi/guoqi/output/WGCNA/")
all_file<-list.files(".")
module_colors<-strsplit(all_file[-9],split="\\.")
colors<-c()
for(i in 1:8){
  temp<-module_colors[[i]][1]
  colors<-c(colors,temp)
}
data<-list()
for(i in 1:8){
  a<-read.table(all_file[i])
  data[[i]]<-a
}
names(data)<-colors
#compute average expression in each module
setwd("/users/PAS1475/liuzeyi/guoqi/output/")
load("/users/PAS1475/liuzeyi/guoqi/output/object_T4857.RData")
load("/users/PAS1475/liuzeyi/guoqi/output/object.2_8.RData")
load("/users/PAS1475/liuzeyi/guoqi/output/object.2_5.RData")
load("/users/PAS1475/liuzeyi/guoqi/output/object.1_1.RData")
load("/users/PAS1475/liuzeyi/guoqi/output/object.2_3.RData")
load("/users/PAS1475/liuzeyi/guoqi/output/object.18_64.RData")
ad_object<-merge(object.2_3,c(object.2_8,object_T4857))
control_object<-merge(object.2_5,c(object.18_64,object.1_1))
count_ad<-ad_object@assays$SCT@data
count_control<-control_object@assays$SCT@data
ad<-c()
control<-c()
for(i in 1:8){
  temp_gene<-data[[i]]$V1
  temp_ad<-count_ad[match(temp_gene,rownames(count_ad)),]
  mean_data_ad<-apply(temp_ad,1,mean)
  mean_ad<-mean(mean_data_ad)
  ad<-c(ad,mean_ad)
  temp_control<-count_control[match(temp_gene,rownames(count_control)),]
  mean_data_control<-apply(temp_control,1,mean)
  mean_control<-mean(mean_data_control)
  control<-c(control,mean_control)
}
set.seed(123)
control_x<-sample(1:9,8,replace=F)
ad_x<-sample(11:20,8,replace=F)
control_data<-data.frame(control_x,control)
ad_data<-data.frame(ad_x,ad)
colnames(ad_data)<-c("x","y")
colnames(control_data)<-c("x","y")
all<-rbind(control_data,ad_data)
all$color<-c(rep("control",8),rep("ad",8))
all$color<-factor(all$color)
#label
temp_list<-strsplit(colors,split="_")
label<-c()
for(i in 1:8){
  temp<-temp_list[[i]][2]
  label<-c(label,temp)
}
label<-rep(label,2)
library("ggplot2")
library(ggrepel)#optimize text of label
g<-ggplot(all,aes(x=lfc,y=y,colour=color))+
  geom_vline(xintercept = 10, col = "black")+
  geom_point(size=5)+
  xlab("Control                                       AD")+
  ylab("Average of expression")+
  labs(title = "WGCNA modules") +
  theme(
    axis.text = element_text(size = 15)
  )+
  geom_text_repel(aes(x, y, label = label),size=7)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size=25,vjust = 0.5, hjust = 0.45),
        plot.title = element_text(size = 30))
g

ggsave(
  plot = g,
  filename ="avg_exp.tiff",
  device = "tiff",
  dpi = 150,
  width = 10,
  height = 10,
  units = "in"
)
ggsave(
  plot = g,
  filename ="avg_exp.png",
  device = "png",
  dpi = 150,
  width = 10,
  height = 10,
  units = "in"
)
all$x/all$y
#MA plot
##mean in 6 sample
setwd("/users/PAS1475/liuzeyi/guoqi/output/")
load("./conservedmarker/specific_marker_final/sample6_merge.RData")
count_6<-sample_6@assays$SCT@data
mean_6<-c()
for(i in 1:8){
  temp_gene<-data[[i]]$V1
  temp_control<-count_6[match(temp_gene,rownames(count_6)),]
  mean_data_control<-apply(temp_control,1,mean)
  mean_control<-mean(mean_data_control)
  mean_6<-c(mean_6,mean_control)
}
logfc<-log(all$y[9:16]/all$y[1:8])
all_re<-data.frame(logfc,mean_6)
all_re$color<-label[1:8]
p<-ggplot(all_re,aes(x=mean_6,y=logfc))+
  geom_point(size=5)+
  xlab("Average of expression")+
  ylab("Logfc")+
  labs(title = "WGCNA modules")+
  theme(
    axis.text = element_text(size = 15)
  )+
  geom_text_repel(aes(mean_6, logfc, label = color),size=7)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size=25,vjust = 0.5, hjust = 0.45),
        plot.title = element_text(size = 30))
ggsave(
  plot = p,
  filename ="WGCNA_mean.tiff",
  device = "tiff",
  dpi = 150,
  width = 10,
  height = 10,
  units = "in"
)
ggsave(
  plot = p,
  filename ="WGCNA_mean.png",
  device = "png",
  dpi = 150,
  width = 10,
  height = 10,
  units = "in"
)
ma_data<-t(data.frame(all_re$mean_6,all$y[9:16],all$y[1:8],all_re$logfc))
colnames(ma_data)<-all_re$color
rownames(ma_data)<-c("Average expression in six samples","Average expression in ad samples","Average expression in control samples","LogFC between average expression between AD v.s. control")
write.csv(ma_data,"ma_data.csv")


#-----------------------intersection of eight module with 20 gene sets 1.17
#-----load data-14
setwd("//users//PAS1475//liuzeyi//guoqi//output//WGCNA//intersection of 8 modules and 14 datasets")
double_col<-list()
for(i in 1:6){
  temp<-read_excel("14 datasets.xlsx", sheet = i)
  double_col[[(2*i-1)]]<-temp$`UP regulated`
  double_col[[2*i]]<-temp$`Down regulated`
  #extract sheet name from xlsx
  temp_name <-excel_sheets(path = "14 datasets.xlsx")[i]
  names(double_col)[[(2*i-1)]]<-paste(temp_name,"Up regulated",sep = "_")
  temp_name <-excel_sheets(path = "14 datasets.xlsx")[i]
  names(double_col)[[2*i]]<-paste(temp_name,"Down regulated",sep = "_")
}
single_col<-list()
for(i in 7:14){
  temp<-read_excel("14 datasets.xlsx", sheet = i)
  colnames(temp)<-c("a")
  single_col[[i-6]]<-temp$a
  #extract sheet name from xlsx
  names(single_col)[i-6] <-excel_sheets(path = "14 datasets.xlsx")[i]
}      
interesting_gene_sets<-append(double_col,single_col)
       
#8 modules
setwd("/users/PAS1475/liuzeyi/guoqi/output/WGCNA/")
all_file<-list.files(".")
all_file<-all_file[grep(pattern = ".txt",all_file)]
module_colors<-strsplit(all_file[-9],split="\\.")
colors<-c()
for(i in 1:8){
  temp<-module_colors[[i]][1]
  colors<-c(colors,temp)
}
eight_modules<-list()
for(i in 1:8){
  temp<-read.table(all_file[i])
  eight_modules[[i]]<-temp$V1
}
names(eight_modules)<-colors

#caculate all gene list-dataset
dataset1<-c()
for(i in 1:length(eight_modules)){
  temp<-eight_modules[[i]]
  dataset1<-union(dataset1,temp)
}
dataset2<-c()
for(i in 1:length(interesting_gene_sets)){
  temp<-interesting_gene_sets[[i]]
  dataset2<-union(dataset2,temp)
}

#Caculation
result_fisher<-data.frame()
for(i in 1:length(eight_modules)){
  for(j in 1:length(interesting_gene_sets)){
    gene1<-eight_modules[[i]]
    gene2<-interesting_gene_sets[[j]]
    result<-Fishertest(gene1,gene2,dataset1,dataset2)
    result_fisher<-rbind(result_fisher,result)
  }
}
#names
library("vctrs")
result_fisher$Module<-vec_rep_each(names(eight_modules),20)
result_fisher$Interesting_geneset<-vec_rep(names(interesting_gene_sets),8)
#adjust_p
adj_p<-c()
for(i in 1:8){
  temp<-p.adjust(result_fisher$P.value[which(result_fisher$Module==names(eight_modules)[i])],method = "fdr",20)
  adj_p<-c(adj_p,temp)
}
result_fisher$Adj_p<-adj_p
a<-result_fisher[,5:6]
result_fisher<-result_fisher[,-c(5:6)]
result_fisher4<-round(result_fisher,4)
result_fisher4<-cbind(a,result_fisher4)
write.csv(result_fisher4,"./intersection of 8 modules and 14 datasets/Result_fisher.test_4.csv")

#---function
Fishertest<-function(gene1,gene2,dataset1,dataset2){
  a<-length(intersect(gene1,gene2))
  b<-length(intersect(setdiff(dataset2,gene2),gene1))
  c<-length(intersect(setdiff(dataset1,gene1),gene2))
  d<-length(intersect(setdiff(dataset1,gene1),setdiff(dataset2,gene2)))
  data_fisher<-data.frame(gene2=c(a,c),dataset2=c(b,d))
  rownames(data_fisher)<-c("gene1","dataset1")
  re<-fisher.test(data_fisher)
  expect<-length(intersect(dataset1,dataset2))*(length(gene1)/length(dataset1))*(length(gene2)/length(dataset2))
  result<-data.frame(Intersection_num=a,Expection=expect,Odd_ratio=re[["estimate"]][["odds ratio"]],P.value=re$p.value)
  return(result)
}

#1.19 do no seperate up-regulated and down regulated genes list
#8 modules
setwd("/users/PAS1475/liuzeyi/guoqi/output/WGCNA/")
all_file<-list.files(".")
all_file<-all_file[grep(pattern = ".txt",all_file)]
module_colors<-strsplit(all_file[-9],split="\\.")
colors<-c()
for(i in 1:8){
  temp<-module_colors[[i]][1]
  colors<-c(colors,temp)
}
eight_modules<-list()
for(i in 1:8){
  temp<-read.table(all_file[i])
  eight_modules[[i]]<-temp$V1
}
names(eight_modules)<-colors
#14 gene list
#-----load data-14
setwd("//users//PAS1475//liuzeyi//guoqi//output//WGCNA//intersection of 8 modules and 14 datasets")
interesting_gene_sets<-list()
double_col<-list()
for(i in 1:6){
  temp<-read_excel("14 datasets.xlsx", sheet = i)
  temp_gene<-c(temp$`UP regulated`,temp$`Down regulated`)
  double_col[[i]]<-temp_gene
}
single_col<-list()
for(i in 7:14){
  temp<-read_excel("14 datasets.xlsx", sheet = i)
  colnames(temp)<-c("V1")
  temp_gene<-c(temp$V1)
  single_col[[(i-6)]]<-temp_gene
}
interesting_gene_sets<-append(double_col,single_col)
#extract sheet name from xlsx
names(interesting_gene_sets) <-excel_sheets(path = "14 datasets.xlsx")
#caculate all gene list-dataset
dataset1<-c()
for(i in 1:length(eight_modules)){
  temp<-eight_modules[[i]]
  dataset1<-union(dataset1,temp)
}
dataset2<-c()
for(i in 1:length(interesting_gene_sets)){
  temp<-interesting_gene_sets[[i]]
  dataset2<-union(dataset2,temp)
}

#Caculation
result_fisher<-data.frame()
for(i in 1:length(eight_modules)){
  for(j in 1:length(interesting_gene_sets)){
    gene1<-eight_modules[[i]]
    gene2<-interesting_gene_sets[[j]]
    result<-Fishertest(gene1,gene2,dataset1,dataset2)
    result_fisher<-rbind(result_fisher,result)
  }
}
#names
library("vctrs")
result_fisher$Module<-vec_rep_each(names(eight_modules),14)
result_fisher$Interesting_geneset<-vec_rep(names(interesting_gene_sets),8)
#adjust_p
adj_p<-c()
for(i in 1:8){
  temp<-p.adjust(result_fisher$P.value[which(result_fisher$Module==names(eight_modules)[i])],method = "fdr",14)
  adj_p<-c(adj_p,temp)
}
result_fisher$Adj_p<-adj_p
a<-result_fisher[,5:6]
result_fisher<-result_fisher[,-c(5:6)]
result_fisher4<-round(result_fisher,4)
result_fisher4<-cbind(a,result_fisher4)
write.csv(result_fisher4,"Result_fisher.test_4_bind.csv")


## calculate closeness centrality for genes from 8 modules.
library(igraph)
# module.cat is a two-column data frame, 
# first column is ID (color name)
# second column is Genes (genes in that module)
#8 modules
setwd("/users/PAS1475/liuzeyi/guoqi/output/WGCNA/")
all_file<-list.files(".")
all_file<-all_file[grep(pattern = ".txt",all_file)]
module_colors<-strsplit(all_file[-9],split="\\.")
colors<-c()
for(i in 1:8){
  temp<-module_colors[[i]][1]
  colors<-c(colors,temp)
}
eight_modules_genes<-c()
eight_modules_id<-c()
for(i in 1:8){
  temp<-read.table(all_file[i])
  eight_modules_genes<-c(eight_modules_genes,temp$V1)
  temp_id<-rep(colors[i],length(temp$V1))
  eight_modules_id<-c(eight_modules_id,temp_id)
}
module.cat<-data.frame(id=eight_modules_id,genes=eight_modules_genes)
#matrix
setwd("/users/PAS1475/liuzeyi/guoqi/output/")
load("./conservedmarker/specific_marker_final/sample6_merge.RData")
counts<-sample_6@assays$SCT@counts
#module category
graph_module <- module.cat
gene_from_8Module <- c()

for (i in 1:8){
genes <- graph_module$genes[graph_module$id ==colors[i]]
all_matrix_DEG <- counts[genes,]
dim(all_matrix_DEG)
all_matrix_DEG <- t(as.matrix(all_matrix_DEG))
similarity.mat <- TOMsimilarityFromExpr(all_matrix_DEG, power = 3)
graph_adj <- adjacency.fromSimilarity( similarity.mat,power = 3,type = "unsigned")
colnames(graph_adj) <- genes
rownames(graph_adj) <- genes
graph_igraph <- graph_from_adjacency_matrix(graph_adj,mode = "undirected", weighted = T)
closebess_centrality <- closeness(graph_igraph)
#data.frame
sorted_centrality <- as.data.frame(cbind(Gene = names(closebess_centrality), 
                                         cloness.centrality = closebess_centrality))
sorted_centrality <- sorted_centrality[order(as.numeric(sorted_centrality$cloness.centrality),decreasing = T),]
write.csv(sorted_centrality,  file = paste0("/users/PAS1475/liuzeyi/guoqi/output/WGCNA/closeness centrality/",colors[i],".csv"),row.names = F)
closebess_centrality_gene_names <- names(sort(closebess_centrality,decreasing = T)[1:15])
#15
gene_from_8Module <- c(gene_from_8Module,closebess_centrality_gene_names)
}


#-----------------------fishertest 1.21 mouse
library(biomaRt)
human <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
mouse <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
geneMn <- rownames(mtx)
geneHs <- getLDS(attributes = "mgi_symbol",filters="mgi_symbol",values= geneMn,mart=mouse,attributesL = "hgnc_symbol",martL=human,uniqueRows = TRUE)
#change mouse name to human name
# interesting_gene_sets_upper<-list()
# for(i in 1:14){
#   temp<-interesting_gene_sets[[i]]
#   temp_upper<-toupper(temp)
#   interesting_gene_sets_upper[[i]]<-temp_upper
# }

interesting_gene_sets_biomart<-list()
for(i in 1:14){
  temp<-interesting_gene_sets[[i]]
  temp_biomart<-getLDS(attributes = "mgi_symbol",filters="mgi_symbol",values= temp,mart=mouse,attributesL = "hgnc_symbol",martL=human,uniqueRows = TRUE)
  interesting_gene_sets_biomart[[i]]<-temp_biomart
}

#check
upper<-c()
lower<-c()
all<-c()
for(i in 1:14){
  #temp<-interesting_gene_sets_upper[[i]]
  tem<-interesting_gene_sets[[i]]
  #a<-c(length(intersect(temp,rownames(sample_6))))
  b<-c(length(intersect(tem,rownames(sample_6))))
  #c<-length(temp)
  # upper<-c(upper,a)
  # lower<-c(lower,b)
  all<-c(all,b)
}


#create new interesting sets
interesting_gene_sets[[9]]<-interesting_gene_sets_biomart[[9]]$HGNC.symbol
interesting_gene_sets[[11]]<-interesting_gene_sets_biomart[[11]]$HGNC.symbol
interesting_gene_sets[[12]]<-interesting_gene_sets_biomart[[12]]$HGNC.symbol
interesting_gene_sets[[13]]<-interesting_gene_sets_biomart[[13]]$HGNC.symbol
interesting_gene_sets[[14]]<-interesting_gene_sets_biomart[[14]]$HGNC.symbol
#dataset
setwd("//users//PAS1475//liuzeyi//guoqi//output//WGCNA//intersection of 8 modules and 14 datasets")
names(interesting_gene_sets) <-excel_sheets(path = "14 datasets.xlsx")
#caculate all gene list-dataset
#8 modules
setwd("/users/PAS1475/liuzeyi/guoqi/output/WGCNA/")
all_file<-list.files(".")
all_file<-all_file[grep(pattern = ".txt",all_file)]
module_colors<-strsplit(all_file[-9],split="\\.")
colors<-c()
for(i in 1:8){
  temp<-module_colors[[i]][1]
  colors<-c(colors,temp)
}
eight_modules<-list()
for(i in 1:8){
  temp<-read.table(all_file[i])
  eight_modules[[i]]<-temp$V1
}
names(eight_modules)<-colors
dataset1<-c()
for(i in 1:length(eight_modules)){
  temp<-eight_modules[[i]]
  dataset1<-union(dataset1,temp)
}
dataset2<-c()
for(i in 1:length(interesting_gene_sets)){
  temp<-interesting_gene_sets[[i]]
  dataset2<-union(dataset2,temp)
}
#updata function about output a dataframe with 14 gene names as colnames
#---function
Fishertest_show<-function(gene1,gene2,dataset1,dataset2){
  a<-length(intersect(gene1,gene2))
  b<-length(intersect(setdiff(dataset2,gene2),gene1))
  c<-length(intersect(setdiff(dataset1,gene1),gene2))
  d<-length(intersect(setdiff(dataset1,gene1),setdiff(dataset2,gene2)))
  data_fisher<-data.frame(gene2=c(a,c),dataset2=c(b,d))
  rownames(data_fisher)<-c("gene1","dataset1")
  re<-fisher.test(data_fisher)
  expect<-length(intersect(dataset1,dataset2))*(length(gene1)/length(dataset1))*(length(gene2)/length(dataset2))
  result<-c(round(a,1),round(re[["estimate"]][["odds ratio"]],4),re$p.value)
  return(result)
}
#Caculation_alldataframe
result_fisher<-data.frame()
for(i in 1:length(eight_modules)){
  for(j in 1:length(interesting_gene_sets)){
    gene1<-eight_modules[[i]]
    gene2<-interesting_gene_sets[[j]]
    result<-Fishertest(gene1,gene2,dataset1,dataset2)
    result_fisher<-rbind(result_fisher,result)
  }
}
#names
library("vctrs")
result_fisher$Module<-vec_rep_each(names(eight_modules),14)
result_fisher$Interesting_geneset<-vec_rep(names(interesting_gene_sets),8)
#adjust_p
adj_p<-c()
for(i in 1:8){
  temp<-p.adjust(result_fisher$P.value[which(result_fisher$Module==names(eight_modules)[i])],method = "fdr",14)
  adj_p<-c(adj_p,temp)
}
result_fisher$Adj_p<-adj_p
a<-result_fisher[,5:6]
result_fisher<-result_fisher[,-c(5:6)]
result_fisher4<-round(result_fisher,4)
result_fisher4<-cbind(a,result_fisher4)

#dataframe show
result_show<-data.frame()
temp_list<-list()
p_row<-c()
all_row<-c()
for(i in 1:8){
  for(j in 1:14){
    temp_row<-Fishertest_show(eight_modules[[i]],interesting_gene_sets[[j]],dataset1,dataset2)
    temp_list[[j]]<-temp_row
  }
  for(p in 1:14){
    a<-temp_list[[p]][3]
    p_row<-c(p_row,a)
  }
  temp_adj<-p.adjust(p_row,method = "fdr",14)
  for(k in 1:14){
    x<-paste(temp_list[[k]][1],temp_list[[k]][2],round(temp_adj[k],4),sep = "/")
    all_row<-c(all_row,x)
  }
  result_show<-rbind(result_show,all_row)
  p_row<-c()
  all_row<-c()
}
colnames(result_show)<-names(interesting_gene_sets)
rownames(result_show)<-names(eight_modules)
write.csv(result_fisher4,"./intersection of 8 modules and 14 datasets/Result_fisher.test_4_human.csv")
write.csv(result_show,"./intersection of 8 modules and 14 datasets/Result_fisher.test_4_human_show.csv")


#GO plot for four given modules
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
#load target data
library(readxl)
setwd("/users/PAS1475/liuzeyi/guoqi/output/WGCNA/4 given modules for GO plot")
module4<-list()
for(i in 1:4){
  temp<-read_excel("No.10, 4 modules GO pathway.xlsx", sheet = i)
  module4[[i]]<- temp
}
names(module4)<-excel_sheets(path = "No.10, 4 modules GO pathway.xlsx")

# library(ggplot2)
# setwd("/users/PAS1475/liuzeyi/guoqi/output/picture/enrichment_dot")
# for(i in 1:4){
#   g<-ggplot(module4[[i]], # you can replace the numbers to the row number of pathway of your interest
#             aes(x = GeneRatio_num, y = Description)) + 
#     geom_point(aes(size = Count, color = p.adjust)) +
#     theme_bw(base_size = 14) +
#     theme(axis.text = element_text(size = 14, face = "bold"),
#     )+
#     scale_colour_gradient(limits=NULL, low="red",high="blue") +
#     ylab(NULL) +
#     ggtitle("GO enrichment")+
#     guides(col = guide_legend(order =0))
#   
#   ggsave(
#     plot = g,
#     filename = paste0(names(module4)[i],"_enrichment.tiff"),
#     device = "tiff",
#     dpi = 150,
#     width = 13,
#     height = 10,
#     units = "in"
#   )
# }
# ggsave(
#   plot = g,
#   filename = paste0(names(module4)[i],"_enrichment.png"),
#   device = "png",
#   dpi = 150,
#   width = 10,
#   height = 10,
#   units = "in"
# )

#load module gene for go object
setwd("//users//PAS1475//liuzeyi//guoqi//output//WGCNA")
all_file <- list.files(pattern = "*txt")
files<-c(all_file[3],all_file[5],all_file[7],all_file[8])

#function input module_name,target_name
plot_list<-list()
for (i in 1:4) {
  temp <- read.table(files[i])
  module_gene <- temp$V1
  target_data<-module4[[i]]
  plot_list[[i]] <- dotplot_WGCNA(module_gene,target_data)
}
setwd("/users/PAS1475/liuzeyi/guoqi/output/picture/enrichment_dot")
for(i in 1:4){
  a<-plot_list[[i]]+  theme(legend.key.size = unit(1.5, 'cm'),
                            legend.title = element_text(size=18), #change legend title font size
                            legend.text = element_text(size=18))
  ggsave(plot = a,
         filename = paste0(names(module4)[i],"_enrichment_NO10.tiff"),
         device = "tiff",
         dpi = 150,
         width = 11,
         height = 12,
         units = "in")
}
#pink problem
i=2
a<-plot_list[[i]]+  theme(legend.key.size = unit(1.5, 'cm'),
                          legend.title = element_text(size=18), #change legend title font size
                          legend.text = element_text(size=18))
ggsave(plot = a,
       filename = paste0(names(module4)[i],"_enrichment_NO10.tiff"),
       device = "tiff",
       dpi = 150,
       width = 17,
       height = 12,
       units = "in")
