# compute logfc for specific gene lists
setwd("/fs/ess/PCON0022/guoqi/NC/output")
load("object.18_64.RData")
load("object.2_5.RData")
load("object.2_8.RData")
load("object_T4857.RData")
load("object.1_1.RData")
load("object.2_3.RData")
#load gene list
setwd("/fs/ess/PCON0022/guoqi/NC/Revision")
deg <- read_excel("Genelist.xlsx", sheet = 1)
update_deg<-data.frame(gene=unique(deg$`GENE list`))
#merge
object_merge <-
  merge(object.18_64,
        y = list(object.2_5, object.2_8, object_T4857,
                 object.1_1, object.2_3))
qsave(object_merge,"/fs/ess/PCON0022/guoqi/NC/Revision/Output/object_merge.qs")

#DEG
object_merge_pre<-qread("/fs/ess/PCON0022/guoqi/NC/Revision/Output/object_merge_pre.qs")
DefaultAssay(object_merge) = "SCT"
object_merge<-PrepSCTFindMarkers(object_merge)
Idents(object_merge)<-object_merge$Layer
all_deg<-data.frame()
layer<-sort(unique(object_merge_pre$Layer))
for(i in 1:length(layer)){
  object_layer<-subset(object_merge_pre,idents=layer[i])
  deg_df <-
    FindMarkers(
      object_layer,
      features = update_deg$gene,
      ident.1 = "AD",
      ident.2 = "control",
      group.by = "category",
      only.pos = F,
      recorrect_umi = FALSE,
      logfc.threshold = 0,
      min.pct = 0,
      min.cells.feature = 0,
      min.cells.group = 0
    )
  deg_df$layer<-layer[i]
  deg_df$gene<-rownames(deg_df)
  all_deg<-rbind(all_deg,deg_df)
}

write.csv(all_deg,"/fs/ess/PCON0022/guoqi/NC/Revision/Output/deg_logfc_layer.csv")