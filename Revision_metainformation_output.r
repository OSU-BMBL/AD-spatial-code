#load packages
library("Seurat")
library("ggplot2")

#load data
load("/fs/ess/PCON0022/guoqi/NC/output/sample6_clusters_n.RData")

#output meta and expression for uploading GEO
meta<-sample6.combined@meta.data
write.csv(meta,"/fs/ess/PCON0022/guoqi/NC/Revision/Output/GEO_cankun/meta.csv",row.names = T,quote = F)
raw_expression<-sample6.combined@assays$Spatial@counts
write.csv(raw_expression,"/fs/ess/PCON0022/guoqi/NC/Revision/Output/GEO_cankun/rawexpression.csv",row.names = T,quote = F)

#output 3000 variable genes
DefaultAssay(sample6.combined)<-"integrated"
length(VariableFeatures(sample6.combined))
fvg<-SpatiallyVariableFeatures(sample6.combined,nfeatures = 3000)
a<-HVFInfo(
  sample6.combined,
  selection.method = c("sct"),
  assay = "integrated"
)
length(intersect(rownames(a),var))
hvg<-a[order(a$variance,decreasing=T)[1:3000],]
identical(VariableFeatures(sample6.combined),rownames(sample6.combined@assays$integrated))
length(intersect(rownames(var_inf),rownames(sample6.combined@assays$integrated)))
write.csv(var_inf,"/fs/ess/PCON0022/guoqi/NC/Revision/Output/integrated3000_gene/3000_integratedgene_scores.csv")
write.table(row.names(var_inf),"/fs/ess/PCON0022/guoqi/NC/Revision/Output/integrated3000_gene/pure_3000selecteintegratefeatures.txt",quote = F,row.names = F)
reads_perspot<-data.frame(spots=colnames(sample6.combined),
                                         reads_perspot=sample6.combined$nCount_Spatial)
write.csv(reads_perspot,"/fs/ess/PCON0022/guoqi/NC/Revision/Output/integrated3000_gene/reads_perspot.csv",row.names = F,quote = F)
