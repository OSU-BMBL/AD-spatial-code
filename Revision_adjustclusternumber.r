#load packages
library("Seurat")
library("ggplot2")

#load data
load("/fs/ess/PCON0022/guoqi/NC/output/sample6_clusters_n.RData")
#ARI
library("mclust")
ari<-function(resolution,object){
  object<-FindClusters(object,resolution = resolution)
  ari<-adjustedRandIndex(object@meta.data[,12], object$Layer)
  return(ari)
}
adjustedRandIndex(sample6.combined$integrated_snn_res.0.17, sample6_adjustclusters$Layer)
adjustedRandIndex(sample6_adjustclusters$integrated_snn_res.0.11, sample6_adjustclusters$Layer)
adjustedRandIndex(sample6_adjustclusters$integrated_snn_res.0.35, sample6_adjustclusters$Layer)
sample6_adjustclusters<-FindClusters(sample6.combined,resolution = 0.5)
adjustedRandIndex(sample6_adjustclusters$integrated_snn_res.0.5, sample6_adjustclusters$Layer)
sample6_adjustclusters<-FindClusters(sample6.combined,resolution = 0.8)
adjustedRandIndex(sample6_adjustclusters$integrated_snn_res.0.8, sample6_adjustclusters$Layer)

#observe distribution
sample_list<-SplitObject(sample6.combined,"patientID")
temp<-sample_list$`2_5`$nCount_Spatial
hist(temp,main=paste0("Reads distribution of",sampleid),
     xlab="nCount")
hist_id<-function(sample,id){
  temp<-sample$nCount_Spatial
  tiff(paste0(id,"_distributiion.tiff"))
  hist(temp,main=paste0("Reads distribution of ",id),
       xlab="nCount")
  dev.off()
}
lapply(sample_list,hist_id )
a<-hist_id(sample_list[[1]],"2_5")
for(i in 1:6){
  hist_id(sample_list[[i]],names(sample_list)[i])
}
setwd("/fs/ess/PCON0022/guoqi/NC/Revision/Output/integrated3000_gene/")

hist(temp,main=paste0("Reads distribution of ",id),
     xlab="nCount")

temp<-sample_list[[1]]$nCount_SCT
hist(temp,main=paste0("Reads distribution of ",id),
     xlab="nCount")