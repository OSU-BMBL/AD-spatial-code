setwd('/fs/ess/PCON0022/guoqi/NC/output')
load("object.2_3.RData")
load("object.2_8.RData")
load("object_T4857.RData")

setwd("/fs/ess/PCON0022/guoqi/NC-snrna/Spatial/pathlogical annotation")
path_all<-list.files(".")
##Aβ
#--------2-3
remove_noise<-function(object,ann_name){
  ann<-read.csv(ann_name)
  print(identical(colnames(object), ann[,1]))
  ann_f<-ann[match(colnames(object), ann[,1]), ]
  print(identical(ann_f[,1],colnames(object)))
  return(ann_f)
}

add_patho<-function(object,id){
  ann_ab<-remove_noise(object,paste0(id," Aβ specific.csv"))
  ann_abp<-ann_ab[which(ann_ab$Aβ=="Aβ"),]
  colnames(ann_abp)<-c("cell","pathology_p")
  ann_at<-remove_noise(object,paste0(id," AT8 specific.csv"))
  ann_atp<-ann_at[which(ann_at[,2]=="AT8"),]
  colnames(ann_atp)<-c("cell","pathology_p")
  ann_abat<-remove_noise(object,paste0(id," AT8,Aβ double positive specific.csv"))
  ann_abatp<-ann_abat[which(ann_abat[,2]=="Aβ/AT8"),]
  colnames(ann_abatp)<-c("cell","pathology_p")
  #rbind
  ann<-rbind(ann_abp,ann_atp)
  ann<-rbind(ann,ann_abatp)
  othercell<-colnames(object)[-match(ann$cell,colnames(object))]
  ann_other<-data.frame(cell=othercell,pathology_p="non")
  ann<-rbind(ann,ann_other)
  ann<-ann[match(colnames(object),ann$cell),]
  print(identical(ann$cell,colnames(object)))
  object$pathologyp<-ann$pathology_p
  return(object)
}
object.2_3<-add_patho(object.2_3,"2-3")

# ann_2_3ab<-remove_noise(object.2_3,"2-3 Aβ specific.csv")
# ann_2_3abp<-ann_2_3ab[which(ann_2_3ab$Aβ=="Aβ"),]
# ann_2_3at<-remove_noise(object.2_3,"2-3 AT8 specific.csv")
# ann_2_3atp<-ann_2_3at[which(ann_2_3at$AT8=="AT8"),]
# ann_2_3atab<-remove_noise(object.2_3,"2-3 AT8,Aβ double positive specific.csv")
# ann_2_3atabp<-ann_2_3atab[which(ann_2_3atab$Aβ.AT8.double.positive=="Aβ/AT8"),]
# colnames(ann_2_3)<-c("cell","Aβ","AT8","AT8Aβ")
# ann_2_3<-rbind(ann_2_3abp,ann_2_3atp)
# ann_2_3<-cbind(ann_2_3,ann_2_3atab[,2])
# colnames(ann_2_3)<-c("cell","Aβ","AT8","AT8Aβ")

#-------2-8
object.2_8<-add_patho(object.2_8,"2-8")

# ann_2_8ab<-remove_noise(object.2_8,"2-8 Aβ specific.csv")
# ann_2_8at<-remove_noise(object.2_8,"2-8 AT8 specific.csv")
# ann_2_8atab<-remove_noise(object.2_8,"2-8 AT8,Aβ double positive specific.csv")
# ann_2_8<-cbind(ann_2_8ab,ann_2_8at[,2])
# ann_2_8<-cbind(ann_2_8,ann_2_8atab[,2])
# colnames(ann_2_8)<-c("cell","Aβ","AT8","AT8Aβ")
# 
# identical(colnames(object.2_8),ann_2_8$cell)
# object.2_8<-AddMetaData(object.2_8,ann_2_8)
#-------T4857
object_T4857<-add_patho(object_T4857,"T4957")

# ann_Tab<-remove_noise(object_T4857,"T4957 Aβ specific.csv")
# ann_Tat<-remove_noise(object_T4857,"T4957 AT8 specific.csv")
# ann_Tatab<-remove_noise(object_T4857,"T4957 AT8,Aβ double positive specific.csv")
# ann_T<-cbind(ann_Tab,ann_Tat[,2])
# ann_T<-cbind(ann_T,ann_Tatab[,2])
# colnames(ann_T)<-c("cell","Aβ","AT8","AT8Aβ")
# 
# identical(colnames(object_T4857),ann_T$cell)
# object_T4857<-AddMetaData(object_T4857,ann_T)

#merge
object<-merge(object.2_3,c(object_T4857, object.2_8))
#findmarker
object<-PrepSCTFindMarkers(object)
pathp_findmarker<-function(object,group1,group2,outputname){
  DefaultAssay(object)<-"SCT"
  markers <-
    FindMarkers(
      object,
      ident.1 = group1,
      ident.2 = group2,
      group.by = "pathologyp",
      verbose = TRUE,
      assay = "SCT",
      only.pos = F
    )
  markers <- cbind(rownames(markers), markers)
  markers <-
    markers[markers$p_val_adj < 0.05, ]
  rownames(markers) <- NULL
  colnames(markers)[1] <- c("Marker")
  marker_up <-
    markers[markers$avg_log2FC > 0, ]
  marker_down <- markers[markers$avg_log2FC < 0, ]
  write.csv(marker_up,paste0("/fs/ess/PCON0022/guoqi/NC/Revision/Output/pathology_deg_1207/",
                             outputname,"_up.csv"),row.names = F)
  write.csv(marker_down,paste0("/fs/ess/PCON0022/guoqi/NC/Revision/Output/pathology_deg_1207/",
                               outputname,"_down.csv"),row.names = F)
}
pathp_findmarker(object = object,"Aβ","AT8")
pathp_findmarker(object = object,"Aβ/AT8","AT8","AβAT8_AT8")
pathp_findmarker(object = object,"Aβ/AT8","Aβ","AβAT8_Aβ")
