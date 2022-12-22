#check data
setwd("/fs/ess/PCON0022/guoqi/NC-snrna/Spatial/Deconvolution_6scRNA/jixin_f_results")
temp1<-read.csv("1-1_Cell2location_results.csv",header=F)
temp2<-read.csv("../Deconvolution_results/1-1_Cell2location_results.csv",header = F)
identical(temp1,temp2)

setwd("/fs/ess/PCON0022/guoqi/NC-snrna/Spatial/Deconvolution_6scRNA/Deconvolution_results_f_jixin/")
samp_file<-list.files(".",pattern = ".csv")
samp_file<-samp_file[-c(2,4)]
ad_file<-samp_file[-c(1,2,4)]
#load ad data (deconvolution)
decon_list<-list()
for(i in 1:3){
  decon_list[[i]]<-read.csv(ad_file[i],header = T)
  rownames(decon_list[[i]])<-decon_list[[i]]$X
  decon_list[[i]]<-decon_list[[i]][,c(2,6)]
}
names(decon_list)<-c("2-3","2-8","T4857")
#load pathlogy data
setwd("/fs/ess/PCON0022/guoqi/NC-snrna/Spatial/pathlogical annotation")
path_all<-list.files(".")
#check if names are equal- equal
#at8
at8_df<-read.csv("2-3 AT8 specific.csv",header = T)
identical(at8_df$Barcode,rownames(decon_list[[1]]))
decon_list[[1]]$AT8<-at8_df$AT8
at8_df_28<-read.csv("2-8 AT8 specific.csv",header = T)
identical(at8_df_28$Barcode,rownames(decon_list[[2]]))
decon_list[[2]]$AT8<-at8_df_28$AT8
at8_df_T4857<-read.csv("./T4957 AT8 specific.csv")
identical(at8_df_T4857$Barcode,rownames(decon_list[[3]]))
decon_list[[3]]$AT8<-at8_df_T4857$AT8
#ab
ab_df<-read.csv("2-3 Aβ specific.csv",header = T)
identical(ab_df$Barcode,rownames(decon_list[[1]]))
decon_list[[1]]$AB<-ab_df$Aβ
ab_df_28<-read.csv("2-8 Aβ specific.csv",header = T)
identical(ab_df_28$Barcode,rownames(decon_list[[2]]))
decon_list[[2]]$AB<-ab_df_28$Aβ
ab_df_T4857<-read.csv("./T4957 Aβ specific.csv")
identical(ab_df_T4857$Barcode,rownames(decon_list[[3]]))
decon_list[[3]]$AB<-ab_df_T4857$Aβ
#abat
abat8_df<-read.csv("2-3 AT8,Aβ double positive specific.csv",header = T)
identical(abat8_df$Barcode,rownames(decon_list[[1]]))
decon_list[[1]]$ABAT<-abat8_df$Aβ.AT8.double.positive
abat8_df_28<-read.csv("2-8 AT8,Aβ double positive specific.csv",header = T)
identical(abat8_df_28$Barcode,rownames(decon_list[[2]]))
decon_list[[2]]$ABAT<-abat8_df_28$Aβ.AT8.double.positive
abat8_df_T4857<-read.csv("./T4957 AT8,Aβ double positive specific.csv")
identical(abat8_df_T4857$Barcode,rownames(decon_list[[3]]))
decon_list[[3]]$ABAT<-abat8_df_T4857$Aβ.AT8.double.positive

#violin
library(ggpubr)
library(ggplot2)
all_dec<-rbind(decon_list[[1]],decon_list[[2]])
all_dec<-rbind(all_dec,decon_list[[3]])
write.csv(all_dec,"all_decon.csv",row.names = T)

#----------read data
all_dec<-read.csv("all_decon.csv")
colnames(all_dec)[1:2]<-c("cell","proportion")

micro_dec<-all_dec[,c(1,3:6)]
micro_at8<-micro_dec[,c(1:3)]
micro_ab<-micro_dec[,c(1:2,4)]
micro_abat<-micro_dec[,c(1:2,5)]

ast_at8<-all_dec[,c(1:2,4)]
ast_ab<-all_dec[,c(1:2,5)]
ast_abat8<-all_dec[,c(1:2,6)]

volin<-function(metadata,outputname){
  pathology<-sort(unique(metadata[,3]))[2]
  level3<-sort(unique(metadata[,3]))[5]
  indexa<-which(metadata[,3]==pathology)
  indexb<-which(metadata[,3]==level3)
  metadata<-metadata[c(indexa,indexb),]
  metadata[,3]<-factor(metadata[,3],levels=c(pathology,level3))
  colnames(metadata)[2]<-"proportion"
  colnames(metadata)[3]<-"category"
  p<-ggplot(metadata, aes(x=category, y=proportion,fill=category)) + 
    geom_violin(trim=FALSE)+
    stat_summary(fun.data=data_summary)+
    theme_classic()+
    #scale_fill_manual(values=c("#56B4E9","#97f7a7"))+
    stat_compare_means( label = "p.signif",method = "wilcox.test")+
    theme(text = element_text(size=25))
  ggsave(
    plot = p,
    filename = paste0("/fs/ess/PCON0022/guoqi/NC-snrna/Spatial/sig_test(astr_mic)/",outputname,"_violin.tiff"),
    device = "tiff",
    dpi = 150,
    width = 12,
    height = 10,
    units = "in"
  )
}
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
#violin plot
volin(micro_at8,"microglia_AT8")
volin(micro_ab,"microglia_Aβ")
volin(micro_abat,"microglia_AβAT8")
volin(ast_ab,"astrocytes_Aβ")
volin(ast_at8,"astrocytes_AT8")
volin(ast_abat8,"astrocytes_AβAT8")

#statistical summary
statistical_sum<-function(metadata,outputname){
  pathology<-sort(unique(metadata[,3]))[2]
  level3<-sort(unique(metadata[,3]))[5]
  indexa<-which(metadata[,3]==pathology)
  indexb<-which(metadata[,3]==level3)
  metadata<-metadata[c(indexa,indexb),]
  metadata[,3]<-factor(metadata[,3],levels=c(pathology,level3))
  colnames(metadata)[2]<-"proportion"
  colnames(metadata)[3]<-"category"
  mean_py<-mean(metadata$proportion[which(metadata$category==pathology)])
  mean_non<-mean(metadata$proportion[which(metadata$category==level3)])
  sd_py<-sd(metadata$proportion[which(metadata$category==pathology)])
  sd_non<-sd(metadata$proportion[which(metadata$category==level3)])
  p<-as.data.frame(compare_means(proportion ~ category,  data = metadata, method = "wilcox.test"))
  py<-c(mean_py,sd_py,as.character(p[1,c(2:5,7:8)]))
  non<-c(mean_non,sd_non,as.character(p[1,c(2:5,7:8)]))
  df<-as.data.frame(t(data.frame(py,non)))
  colnames(df)<-c("mean","sd","group1","group2","p","p.adj","p.signif","method")
  rownames(df)<-c(pathology,level3)
  df$condition=outputname
  write.csv(df,paste0(outputname,"_statistics.csv"))
}
setwd("/fs/ess/PCON0022/guoqi/NC-snrna/Spatial/sig_test(astr_mic)")
statistical_sum(micro_at8,"microglia_AT8")
statistical_sum(micro_ab,"microglia_Aβ")
statistical_sum(micro_abat,"microglia_AβAT8")
statistical_sum(ast_ab,"astrocytes_Aβ")
statistical_sum(ast_at8,"astrocytes_AT8")
statistical_sum(ast_abat8,"astrocytes_AβAT8")

#three
library("Seurat")
load("/fs/ess/PCON0022/guoqi/NC/output/sample6_clusters_n.RData")
setwd("/fs/ess/PCON0022/guoqi/NC-snrna/Spatial/Deconvolution_6scRNA/Deconvolution_results_f_jixin/")
table(sample6.combined$Layer)
#add control's proportion
control_file<-samp_file[c(1,2,4)]
decon_control_list<-list()
pure<-sort(unique(sample6.combined$Layer))[1:6]
Idents(sample6.combined)<-sample6.combined$Layer
pure_combined<-subset(sample6.combined,idents=pure)
sampleid<-c("1_1","18_64","2_5")
control_df<-data.frame()
for(i in 1:3){
  temp<-read.csv(control_file[i],header = T)
  temp$X<-gsub("1",sampleid[i],temp$X)
  b<-temp[match(colnames(pure_combined),temp$X),]
  b<-na.omit(b)
  b<-b[,(c(1,2,6))]
  b$category<-"Control"
  control_df<-rbind(control_df,b)
}

control_micro<-control_df[,c(1,3,4)]
micro_dec<-all_dec[,c(1,3:6)]
micro_at8<-micro_dec[,c(1:3)]
micro_ab<-micro_dec[,c(1:2,4)]
micro_abat<-micro_dec[,c(1:2,5)]

control_ast<-control_df[,c(1,2,4)]
ast_at8<-all_dec[,c(1:2,4)]
ast_ab<-all_dec[,c(1:2,5)]
ast_abat8<-all_dec[,c(1:2,6)]

volin3<-function(metadata,control_metadata,outputname){
  colnames(metadata)[2]<-"proportion"
  colnames(metadata)[3]<-"category"
  colnames(control_metadata)<-c("cell","proportion","category")
  metadata<-rbind(metadata,control_metadata)
  pathology<-sort(unique(metadata[,3]))[2]
  level3<-sort(unique(metadata[,3]))[6]
  control<-sort(unique(metadata[,3]))[3]
  indexa<-which(metadata[,3]==pathology)
  indexb<-which(metadata[,3]==level3)
  indexc<-which(metadata[,3]==control)
  metadata<-metadata[c(indexa,indexb,indexc),]
  metadata[,3]<-factor(metadata[,3],levels=c(pathology,level3,control))
  p<-ggplot(metadata, aes(x=category, y=proportion,fill=category)) + 
    geom_violin(trim=FALSE)+
    stat_summary(fun.data=data_summary)+
    theme_classic()+
    #scale_fill_manual(values=c("#56B4E9","#97f7a7"))+
    stat_compare_means( label = "p.signif",method = "anova")+
    theme(text = element_text(size=25))
  # ggsave(
  #   plot = p,
  #   filename = paste0("/fs/ess/PCON0022/guoqi/NC-snrna/Spatial/sig_test(astr_mic)/three_control/",outputname,"_violin.tiff"),
  #   device = "tiff",
  #   dpi = 150,
  #   width = 12,
  #   height = 10,
  #   units = "in"
  # )
  return(p)
}

volin3(micro_at8,control_micro,"microglia_AT8")
volin3(micro_ab,control_micro,"microglia_Aβ")
volin3(micro_abat,control_micro,"microglia_AβAT8")
volin3(ast_ab,control_ast,"astrocytes_Aβ")
volin3(ast_at8,control_ast,"astrocytes_AT8")
volin3(ast_abat8,control_ast,"astrocytes_AβAT8")

statistical_sum3<-function(metadata,control_metadata,outputname){
  colnames(metadata)<-c("cell","proportion","category")
  colnames(control_metadata)<-c("cell","proportion","category")
  metadata<-rbind(metadata,control_metadata)
  pathology<-sort(unique(metadata[,3]))[2]
  level3<-sort(unique(metadata[,3]))[6]
  control<-sort(unique(metadata[,3]))[3]
  indexa<-which(metadata[,3]==pathology)
  indexb<-which(metadata[,3]==level3)
  indexc<-which(metadata[,3]==control)
  metadata<-metadata[c(indexa,indexb,indexc),]
  mean_py<-mean(metadata$proportion[which(metadata$category==pathology)])
  mean_non<-mean(metadata$proportion[which(metadata$category==level3)])
  mean_control<-mean(metadata$proportion[which(metadata$category==control)])
  sd_py<-sd(metadata$proportion[which(metadata$category==pathology)])
  sd_non<-sd(metadata$proportion[which(metadata$category==level3)])
  sd_control<-sd(metadata$proportion[which(metadata$category==control)])
  p<-as.data.frame(compare_means(proportion ~ category,  data = metadata, method = "anova"))
  py<-c(mean_py,sd_py,as.character(p[c(2:6)]))
  non<-c(mean_non,sd_non,as.character(p[c(2:6)]))
  con<-c(mean_control,sd_control,as.character(p[c(2:6)]))
  df<-as.data.frame(t(data.frame(py,non,con)))
  colnames(df)<-c("mean","sd","p","p.adj","p.format","p.signif","method")
  rownames(df)<-c(pathology,level3,control)
  df$condition=outputname
  df$group<-c(pathology,level3,control)
  write.csv(df,paste0("/fs/ess/PCON0022/guoqi/NC-snrna/Spatial/sig_test(astr_mic)/three_control/",outputname,"_statistics.csv"))
}
statistical_sum3(micro_at8,control_micro,"microglia_AT8")
statistical_sum3(micro_ab,control_micro,"microglia_Aβ")
statistical_sum3(micro_abat,control_micro,"microglia_AβAT8")
statistical_sum3(ast_ab,control_ast,"astrocytes_Aβ")
statistical_sum3(ast_at8,control_ast,"astrocytes_AT8")
statistical_sum3(ast_abat8,control_ast,"astrocytes_AβAT8")

#calculate avg_log2FC
statistical_sum_fc<-function(metadata,outputname,group1,group2){
  indexa<-which(metadata[,3]==group1)
  indexb<-which(metadata[,3]==group2)
  metadata<-metadata[c(indexa,indexb),]
  metadata[,3]<-factor(metadata[,3],levels=c(group1,group2))
  colnames(metadata)[2]<-"proportion"
  colnames(metadata)[3]<-"category"
  mean_py<-mean(metadata$proportion[which(metadata$category==group1)])
  mean_non<-mean(metadata$proportion[which(metadata$category==group2)])
  sd_py<-sd(metadata$proportion[which(metadata$category==group1)])
  sd_non<-sd(metadata$proportion[which(metadata$category==group2)])
  avg_logfc_pc<-log(mean_py/mean_non)
  avg_logfc_cp<-log(mean_non/mean_py)
  p<-as.data.frame(compare_means(proportion ~ category,  data = metadata, method = "wilcox.test"))
  py<-c(group1,group2,avg_logfc_pc,as.character(p[c(4:5,7:8)]))
  non<-c(group2,group1,avg_logfc_cp,as.character(p[c(4:5,7:8)]))
  df<-as.data.frame(t(data.frame(py,non)))
  colnames(df)<-c("group1","group2","avg_logFC","p","p.adj","p.signif","method")
  rownames(df)<-NULL
  write.csv(df,paste0("/fs/ess/PCON0022/guoqi/NC-snrna/Spatial/sig_test(astr_mic)/avg_logfc/",outputname,"_statistics.csv"))
}
#microglia
statistical_sum_fc(micro_ab,"Micro_Aβ_nonlevel3","Aβ","non_Aβ_level3")
colnames(micro_ab)<-c("cell","proportion","category")
colnames(control_micro)<-c("cell","proportion","category")
micro_ab_control<-rbind(micro_ab,control_micro)
statistical_sum_fc(micro_ab_control,"Micro_Aβ_control","Aβ","Control")
statistical_sum_fc(micro_ab_control,"Micro_nonAβlevel3_control2","non_Aβ_level3","Control")

statistical_sum_fc(micro_at8,"Micro_AT8_nonlevel3","AT8","non_AT8_level3")
colnames(micro_at8)<-c("cell","proportion","category")
micro_at8_control<-rbind(micro_at8,control_micro)
statistical_sum_fc(micro_at8_control,"Micro_AT8_control","AT8","Control")
statistical_sum_fc(micro_at8_control,"Micro_nonAT8level3_control2","non_AT8_level3","Control")

statistical_sum_fc(micro_abat,"Micro_AβAT8_nonlevel3","Aβ/AT8","non_Aβ/AT8 level3")
colnames(micro_abat)<-c("cell","proportion","category")
micro_abat_control<-rbind(micro_abat,control_micro)
statistical_sum_fc(micro_abat_control,"Micro_AβAT8_control","Aβ/AT8","Control")
statistical_sum_fc(micro_abat_control,"Micro_nonAβAT8level3_control2","non_Aβ/AT8 level3","Control")

#astrocyte
statistical_sum_fc(ast_ab,"Astrocyte_Aβ_nonlevel3","Aβ","non_Aβ_level3")
colnames(ast_ab)<-c("cell","proportion","category")
colnames(control_ast)<-c("cell","proportion","category")
ast_ab_control<-rbind(ast_ab,control_ast)
statistical_sum_fc(ast_ab_control,"Astrocyte_Aβ_control","Aβ","Control")
statistical_sum_fc(ast_ab_control,"Astrocyte_nonAβlevel3_control2","non_Aβ_level3","Control")

statistical_sum_fc(ast_at8,"Astrocyte_AT8_nonlevel3","AT8","non_AT8_level3")
colnames(ast_at8)<-c("cell","proportion","category")
ast_at8_control<-rbind(ast_at8,control_ast)
statistical_sum_fc(ast_at8_control,"Astrocyte_AT8_control","AT8","Control")
statistical_sum_fc(ast_at8_control,"Astrocyte_nonAT8level3_control2","non_AT8_level3","Control")

statistical_sum_fc(ast_abat8,"Astrocyte_AβAT8_nonlevel3","Aβ/AT8","non_Aβ/AT8 level3")
colnames(ast_abat8)<-c("cell","proportion","category")
ast_abat_control<-rbind(ast_abat8,control_ast)
statistical_sum_fc(ast_abat_control,"Astrocyte_AβAT8_control","Aβ/AT8","Control")
statistical_sum_fc(ast_abat_control,"Astrocyte_nonAβAT8level3_control2","non_Aβ/AT8 level3","Control")
