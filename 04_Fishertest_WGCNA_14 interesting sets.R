library("WGCNA")
library(Seurat)
library(flashClust)
library(dplyr)
library(ggplot2)
library(methylumi)
library("genefilter")
library(biomaRt)
#fisher test
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









#----- load data
#eight modules
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
#load data-14
setwd("//users//PAS1475//liuzeyi//guoqi//output//WGCNA//intersection of 8 modules and 14 datasets")
interesting_gene_sets<-list()
double_col<-list()
for(i in 1:6){
  temp<-read_excel("No.3, 14 datasets update.xlsx", sheet = i)
  temp_gene<-c(temp$`UP regulated`,temp$`Down regulated`)
  double_col[[i]]<-temp_gene
}
single_col<-list()
for(i in 7:14){
  temp<-read_excel("No.3, 14 datasets update.xlsx", sheet = i)
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




#-------Caculation_alldataframe
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
write.csv(result_fisher4,"./intersection of 8 modules and 14 datasets/Result_fisher.test_4_human.csv")



#--------dataframe show
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
write.csv(result_show,"./intersection of 8 modules and 14 datasets/Result_fisher.test_4_human_show.csv")
