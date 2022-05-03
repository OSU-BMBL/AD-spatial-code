setwd("/users/PAS1475/liuzeyi/guoqi/output/pathology")
#Ab
Ab2_3 <- read.csv("2-3 Aβ annotation.csv",
                  header = T,
                  row.names = 1)
T4857 <- read.csv("T4857 Aβ annotation.csv",
                  header = T,
                  row.names = 1)
Ab2_8 <- read.csv("2-8 Aβ annotation.csv",
                  header = T,
                  row.names = 1)
colnames(Ab2_8) <- colnames(T4857)
setwd("/users/PAS1475/liuzeyi/guoqi/output")
load("object.2_3.RData")
load("object.2_8.RData")
load("object_T4857.RData")
load("object.1_1.RData")
load("object.2_5.RData")
load("object.18_64.RData")
object.2_3 <- AddMetaData(object.2_3, metadata = Ab2_3)
object_T4857 <- AddMetaData(object_T4857, metadata = T4857)
object.2_8 <- AddMetaData(object.2_8, metadata = Ab2_8)
object_ab <- merge(object.2_3, c(object_T4857, object.2_8))
DefaultAssay(object_ab) <- "SCT"
table(object_ab$condition)
spot_Ab<-which(object_ab$condition=="Aβ")
spot_Ab_n1<-which(object_ab$condition=="non_Aβ_level1")
spot_Ab_n2<-which(object_ab$condition=="non_Aβ_level2")
spot_Ab_n3<-which(object_ab$condition=="non_Aβ_level3")
counts<-object_ab@assays$SCT@counts
spot_ab_sp<-list(spot_Ab,spot_Ab_n1,spot_Ab_n2,spot_Ab_n3)
#percentage
for(i in 1:4){
  per<-counts[which(rownames(counts)=="SERPINA3"),spot[[i]]]
  print(length(which(per!=0))/length(per))
}
c<-c()
for(i in 1:4){
  temp<-mean(counts[which(rownames(counts)=="C1QB"),spot_ab_sp[[i]]])
  c<-c(c,temp)
}
d<-c()
for(i in 1:4){
  per<-counts[which(rownames(counts)=="C1QB"),spot_ab_sp[[i]]]
  temp<-length(which(per!=0))/length(per)
  d<-c(d,temp)
}
e<-c()
for(i in 1:4){
  temp<-mean(counts[which(rownames(counts)=="GFAP"),spot_ab_sp[[i]]])
  e<-c(e,temp)
}
f<-c()
for(i in 1:4){
  per<-counts[which(rownames(counts)=="GFAP"),spot_ab_sp[[i]]]
  temp<-length(which(per!=0))/length(per)
  f<-c(f,temp)
}
re<-t(data.frame(e,f))
write.csv(re,"Ab.csv")
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
object.2_81 <- AddMetaData(object.2_8, metadata = Ab2_8_sp)
object_ab <- merge(object.2_3, c(object_T4857, object.2_8))
DefaultAssay(object_ab) <- "SCT"
spot_Ab_sp<-which(object_ab$condition=="Aβ")
spot_Ab_sp_n1<-which(object_ab$condition=="non_Aβ_level1")
spot_Ab_sp_n2<-which(object_ab$condition=="non_Aβ_level2")
spot_Ab_sp_n3<-which(object_ab$condition=="non_Aβ_level3")
counts<-object_ab@assays$SCT@counts
spot_ab_sp<-list(spot_Ab_sp,spot_Ab_sp_n1,spot_Ab_sp_n2,spot_Ab_sp_n3)
for(i in 1:4){
  print(exp_Ab_sp<-mean(counts[which(rownames(counts)=="SERPINA3"),spot_ab_sp[[i]]]))
}
for(i in 1:4){
  per<-counts[which(rownames(counts)=="SERPINA3"),spot_ab_sp[[i]]]
  print(length(which(per!=0))/length(per))
}
c<-c()
for(i in 1:4){
  temp<-mean(counts[which(rownames(counts)=="C1QB"),spot_ab_sp[[i]]])
  c<-c(c,temp)
}
d<-c()
for(i in 1:4){
  per<-counts[which(rownames(counts)=="C1QB"),spot_ab_sp[[i]]]
  temp<-length(which(per!=0))/length(per)
  d<-c(d,temp)
}
e<-c()
for(i in 1:4){
  temp<-mean(counts[which(rownames(counts)=="GFAP"),spot_ab_sp[[i]]])
  e<-c(e,temp)
}
f<-c()
for(i in 1:4){
  per<-counts[which(rownames(counts)=="GFAP"),spot_ab_sp[[i]]]
  temp<-length(which(per!=0))/length(per)
  f<-c(f,temp)
}
re<-t(data.frame(e,f))
re<-t(data.frame(c,d))
write.csv(re,"Ab_sp.csv")
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
ab_cat<-sort(unique(object_ab$condition))
spot_Ab_sp<-which(object_ab$condition==ab_cat[2])
spot_Ab_sp_n1<-which(object_ab$condition==ab_cat[3])
spot_Ab_sp_n2<-which(object_ab$condition==ab_cat[4])
spot_Ab_sp_n3<-which(object_ab$condition==ab_cat[5])
counts<-object_ab@assays$SCT@counts
spot_ab_sp<-list(spot_Ab_sp,spot_Ab_sp_n1,spot_Ab_sp_n2,spot_Ab_sp_n3)
a<-c()
for(i in 1:4){
  temp<-mean(counts[which(rownames(counts)=="SERPINA3"),spot_ab_sp[[i]]])
  a<-c(a,temp)
}
b<-c()
for(i in 1:4){
  per<-counts[which(rownames(counts)=="SERPINA3"),spot_ab_sp[[i]]]
  temp<-length(which(per!=0))/length(per)
  b<-c(b,temp)
}
c<-c()
for(i in 1:4){
  temp<-mean(counts[which(rownames(counts)=="C1QB"),spot_ab_sp[[i]]])
  c<-c(c,temp)
}
d<-c()
for(i in 1:4){
  per<-counts[which(rownames(counts)=="C1QB"),spot_ab_sp[[i]]]
  temp<-length(which(per!=0))/length(per)
  d<-c(d,temp)
}
e<-c()
for(i in 1:4){
  temp<-mean(counts[which(rownames(counts)=="GFAP"),spot_ab_sp[[i]]])
  e<-c(e,temp)
}
f<-c()
for(i in 1:4){
  per<-counts[which(rownames(counts)=="GFAP"),spot_ab_sp[[i]]]
  temp<-length(which(per!=0))/length(per)
  f<-c(f,temp)
}
re<-t(data.frame(e,f))
re<-t(data.frame(a,b,c,d))
write.csv(re,"AT8.csv")
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
ab_cat<-sort(unique(object_ab$condition))
spot_Ab_sp<-which(object_ab$condition==ab_cat[2])
spot_Ab_sp_n1<-which(object_ab$condition==ab_cat[3])
spot_Ab_sp_n2<-which(object_ab$condition==ab_cat[4])
spot_Ab_sp_n3<-which(object_ab$condition==ab_cat[5])
counts<-object_ab@assays$SCT@counts
spot_ab_sp<-list(spot_Ab_sp,spot_Ab_sp_n1,spot_Ab_sp_n2,spot_Ab_sp_n3)
a<-c()
for(i in 1:4){
  temp<-mean(counts[which(rownames(counts)=="SERPINA3"),spot_ab_sp[[i]]])
  a<-c(a,temp)
}
b<-c()
for(i in 1:4){
  per<-counts[which(rownames(counts)=="SERPINA3"),spot_ab_sp[[i]]]
  temp<-length(which(per!=0))/length(per)
  b<-c(b,temp)
}
c<-c()
for(i in 1:4){
  temp<-mean(counts[which(rownames(counts)=="C1QB"),spot_ab_sp[[i]]])
  c<-c(c,temp)
}
d<-c()
for(i in 1:4){
  per<-counts[which(rownames(counts)=="C1QB"),spot_ab_sp[[i]]]
  temp<-length(which(per!=0))/length(per)
  d<-c(d,temp)
}
e<-c()
for(i in 1:4){
  temp<-mean(counts[which(rownames(counts)=="GFAP"),spot_ab_sp[[i]]])
  e<-c(e,temp)
}
f<-c()
for(i in 1:4){
  per<-counts[which(rownames(counts)=="GFAP"),spot_ab_sp[[i]]]
  temp<-length(which(per!=0))/length(per)
  f<-c(f,temp)
}
re<-t(data.frame(e,f))
re<-t(data.frame(a,b,c,d))

colnames(re)<-c("AT8_sp","Level2","Leve2","Level3")
write.csv(re,"AT8_sp.csv")
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
spot_Ab_sp<-which(object_ab$condition==ab_cat[2])
spot_Ab_sp_n1<-which(object_ab$condition==ab_cat[3])
spot_Ab_sp_n2<-which(object_ab$condition==ab_cat[4])
spot_Ab_sp_n3<-which(object_ab$condition==ab_cat[5])
counts<-object_ab@assays$SCT@counts
spot_ab_sp<-list(spot_Ab_sp,spot_Ab_sp_n1,spot_Ab_sp_n2,spot_Ab_sp_n3)
a<-c()
for(i in 1:4){
  temp<-mean(counts[which(rownames(counts)=="SERPINA3"),spot_ab_sp[[i]]])
  a<-c(a,temp)
}
b<-c()
for(i in 1:4){
  per<-counts[which(rownames(counts)=="SERPINA3"),spot_ab_sp[[i]]]
  temp<-length(which(per!=0))/length(per)
  b<-c(b,temp)
}
c<-c()
for(i in 1:4){
  temp<-mean(counts[which(rownames(counts)=="C1QB"),spot_ab_sp[[i]]])
  c<-c(c,temp)
}
d<-c()
for(i in 1:4){
  per<-counts[which(rownames(counts)=="C1QB"),spot_ab_sp[[i]]]
  temp<-length(which(per!=0))/length(per)
  d<-c(d,temp)
}
re<-t(data.frame(a,b,c,d))
colnames(re)<-ab_cat[-1]
write.csv(re,"double_sp.csv")
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
spot_Ab_sp<-which(object_ab$condition==ab_cat[2])
spot_Ab_sp_n1<-which(object_ab$condition==ab_cat[3])
spot_Ab_sp_n2<-which(object_ab$condition==ab_cat[4])
spot_Ab_sp_n3<-which(object_ab$condition==ab_cat[5])
counts<-object_ab@assays$SCT@counts
spot_ab_sp<-list(spot_Ab_sp,spot_Ab_sp_n1,spot_Ab_sp_n2,spot_Ab_sp_n3)
a<-c()
for(i in 1:4){
  temp<-mean(counts[which(rownames(counts)=="SERPINA3"),spot_ab_sp[[i]]])
  a<-c(a,temp)
}
b<-c()
for(i in 1:4){
  per<-counts[which(rownames(counts)=="SERPINA3"),spot_ab_sp[[i]]]
  temp<-length(which(per!=0))/length(per)
  b<-c(b,temp)
}
c<-c()
for(i in 1:4){
  temp<-mean(counts[which(rownames(counts)=="C1QB"),spot_ab_sp[[i]]])
  c<-c(c,temp)
}
d<-c()
for(i in 1:4){
  per<-counts[which(rownames(counts)=="C1QB"),spot_ab_sp[[i]]]
  temp<-length(which(per!=0))/length(per)
  d<-c(d,temp)
}
re<-t(data.frame(a,b,c,d))
colnames(re)<-ab_cat[-1]
write.csv(re,"double.csv")
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
spot_Ab_sp<-which(object_ab$condition==ab_cat[2])
spot_Ab_sp_n1<-which(object_ab$condition==ab_cat[3])
spot_Ab_sp_n2<-which(object_ab$condition==ab_cat[4])
spot_Ab_sp_n3<-which(object_ab$condition==ab_cat[5])
counts<-object_ab@assays$SCT@counts
spot_ab_sp<-list(spot_Ab_sp,spot_Ab_sp_n1,spot_Ab_sp_n2,spot_Ab_sp_n3)
a<-c()
for(i in 1:4){
  temp<-mean(counts[which(rownames(counts)=="SERPINA3"),spot_ab_sp[[i]]])
  a<-c(a,temp)
}
b<-c()
for(i in 1:4){
  per<-counts[which(rownames(counts)=="SERPINA3"),spot_ab_sp[[i]]]
  temp<-length(which(per!=0))/length(per)
  b<-c(b,temp)
}
re<-t(data.frame(a,b))
colnames(re)<-ab_cat[-1]
write.csv(re,"w_sp.csv")
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
spot_Ab_sp<-which(object_ab$condition==ab_cat[2])
spot_Ab_sp_n1<-which(object_ab$condition==ab_cat[3])
spot_Ab_sp_n2<-which(object_ab$condition==ab_cat[4])
spot_Ab_sp_n3<-which(object_ab$condition==ab_cat[5])
counts<-object_ab@assays$SCT@counts
spot_ab_sp<-list(spot_Ab_sp,spot_Ab_sp_n1,spot_Ab_sp_n2,spot_Ab_sp_n3)
a<-c()
for(i in 1:4){
  temp<-mean(counts[which(rownames(counts)=="SERPINA3"),spot_ab_sp[[i]]])
  a<-c(a,temp)
}
b<-c()
for(i in 1:4){
  per<-counts[which(rownames(counts)=="SERPINA3"),spot_ab_sp[[i]]]
  temp<-length(which(per!=0))/length(per)
  b<-c(b,temp)
}
c<-c()
for(i in 1:4){
  temp<-mean(counts[which(rownames(counts)=="C1QB"),spot_ab_sp[[i]]])
  c<-c(c,temp)
}
d<-c()
for(i in 1:4){
  per<-counts[which(rownames(counts)=="C1QB"),spot_ab_sp[[i]]]
  temp<-length(which(per!=0))/length(per)
  d<-c(d,temp)
}
re<-t(data.frame(a,b))
colnames(re)<-ab_cat[-1]
write.csv(re,"w.csv")




#############function
expre_per<-function(pathology_files,gene_x){
  #load data
  all_file<-list.files(".")
  file<-all_file[grep(pathology_files,all_file)]
  P2_3 <- read.csv(file[grep("2-3",file)],
                   header = T,
                   row.names = 1)
  PT4857 <- read.csv(file[grep("T4857",file)],
                     header = T,
                     row.names = 1)
  P2_8 <- read.csv(file[grep("2-8",file)],
                   header = T,
                   row.names = 1)
  colnames(P2_8) <- colnames(PT4857)
  #add pathologies information to meta
  pathologies<-list(P2_3,PT4857,P2_8)
  object.2_3 <- AddMetaData(object.2_3, metadata = pathologies[[1]] )
  object_T4857 <- AddMetaData(object_T4857, metadata = pathologies[[2]])
  object.2_8 <- AddMetaData(object.2_8, metadata = pathologies[[3]])
  object_ab <- merge(object.2_3, c(object_T4857, object.2_8))
  DefaultAssay(object_ab) <- "SCT"
  #find specific spots in four conditions
  pathologies_condition<-sort(unique(object_ab$condition))[-1]
  spot_Ab<-which(object_ab$condition==pathologies_condition[1])
  spot_Ab_n1<-which(object_ab$condition==pathologies_condition[2])
  spot_Ab_n2<-which(object_ab$condition==pathologies_condition[3])
  spot_Ab_n3<-which(object_ab$condition==pathologies_condition[4])
  counts<-object_ab@assays$SCT@counts
  spot_ab_sp<-list(spot_Ab,spot_Ab_n1,spot_Ab_n2,spot_Ab_n3)
  #average expression
  expre<-c()
  for(i in 1:4){
    temp<-mean(counts[which(rownames(counts)==gene_x),spot_ab_sp[[i]]])
    print(temp)
    expre<-c(expre,temp)
  }
  #percentage
  percentage<-c()
  for(i in 1:4){
    per<-counts[which(rownames(counts)==gene_x),spot_ab_sp[[i]]]
    temp<-length(which(per!=0))/length(per)
    print(temp)
    percentage<-c(percentage,temp)
  }
  #bind
  re<-t(data.frame(expre,percentage))
  colnames(re)<-pathologies_condition
  return(re)
}

#Load data
setwd("/users/PAS1475/liuzeyi/guoqi/output")
load("object.2_3.RData")
load("object.2_8.RData")
load("object_T4857.RData")
#GFAP
setwd("/users/PAS1475/liuzeyi/guoqi/output/pathology")
re_ab<-expre_per("Aβ annotation.csv","GFAP")
re_ab_specific<-expre_per("Aβ annotation_specific.csv","GFAP")
re_at8<-expre_per("AT8 annotation.csv","GFAP")
re_at8_sp<-expre_per("AT8 annotation_specific.csv","GFAP")
re_GFAP<-gdata::cbindX(re_ab,re_ab_specific,re_at8,re_at8_sp)
setwd("//users//PAS1475//liuzeyi//guoqi//output//pathology")
write.csv(re_GFAP,"GFAP.csv")

#SERPINA3
re_ab<-expre_per("Aβ annotation.csv","SERPINA3")
re_ab_specific<-expre_per("Aβ annotation_specific.csv","SERPINA3")

re_at8<-expre_per("AT8 annotation.csv","SERPINA3")
re_at8_sp<-expre_per("AT8 annotation_specific.csv","SERPINA3")

re_abat8<-expre_per("Double annotation.csv","SERPINA3")
re_abat8_specific<-expre_per("Double annotation_specific.csv","SERPINA3")

re_W<-expre_per("Double_W_AT8 annotation.csv","SERPINA3")
re_W_sp<-expre_per("Double_W_AT8 annotation_specific.csv","SERPINA3")

re_SERPINA3<-gdata::cbindX(re_ab,re_ab_specific,re_at8,re_at8_sp,re_abat8,re_abat8_specific,re_W,re_W_sp)
setwd("//users//PAS1475//liuzeyi//guoqi//output//pathology")
write.csv(re_GFAP,"GFAP.csv")
#C1QB
re_ab<-expre_per("Aβ annotation.csv","C1QB")
re_ab_specific<-expre_per("Aβ annotation_specific.csv","C1QB")

re_at8<-expre_per("AT8 annotation.csv","C1QB")
re_at8_sp<-expre_per("AT8 annotation_specific.csv","C1QB")

re_abat8<-expre_per("Double annotation.csv","C1QB")
re_abat8_specific<-expre_per("Double annotation_specific.csv","C1QB")

#-----------new method-updata for upper code
#############function
expre_per<-function(pathology_files,gene_x){
  #load data
  all_file<-list.files(".")
  file<-all_file[grep(pathology_files,all_file)]
  P2_3 <- read.csv(file[grep("2-3",file)],
                   header = T,
                   row.names = 1)
  PT4857 <- read.csv(file[grep("T4857",file)],
                     header = T,
                     row.names = 1)
  P2_8 <- read.csv(file[grep("2-8",file)],
                   header = T,
                   row.names = 1)
  colnames(P2_8) <- colnames(PT4857)
  #add pathologies information to meta
  pathologies<-list(P2_3,PT4857,P2_8)
  object.2_3 <- AddMetaData(object.2_3, metadata = pathologies[[1]] )
  object_T4857 <- AddMetaData(object_T4857, metadata = pathologies[[2]])
  object.2_8 <- AddMetaData(object.2_8, metadata = pathologies[[3]])
  object_ab <- merge(object.2_3, c(object_T4857, object.2_8))
  DefaultAssay(object_ab) <- "SCT"
  #find specific spots in four conditions
  pathologies_condition<-sort(unique(object_ab$condition))[-1]
  spot_Ab<-which(object_ab$condition==pathologies_condition[1])
  spot_Ab_n1<-which(object_ab$condition==pathologies_condition[2])
  spot_Ab_n2<-which(object_ab$condition==pathologies_condition[3])
  spot_Ab_n3<-which(object_ab$condition==pathologies_condition[4])
  counts<-object_ab@assays$SCT@counts
  spot_ab_sp<-list(spot_Ab,spot_Ab_n1,spot_Ab_n2,spot_Ab_n3)
  #average expression
  expre<-c()
  for(i in 1:4){
    temp<-mean(counts[which(rownames(counts)==gene_x),spot_ab_sp[[i]]])
    print(temp)
    expre<-c(expre,temp)
  }
  #percentage
  percentage<-c()
  for(i in 1:4){
    per<-counts[which(rownames(counts)==gene_x),spot_ab_sp[[i]]]
    temp<-length(which(per!=0))/length(per)
    print(temp)
    percentage<-c(percentage,temp)
  }
  #bind
  re<-t(data.frame(expre,percentage))
  colnames(re)<-pathologies_condition
  return(re)
}

#Load data
setwd("/users/PAS1475/liuzeyi/guoqi/output")
load("object.2_3.RData")
load("object.2_8.RData")
load("object_T4857.RData")
#GFAP
setwd("/users/PAS1475/liuzeyi/guoqi/output/pathology")
re_ab<-expre_per("Aβ annotation.csv","GFAP")
re_ab_specific<-expre_per("Aβ annotation_specific.csv","GFAP")
re_at8<-expre_per("AT8 annotation.csv","GFAP")
re_at8_sp<-expre_per("AT8 annotation_specific.csv","GFAP")
re_GFAP<-gdata::cbindX(re_ab,re_ab_specific,re_at8,re_at8_sp)
setwd("//users//PAS1475//liuzeyi//guoqi//output//pathology")
write.csv(re_GFAP,"GFAP.csv")

#SERPINA3
re_ab<-expre_per("Aβ annotation.csv","SERPINA3")
re_ab_specific<-expre_per("Aβ annotation_specific.csv","SERPINA3")

re_at8<-expre_per("AT8 annotation.csv","SERPINA3")
re_at8_sp<-expre_per("AT8 annotation_specific.csv","SERPINA3")

re_abat8<-expre_per("Double annotation.csv","SERPINA3")
re_abat8_specific<-expre_per("Double annotation_specific.csv","SERPINA3")

re_W<-expre_per("Double_W_AT8 annotation.csv","SERPINA3")
re_W_sp<-expre_per("Double_W_AT8 annotation_specific.csv","SERPINA3")

re_SERPINA3<-gdata::cbindX(re_ab,re_ab_specific,re_at8,re_at8_sp,re_abat8,re_abat8_specific,re_W,re_W_sp)
setwd("//users//PAS1475//liuzeyi//guoqi//output//pathology")
write.csv(re_GFAP,"GFAP.csv")
#C1QB
re_ab<-expre_per("Aβ annotation.csv","C1QB")
re_ab_specific<-expre_per("Aβ annotation_specific.csv","C1QB")

re_at8<-expre_per("AT8 annotation.csv","C1QB")
re_at8_sp<-expre_per("AT8 annotation_specific.csv","C1QB")

re_abat8<-expre_per("Double annotation.csv","C1QB")
re_abat8_specific<-expre_per("Double annotation_specific.csv","C1QB")




#------------------GO for DEGs of WFS1/AT8 layer2    1.17
#load data
library("readxl")
WFS1AT8<-list()
WFS1AT8_specific<-list()
setwd("//users//PAS1475//liuzeyi//guoqi//output//pathology//6 genelists for WFS1AT8 in layer2")
for(i in 1:3){
  a<-read_excel("DEGs of WFS1AT8 v.s. level 123 in layer 2.xlsx", sheet = i)
  a<-a[-1,-1]
  colnames(a)<-a[1,]
  a<-a[-1,]
  rownames(a)<-a$Marker
  WFS1AT8[[i]] <- a
  b<-read_excel("DEGs of WFS1AT8 v.s. level 123 in layer 2-specific.xlsx", sheet = i)
  b<-b[-1,-1]
  colnames(b)<-b[1,]
  b<-b[-1,]
  rownames(b)<-b$Marker
  WFS1AT8_specific[[i]]<-b
}
names(WFS1AT8)<-c("non_level_1","non_level_2","non_level_3")
names(WFS1AT8_specific)<-c("non_level_1_specific","non_level_2_specific","non_level_3_specific")

#GO
RunGO_human <- function(Seurat.DEGs = NULL) {
  library(org.Hs.eg.db)
  library(clusterProfiler)
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
  return(GO_simplied_res)
}
GO_WFS1AT8<-list()
GO_WFS1AT8_specific<-list()
for(i in 1:3){
  GO_WFS1AT8[[i]]<-RunGO_human(WFS1AT8[[i]])
  GO_WFS1AT8_specific[[i]]<-RunGO_human(WFS1AT8_specific[[i]])
}
names(GO_WFS1AT8)<-names(WFS1AT8)
names(GO_WFS1AT8_specific)<-names(WFS1AT8_specific)

saveRDS(GO_WFS1AT8,"GO-DEGs of WFS1AT8 v.s. level 123 in layer 2.RData")
saveRDS(GO_WFS1AT8_specific,"GO-DEGs of WFS1AT8 v.s. level 123 in layer 2-specific.RData")
a<-readRDS("GO-DEGs of WFS1AT8 v.s. level 123 in layer 2.RData")
#export data
# library(xlsx,lib.loc = "/fs/scratch/PCON0022/liyang/lib")
# write.xlsx2(non_specifc[[1]], "GO-DEGs of WFS1AT8 v.s. level 123 in layer 2.xlsx", sheetName = names(non_specifc)[[1]],
#             col.names = TRUE, row.names = F, append = FALSE)
# write.xlsx2(specifc[[1]], "GO-DEGs of WFS1AT8 v.s. level 123 in layer 2-specific.xlsx", sheetName = names(specifc)[[1]],
#             col.names = TRUE, row.names = F, append = FALSE)
# for(i in 2:3){
#   write.xlsx2(non_specifc[[i]], "GO-DEGs of WFS1AT8 v.s. level 123 in layer 2.xlsx", sheetName = names(non_specifc)[[i]],
#               col.names = TRUE, row.names = F, append = TRUE)
#   write.xlsx2(specifc[[i]], "GO-DEGs of WFS1AT8 v.s. level 123 in layer 2-specific.xlsx", sheetName = names(specifc)[[i]],
#               col.names = TRUE, row.names = F, append = TRUE)
# }

#1.18 findmarks for 12 genes
#############function
#merge
findmarkers_merge<-function(pathology_files){
  #load data
  all_file<-list.files(".")
  file<-all_file[grep(pathology_files,all_file)]
  P2_3 <- read.csv(file[grep("2-3",file)],
                   header = T,
                   row.names = 1)
  PT4857 <- read.csv(file[grep("T4857",file)],
                     header = T,
                     row.names = 1)
  P2_8 <- read.csv(file[grep("2-8",file)],
                   header = T,
                   row.names = 1)
  colnames(P2_8) <- colnames(PT4857)
  #add pathologies information to meta
  pathologies<-list(P2_3,PT4857,P2_8)
  object.2_3 <- AddMetaData(object.2_3, metadata = pathologies[[1]] )
  object_T4857 <- AddMetaData(object_T4857, metadata = pathologies[[2]])
  object.2_8 <- AddMetaData(object.2_8, metadata = pathologies[[3]])
  object_ab <- merge(object.2_3, c(object_T4857, object.2_8))
  DefaultAssay(object_ab) <- "SCT"
  Idents(object_ab) <- object_ab$condition
  ab_cat <- sort(unique(object_ab$condition))
  ab_cat <- ab_cat[-1]
  find_markers <- data.frame()
  for (i in ab_cat[-1]) {
    markers <-
      FindMarkers(
        object_ab,
        ident.1 = ab_cat[1],
        ident.2 = i,
        group.by = "condition",
        verbose = TRUE,
        assay = "SCT",
        only.pos = F,
        features = genelist
      )
    markers <- cbind(rownames(markers), markers)
    markers$level <- as.character(rep(i, nrow(markers)))
    find_markers <- rbind(find_markers, markers)
  }
  #object.2_3
  find_markers_object.2_3 <- data.frame()
  for (i in ab_cat[-1]) {
    markers <-
      FindMarkers(
        object.2_3,
        ident.1 = ab_cat[1],
        ident.2 = i,
        group.by = "condition",
        verbose = TRUE,
        assay = "SCT",
        only.pos = F,
        features = genelist
      )
    markers <- cbind(rownames(markers), markers)
    markers$level <- as.character(rep(i, nrow(markers)))
    find_markers_object.2_3 <- rbind(find_markers_object.2_3, markers)
  }
  #object_T4857
  find_markers_object_T4857 <- data.frame()
  for (i in ab_cat[-1]) {
    markers <-
      FindMarkers(
        object_T4857,
        ident.1 = ab_cat[1],
        ident.2 = i,
        group.by = "condition",
        verbose = TRUE,
        assay = "SCT",
        only.pos = F,
        features = genelist
      )
    markers <- cbind(rownames(markers), markers)
    markers$level <- as.character(rep(i, nrow(markers)))
    find_markers_object_T4857 <- rbind(find_markers_object_T4857, markers)
  }
  #object.2_8
  find_markers_object.2_8 <- data.frame()
  for (i in ab_cat[-1]) {
    markers <-
      FindMarkers(
        object.2_8,
        ident.1 = ab_cat[1],
        ident.2 = i,
        group.by = "condition",
        verbose = TRUE,
        assay = "SCT",
        only.pos = F,
        features = genelist
      )
    markers <- cbind(rownames(markers), markers)
    markers$level <- as.character(rep(i, nrow(markers)))
    find_markers_object.2_8 <- rbind(find_markers_object.2_8, markers)
  }
  
  rownames(find_markers) <- NULL
  colnames(find_markers)[1] <- c("Marker")
  list_re<-list(find_markers,find_markers_object.2_3,find_markers_object_T4857,find_markers_object.2_8)
  return(list_re)
}

genelist<-c("SLC1A3","P2RY12","KIF5A","NeuN","SNCG","STMN2","CSRP1","PLP1","GLUL","GFAP","PAQR6","MBP")


# GO for wfs
setwd("/users/PAS1475/liuzeyi/guoqi/output/pathology/6 genelists for WFS1AT8 in layer2/")
library("readxl")
specific<-list()
for(i in 1:3){
  temp<-read_excel("DEGs of WFS1AT8 v.s. level 123 in layer 2-specific-separate level.xlsx", sheet = i)
  colnames(temp)<-temp[2,]
  temp<-temp[-c(1:2),-1]
  rownames(temp)<-temp$Marker
  specific[[i]] <- temp
}
specific<-excel_sheets(path = "DEGs of WFS1AT8 v.s. level 123 in layer 2-specific-separate level.xlsx")

non_specific<-list()
for(i in 1:3){
  temp<-read_excel("DEGs of WFS1AT8 v.s. level 123 in layer 2-spearate level.xlsx", sheet = i)
  colnames(temp)<-temp[2,]
  temp<-temp[-c(1:2),-1]
  rownames(temp)<-temp$Marker
  non_specific[[i]] <- temp
}
non_specific<-excel_sheets(path = "DEGs of WFS1AT8 v.s. level 123 in layer 2-spearate level.xlsx")

#creat merged object with pathological information 1.26
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
#each object
colname<-c("Ab_sp","AT8_sp","AT8_Ab_sp")
pathology<-list()
object_list<-list(object.2_3,object.2_8,object_T4857,object.1_1,object.18_64,object.2_5)
pathology[[1]]<-gdata::cbindX(Ab2_3_sp,AT8_2_3_sp,D_2_3_sp)
pathology[[2]]<-gdata::cbindX(Ab2_8_sp,AT8_2_8_sp,D_2_8_sp)
pathology[[3]]<-gdata::cbindX(Ab4857_sp,AT8_T4857_sp,D_T4857_sp)
pathology[[4]]<-data.frame(rep("control",nrow(object.1_1@meta.data)),rep("control",nrow(object.1_1@meta.data))
                           ,rep("control",nrow(object.1_1@meta.data)))
pathology[[5]]<-data.frame(rep("control",nrow(object.18_64@meta.data)),rep("control",nrow(object.18_64@meta.data))
                           ,rep("control",nrow(object.18_64@meta.data)))
pathology[[6]]<-data.frame(rep("control",nrow(object.2_5@meta.data)),rep("control",nrow(object.2_5@meta.data))
                           ,rep("control",nrow(object.2_5@meta.data)))
object.re<-list()
for(i in 1:6){
  add_meta<-pathology[[i]]
  colnames(add_meta)<-colname
  meta <- cbind(object_list[[i]]@meta.data,add_meta)
  object_list[[i]]@meta.data<-meta
  object.re[[i]]<-object_list[[i]]
}
names(object.re)<-c("object.2_3","object.2_8","object_T4857","object.1_1","object.18_64","object.2_5")
temp<-SplitObject(sample_6,split.by = "group")
sample6_pathology<-merge(object.re[[1]],c(object.re[[2]],object.re[[3]],object.re[[4]],object.re[[5]],object.re[[6]]))
qsave(sample6_pathology, file = './sample6_pathology.qs')
sample6<-qread("sample6_pathology.qs")
#pathology ratio input:pathology output dataframe(row=layer col=sample) 
object<-list(object.2_3,object.2_8,object_T4857)
pathology_ratio<-function(pathology,pathology_list){
  layer<-sort(unique(object.2_3@meta.data$Layer))
  per<-c()
  num_path<-c()
  num_layer<-c()
  path_per<-list()
  for(j in 1:3){
    for(i in 1:7){
      length_layer_spot<-length(which(object[[j]]@meta.data$Layer==layer[i]))
      length_path_spot<-length(intersect(which(object[[j]]@meta.data$Layer==layer[i]),which(pathology_list[[j]]$condition==pathology)))
      percentage<-paste0(round((length_path_spot/length_layer_spot)*100,2),"%")
      per<-c(per,percentage)
      num_path<-c(num_path,length_path_spot)
      num_layer<-c(num_layer,length_layer_spot)
    }
    per_dataframe<-as.data.frame(t(data.frame(Num_spot_layer=num_layer,Num_spot_pathology_layer=num_path,Percentage=per)))
    colnames(per_dataframe)<-layer
    path_per[[j]]<-per_dataframe
    per<-c()
    num_layer<-c()
    num_path<-c()
  }
  return(path_per)
}
#Ab
pathology_list<-list(Ab2_3,Ab2_8,T4857)
path_per_ab<-pathology_ratio("Aβ",pathology_list)
names(path_per_ab)<-c("sample.2_3","sample.2_8","sample_T4857")
#Ab_sp
pathology_list<-list(Ab2_3_sp,Ab2_8_sp,Ab4857_sp)
path_per_ab_sp<-pathology_ratio("Aβ",pathology_list)
names(path_per_ab_sp)<-names(path_per_ab)
#At8
pathology_list_at8<-list(AT8_2_3,AT8_2_8,AT8_T4857)
path_per_AT8<-pathology_ratio("AT8",pathology_list_at8)
names(path_per_AT8)<-names(path_per_ab)
#AT8_sp
pathology_list_at8_sp<-list(AT8_2_3_sp,AT8_2_8_sp,AT8_T4857_sp)
path_per_AT8_sp<-pathology_ratio("AT8",pathology_list_at8_sp)
names(path_per_AT8_sp)<-names(path_per_ab)

#Double
pathology_list_d<-list(D_2_3,D_2_8,D_T4857)
path_per_d<-pathology_ratio("Aβ/AT8",pathology_list_d)
names(path_per_d)<-names(path_per_ab)

#Double_sp
pathology_list_d_sp<-list(D_2_3_sp,D_2_8_sp,D_T4857_sp)
path_per_d_sp<-pathology_ratio("Aβ/AT8",pathology_list_d_sp)
names(path_per_d_sp)<-names(path_per_ab)

#W
pathology_list_w<-list(D_2_3,D_2_8,D_T4857)
path_per_w<-pathology_ratio("WFS1/AT8 double positive",pathology_list_w)
names(path_per_w)<-names(path_per_ab)
#W_sp
pathology_list_w_sp<-list(D_2_3_W_AT8_sp,D_2_8_W_AT8_sp,D_T4857_W_AT8_sp)
path_per_w_sp<-pathology_ratio("WFS1/AT8 double positive",pathology_list_w_sp)
names(path_per_w_sp)<-names(path_per_ab)




#output
dyn.load(
  "/usr/lib/jvm/java-1.6.0-openjdk-1.6.0.41.x86_64/jre/lib/amd64/server/libjvm.so"
)
library(xlsx, lib.loc = "/fs/scratch/PCON0022/liyang/lib")
setwd("/users/PAS1475/liuzeyi/guoqi/output/pathology/percentage_pathology_spot")
list_all<-list(path_per_ab,path_per_ab_sp,path_per_AT8,path_per_AT8_sp,path_per_d,path_per_d_sp,path_per_w,path_per_w_sp)
names(list_all)<-c("Ab","Ab_sp","AT8","AT8_sp","Double_Ab_AT8","Double_Ab_AT8_sp","Double_W_AT8","Double_W_AT8_sp")
for(i in 1:8){
  write.xlsx2(list_all[[i]][[1]], paste0(names(list_all)[i],".xlsx"), sheetName = names(list_all[[1]])[1],
              col.names = TRUE, row.names = T, append = FALSE)
  for(j in 2:3){
    write.xlsx2(list_all[[i]][[j]], paste0(names(list_all)[i],".xlsx"), sheetName = names(list_all[[1]])[j],
                col.names = TRUE, row.names = T, append = TRUE)
  }
}


#heatmap for pathology 2.27
#load data
setwd("/users/PAS1475/liuzeyi/guoqi/output/pathology/heatmap")
data<-read_excel("No.11, Level 3 pathology related heat map genes (Logfc rank) .xlsx",sheet=3)
data<-data[,-4]
colnames(data)<-c("Abeta_up", "AT8_up", "Abeta_AT8_up", "Abeta_down", "AT8_down")
re<-data.frame()
for(i in 1:6){
  temp<-na.omit(data[,i])
  ann<-as.data.frame(rep(colnames(temp),nrow(temp)))
  comb<-cbind(temp,ann)
  colnames(comb)<-c("Marker","Pathology")
  re<-rbind(re,comb)
}
write.csv(re,"interated infor for heatmap.csv")

#heatmap
load("/users/PAS1475/liuzeyi/guoqi/output/sample6_clusters.RData")
expr_data <- as.data.frame(sample6.combined@assays$SCT@data[re$Marker,])
expr_data<-rbind(expr_data,layer=sample6.combined$Layer)
expr_data<-as.data.frame(t(expr_data))
temp_data<-data.frame(rep(0,7))
for(j in 1:length(re$Marker)){
  gene_data<-c()
  for(i in sort(unique(expr_data$layer))){
    temp_gene<-mean(as.numeric(expr_data[which(expr_data$layer==i),j]))
    gene_data<-c(gene_data,temp_gene)
  }
  temp_data<-cbind(temp_data,gene_data)
}
temp_data<-temp_data[,-1]
rownames(temp_data)<-sort(unique(expr_data$layer))
colnames(temp_data)<-re$Marker
heatmap_data<-t(temp_data)
#plot
my_colour = list(
  gene_cat = c( Abeta_up = brewer.pal(n = 7, name = "Dark2")[1],
                AT8_up  = brewer.pal(n = 7, name = "Dark2")[2],
                Abeta_AT8_up = brewer.pal(n = 7, name = "Dark2")[3],
                Abeta_down = brewer.pal(n = 7, name = "Dark2")[4],
                AT8_down = brewer.pal(n = 7, name = "Dark2")[5]))
gene_cat_val <- as.numeric(table(re$Pathology))                
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
  filename = "heatmap_pathology.png",
  device = "png",
  dpi = 150,
  width = 5,
  height = 15,
  units = "in"
)

#Dotplot for pathology nonlevel3
setwd("/users/PAS1475/liuzeyi/guoqi/output/pathology/GO for level3")
#load data
library("readxl")
list_level3<-list()
for(i in 1:6){
  temp<-read_excel("No.11 Go pathway level3 pathology.xlsx", sheet = i)
  temp_sort<-temp[with(temp, order(Count)), ]
  temp_sort$Description<-factor(temp_sort$Description,levels=(temp_sort$Description))
  count_sum<-as.numeric(unlist(strsplit(temp_sort$GeneRatio[1],split="/"))[2])
  temp_sort$GeneRatio_num<-c(temp_sort$Count)/count_sum
  list_level3[[i]]<- temp_sort
}
names(list_level3)<-excel_sheets(path = "No.11 Go pathway level3 pathology.xlsx")


library(ggplot2)
setwd("/users/PAS1475/liuzeyi/guoqi/output/picture/enrichment_dot/")
for(i in 1:6){
  g<-ggplot(list_level3[[i]], # you can replace the numbers to the row number of pathway of your interest
            aes(x = GeneRatio_num, y = Description)) + 
    geom_point(aes(size = Count, color = p.adjust)) +
    theme_bw(base_size = 14) +
    theme(axis.text = element_text(size = 14, face = "bold"),
    )+
    scale_colour_gradient(limits=NULL, low="red",high="blue") +
    ylab(NULL) +
    ggtitle("GO enrichment")
  ggsave(
    plot = g,
    filename = paste0(names(list_level3)[i],"_enrichment_N11.tiff"),
    device = "tiff",
    dpi = 150,
    width = 10,
    height = 10,
    units = "in"
  )
}
for(i in 1:6){
  g<-ggplot(list_level3[[i]], # you can replace the numbers to the row number of pathway of your interest
            aes(x = GeneRatio_num, y = Description)) + 
    geom_point(aes(size = Count, color = p.adjust)) +
    theme_bw(base_size = 14) +
    theme(axis.text = element_text(size = 14, face = "bold"),
    )+
    scale_colour_gradient(limits=NULL,high="blue",low="red") +
    ylab(NULL) +
    ggtitle("GO enrichment")
  ggsave(
    plot = g,
    filename = paste0(names(list_level3)[i],"_enrichment_N11.png"),
    device = "png",
    dpi = 150,
    width = 10,
    height = 10,
    units = "in"
  )
}
ggplot(temp, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio_num, y = Description)) + 
  geom_point(aes(size = Count, color = p.adjust)) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(size = 14, face = "bold"),
  )+
  scale_fill_continuous(guide = guide_legend())
  scale_colour_gradient(limits=NULL,high="blue",low="red")
