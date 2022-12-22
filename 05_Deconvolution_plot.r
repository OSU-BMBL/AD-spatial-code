#spatial plot
dyn.load("/apps/hdf5-serial/gnu/9.1/1.12.0/lib/libhdf5_hl.so.200")
library(hdf5r)
#sptial plot based on cell proportion
setwd("/fs/ess/PCON0022/guoqi/NC/output")
#load spatial object
load("object.2_8_n.RData")
load("object_T4857_n.RData")
load("object.2_3_n.RData")
load("object.2_5_n.RData")
load("object.18_64_n.RData")
load("object.1_1_n.RData")
load("/fs/ess/PCON0022/guoqi/NC/output/object.1_7_n.RData")
load("/fs/ess/PCON0022/guoqi/NC/output/object.2_10_n.RData")

my.all.object.list <- list(
  object.1_1 = object.1_1,
  object.1_7 = object.1_7,
  object.18_64 = object.18_64,
  object.2_10 = object.2_10,
  object.2_3 = object.2_3,
  object.2_5 = object.2_5,
  object.2_8 = object.2_8,
  object_T4857 = object_T4857
)
#load deconvolution results
#load csv
setwd("/fs/ess/PCON0022/guoqi/NC-snrna/Spatial/Deconvolution_6scRNA/Deconvolution_results")
decon_re<-list.files(".",pattern = ".csv")
decon_list<-list()
for(i in 1:8){
  decon_list[[i]]<-read.csv(decon_re[i],row.names = 1)
}
names(decon_list)<-decon_re
#add cell proportion meta
#function
Spatial_pro<-function(object,cell_name,sample){
  xy.cor <- object@images$slice1@coordinates 
  tmp.meta <- object@meta.data 
  tmp.meta <- tmp.meta[rownames(xy.cor),] 
  xy.cor <- as.data.frame(cbind(xy.cor, tmp = tmp.meta[,cell_name])) 
  colnames(xy.cor)[ncol(xy.cor)] <- cell_name 
  p<- ggplot(xy.cor,aes(x = col, y = row,color = as.numeric(xy.cor[,cell_name])))+  
    geom_point(size = 3) + theme_void() + scale_color_gradient(low="#F6F4F2", high="#A44908") + #scale_color_gradientn(colors = colorRampPalette(c("#000000","#0185E1","#06C501","#DD7105","#C30303"))(n =5))+  
    labs(col = cell_name,title = cell_name)+scale_y_reverse() + scale_x_reverse()#+ coord_flip() #+scale_y_reverse()  
  ggsave(p, filename = paste0("/fs/ess/PCON0022/guoqi/NC-snrna/Spatial/Deconvolution_6scRNA/Spatial_plot/",sample,cell_name,".tiff"), 
         device ="tiff", dpi = 150 , width = 10.5, height = 10) 
}
#result
spatial_decon_list<-list()
celltype<-colnames(decon_list[[1]])
for(i in 1:8){
  spatial_decon_list[[i]]<-AddMetaData(my.all.object.list[[i]],decon_list[[i]])
  sample<-paste0(as.character(spatial_decon_list[[i]]$patientID[1]),"_")
  for(j in 1:7){
    Spatial_pro(spatial_decon_list[[i]],celltype[j],sample)
  }
}

#bar plot
#-----------------------------------barplot
library(reshape2)
library(dplyr)
setwd("/fs/ess/PCON0022/guoqi/NC-snrna/Spatial/Deconvolution_6scRNA/Deconvolution_results")
samp_file<-list.files(".",pattern = ".csv")
samp_file<-samp_file[-c(2,4)]
setwd("/fs/project/PAS1475/Yuzhou_Chang/Fu_lab_data/Visium/raw_data/Layer_annotation")
ann_file<-list.files(".",pattern = ".csv")

#load file and delete noise spot
samp_list<-list()
ann_list<-list()
for(i in 1:6){
  setwd("/fs/ess/PCON0022/guoqi/NC-snrna/Spatial/Deconvolution_6scRNA/Deconvolution_results")
  samp_list[[i]]<-read.csv(samp_file[i],row.names = 1)
  setwd("/fs/project/PAS1475/Yuzhou_Chang/Fu_lab_data/Visium/raw_data/Layer_annotation")
  ann_list[[i]]<-read.csv(ann_file[i])
  print(identical(rownames(samp_list[[i]]),ann_list[[i]]$Barcode))
  noise<-which(ann_list[[i]][,2]=="Noise")
  if(length(noise)!=0){
    ann_list[[i]]<-ann_list[[i]][-noise,]
    samp_list[[i]]<-samp_list[[i]][-noise,]
    print(identical(rownames(samp_list[[i]]),ann_list[[i]]$Barcode))
  }
  else{
    ann_list[[i]]<-ann_list[[i]]
    samp_list[[i]]<-samp_list[[i]]
  }
}

#calculate proportion of each cell type
layer<-sort(unique(ann_list[[1]]$Layer.annotation))
data<-list()
data.melt<-list()
proportion<-list()
for(i in 1:6){
  print(i)
  data[[i]]<-cbind(samp_list[[i]],ann_list[[i]][,2])
  colnames(data[[i]])[8]<-c("layer")
  data[[i]]$id<-rownames(data[[i]])
  data.melt[[i]]<-melt(data[[i]],id.var=c("id","layer"))
  proportion[[i]] <- data.melt[[i]] %>%
    group_by(layer,variable) %>%
    summarize(Mean = mean(value, na.rm=TRUE))
}
names(proportion)<-c("1-1","18-64","2-3","2-5","2-8","T4857")

#plot
setwd("/fs/ess/PCON0022/guoqi/NC-snrna/Spatial/Deconvolution_6scRNA/barplot")
for(i in 1:6){
  p<-ggplot(proportion[[i]],aes(layer,Mean,fill=variable))+
    
    geom_bar(stat="identity",position="fill")+
    
    ggtitle("")+
    
    theme_bw()+
    
    theme(axis.ticks.length=unit(0.5,'cm'))
  
  ggsave(
    plot = p,
    filename = paste(names(proportion)[i],"_deconvolution_barplot.tiff",sep = ""),
    device = "tiff",
    dpi = 150,
    width = 10,
    height = 10,
    units = "in"
  )
}
setwd("/fs/ess/PCON0022/guoqi/NC-snrna/Spatial/Deconvolution_6scRNA/cell_proportion")
#___________output_specific format
for(i in 1:6){
  t=tidyr::spread(#better
    data=proportion[[i]],
    key=layer,
    value=Mean
  )
  if(ncol(t)!=8){
    print(i)
    colnames(t)<-c("",'PERLayer2',"PERLayer3","PERLayer4","PERLayer5","PERLayer6","PERWhiteMatter")
  }
  else{
    colnames(t)<-c("","PERLayer1",'PERLayer2',"PERLayer3","PERLayer4","PERLayer5","PERLayer6","PERWhiteMatter")
  }
  write.csv(t,paste(names(proportion)[i],"_cell_proportion.csv"),row.names = F)
}




#----------------beta test
#result
setwd('/fs/ess/PCON0022/guoqi/NC-snrna/Spatial/Deconvolution_6scRNA/cell_proportion')
all_file <- list.files(".")
temp <- read.csv(all_file[1])
beta_regression <- function(pre_data) {
  library("betareg")
  library(dplyr)
  #construct pre_data
  subpopulation <- c()
  layer <- c()
  proportion <- c()
  for (i in 1:nrow(pre_data)) {
    for (j in 2:ncol(pre_data)) {
      sp <- pre_data$X[i]
      subpopulation <- c(subpopulation, sp)
      la <- colnames(pre_data)[j]
      layer <- c(layer, la)
      pro <- pre_data[i, j]
      proportion <- c(proportion, pro)
    }
  }
  data <- data.frame(subpopulation, layer, proportion)
  data[which(data$layer == "PERLayer1"), 2] <- 1
  data[which(data$layer == "PERLayer2"), 2] <- 2
  data[which(data$layer == "PERLayer3"), 2] <- 3
  data[which(data$layer == "PERLayer4"), 2] <- 4
  data[which(data$layer == "PERLayer5"), 2] <- 5
  data[which(data$layer == "PERLayer6"), 2] <- 6
  data[which(data$layer == "PERWhiteMatter"), 2] <- 7
  data$layer <- as.numeric(data$layer)
  #beta regression
  subcelltype <- sort(unique(data$subpopulation))
  br.summary.list <- list()
  br.coefmat.list <- list()
  for (type in subcelltype) {
    br.model <- betareg(
      formula = proportion ~ layer | layer,
      data = data %>%
        dplyr::filter(subpopulation == type),
      link = 'logit',
      type = 'BC'
    )
    br.summary <- summary(br.model)
    br.summary.list[[type]] <- br.summary
    br.coefmat <- as.data.frame(br.summary$coefficients$mean) %>%
      rownames_to_column(var = 'Term') %>%
      dplyr::mutate(Subpopulation = type) %>%
      dplyr::select(Subpopulation, Term, everything())
    br.coefmat.list[[type]] <- br.coefmat
  }
  br.coefmat.combined <- bind_rows(br.coefmat.list) %>%
    dplyr::filter(!grepl('Intercept', Term))
  br.coefmat.combined$p.adj.holm <-
    p.adjust(br.coefmat.combined$`Pr(>|z|)`, method = 'holm')
  return(br.coefmat.combined)
}
#----output
for(i in 1:6){
  temp<-read.csv(all_file[i])
  a<-beta_regression(temp)
  a<-a[,c(1,5,6,7)]
  colnames(a)<-c("subpopulation","z.value","p.value","p.adj.holm")
  write.csv(a,paste("/fs/ess/PCON0022/guoqi/NC-snrna/Spatial/Deconvolution_6scRNA/betatest/",
                    all_file[i],"beta_test.csv",sep = ""))
}
