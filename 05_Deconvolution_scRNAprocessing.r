library("Seurat")
#BiocManager::install("DropletUtils")
#library("DropletUtils")
library("Matrix")
setwd("/fs/ess/PCON0022/guoqi/Data_raw/Hodge_human_MTG")
dyn.load("/apps/hdf5-serial/gnu/9.1/1.12.0/lib/libhdf5_hl.so.200")
library(hdf5r)
exon<-read.csv("human_MTG_2018-06-14_exon-matrix.csv")
gene_refer<-read.csv("human_MTG_2018-06-14_genes-rows.csv")
exon_gene<-gene_refer[match(exon$X,gene_refer$entrez_id),c(1,3)]
identical(exon_gene$entrez_id,exon$X)
rownames(exon)<-exon_gene$gene
exon<-exon[,-1]
nature_mtg <-
  CreateSeuratObject(
    counts = exon,
    project = "hodge",
    min.cells = 0,
    min.features = 0
  )
meta<-nature_mtg@meta.data
write10xCounts(
  '/fs/ess/PCON0022/guoqi/Data_raw/Hodge_human_MTG/mtg_nature.h5',
  nature_mtg@assays$RNA@counts,
  barcodes = colnames(nature_mtg),
  gene.id = rownames(nature_mtg),
  gene.symbol = rownames(nature_mtg),
  gene.type = "Gene Expression",
  overwrite = FALSE,
  type = c("HDF5"),
  genome = "unknown",
  version = c("2"),
  chemistry = "Single Cell 3' v3",
  original.gem.groups = 1L,
  library.ids = "custom"
)
# #other method
# devtools::install_github(repo = "samuel-marsh/scCustomize")
# library("scCustomize")
# exon_matrix<-as.matrix(exon)
# #exon_matrix<-as(exon,"dgCMatrix")
# write.csv(exon_matrix,"gene_expression.csv")
# Create_10X_H5(raw_data_file_path = "gene_expression.csv", source_type = "data.frame", save_file_path = "./",
#               save_name = "mtg_nature.h5")
# #third method
# devtools::install_github('JiekaiLab/RIOH5@HEAD')

#test
test_mtg <- Read10X_h5("mtg_nature.h5")
test_mtg <- CreateSeuratObject(counts = test_mtg)
#meta information
cell_inf<-read.csv("/fs/ess/PCON0022/guoqi/Data_raw/Hodge_human_MTG/human_MTG_2018-06-14_samples-columns.csv")
celltype_res<-c()
celltype<-strsplit(cell_inf$cluster," ")
for(i in 1:length(celltype)){
  celltype_temp<-celltype[[i]][1]
  celltype_res<-c(celltype_res,celltype_temp)
}
cellannotation<-data.frame(barcode=cell_inf$sample_name,celltype=celltype_res)
write.csv(cellannotation,"MTG_annotation.csv",row.names=FALSE,quote=FALSE)
mtg<-read.csv("/fs/ess/PCON0022/guoqi/Data_raw/Hodge_human_MTG/MTG_annotation.csv")
no_cell<-cellannotation$barcode[which(cellannotation$celltype=="no")]
mtg_pure<-mtg[-which(mtg$celltype=="no"),]
write.csv(mtg_pure,"/fs/ess/PCON0022/guoqi/Data_raw/Hodge_human_MTG/MTG_pure_annotation.csv",row.names=FALSE,quote=FALSE)
325
#
nature_mtg$celltype<-cellannotation$celltype
nature_mtg_true<-subset(nature_mtg, subset = celltype != "no")
write10xCounts(
  '/fs/ess/PCON0022/guoqi/Data_raw/Hodge_human_MTG/mtg_nature_pure.h5',
  nature_mtg_true@assays$RNA@counts,
  barcodes = colnames(nature_mtg_true),
  gene.id = rownames(nature_mtg_true),
  gene.symbol = rownames(nature_mtg_true),
  gene.type = "Gene Expression",
  overwrite = FALSE,
  type = c("HDF5"),
  genome = "unknown",
  version = c("2"),
  chemistry = "Single Cell 3' v3",
  original.gem.groups = 1L,
  library.ids = "custom"
)
#original 15928 cells, now 15603
#barcode_information
example<-read.csv("/fs/ess/PCON0022/Yuzhou/NSF_2022/Integrate_barcode_annotation.csv")

example_mtg<-read.csv("/fs/ess/PCON0022/guoqi/Data_raw/Hodge_human_MTG/MTG_pure_annotation.csv")
