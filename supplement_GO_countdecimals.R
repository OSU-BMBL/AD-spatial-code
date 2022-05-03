#problem: count are decimals- This decimals is just an scale, example, legend. Not stand for the real count number.
# If you produce the same length interval, sometimes you'll get a decimal, sometimes you'll get an integer
library(multienrichjam)# It contains the function enrichDF2enrichResult, which can generate a GO object from a specific pathway dataframe
#Generate GO object
#load target data
library(readxl)
setwd("/users/PAS1475/liuzeyi/guoqi/output/WGCNA/4 given modules for GO plot")
module4<-list()
for(i in 1:4){
  temp<-read_excel("No.10, 4 modules GO pathway.xlsx", sheet = i)
  module4[[i]]<- temp
}
names(module4)<-excel_sheets(path = "No.10, 4 modules GO pathway.xlsx")

multi<-enrichDF2enrichResult(module4[["Pink"]])

#plot
p<-dotplot(multi,font.size=20)

#solve count decimals problem
min = min(target_data$Count)
max = max(target_data$Count)
# 10 is length
step = ceiling((max-min)/6)
p <- p + scale_size_continuous(range = c(3,6),breaks=seq(min, max, by=step))
p<-p+ theme(legend.key.size = unit(1.0, 'cm'),
           legend.title = element_text(size=18), #change legend title font size
           legend.text = element_text(size=18))
setwd("/users/PAS1475/liuzeyi/guoqi/output/picture/enrichment_dot")
ggsave(plot = p,
       filename = paste0(names(module4)[i],"_enrichment_NO10.tiff"),
       device = "tiff",
       dpi = 150,
       width = 16,
       height = 12,
       units = "in")

#-----------------------------GO for level3
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


#plot----abeta_down
abeta_down<-enrichDF2enrichResult(list_level3[[2]])
p<-dotplot(abeta_down,font.size=20)

#solve count decimals problem
min = min(list_level3[[2]]$Count)
max = max(list_level3[[2]]$Count)
# 10 is length
step = ceiling((max-min)/3)
p <- p + scale_size_continuous(breaks=seq(min, max, by=step))
p<-p+ theme(legend.key.size = unit(1.0, 'cm'),
            legend.title = element_text(size=18), #change legend title font size
            legend.text = element_text(size=18))
setwd("/users/PAS1475/liuzeyi/guoqi/output/pathology/GO for level3")
ggsave(plot = p,
       filename = paste0(names(list_level3)[2],"_enrichment_N11.tiff"),
       device = "tiff",
       dpi = 150,
       width = 10,
       height = 12,
       units = "in")

#plot----abeta_up
abeta_up<-enrichDF2enrichResult(list_level3[[1]])
#solve count decimals problem
min = min(list_level3[[1]]$Count)
max = max(list_level3[[1]]$Count)
# 10 is length
step = ceiling((max-min)/3)
p <- dotplot(abeta_up,font.size=20) 
p<-p+ theme(legend.key.size = unit(1.0, 'cm'),
            legend.title = element_text(size=18), #change legend title font size
            legend.text = element_text(size=18))+
  scale_size_continuous(breaks=seq(min, max, by=step))
setwd("/users/PAS1475/liuzeyi/guoqi/output/pathology/GO for level3")
ggsave(plot = p,
       filename = paste0(names(list_level3)[1],"_enrichment_NO11.tiff"),
       device = "tiff",
       dpi = 150,
       width = 11,
       height = 12,
       units = "in")

#plot----abeta_at8_up
abat_up<-enrichDF2enrichResult(list_level3[[5]])
#solve count decimals problem
min = min(list_level3[[5]]$Count)
max = max(list_level3[[5]]$Count)
# 10 is length
step = ceiling((max-min)/3)
p <- dotplot(abat_up,font.size=20) 
p<-p+ theme(legend.key.size = unit(1.0, 'cm'),
            legend.title = element_text(size=18), #change legend title font size
            legend.text = element_text(size=18))+
  scale_size_continuous(breaks=seq(min, max, by=step))
setwd("/users/PAS1475/liuzeyi/guoqi/output/pathology/GO for level3")
ggsave(plot = p,
       filename = paste0(names(list_level3)[5],"_enrichment_NO11.tiff"),
       device = "tiff",
       dpi = 150,
       width = 10,
       height = 12,
       units = "in")

#plot----abeta_at8_down
abat_down<-enrichDF2enrichResult(list_level3[[6]])
#solve count decimals problem
min = min(list_level3[[6]]$Count)
max = max(list_level3[[6]]$Count)
# 10 is length
step = ceiling((max-min)/3)
p <- dotplot(abat_down,font.size=20) 
p<-p+ theme(legend.key.size = unit(1.0, 'cm'),
            legend.title = element_text(size=18), #change legend title font size
            legend.text = element_text(size=18))+
  scale_size_continuous(breaks=seq(min, max, by=step))
setwd("/users/PAS1475/liuzeyi/guoqi/output/pathology/GO for level3")
ggsave(plot = p,
       filename = paste0(names(list_level3)[6],"_enrichment_NO11.tiff"),
       device = "tiff",
       dpi = 150,
       width = 11,
       height = 12,
       units = "in")
