#correlation plot
cor<-read.csv("correlation.csv")
ggplotRegression <- function (fit) {
  require(ggplot2)
  ggplot(cor$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    theme_classic2() +
    labs(title = paste("R =",signif(sqrt(summary(fit)$r.squared), 3),
                       #"; R^2.adj = ",signif(summary(fit)$adj.r.squared, 3),
                       "; P =",signif(summary(fit)$coef[2,4], 3)))
}
colnames(cor)<-c("gene","Foldchange_Visium","Foldchange_RNAscope")
#p<-ggplot(cor, aes(Foldchange_ST, Foldchange_RNAscope))
library(ggrepel)
p1<-ggplot(cor, aes(Foldchange_Visium, Foldchange_RNAscope),)+
  geom_point(size=7)+stat_smooth(method = "lm", col = "blue",se=F)+
  geom_text_repel(aes(Foldchange_Visium, Foldchange_RNAscope, label = gene,size=20))+
  theme_classic2()+
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25,face="bold"),legend.position = "none")
setwd("/fs/ess/PCON0022/guoqi/NC/Revision/codes")
ggsave(
  plot = p1,
  filename = "correlation_label.tiff",
  device = "tiff",
  dpi = 150,
  width = 10,
  height = 10,
  units = "in"
)
p2<-ggplot(cor, aes(Foldchange_Visium, Foldchange_RNAscope),)+
  geom_point(size=7)+stat_smooth(method = "lm", col = "blue",se=F)+
 # geom_text_repel(aes(Foldchange_Visium, Foldchange_RNAscope, label = gene,size=20))+
  theme_classic2()+
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25,face="bold"),legend.position = "none")
ggsave(
  plot = p2,
  filename = "correlation_withoutlabel.tiff",
  device = "tiff",
  dpi = 150,
  width = 10,
  height = 10,
  units = "in"
)

log_cor<-cor
log_cor$Foldchange_Visium<-log10(cor$Foldchange_Visium)
log_cor$Foldchange_RNAscope<-log10(cor$Foldchange_RNAscope)
cor(log_cor$Foldchange_Visium, log_cor$Foldchange_RNAscope, method = c("pearson"))
cor(cor$Foldchange_Visium, cor$Foldchange_RNAscope, method = c("pearson"))