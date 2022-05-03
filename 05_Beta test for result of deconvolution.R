#betaregression
install.packages('betareg')
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
#result
setwd('/fs/scratch/PCON0022/LZY/AD/AD_deconvolution_result/')
all_file <- list.files(".")
for (i in 1:length(all_file)) {
  setwd(all_file[i])
  all_files <- list.files(".")
  for (j in all_files) {
    temp <- read.csv(j)
    a <- beta_regression(temp)
    write.csv(
      a,
      paste(
        "/users/PAS1475/liuzeyi/guoqi/output/deconvolution",
        all_file[i],
        j,
        sep = "/"
      )
    )
  }
  setwd("../")
  print(i)
}
