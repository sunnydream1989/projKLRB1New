rm(list=ls())

library(dplyr)
library(stringr)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

clinical_data<-read.delim('../../data/ICGC-LIRI-JP/clean/ICGC-LIRI-JP_clinical_clean_include_normal.csv', sep=',', check.names = F)
clinical_data<-clinical_data %>% 
  filter(str_ends(submitted_sample_id, "_Cancer"))

data<-read.delim('../../data/ICGC-LIRI-JP/clean/ICGC-LIRI-JP_rna_tpm_clean_50_include_normal.csv',row.names=1, sep=',', check.names = F)
data<-data[,clinical_data$icgc_sample_id]
non_na_columns <- !is.na(data["KLRB1", ]) # 有一个KLRB1是NA，需要剔除
data <- data[, non_na_columns]

clinical_data = clinical_data[clinical_data$icgc_sample_id %in% colnames(data), ]
data<-data[,clinical_data$icgc_sample_id]

sample_id = colnames(data)
data = log2(data + 1)###########################################################

gene_list = c("KLRB1", "PDCD1", "CD274", "HAVCR2", "LAG3", "TIGIT", "CD200R1")
gene_list2 = c("KLRB1", "PD1", "PDL1", "TIM3", "LAG3", "TIGIT", "CD200R1")

dat <- data[gene_list, ]

k = c()
R = rownames(dat)
for (k in 1:length(rownames(dat))){
  r = c()
  for (i in 1: length(rownames(dat))){
    a = as.numeric(dat[k, ])
    b = as.numeric(dat[i, ])
    r = c(r, cor(a, b, use="complete.obs"))
  }
  R = cbind(R, r)
}

length2<-length(R[1, ])
colnames(R)[1] = 'gene'
colnames(R)[2:length2] <- gene_list2
R[,1] = gene_list2

write.table(R, '../../data/GO/ICGC/Pearson_matrix_v2.txt', sep = '\t', quot = F, row.names = F)

R = read.delim('../../data/GO/ICGC//Pearson_matrix_v2.txt', row.names = 1)
dat = data.matrix(R)

require(corrgram)

# png('GO/perarson_matrix_plot.png')
pdf("../../data/GO/ICGC//perarson_matrix_plot.pdf")
corrgram(dat, order=F, lower.panel=panel.cor, upper.panel=panel.pie, text.panel=panel.txt,main="ICGC", col.regions=colorRampPalette(c("firebrick1", "white","blue")))
dev.off()

