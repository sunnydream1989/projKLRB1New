rm(list=ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

data<-read.delim('../../data/TCGA-LIHC/clean/TCGA-LIHC_tpm_protein_coding_clean.csv',row.names=1, sep=',', check.names = F)

sample_id = colnames(data)
is_cancer = sapply(sample_id, function(x) substr(x, 14, 14) == '0')
data = data[,is_cancer]

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
    r = c(r, cor(a, b))
  }
  R = cbind(R, r)
}

length2<-length(R[1, ])
colnames(R)[1] = 'gene'
colnames(R)[2:length2] <- gene_list2
R[,1] = gene_list2

write.table(R, '../../data/GO/TCGA-LIHC/Pearson_matrix_v2.txt', sep = '\t', quot = F, row.names = F)

R = read.delim('../../data/GO/TCGA-LIHC/Pearson_matrix_v2.txt', row.names = 1)
dat = data.matrix(R)

require(corrgram)

# png('GO/perarson_matrix_plot.png')
pdf("../../data/GO/TCGA-LIHC/perarson_matrix_plot.pdf")
corrgram(dat, order=F, lower.panel=panel.cor, upper.panel=panel.pie, text.panel=panel.txt,main="TCGA", col.regions=colorRampPalette(c("green1", "white","firebrick1")))
dev.off()
