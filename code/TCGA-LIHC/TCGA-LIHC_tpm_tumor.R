rm(list=ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

data<-read.delim('../../data/TCGA-LIHC/clean/TCGA-LIHC_tpm_protein_coding_clean.csv',row.names=1, sep=',', check.names = F)

sample_id = colnames(data)
is_cancer = sapply(sample_id, function(x) substr(x, 14, 14) == '0')
data = data[,is_cancer]

data$gene_name <- rownames(data)
data <- data[, c(ncol(data), 1:(ncol(data)-1))]

write.csv(data,'../../data/TCGA-LIHC/clean/TCGA-LIHC_tpm_protein_coding_clean_tumor.csv', row.names = F)

write.table(data,'../../data/TCGA-LIHC/clean/TCGA-LIHC_tpm_protein_coding_clean_tumor.txt', sep='\t', row.names = F)
