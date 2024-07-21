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

data <- data[complete.cases(data), ]

data$gene_name <- rownames(data)
data <- data[, c(ncol(data), 1:(ncol(data)-1))]

write.csv(data,'../../data/ICGC-LIRI-JP/clean/ICGC-LIRI-JP_tpm_protein_coding_clean_tumor_no_na.csv', row.names = F)

write.table(data,'../../data/ICGC-LIRI-JP/clean/ICGC-LIRI-JP_tpm_protein_coding_clean_tumor_no_na.txt', sep='\t', row.names = F)
