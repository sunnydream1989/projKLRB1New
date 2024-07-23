rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(dplyr)
library(stringr)
library(MCPcounter)

genes <- data.table::fread("../../data/MCPcounter/genes.txt",data.table = F)
probesets <- data.table::fread("../../data/MCPcounter/probesets.txt",data.table = F,header = F)

# method1
# test <- read.table("../../data/ICGC-LIRI-JP/clean/ICGC-LIRI-JP_tpm_protein_coding_clean_tumor_no_na.txt", header=T, sep="\t", check.names=F,row.names = 1)

# method2
clinical_data<-read.delim('../../data/ICGC-LIRI-JP/clean/ICGC-LIRI-JP_clinical_clean_include_normal.csv', sep=',', check.names = F)
clinical_data<-clinical_data %>% 
  filter(str_ends(submitted_sample_id, "_Cancer"))

data<-read.delim('../../data/ICGC-LIRI-JP/clean/ICGC-LIRI-JP_rna_tpm_clean_50_include_normal.csv',row.names=1, sep=',', check.names = F)
data<-data[,clinical_data$icgc_sample_id]
non_na_columns <- !is.na(data["KLRB1", ]) # 有一个KLRB1是NA，需要剔除
test <- data[, non_na_columns]

results<- MCPcounter.estimate(test,
                              featuresType= "HUGO_symbols", 
                              probesets=probesets,
                              genes=genes)

write.csv(results,file="../../data/MCPcounter/ICGC_MCPCounter_score.csv",row.names =T,quote=F)
