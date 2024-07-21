rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(MCPcounter)

genes <- data.table::fread("../../data/MCPcounter/genes.txt",data.table = F)
probesets <- data.table::fread("../../data/MCPcounter/probesets.txt",data.table = F,header = F)

test <- read.table("../../data/TCGA-LIHC/clean/TCGA-LIHC_tpm_protein_coding_clean_tumor.txt", header=T, sep="\t", check.names=F,row.names = 1)

results<- MCPcounter.estimate(test,
                              featuresType= "HUGO_symbols", 
                              probesets=probesets,
                              genes=genes)

write.csv(results,file="../../data/MCPcounter/TCGA_LIHC_MCPCounter_score.csv",row.names =T,quote=F)
