# install.packages('devtools')
# devtools::install_github('dviraran/xCell', force = TRUE)
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("GSVA","GSEABase"))
# 
# install.packages('Rcpp')

rm(list=ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(xCell)

#########################################################################################################
#TCGA

expr = read.table('../../../data/TCGA-LIHC/clean/TCGA-LIHC_tpm_protein_coding_clean_tumor.txt', header=TRUE, row.names=1, as.is=TRUE, sep='\t', check.names = F)

scores = xCellAnalysis(expr)

result_dir = "../../../data/xCell/TCGA-LIHC/"
if (!dir.exists(result_dir))
{
  dir.create(result_dir, recursive = T)
}
write.csv(scores, paste(result_dir, "TCGA-LIHC_xCell_result.csv", sep=''))
