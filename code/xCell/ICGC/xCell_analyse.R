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
#ICGC

expr = read.table('../../../data/ICGC-LIRI-JP/clean/ICGC-LIRI-JP_tpm_protein_coding_clean_tumor_no_na.txt', header=TRUE, row.names=1, as.is=TRUE, sep='\t')

scores = xCellAnalysis(expr)

result_dir = "../../../data/xCell/ICGC/"
if (!dir.exists(result_dir))
{
  dir.create(result_dir, recursive = T)
}
write.csv(scores, paste(result_dir, "ICGC-LIRI-JP_xCell_result.csv", sep=''))

