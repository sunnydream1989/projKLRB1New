rm(list=ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(dplyr)
library(stringr)
library(GSVA)

clinical_data<-read.delim('../../../data/ICGC-LIRI-JP/clean/ICGC-LIRI-JP_clinical_clean_include_normal.csv', sep=',', check.names = F)
clinical_data<-clinical_data %>% 
  filter(str_ends(submitted_sample_id, "_Cancer"))

data=read.delim('../../../data/ICGC-LIRI-JP/clean/ICGC-LIRI-JP_rna_tpm_clean_50_include_normal.csv',sep=',', row.names=1, header=T)###表达谱数据第一行为ID，第一列为gene symbol
data<-data[,clinical_data$icgc_sample_id]
non_na_columns <- !is.na(data["KLRB1", ])
data <- data[, non_na_columns]

clinical_data = clinical_data[clinical_data$icgc_sample_id %in% colnames(data), ]
data<-data[,clinical_data$icgc_sample_id]

file_name = file.choose()

genelist = read.delim(file_name, row.names=1, header=FALSE, col.names = paste("V", 1:5000, sep = ""))###导入genelist文件，第一列为genelist名，第二列为链接（需要删除），右侧为该功能相关的所有gene名
# 如果不指定最大列数的话，会莫名其妙换行，导致数据读取出错，所以制定一个较大的列数5000

genelist = genelist[, -1]
genelist = as.matrix(genelist) 

# genelist = genelist[grep(pattern="GOBP_ADAPTIVE_IMMUNE_RESPONSE", rownames(genelist)), ]

target = c('GOBP_ADAPTIVE_IMMUNE_RESPONSE',
           'GOBP_IMMUNE_RESPONSE_REGULATING_CELL_SURFACE_RECEPTOR_SIGNALING_PATHWAY',
           'GOBP_NATURAL_KILLER_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL',
           'GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
           'GOBP_T_CELL_ACTIVATION',
           'GOBP_REGULATION_OF_T_CELL_ACTIVATION',
           'GOBP_IMMUNE_RESPONSE',
           'GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_IMMUNE_RESPONSE',
           'GOBP_NATURAL_KILLER_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
           'GOBP_B_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
           'GOMF_C_C_CHEMOKINE_BINDING',
           'GOMF_CYTOKINE_RECEPTOR_ACTIVITY',
           'GOMF_C_X_C_CHEMOKINE_RECEPTOR_ACTIVITY',
           'GOCC_T_CELL_RECEPTOR_COMPLEX',
           'GOCC_EXTERNAL_SIDE_OF_PLASMA_MEMBRANE',
           'KEGG_PRIMARY_IMMUNODEFICIENCY',
           'KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION',
           'KEGG_CHEMOKINE_SIGNALING_PATHWAY',
           'GOBP_IMMUNE_RESPONSE_TO_TUMOR_CELL',
           'GOBP_T_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL',
           'GOBP_LYMPHOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE')

genelist = genelist[intersect(row.names(genelist), target), ]

nrow = nrow(genelist)
data = data.matrix(data)
genelist[genelist == ""] = NA  ##genelist文件空值定义为NA
genesets = c()

for (i in 1:nrow)
{
  a = genelist[i,]
  a = a[!is.na(a)]  ###去除所取行中na
  a = list(as.character(as.matrix(a)))  ###将行转为list格式
  stopifnot(length(a[[1]]) < 4500)
  genesets = c(genesets, a)
}

overlap_num = c()
for(i in 1:nrow(genelist)){
  b = as.character(unlist(genelist[i,]))
  o = intersect(rownames(data), b)
  overlap_num = c(overlap_num, length(o))
}

index = overlap_num > 0
d = rownames(genelist)[index]

result = gsva(data, genesets, mx.diff=FALSE, verbose=TRUE)
result = data.matrix(result)
rownames(result) = tolower(d)
colnames(result) = colnames(data)


result_dir = "../../../data/GSVA/ICGC/"
if (!dir.exists(result_dir))
{
  dir.create(result_dir, recursive = T)
}
write.csv(result, file=paste(result_dir, basename(file_name), ".csv", sep=''), quote=F)


####################################################################################################################
# KLRB1和GSVA的相关性分析

x<-data['KLRB1', ]
x<-as.numeric(x)

x = log2(x+1)###################################################################

pvalue<-c()
cor<-c()

length2<-length(result[,1])
for (i in 1:length2)
{
  y<-result[i, ]
  y<-as.numeric(y)
  p<-cor.test(x,y)$p.value
  pvalue<-c(pvalue,p)
  co<-cor(x,y)
  cor<-c(cor,co)
}

correlation<-cbind(cor,pvalue)
row.names(correlation)<-row.names(result)
correlation <- correlation[order(-correlation[ ,1]), ]

write.csv(correlation, file=paste(result_dir, basename(file_name), "_corr.csv", sep=''), quote=F)

