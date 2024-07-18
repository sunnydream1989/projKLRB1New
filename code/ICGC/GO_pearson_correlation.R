rm(list=ls())
library(dplyr)
library(stringr)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

clinical_data<-read.delim('../../data/ICGC-LIRI-JP/clean/ICGC-LIRI-JP_clinical_clean_include_normal.csv', sep=',', check.names = F)
clinical_data<-clinical_data %>% 
  filter(str_ends(submitted_sample_id, "_Cancer"))

data<-read.delim('../../data/ICGC-LIRI-JP/clean/ICGC-LIRI-JP_rna_tpm_clean_50_include_normal.csv',row.names=1, sep=',', check.names = F)
data<-data[,clinical_data$icgc_sample_id]

sample_id = colnames(data)
# is_cancer = sapply(sample_id, function(x) substr(x, 14, 14) == '0')
# data = data[,is_cancer]

data<-as.matrix(data)
length1<-length(data[1,])
length2<-length(data[,1])
x<-data['KLRB1',1:length1]
x<-as.numeric(x)
pvalue<-c()
cor<-c()
rn = rownames(data)
for (i in 1:length2)
{
  if (rn[i] == 'KLRB1'){
    pvalue<-c(pvalue,NaN)
    cor<-c(cor,NaN)
    next
  }
  y<-data[i,1:length1]
  y<-as.numeric(y)
  if (sum(!is.na(y)) < 200){
    pvalue<-c(pvalue,NaN)
    cor<-c(cor,NaN)
    next
  }
  p<-cor.test(x,y)$p.value
  pvalue<-c(pvalue,p)
  co<-cor(x,y, use="complete.obs")
  cor<-c(cor,co)
}

result_dir = "../../data/GO/ICGC/"
if (!dir.exists(result_dir))
{
  dir.create(result_dir, recursive = T)
}

correlation<-cbind(cor,pvalue)
row.names(correlation)<-row.names(data)
write.csv(correlation,paste(result_dir, "Pearson_corr_with_KLRB1.csv", sep=''))

high_corr <- subset(correlation,correlation[, 1]>=0.5)
write.csv(high_corr,paste(result_dir, "Pearson_corr_with_KLRB1_high.csv", sep=''))
