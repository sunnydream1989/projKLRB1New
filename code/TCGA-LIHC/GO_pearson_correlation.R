rm(list=ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

data<-read.delim('../../data/TCGA-LIHC/clean/TCGA-LIHC_tpm_protein_coding_clean.csv',row.names=1, sep=',', check.names = F)

sample_id = colnames(data)
is_cancer = sapply(sample_id, function(x) substr(x, 14, 14) == '0')
data = data[,is_cancer]

data = log2(data + 1)###########################################################
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
  p<-cor.test(x,y)$p.value
  pvalue<-c(pvalue,p)
  co<-cor(x,y)
  cor<-c(cor,co)
}

result_dir = "../../data/GO/TCGA-LIHC/"
if (!dir.exists(result_dir))
{
  dir.create(result_dir, recursive = T)
}

correlation<-cbind(cor,pvalue)
row.names(correlation)<-row.names(data)
write.csv(correlation,paste(result_dir, "Pearson_corr_with_KLRB1.csv", sep=''))

high_corr <- subset(correlation,correlation[, 1]>=0.5)
high_corr <- subset(high_corr,high_corr[, 2]<0.05)
write.csv(high_corr,paste(result_dir, "Pearson_corr_with_KLRB1_high.csv", sep=''))
