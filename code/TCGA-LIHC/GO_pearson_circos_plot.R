rm(list=ls())

library(corrplot)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

data<-read.delim('../../data/TCGA-LIHC/clean/TCGA-LIHC_tpm_protein_coding_clean.csv',row.names=1, sep=',', check.names = F)

sample_id = colnames(data)
is_cancer = sapply(sample_id, function(x) substr(x, 14, 14) == '0')
data = data[,is_cancer]

data = log2(data + 1)###########################################################

gene_list = c("KLRB1", "PDCD1", "CD274", "HAVCR2", "LAG3", "TIGIT", "CD200R1")
gene_list2 = c("CD161", "PD1", "PDL1", "TIM3", "LAG3", "TIGIT", "CD200R1")

dat <- data[gene_list, ]
dat2 <- data.frame(t(dat))
colnames(dat2) = gene_list2

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
cm = data.matrix(R)

require(corrgram)

# png('GO/perarson_matrix_plot.png')
pdf("../../data/GO/TCGA-LIHC/perarson_matrix_plot.pdf")
corrgram(cm, order=F, lower.panel=panel.cor, upper.panel=panel.pie, text.panel=panel.txt,main="TCGA", col.regions=colorRampPalette(c("green1", "white","firebrick1")))
dev.off()

################################################################################
corr <- cor(dat2)
res <- cor.mtest(dat2, conf.level = .95)
p <- res$p

pdf("../../data/GO/TCGA-LIHC/perarson_matrix_plot2.pdf")

mycol <- colorRampPalette(c("#06a7cd", "white", "#e74a32"), alpha = TRUE)
# mycol <- colorRampPalette(c("#0AA1FF","white", "#F5CC16"), alpha = TRUE)
#上三角图添加显著性水平星号：
corrplot(corr, method = c('pie'), 
         type = c('upper'), 
         col = mycol(100),
         outline = 'grey', 
         order = c('original'), 
         diag = TRUE,
         tl.cex = 1, 
         tl.col = 'black',
         tl.pos = 'd',
         p.mat = p,
         sig.level = c(.001, .01, .05),
         insig = "label_sig", #显著性标注样式："pch", "p-value", "blank", "n", "label_sig"
         pch.cex = 1.2, #显著性标记大小
         pch.col = 'black' #显著性标记颜色
)
#下三角图添加不显著叉号：
corrplot(corr, add = TRUE,
         method = c('number'), 
         type = c('lower'),
         col = mycol(100),
         order = c('original'),
         diag = FALSE, 
         number.cex = 0.9,
         tl.pos = 'n', 
         cl.pos = 'n',
         p.mat = p,
         insig = "pch"
)

# title(main = "TCGA", cex.main = 1.5, font.main = 2)

dev.off()
