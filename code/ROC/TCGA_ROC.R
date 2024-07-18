rm(list=ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(xlsx)
library(pROC)
library(ggpubr)

data = read.xlsx2("../../data/ROC/公共数据库ROC.xlsx", sheetIndex=1)

data$Tumor = as.numeric(data$Group)
data$CD161 = as.numeric(data$CD161)

dfroc1<- roc(data$Tumor, data$CD161)


result_dir = "../../data/ROC/"
if (!dir.exists(result_dir))
{
  dir.create(result_dir, recursive = T)
}
pdf(file=paste(result_dir, "TCGA_ROC_plot.pdf", sep = ""), width = 5,height = 5)

plot(dfroc1,col="red",#颜色
     legacy.axes=T,#y轴格式更改
     print.auc=TRUE,#显示AUC面积
     print.thres=TRUE,#添加截点和95%CI
     print.thres.cex=0.5, #调整截断值的字体大小
     # grid=c(0.1,0.1),
     # grid.col=c("blue","green")
)#网格线设置
dev.off()#输出pdf

auc(dfroc1)
