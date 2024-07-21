rm(list=ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggpubr)
library(tidyr)

xCell_data = as.data.frame(t(read.csv('../../../data/xCell/TCGA-LIHC/TCGA-LIHC_xCell_result.csv', check.names = F, row.names = 1)))
c = colnames(xCell_data)
xCell_data$sample_id = row.names(xCell_data)

xCell_data_L <- pivot_longer(
  xCell_data,
  cols = c,  # 指定要转换的列
  names_to = "cell",  # 新列的名称，用于存储原来的列名
  values_to = "value"  # 新列的名称，用于存储原来的值
)
xCell_data_L = as.data.frame(xCell_data_L)

expr = read.table('../../../data/TCGA-LIHC/clean/TCGA-LIHC_tpm_protein_coding_clean_tumor.txt', header=TRUE, row.names=1, as.is=TRUE, sep='\t', check.names = F)

gene = as.data.frame(t(expr['KLRB1', xCell_data$sample_id]))
median_value <- median(gene$KLRB1)

# 根据中位数将表达值分类为 "high" 或 "low"
gene$KLRB1 <- factor(ifelse(gene$KLRB1 > median_value, "high", "low"))


xCell_data_L$KLRB1 = gene[xCell_data_L$sample_id, 'KLRB1']


####
result_dir = "../../../data/xCell/TCGA-LIHC/"
if (!dir.exists(result_dir))
{
  dir.create(result_dir, recursive = T)
}

for (i in c('ImmuneScore', 'StromaScore', 'MicroenvironmentScore')){
  ImmuneScore = xCell_data_L[xCell_data_L$cell == i, ]
  
  p <- ggboxplot(ImmuneScore, 
                 x = "cell", 
                 y = "value",
                 # order=c("high", "low"),
                 # color = "KLRB1", 
                 fill="KLRB1", 
                 palette = c("orange", "cyan"),
                 add = "boxplot",
                 # width = 0.8,
                 # size=0.5
  ) + 
    stat_compare_means(aes(group = KLRB1), label.y=max(ImmuneScore$value) * 1.1) + 
    guides(fill=guide_legend(title=NULL)) + # 隐藏legend标题
    xlab("KLRB1") + # xlab(indicator_name) +
    ylab(i) +
    labs(title="") +
    theme(plot.title=element_text(hjust=0.45), # title居中
          panel.border = element_rect(color = "black", fill = NA, size=1),  # 增加边框
          legend.title = element_blank(),
          legend.position = c(0.8,0.85),
          legend.background=element_rect(fill = alpha("white", 0)),  # legend背景透明
          legend.text = element_text(size=8),
          axis.text.x = element_blank()
    )
  
  p
  ggsave(file=paste0(result_dir, i, "_KLRB1_boxplot.pdf"),width = 10,height = 10,units = "cm")
  
}

####
xCell_data_L = xCell_data_L[(xCell_data_L$cell %in% rev(c('DC', 'aDC', 'cDC', 'pDC', 'CD8+ naive T-cells', 'CD8+ T-cells', 'CD8+ Tcm', 'CD8+ Tem', 'CD4+ memory T-cells', 'CD4+ naive T-cells', 'CD4+ T-cells', 'CD4+ Tcm', 'CD4+ Tem',  'Eosinophils', 'Fibroblasts', 'Macrophages', 'Macrophages M1', 'Macrophages M2', 'Mast cells', 'Monocytes', 'Neutrophils', 'NK cells', 'NKT', 'Tregs','B-cells', 'Basophils'))), ]
# xCell_data_L = xCell_data_L[(xCell_data_L$cell %in% c("Adipocytes", "Astrocytes", "B-cells", "CD4+ memory T-cells", "CD4+ T-cells", "CD8+ T-cells", "cDC", "DC", "Endothelial cells", "Eosinophils", "Epithelial cells", "Fibroblasts", "GMP", "Hepatocytes", "HSC", "Macrophages", "Macrophages M1", "Macrophages M2", "Mast cells", "Monocytes", "mv Endothelial cells", "Myocytes", "Neutrophils", "NKT", "Smooth muscle", "Th1 cells", "Th2 cells", "Tregs")), ]

xCell_data_L$cell <- factor(xCell_data_L$cell, levels = rev(c('DC', 'aDC', 'cDC', 'pDC', 'CD8+ naive T-cells', 'CD8+ T-cells', 'CD8+ Tcm', 'CD8+ Tem', 'CD4+ memory T-cells', 'CD4+ naive T-cells', 'CD4+ T-cells', 'CD4+ Tcm', 'CD4+ Tem',  'Eosinophils', 'Fibroblasts', 'Macrophages', 'Macrophages M1', 'Macrophages M2', 'Mast cells', 'Monocytes', 'Neutrophils', 'NK cells', 'NKT', 'Tregs','B-cells', 'Basophils')))

p <- ggboxplot(xCell_data_L, x = "cell", y = "value",
               color = "KLRB1", palette = c("#EE0000FF", "#3B4992FF"), 
               orientation = "horizontal") + 
  stat_compare_means(aes(group = KLRB1), label = "p.signif") + 
  guides(fill=guide_legend(title=NULL)) + # 隐藏legend标题
  xlab("") + # xlab(indicator_name) +
  ylab("xCell score") +
  labs(title="") +
  theme(plot.title=element_text(hjust=0.45), # title居中
        panel.border = element_rect(color = "black", fill = NA, size=1),  # 增加边框
        legend.title = element_blank(),
        legend.position = c(0.8,0.05),
        legend.background=element_rect(fill = alpha("white", 0)),  # legend背景透明
        legend.text = element_text(size=8)
  )

p

ggsave(paste(result_dir, "TCGA_xCell_KLRB1_boxplot.pdf", sep=''),width = 20,height = 30,units = "cm")
