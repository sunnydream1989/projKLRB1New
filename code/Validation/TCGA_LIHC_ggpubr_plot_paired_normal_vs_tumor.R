
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggpubr)
library(survival)
library(survminer)
library(haven)
library(tidyr)
library(dplyr)
library(lemon)
library(patchwork)

# 读取基因表达的矩阵
expr_file <- "../../data/TCGA-LIHC/clean/TCGA-LIHC_tpm_protein_coding_clean.csv"
expr=read.csv(expr_file, header=TRUE, stringsAsFactors=FALSE, check.names=F)
row.names(expr) = expr$gene_name
gene_set = c('KLRB1')
expr = expr[gene_set, -1]

dim(expr)
expr = log2(expr + 1)

sample_id = colnames(expr)
is_normal = sapply(sample_id, function(x) substr(x, 14, 14) != '0')
sample_id_normal = sample_id[is_normal]
length(sample_id_normal)

patiendid = sapply(sample_id_normal, function(x) substr(x, 1, 12))

is_cancer = sapply(sample_id, function(x) substr(x, 14, 14) == '0')
sample_id_cancer = sample_id[is_cancer]
length(sample_id_cancer)

sample_id_cancer_df = data.frame(sample_id_cancer, row.names = sapply(sample_id_cancer, function(x) substr(x, 1, 12)))
sample_id_cancer = sample_id_cancer_df[sapply(sample_id_normal, function(x) substr(x, 1, 12)),]

sample_id_combined <- c(sample_id_cancer, sample_id_normal)
expr <- expr[, sample_id_combined]

data = data.frame(t(expr))
data$sampleid = row.names(data)
data$patientid = sapply(row.names(data), function(x) substr(x, 1, 12))
row.names(data) <- NULL
# expr = cbind(gene, expr)

################################################################################
indicator_name = "Tissue type"
data[indicator_name] = ifelse(sapply(data$sampleid, function(x) substr(x, 14, 14) == '0'), "Tumor", "Normal")
# data_clean = data[complete.cases(data[indicator_name]),]

result_dir = "../../data/Validation/"
if (!dir.exists(result_dir))
{
  dir.create(result_dir)
}

write.csv(data, paste(result_dir, 'TCGA_LIHC_ggpubr_plot_paired_normal_vs_tumor.csv', sep=''), row.names = T)

for (gene in gene_set){
  p <- ggpaired(data, 
                 x=indicator_name,
                 y=gene,
                 order=c("Normal", "Tumor"),
                 size=0.5,
                 width = 0.8,
                 palette="aaas", # 颜色板
                 color=indicator_name, 
                 # add="jitter",
                 line.color="gray",
                 line.size=0.2, 
                 add.params=list(color=indicator_name, size=0.4)
                 ) + 
    stat_compare_means(label.x = 1, label.y = max(data[gene]), paired = TRUE) + # 添加p值
    guides(fill=guide_legend(title=NULL)) + # 隐藏legend标题
    xlab("") + # xlab(indicator_name) +
    ylab(paste("Expression of ", gene)) +
    labs(title="TCGA") +
    theme(plot.title=element_text(hjust=0.45), # title居中
          panel.border = element_rect(color = "black", fill = NA, size=1),  # 增加边框
          legend.title = element_blank(),
          legend.position = c(0.85,0.85),
          legend.background=element_rect(fill = alpha("white", 0)),  # legend背景透明
          legend.text = element_text(size=8)
    )
  p
  ggsave(paste(result_dir, 'TCGA_LIHC_', gene, ' ',indicator_name, "_paired.pdf", sep=''),width = 10,height = 10,units = "cm")

}

data2 <- pivot_longer(data = data, cols = gene_set, names_to = "gene",values_to = 'expr')

# Box plot facetted by "dose"
p <- ggpaired(data2, x = indicator_name, y = "expr",
              order=c("Tumor", "Normal"),
              color = indicator_name, palette = "npg", 
              line.color = "gray", line.size = 0.4,
              facet.by = "gene", short.panel.labs = FALSE) + 
  #facet_wrap(~gene, nrow = 1) +
  facet_rep_wrap(. ~gene,scales = 'free',repeat.tick.labels = 'left', nrow = 1, ncol = 4)+
  ylab("Gene Expression") +
  xlab("")+
  theme(legend.position = "right")
p = p + stat_compare_means(label = "p.format", paired = TRUE)

p
ggsave(paste(result_dir, 'TCGA_LIHC_', indicator_name, "_paired.pdf", sep=''),width = 18,height = 10,units = "cm")