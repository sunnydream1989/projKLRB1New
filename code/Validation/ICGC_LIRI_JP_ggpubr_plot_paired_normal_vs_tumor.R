
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggpubr)
library(survival)
library(survminer)
library(haven)
library(tidyr)
library(dplyr)


# 读取基因表达的矩阵
expr_file <- "../../data/ICGC-LIRI-JP/clean/ICGC-LIRI-JP_rna_tpm_clean_50_include_normal.csv"
expr=read.csv(expr_file, header=TRUE, stringsAsFactors=FALSE, check.names=F)
gene = expr['gene_id']

clinical = read.csv('../../data/ICGC-LIRI-JP/clean/ICGC-LIRI-JP_clinical_clean_include_normal.csv', row.names = 1)
clinical$sampleid = row.names(clinical)
clinical = clinical[colnames(expr)[-1],]

clinical_normal = clinical[grepl('Liver', clinical$submitted_sample_id),]
clinical_tumor = clinical[grepl('Cancer', clinical$submitted_sample_id),]
common_donor = intersect(clinical_normal$icgc_donor_id, clinical_tumor$icgc_donor_id)

row.names(clinical_normal) = clinical_normal$icgc_donor_id
row.names(clinical_tumor) = clinical_tumor$icgc_donor_id

clinical_normal = clinical_normal[common_donor, ]
clinical_tumor = clinical_tumor[common_donor, ]

clinical = rbind(clinical_normal, clinical_tumor)

row.names(expr) = expr$gene_id
gene_set = c('KLRB1')
expr = expr[gene_set, -1]
dim(expr)

expr <- expr[, clinical$sampleid]
expr[is.na(expr)] = 0

expr = log2(expr + 1)

data = data.frame(t(expr))
data$icgc_donor_id = clinical$icgc_donor_id
data$submitted_sample_id = clinical$submitted_sample_id
data$sampleid = row.names(data)
row.names(data) <- NULL

################################################################################
indicator_name = "Tissue type"
data[indicator_name] = ifelse(grepl('Cancer', clinical$submitted_sample_id), "Tumor", "Normal")

result_dir = "../../data/Validation/"
if (!dir.exists(result_dir))
{
  dir.create(result_dir)
}

write.csv(data, paste(result_dir, 'ICGC_LIRI_JP_ggpubr_plot_paired_normal_vs_tumor.csv', sep=''), row.names = T)

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
    labs(title="ICGC") +
    theme(plot.title=element_text(hjust=0.45), # title居中
          panel.border = element_rect(color = "black", fill = NA, size=1),  # 增加边框
          legend.title = element_blank(),
          legend.position = c(0.85,0.85),
          legend.background=element_rect(fill = alpha("white", 0)),  # legend背景透明
          legend.text = element_text(size=8)
    )
  p
  ggsave(paste(result_dir, 'ICGC_LIRI_JP_', gene, ' ',indicator_name, "_paired.pdf", sep=''),width = 10,height = 10,units = "cm")

}

data2 <- pivot_longer(data = data,cols = gene_set, names_to = "gene",values_to = 'expr')

# Box plot facetted by "dose"
p <- ggpaired(data2, x = indicator_name, y = "expr",
              order=c("Tumor", "Normal"),
              color = indicator_name, palette = "npg", 
              line.color = "gray", line.size = 0.4,
              facet.by = "gene", short.panel.labs = FALSE) + 
  # facet_wrap(~gene, nrow = 1) +
  facet_rep_wrap(. ~gene,scales = 'free',repeat.tick.labels = 'left', nrow = 1, ncol = 4)+
  ylab("Gene Expression") +
  xlab("")+
  theme(legend.position = "right")
p = p + stat_compare_means(label = "p.format", paired = TRUE)

p
ggsave(paste(result_dir, 'ICGC_LIRI_JP_', indicator_name, "_paired.pdf", sep=''),width = 18,height = 10,units = "cm")

