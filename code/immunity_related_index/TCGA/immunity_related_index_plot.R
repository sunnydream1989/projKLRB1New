rm(list=ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggpubr)
library(dplyr)
library(tidyr)
library(reshape2)

data = read.csv('../../../data/immunity_related_index/TCGA/immunity_related_index.csv', sep=',', check.names = F)

expr<-read.delim('../../../data/TCGA-LIHC/clean/TCGA-LIHC_tpm_protein_coding_clean.csv',row.names=1, sep=',', check.names = F)

sample_id = colnames(expr)
is_cancer = sapply(sample_id, function(x) substr(x, 14, 14) == '0')
expr = expr[,is_cancer]
colnames(expr) = sapply(colnames(expr), function(x) substr(x, 1, 12))

data = data[data$`TCGA Participant Barcode` %in% colnames(expr),]
data$KLRB1 = as.numeric(t(expr['KLRB1', data$`TCGA Participant Barcode`]))

data$KLRB1_Group[data$KLRB1 >= median(data$KLRB1)] = 'High'
data$KLRB1_Group[data$KLRB1 < median(data$KLRB1)] = 'Low'

write.csv(data, '../../../data/immunity_related_index/TCGA/immunity_related_index_new.csv')
# colnames(data) = gsub('[.]', ' ', colnames(data))

index_name = c('Leukocyte Fraction','Stromal Fraction','Macrophage Regulation','Lymphocyte Infiltration Signature Score')
for (i in 1:length(index_name)){
  name = index_name[i]
  print(index_name[i])
  
  p <- ggviolin(data, 
                x="KLRB1_Group",
                y=index_name[i],
                order=c("High", "Low"),
                size=0.5, 
                palette=c("#F5B7B1", "#85C1E9"), #c("#2e9fdf", "#f8766d"), #"aaas", # 颜色板
                # color="KLRB1_Group", 
                fill="KLRB1_Group", 
                add = "boxplot", 
                width = 0.8, 
                # add.params=list(color="KLRB1_Group", size=0.4)
  )+
    stat_compare_means(label.x = 1, label.y = ceiling(max(data[index_name[i]], na.rm=TRUE))) + # 添加p值
    guides(fill=guide_legend(title=NULL)) + # 隐藏legend标题
    xlab("") + # xlab(indicator_name) +
    ylab(index_name[i]) +
    # labs(title=index_name[i]) +
    theme(plot.title=element_text(hjust=0.45), # title居中
          panel.border = element_rect(color = "black", fill = NA, size=1),  # 增加边框
          legend.title = element_blank(),
          legend.position = c(0.9,0.9),
          legend.background=element_rect(fill = alpha("white", 0)),  # legend背景透明
          legend.text = element_text(size=8)
    )
  
  p
  ggsave(paste('../../../data/immunity_related_index/TCGA/', index_name[i], "_plot.pdf", sep=''),width = 10,height = 10,units = "cm")
  
}

################################################################################
# datasub = data[, c('IFN-gamma Response', 'TGF-beta Response', 'KLRB1_Group')]
# datasub = datasub %>% pivot_longer(-c('KLRB1_Group'), names_to = "cell", values_to = "abundances")


################################################################################
datasub = data[, c('BCR Shannon', 'TCR Shannon', 'KLRB1_Group')]
datasub = datasub %>% pivot_longer(-c('KLRB1_Group'), names_to = "cell", values_to = "abundances")


p <- ggviolin(datasub, 
              x="cell",
              y="abundances",
              # color="KLRB1_Group",
              # order=c("High", "Low"),
              size=0.5, 
              palette=c("#F5B7B1", "#85C1E9"), #c("#2e9fdf", "#f8766d"), #"aaas", # 颜色板
              # color="KLRB1_Group", 
              fill="KLRB1_Group", 
              add = "boxplot", 
              width = 0.8, 
              # add.params=list(color="KLRB1_Group", size=0.4)
)+
  stat_compare_means(aes(group = KLRB1_Group), label.x = 1, label.y = ceiling(max(datasub["abundances"], na.rm=TRUE))) + # 添加p值
  guides(fill=guide_legend(title=NULL)) + # 隐藏legend标题
  xlab("") + # xlab(indicator_name) +
  ylab("Z score") +
  # labs(title=index_name[i]) +
  theme(plot.title=element_text(hjust=0.45), # title居中
        panel.border = element_rect(color = "black", fill = NA, size=1),  # 增加边框
        legend.title = element_blank(),
        legend.position = c(0.9,0.7),
        legend.background=element_rect(fill = alpha("white", 0)),  # legend背景透明
        legend.text = element_text(size=8)
  )

p
ggsave(paste('../../../data/immunity_related_index/TCGA/', "Z score", "_plot.pdf", sep=''),width = 15,height = 10,units = "cm")


################################################################################
data = read.csv('../../../data/immunity_related_index/TCGA/immunity_related_index_new.csv', sep=',', check.names = F)

data = data[complete.cases(data["Immune Subtype"]),]

C1 = "Wound Healing"
C2 = "IFN-γ Dominant"
C3 = "Inflammatory"
C4 = "Lymphocyte Depleted"
C6 = "TGF-β Dominant"

data["Immune Subtype"][data["Immune Subtype"] == "C1"] = C1
data["Immune Subtype"][data["Immune Subtype"] == "C2"] = C2
data["Immune Subtype"][data["Immune Subtype"] == "C3"] = C3
data["Immune Subtype"][data["Immune Subtype"] == "C4"] = C4
data["Immune Subtype"][data["Immune Subtype"] == "C6"] = C6

# my_comparisons = list(c("C1", "C2"), c("C2", "C3"), c("C3", "C4"), c("C1", "C3"), c("C2", "C4"), c("C1", "C4"))
my_comparisons = list(c(C1, C2), c(C2, C3), c(C3, C4), c(C1, C3), c(C2, C4), c(C1, C4))
p <- ggviolin(data,
              x="Immune Subtype",
              y="KLRB1",
              order=c(C1, C2, C3, C4, C6),
              size=0.5, 
              palette="Set2", # 颜色板
              # color="KLRB1_Group", 
              fill="Immune Subtype", 
              add = "boxplot", 
              width = 1.2, 
              # add.params=list(color="KLRB1_Group", size=0.4)
)+
  stat_compare_means(label="p.signif", comparisons=my_comparisons) +
  stat_compare_means(label.x = 1, label.y = 80) + # 添加p值
  guides(fill=guide_legend(title=NULL)) + # 隐藏legend标题
  xlab("TCGA Immune Subtype") + # xlab(indicator_name) +
  ylab("Expression of CD161") +
  # labs(title=index_name[i]) +
  theme(plot.title=element_text(hjust=0.45), # title居中
        panel.border = element_rect(color = "black", fill = NA, size=1),  # 增加边框
        legend.title = element_blank(),
        legend.position = c(0.85,0.8),
        legend.background=element_rect(fill = alpha("white", 0)),  # legend背景透明
        legend.text = element_text(size=8),
        axis.text.x = element_text(angle = 25, hjust = 0.5, vjust = 0.6)
  )

p
ggsave(paste('../../../data/immunity_related_index/TCGA/', 'Immune Subtype',"_plot.pdf", sep=''),width = 15,height = 15,units = "cm")








