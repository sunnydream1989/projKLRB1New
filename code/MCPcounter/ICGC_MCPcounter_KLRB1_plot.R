rm(list = ls())

library(ggpubr)
library(dplyr)
library(tidyr)
library(reshape2)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

mcp = read.csv('../../data/MCPcounter/ICGC_MCPCounter_score.csv', row.names = 1, check.names = F)
# mcp = log2(mcp + 1)
mcp = data.frame(t(mcp), check.names = F)
mcp$patient = row.names(mcp)

data <- read.table("../../data/ICGC-LIRI-JP/clean/ICGC-LIRI-JP_tpm_protein_coding_clean_tumor_no_na.txt", header=T, sep="\t", check.names=F,row.names = 1)
data = as.data.frame(t(data['KLRB1',]))
data$patient = row.names(data)
data$KLRB1_Group[data$KLRB1 >= median(data$KLRB1)] = 'high'
data$KLRB1_Group[data$KLRB1 < median(data$KLRB1)] = 'low'
data$KLRB1 = data$KLRB1_Group

mcp = merge(mcp, data[,c('KLRB1', 'patient')], by = "patient", all = FALSE)

mcp = mcp %>% pivot_longer(-c('patient','KLRB1'), names_to = "cell", values_to = "abundances")

mcp$KLRB1 = factor(mcp$KLRB1, levels = c('low', 'high'))

p <- ggboxplot(mcp, x = "cell", y = "abundances",
               color = "KLRB1", palette = c("#F7A24F", "#C6133B")
) + # orientation = "horizontal"
  scale_y_log10() +
  stat_compare_means(aes(group = KLRB1), label = "p.signif") + 
  guides(fill=guide_legend(title=NULL)) + # 隐藏legend标题
  xlab("") + # xlab(indicator_name) +
  ylab("log2(MCPcounter + 1)") +
  # ylim(c(0,2)) + 
  labs(title="ICGC") +
  theme(plot.title=element_text(hjust=0.45), # title居中
        panel.border = element_rect(color = "black", fill = NA, size=1),  # 增加边框
        legend.title = element_blank(),
        legend.position = c(0.95,0.25),
        legend.background=element_rect(fill = alpha("white", 0)),  # legend背景透明
        legend.text = element_text(size=8),
        axis.text.x = element_text(angle=45, hjust=1, vjust=1)
  )

p

result_dir = "../../data/MCPCounter/"
if (!dir.exists(result_dir))
{
  dir.create(result_dir, recursive = T)
}
ggsave(paste(result_dir, "ICGC_MCPCounter_KLRB1_plot.pdf", sep=''),width = 30,height = 20,units = "cm")

