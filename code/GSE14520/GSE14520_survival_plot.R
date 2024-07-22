rm(list=ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# library("xlsx") # https://www.java.com/en/download/manual.jsp
library("survival")
library("survminer")
library(xlsx)
library(haven)

cancer_data <- read.xlsx2("../../data/GSE14520/clean/GSE14520_clinical_KLRB1_merge.xlsx", sheetIndex=1)

cancer_data$KLRB1 = as.numeric(cancer_data$KLRB1)
cancer_data$group[cancer_data$KLRB1 >= median(cancer_data$KLRB1)] = 'high'
cancer_data$group[cancer_data$KLRB1 < median(cancer_data$KLRB1)] = 'low'

cancer_data$RFS = as.numeric(cancer_data$Recurr.months)
cancer_data$RFSstatus = as.numeric(cancer_data$Recurr.status)

fit <- survfit(Surv(RFS, RFSstatus) ~ group, data = cancer_data)

print(fit)

res.sum <- surv_summary(fit)
#res.sum

ggsurv <- ggsurvplot(fit,
                     data = cancer_data,
                     legend.labs = c("high", "low"),
                     pval = TRUE, 
                     conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     # legend.labs = c("low", "high"),
                     risk.table.col = "group", # Change risk table color by groups
                     # linetype = "CD161", # Change line type by groups
                     surv.median.line = "hv", # Specify median survival
                     ggtheme = theme_bw(), # Change ggplot2 theme
                     break.time.by = 6,
                     ylab='Relapse Free Survival',
                     xlab='Time in months',
                     font.title = c(12, "plain"),
                     font.x = c(12, "plain"),
                     font.y = c(12, "plain"),
                     font.tickslab = c(12, "plain"),
                     palette=c("#EE0000FF", "#3B4992FF") # 颜色板 
                     # palette = c("#E7B800", "#2E9FDF")
)

ggsurv

result_dir = "../../data/survival_plot/"
if (!dir.exists(result_dir))
{
  dir.create(result_dir)
}
pdf(paste(result_dir, "GSE14520_RFS.pdf", sep=''), width=7, height=6)
print(ggsurv, newpage = FALSE)
dev.off()

################################################################################
cancer_data$OS = as.numeric(cancer_data$Survival.months)
cancer_data$OSstatus = as.numeric(cancer_data$Survival.status)

fit <- survfit(Surv(OS, OSstatus) ~ group, data = cancer_data)

print(fit)

res.sum <- surv_summary(fit)
#res.sum

ggsurv <- ggsurvplot(fit,
                     data = cancer_data,
                     legend.labs = c("high", "low"),
                     pval = TRUE, 
                     conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "group", # Change risk table color by groups
                     # linetype = "CD161", # Change line type by groups
                     surv.median.line = "hv", # Specify median survival
                     ggtheme = theme_bw(), # Change ggplot2 theme
                     break.time.by = 6,
                     ylab='Overall Survival',
                     xlab='Time in months',
                     font.title = c(12, "plain"),
                     font.x = c(12, "plain"),
                     font.y = c(12, "plain"),
                     font.tickslab = c(12, "plain"),
                     palette=c("#EE0000FF", "#3B4992FF") # 颜色板 
                     # palette = c("#E7B800", "#2E9FDF")
)

ggsurv

pdf(paste(result_dir, "GSE14520_OS.pdf", sep=''), width=7, height=6)
print(ggsurv, newpage = FALSE)
dev.off()

################################################################################
# library(gridExtra)
# library(cowplot)
# library(ggplot2)
# 
# cancer_data = read.csv('clean/GSE14520_clinical_KLRB1_merge.txt', sep='\t')
# clean_data = cancer_data
# clean_data$CD161 = clean_data$KLRB1 >= median(clean_data$KLRB1)
# clean_data$CD161[clean_data$CD161 == TRUE] = 'High'
# clean_data$CD161[clean_data$CD161 == FALSE] = 'Low'
# clean_data$Survival.status[clean_data$Survival.status == 1] = 'Dead'
# clean_data$Survival.status[clean_data$Survival.status == 0] = 'Alive'
# 
# clean_data = clean_data[order(clean_data$KLRB1),]
# 
# p1 = ggplot(clean_data, aes(x=seq(1, nrow(clean_data)), y=KLRB1, color=CD161)) +
#   geom_point() +
#   theme_bw() + 
#   scale_color_manual(values=c('red','blue')) +
#   xlab('Patients(increasing CD161)') +
#   ylab("Expression of CD161") +
#   labs(title="GSE14520") +
#   theme(plot.title=element_text(hjust=0.45), # title居中
#   )
# 
# p2 = ggplot(clean_data, aes(x=seq(1, nrow(clean_data)), y=Survival.months, color=Survival.status)) +
#   geom_point()+
#   theme_bw() + 
#   scale_color_manual(values=c('red','blue')) + 
#   xlab('Patients(increasing CD161)') +
#   ylab("Overall Survival(month)") + 
#   guides(color = guide_legend(title = "Status"))
# 
# # grid.arrange(p1, p2, nrow = 2)
# plot_grid(p1, p2, nrow = 2, align = "v")
# 
# ggsave(paste(result_dir, "GSE14520_OS_v2.pdf", sep=''))
