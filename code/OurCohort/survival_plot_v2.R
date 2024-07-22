rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library("survival")
library("survminer")
library(xlsx)
library(haven)

cancer_data <- read.xlsx2("../../data/OurCohort/COX data.xlsx", sheetIndex=1)

cancer_data$RFS = as.numeric(cancer_data$RFS)
cancer_data$RFS.status = as.numeric(cancer_data$RFSstatus)
cancer_data$OS = as.numeric(cancer_data$OS)
cancer_data$Osstatus = as.numeric(cancer_data$OSstatus)
cancer_data$group = as.numeric(cancer_data$CD161)

cancer_data$g[cancer_data$group == 1] = 'high'
cancer_data$g[cancer_data$group == 0] = 'low'

fit <- survfit(Surv(RFS, RFS.status) ~ g, data = cancer_data)

print(fit)

res.sum <- surv_summary(fit)
# res.sum
par(pin = c(5,3))
ggsurv <- ggsurvplot(fit,
                     data = cancer_data,
                     legend.labs = c("high", "low"),
                     pval = TRUE, 
                     conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "g", # Change risk table color by groups
                     # linetype = "CD161", # Change line type by groups
                     surv.median.line = "hv", # Specify median survival
                     ggtheme = theme_bw(), # Change ggplot2 theme
                     break.time.by = 12,
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
pdf(paste(result_dir, "OurCohort_RFS.pdf", sep=''), width=7, height=6)
print(ggsurv, newpage = FALSE)
dev.off()


################################################################################

fit <- survfit(Surv(OS, Osstatus) ~ g, data = cancer_data)

print(fit)

res.sum <- surv_summary(fit)
#res.sum

ggsurv <- ggsurvplot(fit,
                     data = cancer_data,
                     legend.labs = c("high", "low"),
                     pval = TRUE, 
                     conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "g", # Change risk table color by groups
                     # linetype = "CD161", # Change line type by groups
                     surv.median.line = "hv", # Specify median survival
                     ggtheme = theme_bw(), # Change ggplot2 theme
                     break.time.by = 12,
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

pdf(paste(result_dir, "OurCohort_OS.pdf", sep=''), width=7, height=6)
print(ggsurv, newpage = FALSE)
dev.off()

################################################################################
# library(gridExtra)
# library(cowplot)
# library(ggplot2)
# clean_data$CD161[clean_data$CD161 == 'high'] = 'High'
# clean_data$CD161[clean_data$CD161 == 'low'] = 'Low'
# clean_data$OS_STATUS[clean_data$OS_STATUS == 1] = 'Dead'
# clean_data$OS_STATUS[clean_data$OS_STATUS == 0] = 'Alive'
# 
# clean_data = clean_data[order(clean_data$KLRB1),]
# 
# p1 = ggplot(clean_data, aes(x=seq(1, nrow(clean_data)), y=KLRB1, color=CD161)) +
#   geom_point() +
#   theme_bw() + 
#   scale_color_manual(values=c('red','blue')) +
#   xlab('Patients(increasing CD161)') +
#   ylab("Expression of CD161") +
#   labs(title="TCGA") +
#   theme(plot.title=element_text(hjust=0.45), # title居中
#   )
# 
# p2 = ggplot(clean_data, aes(x=seq(1, nrow(clean_data)), y=OS_MONTHS, color=OS_STATUS)) +
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
# ggsave(paste(result_dir, "TCGA_OS_v2.pdf", sep=''))

