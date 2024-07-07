rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library("survival")
library("survminer")
library(haven)

cancer_data <- read.csv("../../data/TCGA-LIHC/clean/TCGA-LIHC_clinical_cBioPortalData_KLRB1.csv")

clean_data = cancer_data[complete.cases(cancer_data$DFS_MONTHS),]
clean_data = clean_data[clean_data$DFS_MONTHS > 0, ]
clean_data$DFS_STATUS = sub(":Recurred/Progressed",  "", as.matrix(clean_data$DFS_STATUS))
clean_data$DFS_STATUS = sub(":DiseaseFree",  "", as.matrix(clean_data$DFS_STATUS))
clean_data$DFS_STATUS = as.numeric(clean_data$DFS_STATUS)
clean_data$CD161 = clean_data$KLRB1 >= median(clean_data$KLRB1)
clean_data$CD161[clean_data$CD161 == TRUE] = 'high'
clean_data$CD161[clean_data$CD161 == FALSE] = 'low'

fit <- survfit(Surv(DFS_MONTHS, DFS_STATUS) ~ CD161, data = clean_data)

print(fit)

res.sum <- surv_summary(fit)
# res.sum
par(pin = c(5,3))
ggsurv <- ggsurvplot(fit,
                     data = clean_data,
                     legend.labs = c("high", "low"),
                     pval = TRUE, 
                     conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "CD161", # Change risk table color by groups
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
                     palette="aaas" # 颜色板 
                     # palette = c("#E7B800", "#2E9FDF")
)

ggsurv

result_dir = "../../data/survival_plot/"
if (!dir.exists(result_dir))
{
  dir.create(result_dir)
}
pdf(paste(result_dir, "TCGA_RFS.pdf", sep=''), width=7, height=6)
print(ggsurv, newpage = FALSE)
dev.off()


################################################################################
clean_data = cancer_data[complete.cases(cancer_data$OS_MONTHS),]
clean_data = clean_data[clean_data$OS_MONTHS > 0, ]
clean_data$OS_STATUS = sub(":DECEASED",  "", as.matrix(clean_data$OS_STATUS))
clean_data$OS_STATUS = sub(":LIVING",  "", as.matrix(clean_data$OS_STATUS))
clean_data$OS_STATUS = as.numeric(clean_data$OS_STATUS)
clean_data$CD161 = clean_data$KLRB1 >= median(clean_data$KLRB1)
clean_data$CD161[clean_data$CD161 == TRUE] = 'high'
clean_data$CD161[clean_data$CD161 == FALSE] = 'low'

fit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ CD161, data = clean_data)

print(fit)

res.sum <- surv_summary(fit)
#res.sum

ggsurv <- ggsurvplot(fit,
                     data = clean_data,
                     legend.labs = c("high", "low"),
                     pval = TRUE, 
                     conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "CD161", # Change risk table color by groups
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
                     palette="aaas" # 颜色板 
                     # palette = c("#E7B800", "#2E9FDF")
)

ggsurv

pdf(paste(result_dir, "TCGA_OS.pdf", sep=''), width=7, height=6)
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

