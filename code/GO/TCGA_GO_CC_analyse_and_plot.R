rm(list = ls())
library(limma)
library(ggplot2)
library(egg)##控制绘图区大小，以保证标签文字等所占空间不同时，同批图像还是同样大小
library(grid)

filepath = file.choose()

dat = read.delim(filepath)

dat$P = -log10(dat$Benjamini)
dat <- dat[order(dat$Benjamini), ]
n_line <- 5
dat <- dat[1:n_line,]
dat$ID <- seq(n_line,1)

p1 <- ggplot(dat,aes(x=P, y=ID, size=Count)) + 
  theme_bw() + geom_point(colour="#ffa500", shape = 16) +
  ggtitle("CC(TCGA)") + 
  theme(plot.title = element_text(hjust = 0.5)) # legend.position = c(0.8,0.2), 

label = strsplit2(dat$Term, split='~')[, 2]
p2 = p1 + scale_y_continuous(breaks=dat$ID, labels=label, name = '') + 
  theme(axis.text.y = element_text(color = "black"))
p3 = p2 + xlab("-log10(P)") + scale_x_continuous(limits=c(floor(min(dat$P)) - 1, ceiling(max(dat$P)) + 1)) + 
  theme(axis.text.x = element_text(color = "black"))
p3
ggsave(paste(filepath, '.pdf', sep=''), 
       egg::set_panel_size(p3, 
                           width=unit(1, "in"), 
                           height=unit(2, "in")), 
       width = 5, 
       height = 3, 
       units = 'in', 
       dpi = 300)