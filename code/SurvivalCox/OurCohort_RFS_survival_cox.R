rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(survival)
library(xlsx)

# df = read.xlsx("../../data/TCGA-LIHC/clean/TCGA247COX8.xlsx", sheetIndex = 1)
df = read.xlsx("../../data/OurCohort/COX data.xlsx", sheetIndex = 1)
df = df[,c("RFSstatus","RFS","Age","Gender","Hepatitis","Tumor.Size","Tumor.Number","PVTT","MVI","TNM.Stage.AJCC.","CNLC.Stage","BCLC","NLR","PLR","ALBI","AFP","CD161")]
# df <- data.frame(apply(df, 2, as.numeric))

pFilter=0.1 #设一个p值标准，后面用
outResult=data.frame() #建一个空白数据框，后面for循环输出用
sigGenes=c("RFSstatus","RFS") #建一个向量，后面for循环输出用，因为后面还要用到surstat及surtime，所以先放在向量里
for(i in colnames(df[,3:ncol(df)])){ #从第3列开始循环，因为1列2列不是gene，是surstat和surtime
  tdcox <- coxph(Surv(RFS, RFSstatus) ~ df[,i], data = df)#开始逐一循环cox分析
  tdcoxSummary = summary(tdcox) #summary命令对tdcox总结，方面后面提取数据
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"] #提取p值，这个主要是后面提取有意义的gene用
  
  sigGenes=c(sigGenes,i)
  outResult=rbind(outResult,#合并行，实际上是对循环结果的合并，前面设置的空白数据框outResult这里用，循环必须有个开始
                  cbind(id=i,#合并列，是每个基因的统计数据
                        HR=round(tdcoxSummary$conf.int[,"exp(coef)"],3),#提取单个基因的HR
                        L95CI=round(tdcoxSummary$conf.int[,"lower .95"],3),#提取单个基因的HR的95%CI低值
                        H95CI=round(tdcoxSummary$conf.int[,"upper .95"],3),#提取单个基因的HR的95%CI高值
                        pvalue=round(tdcoxSummary$coefficients[,"Pr(>|z|)"],3))#提取单个基因的p值
  )
}

result_dir = "../../data/SurvivalCox/"
if (!dir.exists(result_dir))
{
  dir.create(result_dir, recursive = T)
}

outResult
write.csv(outResult, paste(result_dir, "TCGA_LIHC_RFS_survival_cox_result.csv", sep = ''), row.names = F)

#2.2.1 读入前面输出的分析结果数据UniCoxSurvival.txt
tducs <- read.table("../../data/SurvivalCox/TCGA_LIHC_RFS_survival_cox_result.csv", header=T, sep=",", row.names=1, check.names=F)
row.names(tducs) = c("Age (>60 vs <=60)", "Gender (Male vs Female)", "Hepatitis (Yes vs No)", 
                     "Tumor size (>5cm vs <=5cm)", "Tumor number (Multiple vs Solitary)",
                     "PVTT (Present vs Absent)", "MVI (Present vs Absent)",
                     "TNM stage (III vs I+II)", "CNLC stage (III vs I+II)",
                     "BCLC (C vs B vs A)", "NLR (>3 vs<=3)",
                     "PLR (>130 vs<=130)", "ALBI grade(Grade 2 vs 1)",
                     "AFP level (>400ng/ml vs <=400ng/ml)", "CD161 (High vs Low)")

#2.2.2 提取制图相关数据
gene <- rownames(tducs)#提取基因名
hr <- sprintf("%.3f",tducs$"HR")#从tducs数据中提取HR值，%.3f指数据保留小数点后3位
hrLow  <- sprintf("%.3f",tducs$"L95CI")#提取HR值95%置信区间的低值
hrHigh <- sprintf("%.3f",tducs$"H95CI")#提取HR值95%置信区间的高值
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")#将HR及其置信区间数据组合成一个因子
pValue <- ifelse(tducs$pvalue<0.001, "<0.001", sprintf("%.3f", tducs$pvalue))#提取p值，且当p<0.001时显示为<0.001，而不是确切数值

#2.2.2 开始绘图
# result_dir = "figure/survival_cox/"
# if (!dir.exists(result_dir))
# {
#   dir.create(result_dir, recursive = T)
# }
pdf(file=paste(result_dir, "TCGA_LIHC_RFS_UniCoxSurForestPlot.pdf", sep = ""), width = 12, height = 6)#pdf绘图命令开始打印过程，设置图片宽度6及高度3（自行根据数据设置）
#可以不执行pdf绘图命令，这样执行命令时图像会一一展示出来，执行pdf命令后图像仅在最后输出至pdf文件中
n <- nrow(tducs)#提取出tducs的行数，也就是有多少个基因
nRow <- n+1 #设置一个比基因数多1的数值
ylim <- c(1,nRow) #暂时设置一个y轴长度值，长度比基因数多1个，总要留点空白，总不能跟基因数一样把

layout(matrix(c(1,2),nc=2),width=c(2.5,2)) # 设置页面排版，即1排两图，左图1，右图2。


#开始画图的左边部分，即图1，森林图的文字及数值标注部分
xlim = c(0,2.5)#设置一个x轴长度
par(mar=c(4,2.5,2,1))#设置图形1的边界，即;下4，左2.5,上2，右1
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")#开始绘制图1，此时为空白图，axex=F表示不显示边框，可以设置为T看看
text.cex=0.8 #设置图中文本的大小
text(0,n:1,gene,adj=0,cex=text.cex);text(0,n+1,'Univariate analysis',cex=text.cex,font=2,adj=0)#添加文本，即在图上添加一列基因名，从坐标0开始正向添加，n即前面定义的基因数，gene即基因名
text(1.4,n:1,pValue,adj=1,cex=text.cex);text(1.4,n+1,'P value',cex=text.cex,font=2,adj=1)#添加p值，adj=1从坐标1.4反向添加
text(2.5,n:1,Hazard.ratio,adj=1,cex=text.cex);text(2.5,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)#添加HR值，从坐标2.5反向添加
#图1即画完了，可以不执行pdf绘图命令，这样执行命令时图像会一一展示出来，执行pdf命令后图像仅在最后输出至pdf文件中


#开始画图的右边部分，即图2，森林图的图形
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))#设置边界，mpg设置横标题及坐标轴的位置
xlim = c(0,ceiling(max(as.numeric(hrLow),as.numeric(hrHigh)))) #设置图2的横坐标
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")#开始画图2，axex=F表示不显示边框
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)#HR的横条
abline(v=1,col="black",lty=2,lwd=2)#在横坐标1的地方画一条竖线
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')#设置颜色
points(as.numeric(hr), n:1, pch = 1, col = boxcolor, cex=1.3)#HR处显示为圈，颜色为上面设的颜色。pch可以设置不同的数字，显示不同的图形。
axis(1)#显示横坐标

#画图结束，输出pdf
dev.off()#输出pdf


############################
outResult = outResult[outResult$pvalue < 0.05,]
outResult$id
# cox_result = coxph(formula = Surv(RFS, RFSstatus) ~ Age + Gender + TNM_stage + Vascular_invasion + Risk_Score, data = df, method="breslow")
# cox_result = coxph(formula = Surv(RFS, RFSstatus) ~ TNM_stage + Risk_Score, data = df, method="breslow")

cox_result = step(coxph(Surv(RFS, RFSstatus) ~ Tumor.Size + Tumor.Number + PVTT + TNM.Stage.AJCC. + CNLC.Stage + BCLC + ALBI + AFP + CD161, data = df, method="breslow"), 
     direction = "forward", 
     k = log(nrow(df)), 
     trace = 0)

cox_result = coxph(Surv(RFS, RFSstatus) ~ Tumor.Size + Tumor.Number + PVTT + TNM.Stage.AJCC. + CNLC.Stage + BCLC + ALBI + AFP + CD161, data = df, method="efron")

cox_result = coxph(Surv(RFS, RFSstatus) ~ Tumor.Number + AFP + CD161, data = df, method="efron")

# cox_result = step(coxph(Surv(RFS, RFSstatus) ~ Tumor.Size + Tumor.Number + PVTT + TNM.Stage.AJCC. + CNLC.Stage + BCLC + ALBI + AFP + CD161, data = df, method="breslow"),
#                   direction = "forward")

cox_result_s = summary(cox_result)

df_c = data.frame(cox_result_s$conf.int[,"exp(coef)"], cox_result_s$conf.int[,"lower .95"], cox_result_s$conf.int[,"upper .95"],cox_result_s$coefficients[,"Pr(>|z|)"])
colnames(df_c) = c("HR", "L95CI", "H95CI", "pvalue")

df_c
# result_dir = "../../data/survival_cox/"
write.csv(df_c, paste(result_dir, "TCGA_LIHC_OS_survival_multi_cox_result.csv", sep = ''), row.names = T)

#2.2.1 读入前面输出的分析结果数据UniCoxSurvival.txt
tducs <- read.table(paste(result_dir, "TCGA_LIHC_OS_survival_multi_cox_result.csv", sep = ''), header=T, sep=",", row.names=1, check.names=F)
row.names(tducs) = c("Tumor size (>5cm vs <=5cm)", "Tumor number (Multiple vs Solitary)",
                     "PVTT (Present vs Absent)", "TNM stage (III vs I+II)", "CNLC stage (III vs I+II)",
                     "BCLC (C vs B vs A)", "AFP level (>400ng/ml vs <=400ng/ml)", "CD161 (High vs Low)")

#2.2.2 提取制图相关数据
gene <- rownames(tducs)#提取基因名
hr <- sprintf("%.3f",tducs$"HR")#从tducs数据中提取HR值，%.3f指数据保留小数点后3位
hrLow  <- sprintf("%.3f",tducs$"L95CI")#提取HR值95%置信区间的低值
hrHigh <- sprintf("%.3f",tducs$"H95CI")#提取HR值95%置信区间的高值
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")#将HR及其置信区间数据组合成一个因子
pValue <- ifelse(tducs$pvalue<0.001, "<0.001", sprintf("%.3f", tducs$pvalue))#提取p值，且当p<0.001时显示为<0.001，而不是确切数值

#2.2.2 开始绘图
# result_dir = "figure/survival_cox/"
# if (!dir.exists(result_dir))
# {
#   dir.create(result_dir, recursive = T)
# }
pdf(file=paste(result_dir, "TCGA_LIHC_MultiCoxSurForestPlot.pdf", sep = ""), width = 9,height = 3)#pdf绘图命令开始打印过程，设置图片宽度6及高度3（自行根据数据设置）
#可以不执行pdf绘图命令，这样执行命令时图像会一一展示出来，执行pdf命令后图像仅在最后输出至pdf文件中
n <- nrow(tducs)#提取出tducs的行数，也就是有多少个基因
nRow <- n+1 #设置一个比基因数多1的数值
ylim <- c(0.5,nRow + 1.5) #暂时设置一个y轴长度值，长度比基因数多1个，总要留点空白，总不能跟基因数一样把

layout(matrix(c(1,2),nc=2),width=c(2.5,2)) # 设置页面排版，即1排两图，左图1，右图2。


#开始画图的左边部分，即图1，森林图的文字及数值标注部分
xlim = c(0,2.5)#设置一个x轴长度
par(mar=c(4,2.5,2,1))#设置图形1的边界，即;下4，左2.5,上2，右1
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")#开始绘制图1，此时为空白图，axex=F表示不显示边框，可以设置为T看看
text.cex=0.8 #设置图中文本的大小
text(0,n:1,gene,adj=0,cex=text.cex);text(0,n+1,'Multivariate analysis',cex=text.cex,font=2,adj=0)#添加文本，即在图上添加一列基因名，从坐标0开始正向添加，n即前面定义的基因数，gene即基因名
text(1.4,n:1,pValue,adj=1,cex=text.cex);text(1.4,n+1,'P value',cex=text.cex,font=2,adj=1)#添加p值，adj=1从坐标1.4反向添加
text(2.5,n:1,Hazard.ratio,adj=1,cex=text.cex);text(2.5,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)#添加HR值，从坐标2.5反向添加
#图1即画完了，可以不执行pdf绘图命令，这样执行命令时图像会一一展示出来，执行pdf命令后图像仅在最后输出至pdf文件中


#开始画图的右边部分，即图2，森林图的图形
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))#设置边界，mpg设置横标题及坐标轴的位置
xlim = c(0,ceiling(max(as.numeric(hrLow),as.numeric(hrHigh))) + 1) #设置图2的横坐标
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")#开始画图2，axex=F表示不显示边框
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)#HR的横条
abline(v=1,col="black",lty=2,lwd=2)#在横坐标1的地方画一条竖线
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')#设置颜色
points(as.numeric(hr), n:1, pch = 1, col = boxcolor, cex=1.3)#HR处显示为圈，颜色为上面设的颜色。pch可以设置不同的数字，显示不同的图形。
axis(1)#显示横坐标

#画图结束，输出pdf
dev.off()#输出pdf

