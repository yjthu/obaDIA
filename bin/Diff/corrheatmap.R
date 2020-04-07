##############################################################
args <-commandArgs(TRUE)

library(corrplot)
library(pheatmap)
library(PerformanceAnalytics)


df <- read.csv(args[1],header = T,check.names = F,quote="",sep = "\t")
df <- df[,2:dim(df)[2]]
df <- na.omit(df)
df <- log2(df+1)


### 绘制相关性系数热图
corr<-cor(df,use="complete.obs")
symnum(corr)

pdf(paste(args[1],".","corHeatmap.pdf",sep=""), width=6, height=6, onefile = FALSE)
pheatmap(corr)
dev.off()
png(paste(args[1],".","corHeatmap.png",sep=""), type="cairo-png", width=600, height=600)
pheatmap(corr)
dev.off()


### 绘制相关性系数矩阵图
#color<-colorRampPalette(c("yellow","green4","purple"))
#corrplot(corr,method="number",tl.cex=0.9,number.cex=0.8,col=color(20))
#corrplot(corr=corr,method='color',order="AOE",addCoef.col="green",tl.cex=0.9,number.cex=1.0)
#corrplot(corr=corr,method='pie', tl.cex=1.0,number.cex=1.6)

### To plot correlationplot with 'PerformanceAnalytics' packages.
 
size=dim(df)[2]+1
size2=size*100
pdf(paste(args[1],".","corMatrix.pdf",sep=""), width=size, height=size, onefile = FALSE)
chart.Correlation(df,histogram=TRUE,pch=19)
dev.off()
png(paste(args[1],".","corMatrix.png",sep=""), type="cairo-png", width=size2, height=size2)
chart.Correlation(df,histogram=TRUE,pch=19)
dev.off()


### 绘制表达丰度热图

pdf(paste(args[1],".","abHeatmap.pdf",sep=""), width=6, height=8, onefile = FALSE)
pheatmap(df,color=colorRampPalette(rev(c("red","white","blue")))(100))
dev.off()
png(paste(args[1],".","abHeatmap.png",sep=""), type="cairo-png", width=600, height=800)
pheatmap(df,color=colorRampPalette(rev(c("red","white","blue")))(100))
dev.off()
