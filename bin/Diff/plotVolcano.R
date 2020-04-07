args <-commandArgs(TRUE)
library("ggplot2")

input <- args[1]
fdr_cutoff <- as.numeric(args[2])
fc_cutoff <- as.numeric(args[3])
fc_cutoff2 <- 0-fc_cutoff

#library(ggthemes) #requires R >= 3.3.0
temp <- read.csv(input,header = T,check.names = F,quote="",sep = "\t")
temp <- na.omit(temp)
temp$threshold[temp$FDR >=0 ] = "non"
temp$threshold[temp$FDR <= fdr_cutoff & temp$log2FC>=fc_cutoff ] = "up"
temp$threshold[temp$FDR <= fdr_cutoff & temp$log2FC<=fc_cutoff2 ] = "down"


p<-ggplot(temp,aes(x=temp$log2FC,y=-log10(temp$FDR),colour=threshold))+xlab("log2 Fold Change")+ylab("-log10 FDR")+
	geom_point(size=1.5,alpha=0.8)+
	scale_color_manual(values =c("#0072B5","darkgrey","#BC3C28")) +
	geom_hline(aes(yintercept=-log10(fdr_cutoff)),colour="grey",size=1 ,linetype=2) + 
	geom_vline(aes(xintercept=fc_cutoff), colour="grey",size=1 ,linetype=2)+ 
	geom_vline(aes(xintercept=fc_cutoff2), colour="grey",size=1 ,linetype=2) 
#	theme_few()+theme(legend.title = element_blank()) 
p <- p + theme(panel.background = element_rect(fill='white', colour='black'),panel.grid.minor=element_blank(),panel.grid.major=element_blank())


ggsave(paste(args[1],".","volcano.png",sep=""), plot=p,type="cairo-png", width=6, height=6, dpi=700)
ggsave(paste(args[1],".","volcano.pdf",sep=""), plot=p, width=6, height=6)
