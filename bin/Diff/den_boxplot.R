args <-commandArgs(TRUE)

library('reshape2')
library("ggplot2")
library("pheatmap")

rp<-read.table(args[1],head=T)

#==============================================================================
# prepare data for density and box plotting 
#
df<-log10(rp[-1]+1)
df<-melt(df)
colnames(df)<-c("Group","value")
#------------------------------------------------------------------------------
p<-ggplot(df, aes(x=value, colour=Group, group=Group, fill=Group)) +  geom_density(alpha=0.4) + xlab("relative abundance") + ylab("density") + labs(title="")
p <- p + theme(panel.background = element_rect(fill='white', colour='black'),panel.grid.minor=element_blank(),panel.grid.major=element_blank())

ggsave(paste(args[1],".","density.png",sep=""), plot=p,type="cairo-png", width=6, height=5, dpi=700)
ggsave(paste(args[1],".","density.pdf",sep=""), plot=p, width=6, height=5)

#-1----------------------------------------------------------------------------
p<-ggplot(df, aes(Group, value)) + geom_boxplot(aes(fill = Group),alpha=0.8,outlier.size=1) + xlab("samples") + ylab("relative abundance") + labs(title="")
p <- p + theme(panel.background = element_rect(fill='white', colour='black'),panel.grid.minor=element_blank(),panel.grid.major=element_blank())

ggsave(paste(args[1],".","boxplot.png",sep=""), plot=p,type="cairo-png", width=6, height=5, dpi=700)
ggsave(paste(args[1],".","boxplot.pdf",sep=""), plot=p, width=6, height=5)

#-1----------------------------------------------------------------------------
p<-ggplot(df, aes(Group, value,fill=Group,colour=Group)) +geom_violin(alpha=0.4, width=1)+geom_boxplot(alpha=0.8,width=0.1,outlier.colour=NA) + xlab("samples") + ylab("relative abundance")+ labs(title="")
p <- p + theme(panel.background = element_rect(fill='white', colour='black'),panel.grid.minor=element_blank(),panel.grid.major=element_blank())

ggsave(paste(args[1],".","violin.png",sep=""), plot=p,type="cairo-png", width=6, height=5, dpi=700)
ggsave(paste(args[1],".","violon.pdf",sep=""), plot=p, width=6, height=5)



