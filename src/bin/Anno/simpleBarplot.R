args<-commandArgs(TRUE)
input=args[1]
h=as.numeric(args[2])
w=as.numeric(args[3])
#-------------------------------------------------

library("ggplot2")
require("RColorBrewer")

df<-read.table(input)
colnames(df)<-c("id","class","cate")
colourCount = length(unique(df$class))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

p <- ggplot(df, aes(class,fill=class))
p <- p+geom_bar()
#p <- p + theme(axis.text.x=element_text(hjust=1,angle=45))
p <- p + theme(panel.background = element_rect(fill='white', colour='grey'),panel.grid.minor=element_blank(),panel.grid.major=element_blank())
p <- p + theme(legend.position="right")
p <- p + scale_fill_manual(values = getPalette(colourCount),name="categories",breaks=sort(as.vector(df$class)),labels=sort(as.vector(df$cate)))
p <- p + labs(title="") + xlab("categories")
ggsave(filename=paste(input,".bar.pdf",sep=''), height=h, width=w,plot=p)
ggsave(filename=paste(input,".bar.png",sep=''), height=h, width=w,type="cairo-png", plot=p)
