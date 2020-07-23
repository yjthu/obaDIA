args <-commandArgs(TRUE);##sampleName ##output dir ##identify file
library(ggplot2)
path=read.table(args[3],sep="\t",header=T)

#------------------------------------------------------------
trim_string<-function(str,num){
    strTrim<-str
        if(nchar(str)>num){
            strTrim<-paste(substr(str,1,num),"...",sep="");
        }
    strTrim
}
change_name<-function(old_name){
        change_nu<-length(old_name)
        tmp0<-c(rep=0,times=change_nu)
        j=50
        for(i in 1:change_nu){
                tmp<-trim_string(as.character(old_name[i]),j)
                if(tmp %in% tmp0){
                        j=j+1
                }
                tmp0[i]<-trim_string(as.character(old_name[i]),j)
        }
        tmp0
}
#-------------------------------------------------------------


if(nrow(path)==0)
{
	cat("This compare has no DGE annotation, please check!\n")
	quit()
}
colnames(path)<-c("Pathway_term","Rich_factor","Qvalue","Number")
path$Pathway_term<-change_name(path$Pathway_term)
p<-ggplot(path, aes(Pathway_term,Rich_factor))
p<-p+geom_point(aes(colour=Qvalue,size=Number))+scale_colour_gradientn(colours=rainbow(2),guide = "colourbar") +expand_limits(color=seq(0, 1, by=0.25))
p<-p+ggtitle("Statistics of Enrichment") + ylab("Rich factor") +xlab("Term")
#p<-p+theme_bw()+theme(axis.text=element_text( color="black", size=10))
#p<-p+theme(panel.border=element_rect(colour = "black"))
p<-p+theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position="left")
p<-p+theme(plot.title=element_text(vjust=1), legend.key=element_blank())
p <- p + theme(panel.background = element_rect(fill='white', colour='grey'),panel.grid.minor=element_blank(),panel.grid.major=element_blank())

p
ggsave(paste(args[2],"/",args[1],".","Enrich.scatterplot.png",sep=""), plot=p,type="cairo-png", width=8, height=6, dpi=700)
ggsave(paste(args[2],"/",args[1],".","Enrich.scatterplot.pdf",sep=""), plot=p, width=8, height=6)
