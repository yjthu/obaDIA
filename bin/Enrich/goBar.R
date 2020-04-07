### ARGS: [enrichment result] [output dir] [title]

#==========================================================
# the input enrichment result like :
# GO_accession	Description	Term_type	Over_represented_pValue	Corrected_pValue	DEG_item	DEG_list	Bg_item	Bg_list	up	down
# GO:0008152	metabolic process	biological_process	4.60E-34	1.52E-30	968	1512	12742	26634	494	474
#
# Author : yj
#==========================================================

args<-commandArgs(TRUE)
input=args[1]
dir=args[2]
title=args[3]

library(ggplot2)
library(reshape2)
setwd(dir)
df<-read.csv(input,header=TRUE,sep="\t",fill=TRUE)
padj=ifelse(df[8]<1,'*','')
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
# bar plot
data<-cbind(df[1],df[4],df[5],padj)
data<-head(data,30)
colnames(data)<-c("name","type","number","padj")

data$name<-change_name(data$name)

bp<-subset(data,data$type=="biological_process")
cc<-subset(data,data$type=="cellular_component")
mf<-subset(data,data$type=="molecular_function")

mt<-rbind(bp,cc,mf)
colnames(mt)<-c("name","type","number","padj")
mt$name <- factor(mt$name, levels = mt$name)

p <- ggplot(data=rev(mt), aes(x=name,y=number,fill=type))+
geom_bar(position=position_dodge(),stat="identity")+
#coord_flip()+
xlab("GO term") + ylab("Number") + 
labs(title=paste("The Most Enriched GO Terms (",title,")",sep=""))+
scale_fill_brewer(palette="Accent")+
theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position="left")+
geom_text(aes(label=padj), vjust=0.5,hjust=0.5,color=1)
p <- p + theme(panel.background = element_rect(fill='white', colour='grey'),panel.grid.minor=element_blank(),panel.grid.major=element_blank())

ggsave(filename=paste(title,".go_bar.pdf",sep=""), plot=p, height=6, width=8)
ggsave(filename=paste(title,".go_bar.png",sep=""), type="cairo-png", plot=p, height=6, width=8)
