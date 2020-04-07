args<-commandArgs(TRUE)
workdir=args[1]
goResource=args[2]
goEnrich=args[3]
diff=args[4]
output=args[5]

##===================================================================================
setwd(workdir)
library(GO.db)
library(goTools)
library(topGO)

goseq<-read.csv(args[3],header=TRUE,sep="\t")

source(args[2]);  ##GO annotation data
universe<-names(gene2GO);
genelist <- scan(args[4], what="character",quiet= TRUE);
gene.vector=factor(as.integer(universe %in% genelist));
names(gene.vector)=universe;

#goseq<-subset(goseq,goseq$Corrected_P_Value<0.05)
goseq<-subset(goseq,goseq$P_Value<0.05)
rownames(goseq)<-goseq$ID;

if(sum(goseq$Term_type=="biological_process")>0){
	GObpdata<-new("topGOdata",description="BP",ontology="BP",allGenes=gene.vector,annot=annFUN.gene2GO,gene2GO=gene2GO,nodeSize=10);
	bpNodes<-nodes(graph(GObpdata));
	bpScores<-goseq[bpNodes,5];
	names(bpScores)<-bpNodes;
	if(sum(!is.na(bpScores))>10){
		bps<-10;
	}else{
		bps<-sum(!is.na(bpScores));
	}
	if(bps>=1){
		pdf(paste(output,".GO_BP_DAG.pdf",sep=""));
		showSigOfNodes(GObpdata,useInfo="all",bpScores,firstSigNodes=bps);
		dev.off();

		png(paste(output,".GO_BP_DAG.png",sep=""),type="cairo-png",width=480*4,height=480*4,res=72*4);
		showSigOfNodes(GObpdata,useInfo="all",bpScores,firstSigNodes=bps);
		dev.off();
	}
}

if(sum(goseq$Term_type=="cellular_component")>0){
	GOccdata<-new("topGOdata",description="CC",ontology="CC",allGenes=gene.vector,annot=annFUN.gene2GO,gene2GO=gene2GO,nodeSize=10);
	ccNodes<-nodes(graph(GOccdata));
	ccScores<-goseq[ccNodes,5];
	names(ccScores)<-ccNodes;
	if(sum(!is.na(ccScores))>10){
		ccs<-10;
	}else{
		ccs<-sum(!is.na(ccScores));
	}
	if(ccs>=1){
		pdf(paste(output,".GO_CC_DAG.pdf",sep=""));
		showSigOfNodes(GOccdata,useInfo="all",ccScores,firstSigNodes=ccs);
		dev.off();

		png(paste(output,".GO_CC_DAG.png",sep=""),type="cairo-png",width=480*4,height=480*4,res=72*4);
		showSigOfNodes(GOccdata,useInfo="all",ccScores,firstSigNodes=ccs);
		dev.off();
	}
}
if(sum(goseq$Term_type=="molecular_function")>0){
	GOmfdata<-new("topGOdata",description="MF",ontology="MF",allGenes=gene.vector,annot=annFUN.gene2GO,gene2GO=gene2GO,nodeSize=10);
	mfNodes<-nodes(graph(GOmfdata));
	mfScores<-goseq[mfNodes,5];
	names(mfScores)<-mfNodes;
	if(sum(!is.na(mfScores))>10){
		mfs<-10;
	}else{
		mfs<-sum(!is.na(mfScores));
	}
	if(mfs>=1){
		pdf(paste(output,".GO_MF_DAG.pdf",sep=""));
		showSigOfNodes(GOmfdata,useInfo="all",mfScores,firstSigNodes=mfs);
		dev.off();
	
		png(paste(output,".GO_MF_DAG.png",sep=""),type="cairo-png",width=480*4,height=480*4,res=72*4);
		showSigOfNodes(GOmfdata,useInfo="all",mfScores,firstSigNodes=mfs);
		dev.off();
	}
}

#======================================================================================
