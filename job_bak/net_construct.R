args <- commandArgs(T)
library(GeneNet)
library(DMwR)
intensity_df <- read.csv(args[1],sep="\t",header=T,row.names=1)
sig_gene <- rownames(intensity_df)[which(intensity_df$Bonferroni.corrected.P.value<0.05)]
intensity_df <- intensity_df[,grep("z.statistic",colnames(intensity_df))]
if(length(which(is.na(intensity_df)))>0){
	intensity_df <- knnImputation(intensity_df,k=10,scale=F,meth="weighAvg",distData = NULL)
}
pcor<-ggm.estimate.pcor(t(intensity_df))
pcor_test<-network.test.edges(pcor,plot=FALSE)
pcor_gene_p <- pcor_test[which(pcor_test$prob>0.8 & pcor_test$pcor>0),]
pcor_gene_p$node1<-as.character(rownames(intensity_df)[pcor_gene_p$node1])
pcor_gene_p$node2<-as.character(rownames(intensity_df)[pcor_gene_p$node2])
n <- intersect(which(pcor_gene_p$node1 %in% sig_gene),which(pcor_gene_p$node2 %in% sig_gene))
pcor_gene_p<-pcor_gene_p[n,]
write.table(pcor_gene_p,paste0(c(unlist(strsplit(args[1],"/"))[1:(length(unlist(strsplit(args[1],"/")))-1)],"Protein_partial_correlation.txt"),collapse="/"),sep="\t",row.names=F,col.names=T,quote=F)


