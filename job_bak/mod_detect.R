
library(dynamicTreeCut)
library(igraph)
library(RColorBrewer)
args <- commandArgs(T)
test_pcor<-read.table(args[1],sep="\t",header=T)
# raw_data <- read.table(args[3],sep="\t",header=T)
# raw_data <- raw_data[,grep("z.statistic",colnames(raw_data))]
method<-args[2]
test_pcor<-rbind(test_pcor[,c("node1","node2")],data.frame(node1=test_pcor$node2,node2=test_pcor$node1))
test_pcor$prob<-1
pcor_mat<-reshape2::acast(test_pcor,node1~node2,value.var="prob")
pcor_mat[is.na(pcor_mat)]<-0
DT<-as.matrix(pcor_mat) %*% t(as.matrix(pcor_mat))
DT_k<-diag(DT)
DT_k_m<-sapply(DT_k,function(x){
	ifelse(DT_k>x,x,DT_k)
	})
DT_w<-(DT+pcor_mat)/(DT_k_m+1-pcor_mat)
diag(DT_w)<-1
if(method=="static"){
	mod_score<-c()
	for(h in seq(0.3,0.9,0.05)){
		# mergeCloseModules(t(raw_data[rownames(DT_w),]),cutreeDynamic(hclust(as.dist(1-DT_w)),distM=1-DT_w))$colors
		sample_clust<-cutree(hclust(as.dist(1-DT_w),"ave"),h=h)
		clust<-names(table(sample_clust))[which(table(sample_clust)>5)]
		if(length(clust) > 28){
			mod_score<-c(mod_score,0)
		}else{
			clust_gene<-rownames(pcor_mat)[which(sample_clust %in% clust)]
			gene_net<-test_pcor[intersect(which(test_pcor$node1 %in% clust_gene),which(test_pcor$node2 %in% clust_gene)),]
			gene_net<-subset(gene_net,as.character(node1)>as.character(node2))
			edges<-c()
			labels<-c()
			for(i in 1:nrow(gene_net)){
				edges<-c(edges,as.character(gene_net[i,1]),as.character(gene_net[i,2]))
				labels<-c(labels,sample_clust[which(names(sample_clust)==as.character(gene_net[i,1]))],sample_clust[which(names(sample_clust)==as.character(gene_net[i,2]))])
			}
			membership<-data.frame(edges=edges,labels=match(labels,unique(labels)))
			membership_label<-membership$labels[which(!duplicated(membership$edges))]
			g<-make_graph(edges,directed=FALSE)
			mod_score<-c(mod_score, modularity(g,membership_label))
		}	
	}
	h<-seq(0.3,0.9,0.05)[which.max(mod_score)]
	sample_clust<-cutree(hclust(as.dist(1-DT_w),"ave"),h=h)
	clust<-names(table(sample_clust))[which(table(sample_clust)>5)]
	clust_gene<-rownames(pcor_mat)[which(sample_clust %in% clust)]
	gene_net<-test_pcor[intersect(which(test_pcor$node1 %in% clust_gene),which(test_pcor$node2 %in% clust_gene)),]
	gene_net$prob<-NULL
	gene_net<-subset(gene_net,as.character(node1)>as.character(node2))
	gene_label<-data.frame(Protein=clust_gene,label=as.character(sample_clust)[match(clust_gene,rownames(pcor_mat))])
	gene_label$label<-match(gene_label$label,unique(gene_label$label))
	
}else{
	sample_clust<-cutreeDynamic(hclust(as.dist(1-DT_w)),distM=1-DT_w,minClusterSize = 8)
	clust<-names(table(sample_clust))[which(table(sample_clust)>5)]
	clust<-clust[clust!="0"]
  mod_score<-c()
	if(length(clust) > 30){
		mod_score<-c(mod_score,0)
	}else{
		clust_gene<-rownames(pcor_mat)[which(sample_clust %in% clust)]
		gene_net<-test_pcor[intersect(which(test_pcor$node1 %in% clust_gene),which(test_pcor$node2 %in% clust_gene)),]
		gene_net<-subset(gene_net,as.character(node1)>as.character(node2))
		edges<-c()
		labels<-c()
		for(i in 1:nrow(gene_net)){
			edges<-c(edges,as.character(gene_net[i,1]),as.character(gene_net[i,2]))
			labels<-c(labels,sample_clust[which(rownames(pcor_mat)==as.character(gene_net[i,1]))],sample_clust[which(rownames(pcor_mat)==as.character(gene_net[i,2]))])
		}
		membership<-data.frame(edges=edges,labels=match(labels,unique(labels)))
		membership_label<-membership$labels[which(!duplicated(membership$edges))]
		g<-make_graph(edges,directed=FALSE)
		mod_score<-modularity(g,membership_label)
		
		
	}
}
gene_label<-data.frame(Protein=clust_gene,label=as.character(sample_clust)[match(clust_gene,rownames(pcor_mat))])
gene_label$label<-match(gene_label$label,unique(gene_label$label))
gene_net$prob<-NULL
write.table(gene_net,paste0(c(unlist(strsplit(args[1],"/"))[1:(length(unlist(strsplit(args[1],"/")))-1)],"final_pcor_protein_net.txt"),collapse="/"),sep="\t",row.names=F,col.names=T,quote=F)
write.table(gene_label,paste0(c(unlist(strsplit(args[1],"/"))[1:(length(unlist(strsplit(args[1],"/")))-1)],"final_pcor_protein_label.txt"),collapse="/"),sep="\t",row.names=F,col.names=T,quote=F)
color_opt<-c(brewer.pal(8,"Dark2"),brewer.pal(8,"Set2"),brewer.pal(12,"Paired"))
label_color<-data.frame(label=unique(gene_label$label),color=color_opt[unique(gene_label$label)])
write.table(label_color,paste0(c(unlist(strsplit(args[1],"/"))[1:(length(unlist(strsplit(args[1],"/")))-1)],"cluster_mapped_color.txt"),collapse="/"),col.names=F,row.names=F,quote=F,sep="\t")
files <- list.files(paste0(unlist(strsplit(args[1],"/"))[1:(length(unlist(strsplit(args[1],"/")))-1)],collapse="/"))
file <- paste0(paste0(unlist(strsplit(args[1],"/"))[1:(length(unlist(strsplit(args[1],"/")))-1)],collapse="/"),"/",files[grep(".xls",files)])

all_gene_number <- as.numeric(unlist(strsplit(system(paste0("wc -l ",file),intern=TRUE)," "))[1])-1

module_number<-c()
module_anno<-c()
module_padj<-c()
module_cover<-c()
gene_ontology<-strsplit(readLines("~/zMap/jobs/gmt_files/c5.all.v7.1.symbols.gmt"),"\t")
for(i in unique(gene_label$label)){
	gene_set<-gene_label$Protein[which(gene_label$label==i)]
	e_1<-c()
	e_2<-c()
	e_3<-c()
	for(j in 1:length(gene_ontology)){
	gene_intersect<-intersect(gene_set,gene_ontology[[j]][3:length(gene_ontology[[j]])])
	num1<-length(gene_intersect)
	if(num1!=0){
	  num2<-length(gene_set)-num1
	  num3<-length(gene_ontology[[j]])-num1-2
	  num4<-all_gene_number-num1-num2-num3
	  four_test<-fisher.test(matrix(c(num1,num2,num3,num4),nrow=2),alternative="greater")

	  e_1<-c(e_1,gsub("_"," ",tolower(substring(as.character(gene_ontology[[j]][1]),4))))
	  e_2<-c(e_2,four_test$p.value)
	  e_3<-c(e_3,paste(num1,"/",num1+num2))
	}else{
	  e_1<-c(e_1,as.character(gene_ontology[[j]][1]))
	  e_2<-c(e_2,1)
	  e_3<-c(e_3,"")
	}

	}
	enrichment_result<-data.frame(e_1,e_2,e_3)
	enrichment_result$q<-p.adjust(e_2,method="BH",n=length(e_2))
	enrichment_result<-enrichment_result[which(enrichment_result$q<0.05),]
	enrichment_result<-enrichment_result[order(enrichment_result$q),]
	if(nrow(enrichment_result)>5){
		module_number<-c(module_number,rep(paste0("clust",i),5))
		module_anno<-c(module_anno,as.character(enrichment_result$e_1[1:5]))
		module_cover<-c(module_cover,as.character(enrichment_result$e_3[1:5]))
		module_padj<-c(module_padj,enrichment_result$q[1:5])
	}else{
		module_number<-c(module_number,rep(paste0("clust",i),nrow(enrichment_result)))
		module_anno<-c(module_anno,as.character(enrichment_result$e_1))
		module_cover<-c(module_cover,as.character(enrichment_result$e_3))
		module_padj<-c(module_padj,enrichment_result$q)
	}
}

write.table(data.frame(module_number,module_anno,module_padj,module_cover),paste0(c(unlist(strsplit(args[1],"/"))[1:(length(unlist(strsplit(args[1],"/")))-1)],"cluster_go_enrichment.txt"),collapse="/"),col.names=F,row.names=F,quote=F,sep="\t")