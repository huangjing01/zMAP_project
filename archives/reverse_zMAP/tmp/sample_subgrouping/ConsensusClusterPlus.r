setwd('/mnt/storage/user/zmap_proteomics/script/guixiuqi/reverse_zMAP/sample_subgrouping')
library(ConsensusClusterPlus)
dc = read.csv('/mnt/storage/user/zmap_proteomics/script/guixiuqi/reverse_zMAP/sample_subgrouping/top_variance_tumor_z_df_dropna_3000_1673.txt',sep='	',row.names = 1,check.names=FALSE)
dc = as.matrix(dc)
rcc = ConsensusClusterPlus(dc,maxK=5,reps=1000,pItem=0.8,pFeature=1,title='euclidean_km',                            distance='euclidean',clusterAlg='km',plot='pdf',seed=1262118322)
cluster <- rcc[[3]]$consensusClass
write.csv(cluster,file='/mnt/storage/user/zmap_proteomics/script/guixiuqi/reverse_zMAP/sample_subgrouping/cluster_3.csv', quote = FALSE)
cluster <- rcc[[4]]$consensusClass
write.csv(cluster,file='/mnt/storage/user/zmap_proteomics/script/guixiuqi/reverse_zMAP/sample_subgrouping/cluster_4.csv', quote = FALSE)
cluster <- rcc[[5]]$consensusClass
write.csv(cluster,file='/mnt/storage/user/zmap_proteomics/script/guixiuqi/reverse_zMAP/sample_subgrouping/cluster_5.csv', quote = FALSE)
