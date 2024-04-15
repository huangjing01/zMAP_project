#zMAP
python zmap_step1_quantity_anaysis_2023_9_14.py --protein_intensity_file raw_protein_intensity_in_gene_level_for_web.txt --sample_info zmap_sample_info_for_web.txt --outdir /mnt/storage/user/zmap_proteomics/script/guixiuqi/zMAP/zMAP_results --window_size 400 --step_size 100 --percent 30 --method exponential_function

######## downstream analysis ###########


#HVPs identification and clustering
python zmap_hypervariable_proteins_cluster.py --pvalue_results /mnt/storage/user/zmap_proteomics/script/guixiuqi/zMAP/zMAP_results/zmap_chi_square_pvalue.txt --z_statistic_matrix /mnt/storage/user/zmap_proteomics/script/guixiuqi/zMAP/zMAP_results/z_statistic_table.txt --sample_info zmap_sample_info_for_web.txt --outdir /mnt/storage/user/zmap_proteomics/script/guixiuqi/zMAP/zmap_hypervariable_proteins_cluster --cluster_number_for_hypervariable 15 --minclustersize 20 --top 100 --cluster_number_for_top_proteins 8 --bh_pvalue_cutoff 0.05


#Gene set variation analysis
python gsva.py --z_statistic_matrix /mnt/storage/user/zmap_proteomics/script/guixiuqi/zMAP/zMAP_results/z_statistic_table.txt --sample_info /mnt/storage/user/zmap_proteomics/script/guixiuqi/zMAP/gsva_sample_info.txt --outdir /mnt/storage/user/zmap_proteomics/script/guixiuqi/zMAP/gsva --top_n 50 --adj_P_Val 0.05
