#reverse-zMAP
#python reverse_zmap_step1_for_linux_for_xlsx_1.5.py --intensity_data small_data_protein_level_intensity.txt --sample_information small_data_sample_info_add_internal_ref.txt --outdir /mnt/storage/user/zmap_proteomics/script/guixiuqi/reverse_zMAP/reverse_zMAP_results --window_size 400 --step_size 100 --percent 50 --method natural_cubic_spline

############# downstream analysis #############


#sample quality control
#python sample_quality_control.py --z_statistic_matrix /mnt/storage/user/zmap_proteomics/script/guixiuqi/reverse_zMAP/reverse_zMAP_results/z_statistic_table.txt --sample_info small_data_sample_info_add_internal_ref.txt --outdir /mnt/storage/user/zmap_proteomics/script/guixiuqi/reverse_zMAP//sample_quality_control


#subgrouping
#python sample_subgrouping_ConsensusClusterPlus.py --z_statistic_matrix cervical_cancer_z_statistic.txt --sample_info cervical_cancer_sample_info.txt --sample_condition Tumor_tissue --outdir /mnt/storage/user/zmap_proteomics/script/guixiuqi/web_html_mannual/input_file/reverse_zMAP/sample_subgrouping --top_n 3000


#Association with clinical features
#python association_with_clinical_and_molecule_feature.py --z_statistic_matrix cervical_cancer_z_statistic.txt --cluster_f /mnt/storage/user/zmap_proteomics/script/guixiuqi/web_html_mannual/input_file/reverse_zMAP/sample_subgrouping/cluster_3.csv --clinical_info cervical_cancer_clinical_info.txt --discrete 'Histology,Degree,Stage,Lymph_node,SCC_ng_perl_ml_higher_than_1_point_5' --continuous 'Tumor_size_cm,Age' --color_f discrete_feature_color_RGB.txt --colorbar_f continuous_feature_colorbar.txt --outdir /mnt/storage/user/zmap_proteomics/script/guixiuqi/web_html_mannual/input_file/reverse_zMAP/association_with_clinical_and_molecule_feature  --fdr 0.05 --cluster_n 4 

#Survival analysis
#python sample_group_survival_analysis.py --input_file survival_analysis_data_for_web.txt --outdir /mnt/storage/user/zmap_proteomics/script/guixiuqi/reverse_zMAP/survival_analysis

#Association with survival data
#python association_z_statistic_with_survival.py --survival_f survival_analysis_data_for_web.txt --z_statistic_matrix HCC_z_statistic_df.txt --outdir /mnt/storage/user/zmap_proteomics/script/guixiuqi/reverse_zMAP/association_z_statistic_with_survival

#Association with mutation data
#python association_z_statistic_with_mutation_data.py --mutation_f /mnt/storage/user/zmap_proteomics/script/guixiuqi/reverse_zMAP/snp_df_10.txt --z_statistic_matrix /mnt/storage/user/zmap_proteomics/script/guixiuqi/reverse_zMAP/HCC_z_statistic_df.txt --covariates_f /mnt/storage/user/zmap_proteomics/script/guixiuqi/reverse_zMAP/covariates_df.txt --gene_tss_location /mnt/storage/user/zmap_proteomics/script/guixiuqi/reverse_zMAP/gene_tss_location.txt --chr_length /mnt/storage/user/zmap_proteomics/script/guixiuqi/reverse_zMAP/hg19_chrom_sizes_remove_xy.txt --outdir /mnt/storage/user/zmap_proteomics/script/guixiuqi/reverse_zMAP/association_z_statistic_with_mutation_data --fdr 0.05


python reverse_zmap_step1_for_linux_for_xlsx_1.5_2023_11_15.py --intensity_data ~/zMap/app/templates/input_file/reverse_zMAP/small_data_protein_level_intensity.txt --sample_information ~/zMap/app/templates/input_file/reverse_zMAP/small_data_sample_info_add_internal_ref.txt --outdir ~/zMap/app/templates/input_file/reverse_zMAP/test --window_size 400 --step_size 100 --percent 50 --method natural_cubic_spline

python /home/shao/zMap/app/templates/input_file/reverse_zMAP/sample_quality_control.py --z_statistic_matrix ~/zMap/app/templates/input_file/reverse_zMAP/test/z_statistic_table.txt --sample_info ~/zMap/app/templates/input_file/reverse_zMAP/small_data_sample_info_add_internal_ref.txt --outdir ~/zMap/app/templates/input_file/reverse_zMAP/sample_quality_control



python /home/shao/zMap/app/templates/input_file/reverse_zMAP/sample_subgrouping_ConsensusClusterPlus.py --z_statistic_matrix /home/shao/zMap/app/templates/input_file/reverse_zMAP/cervical_cancer_z_statistic.txt --sample_info /home/shao/zMap/app/templates/input_file/reverse_zMAP/cervical_cancer_sample_info.txt --sample_condition Tumor_tissue --outdir /home/shao/zMap/app/templates/input_file/reverse_zMAP/sample_subgrouping --top_n 3000

python /home/shao/zMap/app/templates/input_file/reverse_zMAP/association_with_clinical_and_molecule_feature.py --z_statistic_matrix /home/shao/zMap/app/templates/input_file/reverse_zMAP/cervical_cancer_z_statistic.txt --cluster_f /home/shao/zMap/app/templates/input_file/reverse_zMAP/sample_subgrouping/cluster_3.csv --clinical_info /home/shao/zMap/app/templates/input_file/reverse_zMAP/cervical_cancer_clinical_info.txt --discrete 'Histology,Degree,Stage,Lymph_node,Abnormal_SCC' --continuous 'Tumor_size_cm,Age' --color_f /home/shao/zMap/app/templates/input_file/reverse_zMAP/discrete_feature_color_RGB.txt --colorbar_f  /home/shao/zMap/app/templates/input_file/reverse_zMAP/continuous_feature_colorbar.txt --outdir /home/shao/zMap/app/templates/input_file/reverse_zMAP/association_with_clinical_and_molecule_feature --fdr 0.05 --cluster_n 4


python /home/shao/zMap/app/templates/input_file/reverse_zMAP/sample_group_survival_analysis.py --input_file /home/shao/zMap/app/templates/input_file/reverse_zMAP/survival_analysis_data_for_web.txt --outdir /home/shao/zMap/app/templates/input_file/reverse_zMAP/survival_analysis

python /home/shao/zMap/app/templates/input_file/reverse_zMAP/association_z_statistic_with_survival.py --survival_f /home/shao/zMap/app/templates/input_file/reverse_zMAP/survival_analysis_data_for_web.txt --z_statistic_matrix /home/shao/zMap/app/templates/input_file/reverse_zMAP/HCC_z_statistic_df.txt --outdir /home/shao/zMap/app/templates/input_file/reverse_zMAP/association_z_statistic_with_mutation_data


python /home/shao/zMap/app/templates/input_file/reverse_zMAP/association_z_statistic_with_mutation_data.py --mutation_f /home/shao/zMap/app/templates/input_file/reverse_zMAP/snp_df_10.txt --z_statistic_matrix /home/shao/zMap/app/templates/input_file/reverse_zMAP/HCC_z_statistic_df.txt --covariates_f /home/shao/zMap/app/templates/input_file/reverse_zMAP/covariates_df.txt --gene_tss_location /home/shao/zMap/app/templates/input_file/reverse_zMAP/gene_tss_location.txt --chr_length /home/shao/zMap/app/templates/input_file/reverse_zMAP/hg19_chrom_sizes_remove_xy.txt --outdir /home/shao/zMap/app/templates/input_file/reverse_zMAP/association_z_statistic_with_mutation_data --fdr 0.05
























