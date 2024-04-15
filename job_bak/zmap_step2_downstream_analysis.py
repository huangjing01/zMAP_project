# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 10:12:48 2020

@author: guixiuqi
"""

import pandas as pd
import numpy as np
#import Cluster_Ensembles as CE
from sklearn.cluster import KMeans
from collections import defaultdict
import os
import seaborn as sns
import argparse
from sklearn.decomposition import PCA 
from matplotlib.lines import Line2D
import scipy

import random
import matplotlib.pyplot as plt
from collections import Counter

def read_file(filein):
    filein_df = pd.read_csv(filein,sep = "\t",header = [0,1],index_col = 0)
    return filein_df

def get_hypervariable_df(filein_df,pvalue_cutoff = 0.01):
    
    hypervariable_df = filein_df.loc[filein_df["-"]["BH-corrected average pvalue"] < 0.01]
    
    return hypervariable_df

def get_z_trans_df(hypervariable_df):
    columns_set = list(hypervariable_df.columns)
    header1 = []
    header2 = []
    sample_num = 0
    z_trans = []
    for i in columns_set:
        header1.append(i[0])
        header2.append(i[1])
        if i[1].startswith("mean"):
            sample_num += 1
        if i[1].startswith("z-transformed"):
            if i[0] in z_trans:
                pass
            else:
                z_trans.append(i[1])            
    rep_num = len(set(header1))-1    
    z_trans_df = hypervariable_df.loc(axis=1)[:,z_trans]
    return z_trans_df

def rename_columns(z_trans_df):
    recolumn_z_trans_df = z_trans_df.copy()
    new_column = []
    for column in z_trans_df.columns:
        rep = column[0]
        sample = column[1].split(" ")[-1]        
        new_column.append(rep + " " + sample)
    recolumn_z_trans_df.columns = new_column
    return recolumn_z_trans_df

def rearrange_column(recolumn_z_trans_df):
    condition_dict = {}
    for column in recolumn_z_trans_df.columns:
        condition = column.split(" ")[1]
        if condition in condition_dict.keys():
            condition_dict[condition].append(column)
        else:
            condition_dict[condition] = [column]
    new_column = []
    for key in condition_dict.keys():
        new_column += condition_dict[key]
    rearrange_z_trans_df = recolumn_z_trans_df[new_column]
    return rearrange_z_trans_df
        
def pcc(recolumn_z_trans_df,outdir):
    pearsonr_array=[]
    for i in recolumn_z_trans_df.columns:
        pcc_list=[]
        for j in recolumn_z_trans_df.columns:
            if i==j:
                pcc_list.append(1)
            else:
                two_columns_df = recolumn_z_trans_df.loc[:,[i,j]]
                new_two_columns_df = two_columns_df.loc[(two_columns_df[i].notnull()) & (two_columns_df[j].notnull())]
                print (i +"\t"+j+"\t"+"Number of genes detected simultaneously:"+str(len(new_two_columns_df)))
                corr = scipy.stats.pearsonr(new_two_columns_df.loc[:,i],new_two_columns_df.loc[:,j])[0]
                pcc_list.append(corr)
        pearsonr_array.append(pcc_list)
    pearsonr_array=pd.DataFrame(pearsonr_array,index=recolumn_z_trans_df.columns,columns = recolumn_z_trans_df.columns)
    pearsonr_array.to_csv(outdir + "/pearsonr_correlation.xlsx",sep="\t")
    return pearsonr_array

def get_color(columns):
    rep_list = []
    condition_list = []
    for column in columns:
        rep = column.split(" ")[0]
        condition = column.split(" ")[1]
        if rep in rep_list:
            pass
        else:
            rep_list.append(rep)
        if condition in condition_list:
            pass
        else:
            condition_list.append(condition)
    rep_num = len(rep_list)
    condition_num = len(condition_list)
    
    rep_color_list = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])for i in range(rep_num)]
    
    condition_color_list = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(condition_num)]
    
    rep_marker_list = ["s","^","o","d","P","<",">","v"]    
    rep_color_dict = dict(zip(rep_list,rep_color_list))
    rep_marker_dict = dict(zip(rep_list,rep_marker_list[0:rep_num]))
    condition_color_dict = dict(zip(condition_list,condition_color_list))
    
    batch_color = []
    condition_color = []
    batch_marker = []
    for column in columns:
        rep = column.split(" ")[0]
        condition = column.split(" ")[1]
        batch_color.append(rep_color_dict[rep])
        condition_color.append(condition_color_dict[condition])
        
    df = pd.DataFrame([batch_color,condition_color],index=["batch","condition"],columns=columns)
    df_T = df.T
    print(df_T)
    return rep_color_dict,rep_marker_dict,condition_color_dict,df_T
    
def hierarchical_clustering(pearsonr_array,df_T,outdir):
    g = sns.clustermap(pearsonr_array,cmap = "bwr",linewidths = 0.5,
                   cbar_kws = {"shrink": .5},vmin = -1,vmax = 1,col_colors = df_T,xticklabels= False,yticklabels=pearsonr_array.columns)
    g.fig.suptitle("     Pearsonr correlation coefficient of z-transformed",size=20)
    plt.savefig(outdir + "/pearsonr_correlation_coefficient_of_z_transformed_hierarchical_clustering.png",dpi = 300,bbox_inches = "tight")
    return g.data2d.columns


   
def pca_analysis(pearsonr_array,outdir):
    X_new =  pearsonr_array
    pca = PCA(n_components=2)
    X_r = pca.fit(X_new).transform(X_new)
    print('explained variance ratio (first three components): %s'
      % str(pca.explained_variance_ratio_))
    pc_df = pd.DataFrame(X_r,index = X_new.index,columns = ["PC1","PC2"])
    pc_df.to_csv(outdir + "/pca_analysis.xlsx",sep="\t")
    return pca,pc_df


def pca_scatter(pc_df,pca,rep_marker_dict,condition_color_dict,outdir):
    
    pc1_ratio = "%.2f"% (pca.explained_variance_ratio_[0]*100)
    pc2_ratio = "%.2f"% (pca.explained_variance_ratio_[1]*100)
    
    plt.figure(figsize=(5,5))
    for index in pc_df.index:
        rep = index.split(" ")[0]
        condition = index.split(" ")[1]
        plt.scatter(pc_df.loc[index,"PC1"],pc_df.loc[index,"PC2"],
                    edgecolors = condition_color_dict[condition],
                    linewidths=2,facecolors="white",marker=rep_marker_dict[rep],s=70)
    plt.xlabel("PC1(%s"%(pc1_ratio)+"%)",size=15)
    plt.ylabel("PC2(%s"%(pc2_ratio)+"%)",size=15)
    
    legend_elements = []
    
    for key in rep_marker_dict.keys():
        
        legend_elements.append(Line2D([0], [0],color="black", 
                                      markeredgewidth = 2,
                                      markerfacecolor="white",
                                      marker = rep_marker_dict[key],label=key,
                                      linestyle='None',
                                      markersize=8))
    for key in condition_color_dict.keys():
        color = condition_color_dict[key]
        legend_elements.append(Line2D([0], [0], color = color, lw=3, label=key))
    plt.legend(handles=legend_elements, loc="best")
    
    plt.savefig(outdir + "/sample_cluster_pca.png",dpi=300,bbox_inches="tight")

def get_hyper_z_trans_df_mean(hypervariable_df):
    column_sets = hypervariable_df.columns.tolist()
    mean_list = []
    for column in column_sets:
        column_2 = column[1]
        if column_2.startswith("mean"):
            mean_list.append(column_2)
    mean_hypervariable_df = hypervariable_df.loc(axis=1)[:,mean_list]
    mean_hypervariable_df.columns = mean_list
    mean_hypervariable_df.to_csv(outdir +"/mean_hypervariable_df.xlsx",sep="\t")
    
    return mean_hypervariable_df
    
    
def gene_cluster_by_mean_value(mean_hypervariable_df,rearrange_z_trans_df,filein_df,vmin_value = -6,vmax_value = 6,
                               cluster_number = 8,cutoff = 50):
    
    ## cluster
    map_color = "bwr"
    g1=sns.clustermap(mean_hypervariable_df,cmap=map_color,metric='correlation',col_cluster = False,
                 figsize=(3,12),vmin = vmin_value,vmax = vmax_value,xticklabels = mean_hypervariable_df.columns,
                 cbar_pos = (0,0.75,0.05,0.05))
    
    proteins_indexs=scipy.cluster.hierarchy.fcluster(g1.dendrogram_row.calculated_linkage,
                                                     t = cluster_number,criterion='maxclust')
    
    results=Counter(proteins_indexs)
    sub_protein_clusters={}
    
    for protein ,cluster in zip(mean_hypervariable_df.index,proteins_indexs):
        if cluster in sub_protein_clusters.keys():
            sub_protein_clusters[cluster].append(protein)
        else:
            sub_protein_clusters[cluster] = [protein]
            
    cluster_list = []
    combine_cluster = []
    for cluster in sub_protein_clusters.keys():
        if len(sub_protein_clusters[cluster]) > cutoff:
            cluster_list.append(cluster)
        else:
            combine_cluster.append(cluster)
    
    cluster_list.append(combine_cluster)
    gene_list = []
    for cluster in cluster_list:
        if isinstance(cluster,list):
            for small_cluster in cluster:
                gene_list += sub_protein_clusters[small_cluster]
        else:
            gene_list += sub_protein_clusters[cluster]
            
    ##cluster_heatmap
    plt.figure(figsize=(5,15))
    plt.box(False)
    c_map = "bwr"        
    vmin = -6
    vmax=6
    heat = plt.imshow(rearrange_z_trans_df.loc[gene_list],aspect="auto",cmap=c_map,vmin=vmin, vmax=vmax)
    a = -0.5    
    extend = 10
    x_value = len(rearrange_z_trans_df.columns) - 0.2    
    text_value = len(rearrange_z_trans_df.columns)
    
    i = 1
    for key in cluster_list:
        if isinstance(key,list):
            b=0
            for j in key:
                b+=len(sub_protein_clusters[j])
            plt.plot([x_value,x_value],[a+extend,a+b-extend],linewidth=5)
            plt.text(text_value,a+b/2,"cluster_0")
            a=a+b          
        else:
            b=len(sub_protein_clusters[key])
            plt.plot([x_value,x_value],[a+extend,a+len(sub_protein_clusters[key])-extend],
                     linewidth=5)
            plt.text(text_value,a+b/2,"cluster_%s"%(i))
            
            i +=1
            a = a+len(sub_protein_clusters[key])
    plt.xticks(range(0,len(rearrange_z_trans_df.columns)),rearrange_z_trans_df.columns,rotation=90)
    plt.yticks([])
    cbar = plt.colorbar(heat, shrink=0.5,orientation="horizontal")
    cbar.set_ticks([-6,-3,0,3,6],update_ticks=True)
    cbar.set_ticklabels([-6,-3,0,3,6])    
    plt.savefig(outdir + "/gene_cluster_heatmap.png",dpi=300,bbox_inches="tight")
    
   # mkdir(outdir+"/clusters")
    new_dir = outdir+"/clusters"
    mkdir(new_dir)
    ## cluster results
    i=1
    for cluster in cluster_list:
        
        if isinstance(cluster, list):
            sub_gene_list = []
            
            for small_cluster in cluster:
                sub_gene_list += sub_protein_clusters[small_cluster]
                
            cluster_df = filein_df.loc[sub_gene_list]
            cluster_df.to_csv(new_dir+ "/" +"cluster_0.xlsx",sep="\t")
            
        else:
            cluster_df = filein_df.loc[sub_protein_clusters[cluster]]
            cluster_df.to_csv(new_dir + "/" + "cluster_%s.xlsx"%(i),sep="\t")
            i+=1
    return sub_protein_clusters,cluster_list

from . import enrichment_protein
#import enrichment_protein
def pathway_background():
    kegg_pathway_file = "/home/shao/zMap/jobs/gmt_files/enrichr.KEGG_2016.gmt"
    go_bp_pathway_file = "/home/shao/zMap/jobs/gmt_files/human_GO_Biological_Process_2015.gmt"
    go_cc_pathway_file = "/home/shao/zMap/jobs/gmt_files/human_GO_Cellular_Component_2015.gmt"
    go_mf_pathway_file = "/home/shao/zMap/jobs/gmt_files/human_GO_Molecular_Function_2015.gmt" 
    kegg_pathway_dict = enrichment_protein.get_kegg_pathway_information(kegg_pathway_file,filein_df.index)
    go_pathway_bp_dict = enrichment_protein.get_go_pathway_information(go_bp_pathway_file,filein_df.index)
    go_pathway_cc_dict = enrichment_protein.get_go_pathway_information(go_cc_pathway_file,filein_df.index)
    go_pathway_mf_dict = enrichment_protein.get_go_pathway_information(go_mf_pathway_file,filein_df.index)
    return [kegg_pathway_dict,go_pathway_bp_dict,go_pathway_cc_dict,go_pathway_mf_dict]


def pathway_analysis(sub_protein_clusters,cluster_list,outdir,pathway_dict_list):
    
    enrich_dir = outdir + "/enrichment_results/"
    mkdir(enrich_dir)
    str_list = ["kegg","go_bp","go_cc","go_mf"]
    str_list1 = ["KEGG","GO_Biological_Process","GO_Cellular_Component","GO_Molecular_Function"]
    
    for cluster in cluster_list:
        if isinstance(cluster,list):
            pass
        else:
            sub_gene_list = sub_protein_clusters[cluster]
            i = 0
            for pathway_dict in pathway_dict_list:
                enrich_result = enrichment_protein.enrich(pathway_dict,sub_gene_list,filein_df.index)
                enrich_result_sort = enrich_result.sort_values(by=["padjust_by_BH"])
        
                #enrich_file = enrich_dir + "cluster%s_%s_enrich_result.pdf"%(cluster,str_list[i])
                enrich_file1 = enrich_dir + "cluster%s_%s_enrich_result.png"%(cluster,str_list[i])
        
                #enrichment_protein.enrich_plot(enrich_result_sort,title= str(cluster) + "_" + str_list1[i] ,pvalue=0.05,fig_file = enrich_file)
                enrichment_protein.enrich_plot(enrich_result_sort,title= str(cluster) + "_" + str_list1[i],pvalue=0.05,fig_file = enrich_file1)
        
                enrich_result.to_csv(enrich_dir+"cluster%s_%s.txt"%(cluster,str_list[i]),sep = "\t")
                i+=1



def top_gene(rearrange_z_trans_df,hypervariable_df,mean_hypervariable_df,top_gene = 100,cluster_num = 8):
    top_hypervariable_dir = outdir +"/top_hyper_%s"%(top_gene)
    mkdir(top_hypervariable_dir)
    hypervariable_df_sort =  hypervariable_df.sort_values([('-','Second best pvalue')])
    top_gene_list = hypervariable_df_sort.iloc[0:top_gene].index
    top_gene_df = rearrange_z_trans_df.loc[top_gene_list]
    top_mean_df = mean_hypervariable_df.loc[top_gene_list]
    top_mean_df.to_csv(top_hypervariable_dir+"/top_mean_df.xlsx",sep="\t")
    #sns.heatmap(top_gene_df)
    
    estimator = KMeans(n_clusters=cluster_num,random_state =222)
    estimator.fit(top_mean_df)
    label_pred = estimator.labels_
    #color=["lightcoral","darkorange","yellowgreen","lime","lightseagreen","grey","lightskyblue","blue","black"]
    #m=331
    #plt.figure(figsize=(4,15))
    protein = []
    for i in range(0,cluster_num):
        cluster_df = top_mean_df[label_pred == i]
        if len(cluster_df)>2:
            sort_rank_dict = defaultdict(list)
            for gene in cluster_df.index:
                expression_value = list(cluster_df.loc[gene].values)
                sorted_expression_value = sorted(expression_value)
                str_rank=""
                for j in sorted_expression_value[-2:]:
                    str_rank += str(expression_value.index(j))
                sort_rank_dict[str_rank].append(gene)
            for key in sort_rank_dict.keys():
                protein +=sort_rank_dict[key]
        else:                    
            protein += list(top_mean_df[label_pred == i].index)
    
    new_top_df = top_gene_df.loc[protein]
    new_top_df.to_csv(top_hypervariable_dir + "/top_hyper.xlsx",sep="\t")
    plt.figure(figsize=(4,15))      
    sns.heatmap(new_top_df,cmap="bwr",vmin=-10,vmax=10,square=True,cbar_kws = {"shrink":0.2})
    plt.savefig(top_hypervariable_dir+"/top_hyper_heatmap.png",dpi=300,bbox_inches="tight")    
    
   
def mkdir(path_dir):
    isexists = os.path.exists(path_dir)
    if not isexists:
        os.makedirs(path_dir)
    else:
        pass    
    
        
    
    
if __name__=="__main__":
    parser=argparse.ArgumentParser(description="")
    parser.add_argument("--filein",help="step1 output results",required=True)
    parser.add_argument("--outdir",help="outdir",required=True)
    parser.add_argument("--pvalue",help="BH-corrected pvalue cutoff used to define hypervariable proteins,default is 0.01",required=True)
    parser.add_argument("--cluster_number_for_hypervariable",help="cluster number of hypervariable proteins defined by user",required=True)
    parser.add_argument("--minclustersize",help="minimum protein number for each cluster",required=True)
    parser.add_argument("--top",help="number of top hypervariable proteins",required=True)
    parser.add_argument("--cluster_number_for_top_proteins",help="cluster number for top hypervariable proteins",required=True)
    argv = vars(parser.parse_args())
    filein = argv['filein']
    outdir = argv['outdir']
    mkdir(outdir)
    pvalue = float(argv['pvalue'])
    cluster_num = int(argv['cluster_number_for_hypervariable'])
    minclustersize = int(argv['minclustersize'])
    top = int(argv['top'])
    cluster_number2 = int(argv["cluster_number_for_top_proteins"])
  #  filein = r"H:\project\protein\web_server\zmap\all_batch_final_results.xls"
    
    filein_df = read_file(filein) 
    hypervariable_df = get_hypervariable_df(filein_df,pvalue_cutoff=pvalue)
    z_trans_df = get_z_trans_df(hypervariable_df)
    recolumn_z_trans_df = rename_columns(z_trans_df)
    rearrange_z_trans_df = rearrange_column(recolumn_z_trans_df)
    pearsonr_array = pcc(recolumn_z_trans_df,outdir)
    
    rep_color_dict,rep_marker_dict,condition_color_dict,df_T = get_color(pearsonr_array.columns)
    
    hierarchical_clustering(pearsonr_array,df_T,outdir)      
    
    pca,pc_df = pca_analysis(pearsonr_array,outdir)
    
    pca_scatter(pc_df,pca,rep_marker_dict,condition_color_dict,outdir)
    
    mean_hypervariable_df  = get_hyper_z_trans_df_mean(hypervariable_df)
    
    sub_protein_clusters,cluster_list = gene_cluster_by_mean_value(mean_hypervariable_df,rearrange_z_trans_df,filein_df,vmin_value = -6,vmax_value = 6,
                               cluster_number = cluster_num,cutoff = minclustersize)
    
    pathway_dict_list = pathway_background()
    
    pathway_analysis(sub_protein_clusters,cluster_list,outdir,pathway_dict_list)
    
    top_gene(rearrange_z_trans_df,hypervariable_df,mean_hypervariable_df,top_gene = top,cluster_num = cluster_number2)
    
    
    
    
    
    
    
    
    
    
    
    












