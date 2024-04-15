from app import celery
import os,zipfile,json,random
from flask import current_app
import numpy as np

from celery.task.control import revoke
from .zmap_step2_downstream_analysis import *


def zip_ya(start_dir):
    start_dir = start_dir #compress directory path
    file_news = start_dir + '.zip'

    z = zipfile.ZipFile(file_news, 'w', zipfile.ZIP_DEFLATED)
    for dir_path, dir_names, file_names in os.walk(start_dir):
        f_path = dir_path.replace(start_dir, '')  # from current directory start to copy
        f_path = f_path and f_path + os.sep or '' # 
        for filename in file_names:
            z.write(os.path.join(dir_path, filename), f_path + filename)
    z.close()
    return file_news


def gene_cluster(intensity_df,sample_order,cluster_num,cutoff,vmin_value,vmax_value,outdir):
    map_color = "bwr"
    intensity_df = intensity_df[sample_order]
    g1=sns.clustermap(intensity_df,cmap=map_color,metric='correlation',col_cluster = False,
                 figsize=(3,12),vmin = vmin_value,vmax = vmax_value,xticklabels = intensity_df.columns,
                 cbar_pos = (0,0.75,0.05,0.05))
    proteins_indexs=scipy.cluster.hierarchy.fcluster(g1.dendrogram_row.calculated_linkage,
                                                     t = cluster_num,criterion='maxclust')
    results=Counter(proteins_indexs)
    sub_protein_clusters={}
    
    for protein ,cluster in zip(intensity_df.index,proteins_indexs):
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
    heat = plt.imshow(intensity_df.loc[gene_list],aspect="auto",cmap=c_map,vmin=vmin, vmax=vmax)
    a = -0.5    
    extend = 10
    x_value = len(intensity_df.columns) - 0.2    
    text_value = len(intensity_df.columns)
    
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
    plt.xticks(range(0,len(intensity_df.columns)),intensity_df.columns,rotation=90)
    plt.yticks([])
    cbar = plt.colorbar(heat, shrink=0.5,orientation="horizontal")
    cbar.set_ticks([-6,-3,0,3,6],update_ticks=True)
    cbar.set_ticklabels([-6,-3,0,3,6])    
    plt.savefig(outdir + "/gene_cluster_heatmap.png",dpi=300,bbox_inches="tight")
    
   # mkdir(outdir+"/clusters")
    new_dir = outdir+"/clusters"
    os.mkdir(new_dir)
    ## cluster results
    i=1
    for cluster in cluster_list:
        
        if isinstance(cluster, list):
            sub_gene_list = []
            
            for small_cluster in cluster:
                sub_gene_list += sub_protein_clusters[small_cluster]
                
            cluster_df = intensity_df.loc[sub_gene_list]
            cluster_df.to_csv(new_dir+ "/" +"cluster_0.xlsx",sep="\t")
            
        else:
            cluster_df = intensity_df.loc[sub_protein_clusters[cluster]]
            cluster_df.to_csv(new_dir + "/" + "cluster_%s.xlsx"%(i),sep="\t")
            i+=1
    return sub_protein_clusters,cluster_list


@celery.task(name='async_cluster',bind=True)
def async_cluster(self,filepath,cluster_num,task_type):
    if task_type == 'upload':
        self.update_state(state='PROGRESS',meta={'status':'start pca processing'})
        raw_intensity_data = pd.ExcelFile(filepath)
        sheet_names = raw_intensity_data.sheet_names
        intensity_df = raw_intensity_data.parse(sheet_names[0],index_col=0)
        condition_df = raw_intensity_data.parse(sheet_names[1],index_col=0)
        if len(set(intensity_df.index))!=len(intensity_df.index):
            raise Exception("Gene IDs were not uniq")
        elif len([i for i in list(intensity_df.index) if isinstance(i,float)]) > 0:
            raise Exception("Nans are included in the gene ID")
        elif len(set(intensity_df.columns))!=len(intensity_df.columns):
            raise Exception("Sample IDs were not uniq")

        elif len([i for i in list(intensity_df.index) if str(i).startswith('Unnamed')]) > 0:
            raise Exception("Nans are included in the sample ID")

        elif set(intensity_df.columns) != set(condition_df.index):
            raise Exception("Sample list not exactly the same")
        elif condition_df.isna().any().any():
            raise Exception("Nans are included in condition info.")
        else:
            self.update_state(state='PROGRESS',meta={'status':'calculate sample correlation'})
            os.chdir(os.path.dirname(filepath))
            os.mkdir(os.path.join(os.getcwd(),"cluster_num_" + str(cluster_num)))
            res_dir = os.path.join(os.getcwd(),"cluster_num_" + str(cluster_num))
            sample_pcc = pcc(intensity_df,res_dir)
            self.update_state(state='PROGRESS',meta={'status':'start cluster'})
            condition_color_list = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(len(set(condition_df[condition_df.columns[0]])))]
            rep_color_list = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(len(set(condition_df[condition_df.columns[1]])))]
            condition_color_dict = dict(zip(list(set(condition_df[condition_df.columns[0]])),condition_color_list))
            rep_color_dict = dict(zip(list(set(condition_df[condition_df.columns[1]])),rep_color_list))
            
            condition_color = map(condition_color_dict.get,list(condition_df[condition_df.columns[0]]))
            rep_color = map(rep_color_dict.get,list(condition_df[condition_df.columns[1]]))
            df = pd.DataFrame([condition_color,rep_color],index=["condition","batch"],columns=condition_df.index)
            df = df.T
            sample_order = hierarchical_clustering(sample_pcc,df,res_dir)
            sub_protein_clusters,cluster_list = gene_cluster(intensity_df,sample_order,cluster_num,50,-6,6,res_dir)
            #condition order
            zip_ya(os.path.dirname(filepath))
            return {'status': 'cluster finished'}                   

