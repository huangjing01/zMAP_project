from __future__ import with_statement
import signal
from contextlib import contextmanager

class TimeoutException(Exception): pass

@contextmanager
def time_limit(seconds):
    def signal_handler(signum, frame):
        raise TimeoutException
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)

from app import celery
import os
from flask import current_app
from celery.task.control import revoke
import zipfile
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

import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
#import scipy.stats as stats
#import seaborn as sns
from scipy.stats import spearmanr
from scipy.stats import pearsonr
#from matplotlib_venn import venn3
#import scipy.stats as stats
#from scipy.stats import chi2_contingency
#import scipy
import matplotlib
import argparse
import os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import scipy.stats as stats
import seaborn as sns
#matplotlib.rcParams['font.family'] ='Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import seaborn as sns
from . import enrichment_protein
import collections
import scipy
import matplotlib as mpl
import matplotlib.cm as cm
from sklearn import preprocessing
from scipy.stats import chi2_contingency
import statsmodels.api as sm
from statsmodels.formula.api import ols


def get_color_list(cmap=cm.bwr,value_list=[],vmin=-3,vmax=8):
    pc1_norm = mpl.colors.Normalize(vmin = vmin,vmax = vmax)
    pc1_cmap = cmap
    pc1_m=cm.ScalarMappable(norm=pc1_norm,cmap=pc1_cmap)
    pc1_color_list = []
    for pc in value_list:
        pc1_color_list.append(pc1_m.to_rgba(pc))
    return pc1_m,pc1_color_list



def data_for_feature(df):

    subgroup_list = set(df.iloc[:,1].values)
    color_list = list(set(df.iloc[:,0].values))
    all_table_list = []
    for subgroup in subgroup_list:
        new_df = df.loc[df.iloc[:,1]==subgroup]
        sub_color_list = list(new_df.iloc[:,0].values)
        sub_color_dict = collections.Counter(sub_color_list)
        table_list = []
        for color in color_list:
            table_list.append(sub_color_dict[color])
        all_table_list.append(table_list)

    return all_table_list


def chi_pvalue(discrete,continuous,cluster_df,clinical_df):
    
    pvalue_list = []
    clinical_df = clinical_df.loc[cluster_df.index]
    add_feature_cluster_df = pd.concat([clinical_df,cluster_df],axis=1)

    for column in discrete:
        df = add_feature_cluster_df[[column,"x"]]
        df = df.dropna(how='any')
        a = data_for_feature(df)
        g, pvalue, dof, expctd = chi2_contingency(a, lambda_="log-likelihood")
        pvalue_list.append(pvalue)
        
        
    for column in continuous:
        
        
        
        a='%s ~ C(x)'%(column)
        model = ols(a, data=add_feature_cluster_df).fit()
        anova_table = sm.stats.anova_lm(model,typ=2)
        pc1_anova_pvalue = anova_table.loc['C(x)','PR(>F)']
        
        pvalue_list.append(pc1_anova_pvalue)
      
        
    return pvalue_list
        
    
    

def cluster_map_sample_type(cluster_df,ax,clinical_df,discrete_color_dict,continuous_color_dict,discrete,continuous):  

    clinical_df = clinical_df.loc[cluster_df.index]

    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    
    hight=-4
    width=1
    hight_interval = -1
    bottom=-1
    y_position_list = []    
    cluster_color_dict = {1:"tab:green",2:"tab:blue",3:"tab:red",4:"tab:orange",5:'tab:purple',6:'tab:brown',7:'tab:pink',8:'tab:gray',9:'tab:olive',10:'tab:cyan'}
    sample_cluster_color_list = []
    for index in cluster_df.index:
        color = cluster_color_dict[cluster_df.loc[index,"x"]]
        sample_cluster_color_list.append(color)   
#    print cluster_df
    #cluster_bar
    ax.bar(range(0,len(cluster_df)),[hight]*len(cluster_df),width=1,align="edge",
         bottom =[bottom]*len(cluster_df),color = sample_cluster_color_list)
    pos= bottom+0.5*hight
    bottom =bottom+(hight+hight_interval)
    y_position_list.append(pos)   

    for feature in discrete:
        feature_list = clinical_df.loc[cluster_df.index,feature]

        ax.bar(range(0,len(cluster_df)),[hight]*len(cluster_df),width=1,align="edge",
            bottom =[bottom]*len(cluster_df),color = [discrete_color_dict[i] for i in feature_list])
        pos= bottom+0.5*hight
        bottom =bottom+(hight+hight_interval)

        y_position_list.append(pos)
    
    for feature in continuous:
           
        min_ = min(clinical_df[feature].values)
        max_ = max(clinical_df[feature].values)

        pc1_m,pc1_color_list = get_color_list(cmap=continuous_color_dict[feature],value_list=clinical_df[feature].values,vmin=min_,vmax=max_)  

        ax.bar(range(0,len(clinical_df)),[hight]*len(clinical_df),width=1,
         bottom =[bottom]*len(clinical_df),color = pc1_color_list,align='edge')
        pos= bottom+0.5*hight
        bottom =bottom+(hight+hight_interval)
        y_position_list.append(pos)
  
    pvalue_list = chi_pvalue(discrete,continuous,cluster_df,clinical_df)
#    print(pvalue_list)

    yticks_list = discrete+continuous

    for i in range(0,len(yticks_list)):
        ax.text(-1,y_position_list[i+1],yticks_list[i],horizontalalignment="right",verticalalignment="center",size=8)    
    #    if i < len(yticks_list)-1:        
        ax.text(len(clinical_df),y_position_list[i+1],"%.3e"%(pvalue_list[i]),horizontalalignment="left",verticalalignment="center",size=8) 

def pcc_heatmap(sub_protein_clusters,cluster_df,final_imputer_df,ax,vmin_=-0.5,vmax_=1):
    cluster_list = sorted(list(sub_protein_clusters.keys()))
    all_protein_list = []
    for cluster in cluster_list:
        
        sub_gene_list = sub_protein_clusters[cluster]
        all_protein_list += sub_gene_list
        
    
    sns.heatmap(final_imputer_df.loc[all_protein_list,list(cluster_df.index)],vmin=vmin_,vmax=vmax_,cmap='bwr',cbar=False,yticklabels=False,xticklabels=False)
 

def other_square(ax,discrete,discrete_color_dict,cluster_df,clinical_df):
     
    cluster_color_dict = {1:"tab:green",2:"tab:blue",3:"tab:red",4:"tab:orange",5:'tab:purple',6:'tab:brown',7:'tab:pink',8:'tab:gray',9:'tab:olive',10:'tab:cyan'}
    
    our_cluster_color_dict = {}
    cluster_l = sorted(list(set(cluster_df['x'].values)))
    for cluster in cluster_l:
        our_cluster_color_dict[cluster]=cluster_color_dict[cluster]


    feature_color_dict_list = [our_cluster_color_dict]

    for feature in discrete:
        dict_ = {}
        value_ = list(set(clinical_df[feature].values))
        for v_ in value_:
     #       print(v_)
            if pd.isna(v_):
                pass
            else:
                dict_[v_] = discrete_color_dict[v_]

        feature_color_dict_list.append(dict_)

        
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    
    
    i=0
    start=0
    interval = 1
    height=-4
    bar_width=0.1
    bottom =-1
    
    text_size=8
    
    i=0

    m_type_list = []

    for dict_ in feature_color_dict_list:

        for m_type in dict_.keys():
            if m_type in m_type_list:
                pass
            else:
                ax.bar([i*interval+start],[height],bottom=bottom,width=bar_width,color=dict_[m_type],align='edge')
                ax.text(i*interval+start+bar_width,bottom+height/2,m_type,size=text_size,verticalalignment="center")
                bottom = bottom+(-interval+height)
        
        bottom = bottom+(-interval+height)
        
    
def drew_bar_continuous(ax,feature,continuous_color_dict,clinical_df):

        
    min_ = min(clinical_df[feature].values)
    max_ = max(clinical_df[feature].values)

    
    pc1_norm = mpl.colors.Normalize(vmin = min_,vmax = max_)
    cmap = cm.get_cmap(continuous_color_dict[feature])

    cb1=mpl.colorbar.ColorbarBase(ax,cmap=cmap,norm=pc1_norm,ticks=[min_,max_],orientation='horizontal')

    ax.set_title(feature,size=6)
    cb1.ax.set_xticklabels(['Min','Max'],rotation=0,fontdict={'fontsize':5})

def drew_bar_z(ax,vmin = -3,vmax = 3):
    
    pc1_norm = mpl.colors.Normalize(vmin = vmin,vmax = vmax)
    cb1=mpl.colorbar.ColorbarBase(ax,cmap=cm.bwr,norm=pc1_norm,ticks=[vmin,vmax],orientation='horizontal')
   
    ax.set_title("z-statistic",size=6)
    cb1.ax.set_xticklabels([str(vmin),str(vmax)],rotation=0,fontdict={'fontsize':5})

        
def mkdir(path_dir):
    isexists = os.path.exists(path_dir)
    if not isexists:
        os.makedirs(path_dir)
    else:
        pass

def knn_imputation(z_df):
    import numpy as np
    from sklearn.impute import KNNImputer
    imputer = KNNImputer(n_neighbors=5)
    imputer_df = imputer.fit_transform(z_df)
    final_imputer_df = pd.DataFrame(imputer_df,index = z_df.index,columns=z_df.columns)
    return final_imputer_df

def annova_test(final_imputer_df,cluster_df):
    pvalue_list = []
    protein_list = []
    for index in final_imputer_df.index:
        protein_list.append(index)
        gene_df = final_imputer_df.loc[index].to_frame()
        concat_df = pd.concat([gene_df,cluster_df],axis=1)
        concat_df = concat_df.dropna(how="any")
        concat_df.columns = ['protein','group']
        #print(concat_df)
        
        a='protein ~ C(group)'
        
        model = ols(a, data= concat_df).fit()
        anova_table = sm.stats.anova_lm(model,typ=2)

        anova_pvalue = anova_table.loc['C(group)','PR(>F)']
        
        pvalue_list.append(anova_pvalue)
        
    final_df = pd.DataFrame([protein_list,pvalue_list]).T
       
    final_df.columns = ["protein","anova_pvalue"]
    final_df["anova_fdr"] =[i if i<1 else 1 for i in final_df["anova_pvalue"]/final_df["anova_pvalue"].rank(axis=0)*len(final_df)]
    return final_df

def pathway_background(filein_df):
    kegg_pathway_file = "/home/shao/zMap/app/templates/input_file/reverse_zMAP/gmt/enrichr.KEGG_2016.gmt"
    go_bp_pathway_file = "/home/shao/zMap/app/templates/input_file/reverse_zMAP/gmt/human_GO_Biological_Process_2015.gmt"
    go_cc_pathway_file = "/home/shao/zMap/app/templates/input_file/reverse_zMAP/gmt/human_GO_Cellular_Component_2015.gmt"
    go_mf_pathway_file = "/home/shao/zMap/app/templates/input_file/reverse_zMAP/gmt/human_GO_Molecular_Function_2015.gmt"
    kegg_pathway_dict = enrichment_protein.get_kegg_pathway_information(kegg_pathway_file,filein_df.index)
    go_pathway_bp_dict = enrichment_protein.get_go_pathway_information(go_bp_pathway_file,filein_df.index)
    go_pathway_cc_dict = enrichment_protein.get_go_pathway_information(go_cc_pathway_file,filein_df.index)
    go_pathway_mf_dict = enrichment_protein.get_go_pathway_information(go_mf_pathway_file,filein_df.index)
    return [kegg_pathway_dict,go_pathway_bp_dict,go_pathway_cc_dict,go_pathway_mf_dict]

def clusterbar(ax,sub_protein_clusters):
    bottom=0
    for key in range(1,len(sub_protein_clusters)+1):
        cluster_gene = sub_protein_clusters[key]
        ax.bar(0,-len(cluster_gene),bottom=bottom,width=1,align="edge")
        bottom -= len(cluster_gene)
    ax.set_xlim([0,3])
    ax.set_ylim([bottom,0])
    ax.set_xticks([])
    ax.set_yticks([])
    
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)


def dep_cluster_enrich(final_imputer_df,anova_df,cluster_df,outdir,pathway_dict_list,fdr=0.05,cluster_n=4):

    #pathway_dict_list = pathway_background(filein_df) 
    sig_df = anova_df.loc[anova_df['anova_fdr']<fdr]
    
    sig_z_df_raw = final_imputer_df.loc[list(sig_df["protein"])]
    
    sig_z_df = sig_z_df_raw.sub(sig_z_df_raw.median(axis=1), axis=0)
    
    
    
    g = sns.clustermap(sig_z_df,cmap="bwr",col_cluster=False,row_cluster=True,vmin=-2,vmax=2,metric="correlation")
    
    
    proteins_indexs=scipy.cluster.hierarchy.fcluster(g.dendrogram_row.calculated_linkage,t=cluster_n,criterion='maxclust')
    
    
    sub_protein_clusters={}
    for protein ,cluster in zip(sig_z_df.index,proteins_indexs):
        if cluster in sub_protein_clusters.keys():
            sub_protein_clusters[cluster].append(protein)
        else:
            sub_protein_clusters[cluster] = [protein]

    enrich_dir = outdir + "/enrichment_results/"
    mkdir(enrich_dir)
    str_list = ["kegg","go_bp","go_cc","go_mf"]
    str_list1 = ["KEGG","GO_Biological_Process","GO_Cellular_Component","GO_Molecular_Function"]
    
    cluster_list = sorted(list(sub_protein_clusters.keys()))
    for cluster in cluster_list:
        
        sub_gene_list = sub_protein_clusters[cluster]
        df_ = pd.DataFrame(sub_gene_list)
        df_.to_csv(enrich_dir + "cluster%s_protein_list.csv"%cluster)
        i = 0
        for pathway_dict in pathway_dict_list:
            enrich_result = enrichment_protein.enrich(pathway_dict,sub_gene_list,final_imputer_df.index)
            enrich_result_sort = enrich_result.sort_values(by=["padjust_by_BH"])
     
            enrich_file = enrich_dir + "cluster%s_%s_enrich_result.pdf"%(cluster,str_list[i])
            enrich_file1 = enrich_dir + "cluster%s_%s_enrich_result.png"%(cluster,str_list[i])
        
            enrichment_protein.enrich_plot(enrich_result_sort,title= str(cluster) + "_" + str_list1[i] ,pvalue=0.05,fig_file = enrich_file)
            enrichment_protein.enrich_plot(enrich_result_sort,title= str(cluster) + "_" + str_list1[i],pvalue=0.05,fig_file = enrich_file1)
        
            enrich_result.to_csv(enrich_dir+"cluster%s_%s.txt"%(cluster,str_list[i]),sep = "\t")
            i+=1

    return sub_protein_clusters,sig_z_df
   
def prepare_feature_dict(color_f):

    feature_color_dict = {}
    df_ = pd.read_csv(color_f,sep="\t",header=None,index_col=0)
#    print (df_)
    for index in df_.index:
        value_ = df_.loc[index,1]
#        print(index,value_)
        feature_color_dict[index] = value_
    return feature_color_dict
    



"""
if __name__=="__main__":

    parser=argparse.ArgumentParser(description="")
    parser.add_argument("--z_statistic_matrix",help="z-transformed log2-intensity",required=True)
    parser.add_argument("--cluster_f",help="cluster",required=True)

    parser.add_argument("--clinical_info",help="sample_information",required=True)
    
    parser.add_argument("--discrete",help="discrete clinical feature",required=True)
    parser.add_argument("--continuous",help="continuous clinical feature",required=True)

    parser.add_argument("--color_f",help="color for discrete clinical feature",required=True)
    parser.add_argument("--colorbar_f",help="colorbar for continuous clinical feature",required=True)
 
    parser.add_argument("--outdir",help="outdir for result",required=True)
    parser.add_argument("--fdr",help="FDR cutoff for anova",type=float, default=0.05)
    parser.add_argument("--cluster_n",help="number of signiture",required=True)

    argv = vars(parser.parse_args())
    z_statistic_matrix = argv['z_statistic_matrix']
    cluster = argv['cluster_f']
    clinical_info = argv['clinical_info']
    discrete = argv['discrete']
    continuous = argv['continuous']
    discrete = discrete.strip().split(',')
    continuous = continuous.strip().split(',')

    color_f=argv['color_f']
    colorbar_f=argv['colorbar_f']
    
    outdir = argv['outdir']

    fdr = float(argv['fdr'])
    cluster_n = int(argv['cluster_n'])
"""





@celery.task(name='association_with_clinical_feature_',bind=True)
def association_with_clinical_feature_(self,z_statistic_matrix,cluster,clinical_info,discrete,continuous,color_f,colorbar_f,outdir,fdr,cluster_n):
    try:
        with time_limit(60*60*6):
            status = ["start processing"]
            self.update_state(state='PROGRESS',meta={'status': status})
            z_df = pd.read_csv(z_statistic_matrix,sep="\t",index_col=0)
            sample_n = int(len(z_df.columns)/2)
            z_df  = z_df.dropna(thresh=sample_n)
            final_imputer_df = knn_imputation(z_df)

            cluster_df = pd.read_csv(cluster,sep=",",index_col=0)
            cluster_df.columns=['x']
            clinical_df = pd.read_csv(clinical_info,sep="\t",index_col=0)
            
            overlap_sample = list(set(final_imputer_df.columns)&set(cluster_df.index)&set(clinical_df.index))

            status.append("%s common samples in input files were included in subsequant analysis."%(len(overlap_sample)))
            self.update_state(state='PROGRESS',meta={'status': status})

            final_imputer_df = final_imputer_df[overlap_sample]
            cluster_df = cluster_df.loc[overlap_sample]
            clinical_df = clinical_df.loc[overlap_sample]
         
            
            cluster_df = cluster_df.sort_values(by='x')

            status.append("Identifying differentially expressed proteins between groups.")
            self.update_state(state='PROGRESS',meta={'status': status})

            anova_df = annova_test(final_imputer_df,cluster_df)
            discrete_color_dict = prepare_feature_dict(color_f)
            discrete_color_dict[np.nan] = 'white'
            continuous_color_dict = prepare_feature_dict(colorbar_f)
            pathway_dict_list = pathway_background(z_df)
            
            status.append("performing clustering on differentially expressed proteins(DEPs) to recognize expression signatures. DEPs will be clustered into %s clusters"%(str(cluster_n)))
            status.append("performing pathway enrichment analysis on protein sets coressponging each cluster.")
            self.update_state(state='PROGRESS',meta={'status': status})

            sub_protein_clusters,sig_z_df = dep_cluster_enrich(final_imputer_df,anova_df,cluster_df,outdir,pathway_dict_list,fdr=fdr,cluster_n=cluster_n)   
            vmin=-2
            vmax=2
            method='ward'
            
            status.append("Generating the result picture.")
            self.update_state(state='PROGRESS',meta={'status': status})
            import argparse
            import matplotlib.gridspec as gridspec
            from matplotlib.gridspec import GridSpec
            
            fig = plt.figure(figsize=(10,10),dpi=300)
            gs = GridSpec(60,30,figure=fig)
            
            
            ax2 = fig.add_subplot(gs[0:20,8:25])
            
            cluster_map_sample_type(cluster_df,ax2,clinical_df,discrete_color_dict,continuous_color_dict,discrete,continuous)

            
            ax3 = fig.add_subplot(gs[20:60,8:25])
            vmin=-2
            vmax=2
            
            pcc_heatmap(sub_protein_clusters,cluster_df,sig_z_df,ax3,vmin_=vmin,vmax_=vmax)
            
            
            ax5 = fig.add_subplot(gs[0:30,0:1])
            other_square(ax5,discrete,discrete_color_dict,cluster_df,clinical_df)

            ax9 = fig.add_subplot(gs[20:60,25:26])
            clusterbar(ax9,sub_protein_clusters)

             
           
            x =0
            y=2
            
            m = 32
            for feature in continuous:
                n = m+1
                ax6 = fig.add_subplot(gs[m:n,x:y])
                drew_bar_continuous(ax6,feature,continuous_color_dict,clinical_df)
                m+=4
            
            

            ax8 = fig.add_subplot(gs[58:59,x:y])
            drew_bar_z(ax8,vmin=vmin,vmax=vmax)
            
            plt.savefig(outdir+"/sample_clustering_association_with_clinical_and_molecule_feature.pdf",dpi=300,bbox_inches="tight")

            zip_ya(outdir)    
            status.append("All finished.")
            self.update_state(state='PROGRESS',meta={'status': status})

            return {'status': status}
    except TimeoutException:
        status.append("Error : Excution time out of limit. please check your file and feel free to contact us if you have any question.")
        return {'status':status}
    except Exception as e:
        status.append("Error : " + str(e) + ". please check your file and feel free to contact us if you have any question.")
        return {'status':status}





