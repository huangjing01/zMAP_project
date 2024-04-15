# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 11:06:48 2020

@author: guixiuqi
"""

import numpy as np
import pandas as pd

import random
from scipy.stats import spearmanr
from scipy.stats import pearsonr

from scipy.stats import mannwhitneyu
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib.patches import Circle
from matplotlib.patheffects import withStroke


def get_batch_color_dict(batch_df):
    batch_color_dict = {}
    color_list = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])for i in range(len(batch_df))]
    i=0
    for index in batch_df.index:
        for sample in batch_df.loc[index]:
            batch_color_dict[sample] = color_list[i]
        i+=1
    return batch_color_dict




def get_color_list(columns,batch_color_dict):
    condition_color_list = []
    batch_color_list = []
    for i in columns:
        if i.startswith("T"):
            condition_color_list.append("lightcoral")
        else:
            condition_color_list.append("lightseagreen")
            
        batch_color_list.append(batch_color_dict[i])
        
    color_df = pd.DataFrame([batch_color_list,condition_color_list],index=["batch","condition"],columns=columns)
    return color_df.T



#####################################sample_cluster##################################
def correlation(df,method="pearsonr"):
    ### method can be pearsonr or spearmanr
    correlation_array=[]
    for i in df.columns:
        cc_list=[]
        for j in df.columns:
            if i==j:
                cc_list.append(1)
            else:
                two_columns_df = df.loc[:,[i,j]]
                detected_index =two_columns_df.index[(1-(np.isnan(two_columns_df[i]) | np.isnan(two_columns_df[j])).values).astype(np.bool)]
              #  print (i +"\t"+j+"\t"+"Number of genes detected simultaneously:"+str(len(detected_index)))
                new_two_columns_df = two_columns_df.loc[detected_index,:]
                
                if method =="spearmanr":
                    corr = spearmanr(new_two_columns_df.loc[:,i],new_two_columns_df.loc[:,j])[0]
                else:
                    corr = pearsonr(new_two_columns_df.loc[:,i],new_two_columns_df.loc[:,j])[0]
                cc_list.append(corr)
        correlation_array.append(cc_list)
    correlation_array = pd.DataFrame(correlation_array,index=df.columns,columns=df.columns)
    return correlation_array




###cluster_map##
def cluster_map(df,pearsonr_cor,color_list_df,vmin=-0.5,vmax=1,fig = ""):   
    tumor_patch=mpatches.Patch(color='lightcoral',label='Tumor')
    normal_patch=mpatches.Patch(color='lightseagreen',label='Normal')
    g = sns.clustermap(pearsonr_cor,cmap="bwr",col_colors = color_list_df,
                       vmin=vmin,vmax=vmax,cbar_kws={"ticks":[-0.5,0,0.5,1],'label': 'PCC'},
                       xticklabels=False,yticklabels=False)
    g.ax_row_dendrogram.set_visible(False)
    legends=[tumor_patch,normal_patch]
    l2=g.ax_heatmap.legend(loc='center left',bbox_to_anchor=(0.7,1.25),handles=legends,frameon=True)
    l2.set_title(title='Tissue type',prop={'size':15})
    plt.setp(g.ax_heatmap.get_legend().get_texts(),fontsize='15')
    plt.savefig(fig+".pdf",dpi=300,bbox_inches="tight")
    plt.savefig(fig+".png",dpi=300,bbox_inches="tight")
    
#batch_color_dict = get_batch_color_dict(batch_df)
#color_list_df = get_color_list(final_df.columns)
#pearsonr_cor = correlation(final_df,method="pearsonr")
       
#outdir=r"H:\project\protein\web_server\fig_for_reversezmap_workflow\\"
#cluster_map(final_df,pearsonr_cor,vmin=-0.5,vmax=1,fig=outdir+"clustermap")
