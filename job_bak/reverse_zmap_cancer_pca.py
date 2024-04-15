# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 13:29:02 2020

@author: guixiuqi
"""

import numpy as np
import pandas as pd
import os
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
import matplotlib.gridspec as gridspec


def pca_analysis(pearsonr_cor):
    from sklearn.decomposition import PCA
    X_new =  pearsonr_cor.T
    pca = PCA(n_components=2)
    X_r = pca.fit(X_new).transform(X_new)
    pc_df = pd.DataFrame(X_r,index = X_new.index,columns = ["PC1","PC2"])
    return pca,pc_df




def pc_scatter(pc_df,pca,batch_color_dict,outdir):
    plt.rcParams["figure.subplot.right"] = 0.88
    plt.rcParams["figure.subplot.top"] = 0.88
    batch_color_list = []
    condition_shape_list =[]
    for index in batch_color_dict.keys():       
        if index in pc_df.index:
            if "T" in index:
                condition_shape_list.append("^")
                batch_color_list.append(batch_color_dict[index])
            else:
                condition_shape_list.append("o")
                batch_color_list.append(batch_color_dict[index])

    
    pc_number = len(pc_df.columns)
    for i in range(pc_number-1):
        for j in range(i+1,pc_number):

            x_pc = pc_df.columns[i]
            y_pc = pc_df.columns[j]
            batch_number = 0
            condition_number = 0

            fig = plt.figure(figsize=(9,6))
            gs=gridspec.GridSpec(nrows=60,ncols= 90)
            ax = fig.add_subplot(gs[0:60,0:60])
            
            
            m=1
            batch_g_list = []
            color_list=[]
            for index in batch_color_dict.keys():
                if index in pc_df.index:
                    if batch_color_list[batch_number] in color_list:
                        ax.scatter(pc_df.loc[index,x_pc],pc_df.loc[index,y_pc],edgecolors=batch_color_list[batch_number],
                            facecolors = "white",linewidths=2,
                            marker=condition_shape_list[condition_number])
                    else:
                        if index.startswith("N"):
                            g1=ax.scatter(pc_df.loc[index,x_pc],pc_df.loc[index,y_pc],edgecolors="black",
                            facecolors = "white",linewidths=2,
                            marker=condition_shape_list[condition_number])
                            
                            
                            g=ax.scatter(pc_df.loc[index,x_pc],pc_df.loc[index,y_pc],edgecolors=batch_color_list[batch_number],
                            facecolors = "white",linewidths=2,
                            marker=condition_shape_list[condition_number],label= str(m))
                            m+=1
                            color_list.append(batch_color_list[batch_number])
                            batch_g_list.append(g)
                        else:
                            
                            g2=ax.scatter(pc_df.loc[index,x_pc],pc_df.loc[index,y_pc],edgecolors="black",
                            facecolors = "white",linewidths=2,
                            marker=condition_shape_list[condition_number])
                            ax.scatter(pc_df.loc[index,x_pc],pc_df.loc[index,y_pc],edgecolors=batch_color_list[batch_number],
                            facecolors = "white",linewidths=2,
                            marker=condition_shape_list[condition_number])  
                    batch_number+=1
                    condition_number+=1
                else:
                    pass
            
            ax.set_xlabel(x_pc+"(%s%%)"%("%.1f"%(pca.explained_variance_ratio_[i]*100)),size=15)
            ax.set_ylabel(y_pc+"(%s%%)"%("%.1f"%(pca.explained_variance_ratio_[j]*100)),size=15)
            ax.legend(handles=batch_g_list,labels=range(1,m+1),ncol=2, title="Batch", title_fontsize=15,loc=1, bbox_to_anchor=(1.38,0.82))
            fig.legend(handles=[g2,g1],labels=["Tumor","Normal"],bbox_to_anchor=(0.72,0.84),ncol=1,title="Class",title_fontsize=15,loc=1)
#            ax.scatter([8.5],[3],marker="o",c="scatter")
            #circle(8.5,3.5,"black")
            plt.savefig(outdir+"/%s_%s_scatterplot.png"%(x_pc,y_pc),dpi=300,bbox_inches="tight")
    plt.close("all")

#pca,pc_df = pca_analysis(pearsonr_cor)            
#pc_scatter(pc_df,pca,batch_color_dict,outdir)
