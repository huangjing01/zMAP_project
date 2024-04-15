# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 15:59:34 2022

@author: guixiuqi
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import matplotlib
#matplotlib.use('Agg') 
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
#matplotlib.rcParams['font.family'] ='Arial'
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
import os
import math
from scipy import stats
import seaborn as sns
import argparse
import matplotlib.gridspec as gridspec

def mkdir(path_dir):
    isexists = os.path.exists(path_dir)
    if not isexists:
        os.makedirs(path_dir)
    else:
        pass



def samplelist_for_cluster(sample_info_df,sample_type):
    
    sample_type=sample_type.strip().split(',')
    sample_list = []
    for index in sample_info_df.index:
        sample_condition = sample_info_df.loc[index,"Sample_condition"]
        if sample_condition in sample_type:
            sample_list.append(index)
        else:
            pass
     
    return sample_list


def tumor_mean_variance(tumor_z_df,outdir):
    
    mean_ = tumor_z_df.mean(axis=1)
    
    variance_ = tumor_z_df.var(axis=1)
    x=mean_
    y=variance_
    
    sns.scatterplot(x=x, y=y, s=5, color=".15")
    sns.histplot(x=x, y=y, bins=50, pthresh=.1, cmap="mako")
    sns.kdeplot(x=x, y=y, levels=5, color="w", linewidths=1)
    plt.xlabel("Average of z-statistic")
    plt.ylabel("Variance of z-statistic")
    plt.savefig(outdir+"/mean_variance_scatterplot.png",dpi=300)
    plt.savefig(outdir+"/mean_variance_scatterplot.pdf",dpi=300)
    

def top_tumor_variance_table(z_df,tumor_z_df,outdir,top=3000):
    
    mean_ = tumor_z_df.mean(axis=1)
    
    variance_ = tumor_z_df.var(axis=1).to_frame()
    variance_sort = variance_.sort_values(by=0,ascending=False)
    
    top_index = variance_sort.index[0:top]
        
#    top_variance_all_z_df = tumor_z_df.loc[top_index]

    top_variance_z_df_dropna = z_df.loc[top_index].dropna(how="any")

    top_variance_tumor_z_df_dropna = top_variance_z_df_dropna[tumor_z_df.columns]
    
   

    file_ = outdir+"/top_variance_tumor_z_df_dropna_%s_%s.txt"%(top,len(top_variance_tumor_z_df_dropna))

    top_variance_tumor_z_df_dropna.to_csv(file_,sep="\t")
        
    return file_

def prepare_ConsensusClusterPlus(file_,outdir): 
    out = open(outdir+"/ConsensusClusterPlus.r","w")   
    out.write("setwd('%s')\n"%(outdir))
    out.write("library(ConsensusClusterPlus)\n")
    out.write("dc = read.csv('%s',sep='\t',row.names = 1,check.names=FALSE)\n"%(file_))
    out.write("dc = as.matrix(dc)\n")
    out.write("rcc = ConsensusClusterPlus(dc,maxK=5,reps=1000,pItem=0.8,pFeature=1,title='euclidean_km', \
                           distance='euclidean',clusterAlg='km',plot='pdf',seed=1262118322)\n")
    out.write("cluster <- rcc[[3]]$consensusClass\n")
    out.write("write.csv(cluster,file='%s', quote = FALSE)\n"%(outdir+'/cluster_3.csv'))
    out.write("cluster <- rcc[[4]]$consensusClass\n")
    out.write("write.csv(cluster,file='%s', quote = FALSE)\n"%(outdir+'/cluster_4.csv'))
    out.write("cluster <- rcc[[5]]$consensusClass\n")
    out.write("write.csv(cluster,file='%s', quote = FALSE)\n"%(outdir+'/cluster_5.csv'))
    out.close()
    os.system("Rscript %s"%(outdir+"/ConsensusClusterPlus.r"))



 
if __name__=="__main__":
    parser=argparse.ArgumentParser(description="")
    parser.add_argument("--z_statistic_matrix",help="z-transformed log2-intensity",required=True)
    parser.add_argument("--sample_info",help="sample_information",required=True)
    parser.add_argument("--sample_condition",help="perform clustering on samples under specific condition",required=True)
    parser.add_argument("--outdir",help="results dictionary",required=True)
    parser.add_argument("--top_n",help="proteins were ranked based on their variance across samples, and the z-statistic matrices of the top_n proteins were used for clustering. By default, the top_n is set to 3000",default = 3000)
   
    argv = vars(parser.parse_args())
    z_statistic_matrix = argv['z_statistic_matrix']
    sample_info = argv['sample_info']
    sample_type = argv['sample_condition']
    outdir = argv['outdir']
    top_n = int(argv['top_n'])
    mkdir(outdir)

    z_df = pd.read_csv(z_statistic_matrix,sep="\t",index_col=0)

    sample_info_df = pd.read_csv(sample_info,sep="\t",index_col=0)
     
    sample_list = samplelist_for_cluster(sample_info_df,sample_type)

    overlap_sample = list(set(sample_list)&set(z_df.columns))
    
    tumor_z_df = z_df[overlap_sample]

    tumor_mean_variance(tumor_z_df,outdir)
    print ("Generating mean-variance scatterplot.")
    top_variance_z_file = top_tumor_variance_table(z_df,tumor_z_df,outdir,top=top_n)
    print("Perform clustering on %s %s samples."%(len(overlap_sample),sample_type))
    prepare_ConsensusClusterPlus(top_variance_z_file,outdir)
    print("All finished.")






    
   
    
    
    
    
    
    
    
    
    
    
    
        
        
        
        
    
