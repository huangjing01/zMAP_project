#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 17:37:58 2020

@author: guixiuqi
"""

import os
import glob
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy
from scipy import stats
import argparse
import random
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

def mkdir(path_dir):
    isexists = os.path.exists(path_dir)
    if not isexists:
        os.makedirs(path_dir)
    else:
        pass

def get_batch_color_dict(batch_df,col_=""):
    
    batch_df = batch_df.loc[batch_df["internal_ref"]!="Yes"]
    num = len(set(batch_df[col_].values))
    #print(num)
    batch_color_dict = {}
    color_list = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])for i in range(num)]

    i=0
    list_ = list(set(batch_df[col_].values))
    list_ = sorted(list_)
    for index in list_:

        batch_color_dict[index] = color_list[i]
        i+=1
    return batch_color_dict

def combine_all_sample_results(indir,outdir,batch_df):    

    condition_color_dict = get_batch_color_dict(batch_df,col_="Sample_condition")
    out_files = indir
    df_list = []
    fig = plt.figure(figsize=(8,10))
    ax = fig.add_subplot(211)
    ax1 = fig.add_subplot(212)
    i=0
    j=0
    for out_file in out_files:
        sample=out_file.split("/")[-2]
        sample_type = batch_df.loc[sample,'Sample_condition']

        color = condition_color_dict[sample_type]
         
        out_file_df = pd.read_csv(out_file,sep="\t",index_col=0)
        out_file_df.columns = [sample+" " +column for column in out_file_df.columns]
        df_list.append(out_file_df)
        m_value = out_file_df["%s adjust mvalue"%(sample)].values
        data = m_value
        data_sorted = np.sort(data)
###CDF picture
# calculate the proportional values of samples
        p = 1. * np.arange(len(data)) / (len(data) - 1)    
        z_sta = out_file_df["%s z statistic"%(sample)].values
        z_sta_sorted = np.sort(z_sta)
###CDF picture
# calculate the proportional values of samples
        p1 = 1. * np.arange(len(z_sta)) / (len(z_sta) - 1)

        ax.plot(data_sorted,p,lw=0.2,c=color)    
        ax1.plot(z_sta_sorted,p1,lw=0.2,c=color)
    color_handles = [mpatches.Patch(color=list(condition_color_dict.values())[i]) for i in range(len(condition_color_dict))]
    color_legend = ax.legend(color_handles, list(condition_color_dict.keys()), loc='upper left', title='Sample_condition',bbox_to_anchor=(1.05, 0.5),ncol=1)
    ax.add_artist(color_legend)
#    ax1.add_artist(color_legend)
       
    ax.set_xlim(-10,10)
    ax1.set_xlim(-10,10)
   # ax.legend(loc="upper left")
   # ax1.legend(loc="upper left")
    ax.set_xlabel("M value",size=15)
    ax.set_ylabel("CDF",size=15)
    ax.plot([-10,10],[0.5,0.5],c="grey",linestyle="-.")
    ax.plot([0,0],[0,1],c="grey",linestyle="-.")
    ax1.set_xlabel("z statistic",size =15)
    ax1.set_ylabel("CDF",size=15)
    ax1.plot([-10,10],[0.5,0.5],c="grey",linestyle="-.")
    ax1.plot([0,0],[0,1],c="grey",linestyle="-.")
    plt.savefig(outdir+"/cumulative_distribution.png",dpi=300,bbox_inches="tight")
    plt.savefig(outdir+"/cumulative_distribution.pdf",dpi=300,bbox_inches="tight")
    plt.close("all")
    all_sample_df = pd.concat(df_list,axis=1,sort=False)
    return  all_sample_df

def get_z_sta_df(all_sample_df,outdir):
    all_z_sta_columns = []
    for column in all_sample_df.columns:
        if "statistic" in column:
            all_z_sta_columns.append(column)
        else:
            pass
     
    all_z_sta_df = all_sample_df[all_z_sta_columns]
    return all_z_sta_df


def calculate_chi2_pvalue(z_sta_table,outdir):
    
    number_of_nan_value_each_protein = list(np.sum(pd.isna(z_sta_table),axis=1).values)
    square_z_sta_df = z_sta_table**2
    chi2_statistic = list(square_z_sta_df.sum(1).values)
    
    pvalue_list = [stats.chi2.sf(chi2_statistic[i],
                                 len(z_sta_table.columns)-number_of_nan_value_each_protein[i])\
    for i in range(len(z_sta_table))]
    z_sta_table["number_of_undetected_samples"] = number_of_nan_value_each_protein
    z_sta_table["chi^2"] = chi2_statistic
    z_sta_table["pvalue"] = pvalue_list
    
    sorted_by_pvalue_table = z_sta_table.sort_values(by="pvalue")
    sorted_by_pvalue_table["BH-corrected P-value"] =[i if i<1 else 1 for i 
                          in sorted_by_pvalue_table["pvalue"]/sorted_by_pvalue_table["pvalue"].rank(axis=0)*len(sorted_by_pvalue_table)]
    
    sorted_by_pvalue_table["Bonferroni-corrected P-value"] =[i if i<1 else 1 for i 
                          in sorted_by_pvalue_table["pvalue"]*len(sorted_by_pvalue_table)]
    sub_sorted_by_pvalue_table = sorted_by_pvalue_table.iloc[:,-5:]
    sub_sorted_by_pvalue_table.to_csv(outdir+"/reverse_zmap_chi_square_pvalue.txt",sep="\t")
    return sorted_by_pvalue_table
    
  


   
