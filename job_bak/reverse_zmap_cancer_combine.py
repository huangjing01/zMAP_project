#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 17:37:58 2020

@author: guixiuqi
"""

import os
import glob
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy
from scipy import stats
import argparse

def mkdir(path_dir):
    isexists = os.path.exists(path_dir)
    if not isexists:
        os.makedirs(path_dir)
    else:
        pass


def combine_all_sample_results(indir,outdir):    
  
    out_files = indir
    df_list = []
    fig = plt.figure(figsize=(8,10))
    ax = fig.add_subplot(211)
    ax1 = fig.add_subplot(212)
    i=0
    j=0
    for out_file in out_files:
        sample=out_file.split("/")[-2]
        if out_file.split("/")[-2][0] == "T":
            i+=1
            sample_type="Tumor"
            color = "red"
        else:
            sample_type="Normal"
            color = "green"
            j+=1
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

        if (sample_type=="Tumor" and i==1) or (sample_type=="Normal" and j==1):
            ax.plot(data_sorted,p,lw=0.2,c=color,label=sample_type)
        #    ax.set_xlim(-10,10)
            ax1.plot(z_sta_sorted,p1,lw=0.2,c=color,label=sample_type)
        #    ax1.set_xlim(-10,10)
        else:

            ax.plot(data_sorted,p,lw=0.2,c=color)
        #ax.set_xlim(-10,10)
            ax1.plot(z_sta_sorted,p1,lw=0.2,c=color)
        #ax1.set_xlim(-10,10)
    ax.set_xlim(-10,10)
    ax1.set_xlim(-10,10)
    ax.legend()
    ax1.legend()
    ax.set_xlabel("M value",size=15)
    ax.set_ylabel("Summation value",size=15)
    ax.plot([-10,10],[0.5,0.5],c="grey",linestyle="-.")
    ax.plot([0,0],[0,1],c="grey",linestyle="-.")
    ax1.set_xlabel("z statistic",size =15)
    ax1.set_ylabel("Summation value",size=15)
    ax1.plot([-10,10],[0.5,0.5],c="grey",linestyle="-.")
    ax1.plot([0,0],[0,1],c="grey",linestyle="-.")
    plt.savefig(outdir+"/cumulation.png",dpi=300)
    plt.savefig(outdir+"/cumulation.pdf",dpi=300)
    all_sample_df = pd.concat(df_list,axis=1)
    plt.close("all")
    return  all_sample_df

def get_z_sta_df(all_sample_df,outdir):
    all_z_sta_columns = []
    all_m_value_columns = []
    for column in all_sample_df.columns:
        if "statistic" in column:
            all_z_sta_columns.append(column)
        if "mvalue" in column:
            all_m_value_columns.append(column)
        
    all_z_sta_df = all_sample_df[all_z_sta_columns]
#    final_z_sta_df = all_z_sta_df.dropna(axis=0,how='any')
    all_m_value_df = all_sample_df[all_m_value_columns]
#    final_mvalue_df = all_m_value_df.dropna(axis=0,how="any")

    all_tumor_z_sta_columns = []
    all_normal_z_sta_columns = []
    all_tumor_m_value_columns=[]
    all_normal_m_value_columns=[]
    for column in all_z_sta_df.columns:
        if column.startswith("N"):
            all_normal_z_sta_columns.append(column)
        else:
            all_tumor_z_sta_columns.append(column)

    for column in all_m_value_df.columns:
        if column.startswith("N"):
            all_normal_m_value_columns.append(column)
        else:
            all_tumor_m_value_columns.append(column)
        
    all_normal_z_sta = list(all_z_sta_df[all_normal_z_sta_columns].median(axis=0))
    all_tumor_z_sta = list(all_z_sta_df[all_tumor_z_sta_columns].median(axis=0))
    all_normal_m_value = list(all_m_value_df[all_normal_m_value_columns].median(axis=0))
    all_tumor_m_value = list(all_m_value_df[all_tumor_m_value_columns].median(axis=0))
    plt.figure(figsize=(8,5))
    plt.boxplot([all_normal_z_sta,all_tumor_z_sta,all_normal_m_value,all_tumor_m_value])
    plt.plot([0.5,4.5],[0,0])
    plt.ylim(-0.3,0.3)
    plt.xticks([1,2,3,4],["z statistic \nof normal","z statistic\nof tumor","M value\nof normal",
           "M value\nof tumor"],rotation=60,size=15)
    plt.savefig(outdir +"/distribution_of_m_z.png",dpi=300,bbox_inches="tight")
    plt.savefig(outdir +"/distribution_of_m_z.pdf",dpi=300,bbox_inches="tight")
    plt.close("all")
    return all_z_sta_df

#out_files = glob.glob("/Users/guixiuqi/Downloads/reverse_zmap_step1/map_results_2/*/map_output_results.xls") 

def calculate_chi2_pvalue(z_sta_table,outdir):
    
    number_of_nan_value_each_protein = list(np.sum(pd.isna(z_sta_table),axis=1).values)
    square_z_sta_df = z_sta_table**2
    chi2_statistic = list(square_z_sta_df.sum(1).values)
    
    pvalue_list = [stats.chi2.sf(chi2_statistic[i],
                                 len(z_sta_table.columns)-number_of_nan_value_each_protein[i])\
    for i in range(len(z_sta_table))]
    z_sta_table["number_of_undetected_sample"] = number_of_nan_value_each_protein
    z_sta_table["chi2_statistics"] = chi2_statistic
    z_sta_table["pvalue"] = pvalue_list
    
    sorted_by_pvalue_table = z_sta_table.sort_values(by="pvalue")
    sorted_by_pvalue_table["BH-corrected P-value"] =[i if i<1 else 1 for i 
                          in sorted_by_pvalue_table["pvalue"]/sorted_by_pvalue_table["pvalue"].rank(axis=0)*len(sorted_by_pvalue_table)]
    
    sorted_by_pvalue_table["Bonferroni-corrected P-value"] =[i if i<1 else 1 for i 
                          in sorted_by_pvalue_table["pvalue"]*len(sorted_by_pvalue_table)]
    #writer = pd.ExcelWriter(outdir+"/all_sample_z_statistic_df.xls") 
    sorted_by_pvalue_table.to_csv(outdir+"/all_sample_z_statistic_df.xls",sep="\t")
    return sorted_by_pvalue_table


   
