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


def combine_batch_results(indir):    
  
    out_files = indir
    df_list = []
    all_z_sta_columns = []
    for out_file in out_files:
        out_file_df = pd.read_csv(out_file,sep="\t",index_col=0)
        del out_file_df['mean']
        column_name = list(out_file_df)
        column_name[1:] = ['{} {}'.format(column_name[0],a) for a in column_name[1:]]        
        out_file_df.columns = column_name
        df_list.append(out_file_df)
        for column in out_file_df.columns:
            if "statistic" in column:
                all_z_sta_columns.append(column)
    
        
            
    all_batch_df = pd.concat(df_list,axis=1)
    z_sta_table = all_batch_df[all_z_sta_columns]
    
    number_of_nan_value_each_protein = list(np.sum(pd.isna(z_sta_table),axis=1).values)
    square_z_sta_df = z_sta_table**2
    chi2_statistic = list(square_z_sta_df.sum(1).values)
    pvalue_list = [stats.chi2.sf(chi2_statistic[i],len(z_sta_table.columns)-number_of_nan_value_each_protein[i]) for i in range(len(z_sta_table))]   
    all_batch_df["Chi^2-statistic"] = chi2_statistic
    all_batch_df["P value"] = pvalue_list    
    sorted_by_pvalue_table = all_batch_df.sort_values(by="P value")
    sorted_by_pvalue_table["BH-corrected P-value"] =[i if i<1 else 1 for i in sorted_by_pvalue_table["P value"]/sorted_by_pvalue_table["P value"].rank(axis=0)*len(sorted_by_pvalue_table)]
    sorted_by_pvalue_table["Bonferroni-corrected P-value"] =[i if i<1 else 1 for i in sorted_by_pvalue_table["P value"]*len(sorted_by_pvalue_table)]
    #writer = pd.ExcelWriter(outdir+"/all_sample_z_statistic_df.xls") 
    sorted_by_pvalue_table.to_csv(os.path.join(os.path.dirname(os.path.dirname(indir[0])) , indir[0].split("/")[-3] + "_reverse_map_output_results.xls"),sep="\t")


    return sorted_by_pvalue_table

def combine_result(map_output_results,condition,out_dir):
    results_df_list=[]
    batch_id_list = []
    for result in map_output_results:
        result_df = pd.read_csv(result,sep="\t",index_col = 0)
        batch_id = os.path.dirname(result).split("/")[-1]
        batch_id_list.append(batch_id)
        result_df.columns = [np.array([batch_id]*len(result_df.columns)),result_df.columns]
        results_df_list.append(result_df)
    all_df = pd.concat(results_df_list,axis=1,sort=False)

    summary_df = pd.DataFrame(index = all_df.index)
    for sample_id in  condition:
        summary_df["mean " + sample_id] = all_df.loc(axis=1)[:,sample_id + " z statistic"].mean(axis=1)
    summary_df["Best pvalue"] = all_df.filter(like="P value").min(axis=1)

    second_pvalue_list=[]
    for index in all_df.index:
        second_pvalue_list.append( all_df.filter(like="P value").loc[index].sort_values()[1])

    summary_df["Second best pvalue"] = second_pvalue_list
    summary_df["Average Chi^2-statistic"] = all_df.filter(like="Chi^2-statistic").mean(axis=1)

    N = len(condition)

    summary_df["Number of detected"] = all_df.filter(like="P value").count(axis=1)
    summary_df["Average pvalue"]= [scipy.stats.gamma.sf(summary_df["Average Chi^2-statistic"].loc[index],
                       summary_df["Number of detected"].loc[index]*(N-1)/2,
                       scale = 2/summary_df["Number of detected"].loc[index]) for index in summary_df.index]
    all_sort_df = summary_df.sort_values(by = "Average pvalue")
    all_sort_df["BH-corrected average pvalue"] =[i if i<1 else 1 for i in \
               all_sort_df["Average pvalue"]/all_sort_df["Average pvalue"].rank(axis=0)*len(all_sort_df)]

    all_sort_df.columns = [np.array(["-"]*len(all_sort_df.columns)),all_sort_df.columns]
    all_batch_final_results = pd.concat([all_sort_df,all_df],axis=1,sort = False)

    all_batch_final_results.to_csv(out_dir+"/"+"all_batch_final_results.xls",sep="\t")
