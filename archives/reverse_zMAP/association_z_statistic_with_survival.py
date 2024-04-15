# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 17:14:58 2020

@author: guixiuqi
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt

import argparse
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12

from lifelines import CoxPHFitter
import os

def cox_regression(z_statistic_f,subgroup_df,outdir):
    
    half_n = int(len(z_statistic_f.columns)/2)

    z_statistic_f = z_statistic_f.dropna(thresh=half_n)

    protein_list = []
    results_list = []

    for index in z_statistic_f.index:
        gene_df = z_statistic_f.loc[index].to_frame()
        gene_df["survival_time"] = subgroup_df.loc[gene_df.index,"survival_time"]
        gene_df["death_or_not"] = subgroup_df.loc[gene_df.index,"death_or_not"]
        gene_df = gene_df.dropna(how="any")
    
        cph = CoxPHFitter()
        cph.fit(gene_df, duration_col='survival_time', event_col='death_or_not')
        results = cph.summary
        results_list.append(results)

    all_results_df = pd.concat(results_list)
    all_results_df["BH-corrected pvalue"] =all_results_df["p"]/all_results_df["p"].rank(axis=0) * all_results_df["p"].count()
    
    sub_all_results_df = all_results_df[['coef','coef lower 95%','coef upper 95%','p','BH-corrected pvalue']]
    sub_all_results_df.columns = ['log(HR)','log(HR) lower 95%','log(HR) upper 95%','pvalue','BH-corrected pvalue']
    protein_type = {}
    for index in sub_all_results_df.index:
        if sub_all_results_df.loc[index,'BH-corrected pvalue']<0.05:
            if sub_all_results_df.loc[index,'log(HR) upper 95%'] > 0:
                protein_type[index] = "unfavorable"
            else:
                protein_type[index] = "favorable"
        else:
            protein_type[index] = "not prognostic"
    new_df = pd.DataFrame.from_dict(protein_type,orient="index")
    new_df.columns=["prognostic association"]
#    print(new_df)
    all_df = pd.concat([sub_all_results_df,new_df],axis=1)

    all_df.to_csv(outdir+"/results_cox_regression.txt",sep = "\t")    

def mkdir(path_dir):
    isexists = os.path.exists(path_dir)
    if not isexists:
        os.makedirs(path_dir)
    else:
        pass



if __name__=="__main__":

    parser=argparse.ArgumentParser(description="")
    parser.add_argument("--survival_f",help="patients survival data",required=True)
    parser.add_argument("--z_statistic_matrix",help="z-transformed log2-intensity",required=True)
    parser.add_argument("--outdir",help="outdir for result",required=True)
    argv = vars(parser.parse_args())
    survival_f = argv['survival_f']
    z_statistic = argv['z_statistic_matrix']
    outdir = argv['outdir']
    z_statistic_f = pd.read_csv(z_statistic,sep="\t",index_col=0)
    survival_df = pd.read_csv(survival_f,sep="\t",index_col=0)
    overlap_sample = list(set(z_statistic_f.columns)&set(survival_df.index))
    n=len(overlap_sample) 
    print("Only %s common samples in input files were included in cox proportional hazards regression analysis."%str(n))
    mkdir(outdir)
    z_statistic_f = z_statistic_f[overlap_sample]
    survival_df = survival_df.loc[overlap_sample]
    print("Starting cox proportional hazards regression analysis.This may take a few minutes.")
    cox_regression(z_statistic_f,survival_df,outdir)
    print("All finished.")


