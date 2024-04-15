# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 15:32:47 2020

@author: guixiuqi
"""

import pandas as pd
import numpy as np
import matplotlib
import os
matplotlib.use('Agg') 
import matplotlib.pyplot as plt

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
#matplotlib.rcParams['font.family'] ='Arial'
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.statistics import multivariate_logrank_test
import argparse
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec


def survival_analysis(survival_df,outdir):
    survival_df = survival_df.dropna(how="any")
    kmf = KaplanMeierFitter()
    fig = plt.figure(figsize=(5,5),dpi=300)

    color_dict = {1:"tab:green",2:"tab:blue",3:"tab:red",4:"tab:orange",5:'tab:purple',6:'tab:brown',7:'tab:pink',8:'tab:gray',9:'tab:olive',10:'tab:cyan'}

    ax = fig.gca()

    for name, grouped_df in survival_df.groupby("group"):
        #print (group_df)
        kmf.fit(grouped_df["survival_time"], grouped_df["death_or_not"], label="group"+str(name)+"(%s)"%(len(grouped_df)),
                alpha =0.1)

        kmf.plot(ax=ax,show_censors=True,ci_show=False,color=color_dict[name],
               censor_styles={'ms': 6},linewidth=3)
    from lifelines.statistics import multivariate_logrank_test

    results = multivariate_logrank_test(survival_df['survival_time'],
                                   survival_df['group'], 
                                   survival_df['death_or_not'])
    pvalue = results.p_value
    plt.rcParams['xtick.labelsize']=12
    plt.rcParams['ytick.labelsize']=12

    plt.ylim([0,1.1])
    plt.text(5,0.5,"P-value:" +'{:.2e}'.format(pvalue),size=12)
    ax.set_xlabel("Survival time",size=15)
    ax.set_ylabel("Survival ratio",size=15)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    plt.savefig(outdir+"/survival_analysis.pdf",dpi=300,bbox_inches="tight")
def mkdir(path_dir):
    isexists = os.path.exists(path_dir)
    if not isexists:
        os.makedirs(path_dir)
    else:
        pass


if __name__=="__main__":

    parser=argparse.ArgumentParser(description="")
    parser.add_argument("--input_file",help="survival and group data of patients",required=True)
    parser.add_argument("--outdir",help="outdir for result",required=True)
    argv = vars(parser.parse_args())
    input_file = argv['input_file']
    outdir = argv['outdir']

    survival_df = pd.read_csv(input_file,sep="\t",index_col=0)
    mkdir(outdir)
    print("Plots a figure of the Kaplan–Meier estimate model")
    survival_analysis(survival_df,outdir)
    print("All finished.")




