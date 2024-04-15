# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 10:11:16 2020

@author: guixiuqi
"""

import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
matplotlib.use('Agg') 
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
#matplotlib.rcParams['font.family'] ='Arial'


def barplot_for_u_fitting(outdir):
	
#r2_nonlinear = r"H:\project\protein\web_server\reverse_zmap_script\workflow\r2_of_u_natural_cubic_model.csv"
#import os
#os.chdir(r"H:\project\protein\web_server\reverse_zmap_script\workflow")

	r2_df = pd.read_csv(outdir+"/r2_of_natural_cubic_model_fitting_for_intercept_u.csv",sep="\t",index_col=0)

	r2_list = list(r2_df["r^2_of_u_natural_cubic_model"].values)

	r2_list_sort = sorted(r2_list)
	q1 = r2_list_sort[int(len(r2_list)*0.25)]
	q2 = r2_list_sort[int(len(r2_list)*0.5)]
	q3 = r2_list_sort[int(len(r2_list)*0.75)]

#sns.set(style="ticks", palette="muted", color_codes=True)

    
	fig=plt.figure(figsize=(6,3),dpi=300)
	ax=fig.gca()
	plt.hist(r2_list,bins=50,color="#2980b9")

	plt.plot([q1,q1],[0,1.5],color="red")
	plt.plot([q2,q2],[0,1.5],color="orange")
	plt.plot([q3,q3],[0,1.5],color="green")
	plt.text(q1-0.025,3,"Q1")
	plt.text(q2-0.025,3,"Q2")
	plt.text(q3-0.025,3,"Q3")
	plt.xticks([0.2,0.4,0.6,0.8])
	plt.rcParams['xtick.labelsize']=12
	plt.rcParams['ytick.labelsize']=12
#	plt.yticks(range(0,31,10))

	ax.xaxis.set_minor_locator(AutoMinorLocator())
	ax.yaxis.set_minor_locator(AutoMinorLocator())


	plt.ylabel("Sample number",size=15)
	plt.xlabel("$R^2$ of nonlinear model for $\mu$",size=15)
	plt.savefig(outdir + "/r2_distribution_of_nonlinear_fitting_for_intercept_u.pdf",dpi=300,bbox_inches="tight")

#barplot_for_u_fitting("/picb/rsgeno/guixiuqi/protein/reverse_zmap/reverse_zmap_script/test")
