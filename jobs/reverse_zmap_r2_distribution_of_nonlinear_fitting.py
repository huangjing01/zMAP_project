# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 16:45:38 2020

@author: guixiuqi
"""

import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
#matplotlib.rcParams['font.family'] ='Arial'

def barplot_r2_estimated_variance(outdir):

#r2_nonlinear = r"H:\project\protein\web_server\reverse_zmap_script\workflow\r2_of_nonlinear_model.csv"
#import os
#os.chdir(r"H:\project\protein\web_server\reverse_zmap_script\workflow")

	r2_df = pd.read_csv(outdir+"/r2_of_nonlinear_model_fitting_for_estimated_variance.csv",sep="\t",index_col=0)

	r2_list = list(r2_df["r^2_of_nonlinear_model"].values)

	r2_list_sort = sorted(r2_list)
	q1 = r2_list_sort[int(len(r2_list)*0.25)]
	q2 = r2_list_sort[int(len(r2_list)*0.5)]
	q3 = r2_list_sort[int(len(r2_list)*0.75)]

#sns.set(style="ticks", palette="muted", color_codes=True)
    	#matplotlib.rcdefaults()
	fig = plt.figure(figsize=(6,3),dpi=300)
	ax=fig.gca()
# sns.despine(left=True)
# # sns.distplot(r2_list,bins=50,kde=False,axlabel="$R^2$ of nonlinear model for $\sigma^2$",color="#0984e3",
# #              hist_kws=dict(edgecolor="#0984e3"))

# #74b9ff
# sns.distplot(r2_list,bins=50,kde=False,axlabel="$R^2$ of nonlinear model for $\sigma^2$",color="#2980b9",
#              hist_kws=dict(edgecolor="#2980b9"))
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
	plt.xlabel("$R^2$ of nonlinear model for $\sigma^2$",size=15)
	plt.ylabel("Sample number",size=15)

	plt.savefig(outdir+"/r2_distribution_of_nonlinear_model_fitting_for_estimated_variance.pdf.pdf",dpi=300,bbox_inches="tight")

#barplot_r2_estimated_variance("/picb/rsgeno/guixiuqi/protein/reverse_zmap/reverse_zmap_script/test_new")

