# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 13:54:18 2020

@author: guixiuqi
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import scipy
import scipy.stats as stats
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import argparse
import sys
import re

from . import reverse_zmap_combine


import matplotlib
matplotlib.use('Agg') 

def mkdir(path_dir):
    isexists = os.path.exists(path_dir)
    if not isexists:
        os.makedirs(path_dir)
    else:
        pass

def read_protein_intensity_file(file_for_map):
    #plt.figure(figsize=(4,10))
    intensity_df = pd.read_csv(file_for_map,sep="\t",index_col=0)
    #log2_intensity_df = np.log2(intensity_df+0.5)
    #boxplot1 = intensity_df.boxplot(grid=False,fontsize=15)
    #boxplot2 = log2_intensity_df.boxplot(fontsize=15)
    return intensity_df

def linear_regression(x,y):
    p = np.polyfit(x,y,1)
    regr = LinearRegression()
    regr.fit(x.reshape(-1,1),y.reshape(-1,1))
    r2 = regr.score(x.reshape(-1,1),y.reshape(-1,1))            
    if p[1] > 0:
        p_function = "y= %s x + %s, r-square = %s" %("%.2f"%p[0], "%.2f"%p[1], "%.2f"%r2)
    elif p[1] < 0:
        p_function = "y= %s x - %s, r-square = %s" %("%.2f"%p[0], "%.2f"%p[1], "%.2f"%r2)
    else:
        p_function = "y= %s x, r-square = %s" %("%.2f"%p[0], "%.2f"%r2)
        
    return p_function

def getIndexes(y_predict, y_data):
    
    n = y_data.size
    

    SSE = ((y_data-y_predict)**2).sum()
    
    MSE = SSE/n
    
    RMSE = np.sqrt(MSE)
    
    u = y_data.mean()
    SST = ((y_data -u)**2).sum()
    SSR = SST -SSE
    R_square = SSR/SST
    
    return R_square

       
def normalization(raw_intensity_df,output_dir,l=1.5,batch=""):
    #raw_intensity_df = pd.read_csv(raw_intensity_table,sep="\t",index_col = "Gene Name")
    good_index = []
    for index in raw_intensity_df.index:
        index_df = raw_intensity_df.loc[index] > 10
        num = sum(index_df.values)
        if num > 0:
            good_index.append(index)
        else:
            pass
    raw_intensity_df = raw_intensity_df.loc[good_index]
    sum_trimmed_intensity_list = []
    for column in raw_intensity_df.columns:
        q3 = np.percentile(raw_intensity_df[column],75)
        q1 = np.percentile(raw_intensity_df[column],25)
        top_threshold_intensity = q3 + l*(q3 - q1)
        trimmed_intensity_df = raw_intensity_df[column].loc[raw_intensity_df[column]<=top_threshold_intensity]
        sum_trimmed_intensity = trimmed_intensity_df.sum(axis=0)
        sum_trimmed_intensity_list.append(sum_trimmed_intensity)
    mean_of_sum_tirmmed_intensity = np.mean(sum_trimmed_intensity_list)
    size_factor = np.array(sum_trimmed_intensity_list)/mean_of_sum_tirmmed_intensity
    normalized_intensity_df = raw_intensity_df.div(size_factor)
    log2_normalized_intensity_df = np.log2(normalized_intensity_df + 1)
    log2_raw_intensity_df = np.log2(raw_intensity_df + 1)
    
    plt.figure(figsize=(8,8))
    
    ax = plt.subplot(211)
    ax1 = plt.subplot(212)

    boxprops = dict(linestyle = '-', linewidth = 2, color = 'darkorange')
    medianprops = dict(linestyle = '-', linewidth = 1.5, color = 'firebrick')
    flierprops = dict(marker = 'o', markerfacecolor = 'grey',
                  linestyle = 'none')
    whiskerprops = dict(linestyle = '-', linewidth = 2, color = 'grey')

    capprops = dict(linestyle='-', linewidth=2, color='grey')
    
    ax.boxplot(log2_raw_intensity_df.values,labels= log2_raw_intensity_df.columns,
                vert = False,widths=0.7,
            boxprops = boxprops,medianprops=medianprops,whiskerprops=whiskerprops,
            capprops=capprops,flierprops=flierprops)
    ax.set_xlabel("Raw log$_2$-intensity",size=15)
    ax.set_ylabel("Sample",size=15)
	
    
    ax1.boxplot(log2_normalized_intensity_df.values,labels= log2_normalized_intensity_df.columns,
                vert = False,widths=0.7,
            boxprops = boxprops,medianprops=medianprops,whiskerprops=whiskerprops,
            capprops=capprops,flierprops=flierprops)
    ax1.set_xlabel("Normalized log$_2$-intensity",size=15)
    ax1.set_ylabel("Sample",size=15)
   # plt.title(batch,size=20)
    plt.suptitle(batch,size=20)
    plt.savefig(output_dir+r"/"+"%s_protein_intensity_boxplot.png"%(batch),dpi=300,bbox_inches="tight")
    plt.savefig(output_dir+r"/"+"%s_protein_intensity_boxplot.pdf"%(batch),dpi=300,bbox_inches="tight")
    plt.close("all")
    return log2_normalized_intensity_df


def MS_plotting_position(start_pos, end_pos, total_pos, aa):
    ### start_Pos=1 end_pos = window size+1 total_pos = window size aa =0.375
    from functools import reduce
    constant = (total_pos - aa + 1) / (total_pos - 2*aa + 1)
    
    def get_inner(x):
        r_digital = (x - aa) / (x - aa + 1)
        return r_digital

    pp_list = []
    for j in range(start_pos, end_pos):
        p_i_list = map(get_inner, range(j, total_pos+1))
        p_i = constant * reduce(lambda m,n: m*n, p_i_list)
        pp_list.append(p_i)
    return pp_list

def MA_linear_regression(sample_outdir,log2_normalized_intensity_df,wsize=400,step=100,middle_percent = 0.5):
    sort_normed_intensity_df = log2_normalized_intensity_df.sort_values(by="A value",ascending = True)

    
    sliding_region_left = 0
    sliding_region_right = len(sort_normed_intensity_df)
   # maplot_dir = sample_outdir +r"/maplot"
   # linear_dir = sample_outdir+r"/linear_regression"
   # mkdir(maplot_dir)
   # mkdir(linear_dir)
    all_window_mean_value_list = []
    all_window_variance_list = []
    all_window_intercept_list = []
    r2_list = []
    left = sliding_region_left
    sorted_table = sort_normed_intensity_df
    while left < sliding_region_right:
        if left + wsize > sliding_region_right:
            non_target_window_df = sorted_table.iloc[0:left,:]
        else:
            non_target_window_df = pd.concat([sorted_table.iloc[0:left,:],sorted_table.iloc[left+wsize:,:]],axis=0)
        target_df = sorted_table.iloc[left:min(left+wsize,sliding_region_right),:]
        target_window_df = target_df.sort_values(by="M value",ascending = True)
        
        #target_window_df_M_value = target_window_df["M value"]
        gene_number_in_window = len(target_window_df)
        gene_number_in_window = len(target_window_df)
        if gene_number_in_window > 100:

        
            middle_protein_index = np.arange(gene_number_in_window*(0.5 - middle_percent/2),gene_number_in_window*(0.5 + middle_percent/2))
            edge_protein_index = np.hstack((np.arange(gene_number_in_window*(0.5 - middle_percent/2)),
                                       np.arange(gene_number_in_window*(0.5 + middle_percent/2),gene_number_in_window)))
        
            target_middle_df = target_window_df.iloc[middle_protein_index,:]
            target_edge_df = target_window_df.iloc[edge_protein_index,:]
            window_mean_A = np.mean(target_window_df["A value"])
            all_window_mean_value_list.append(window_mean_A)
        
            ### maplot for each window
            """
            plt.scatter(non_target_window_df["A value"],non_target_window_df["M value"],s=20,facecolors='none', edgecolors="#D3D3D3",
            label = "Proteins out of sliding window")
            plt.scatter(target_edge_df["A value"],target_edge_df["M value"],s=20,
            facecolor = "grey",edgecolors= "black",linewidths=0.2,alpha = 0.9,
            label = "The other proteins in sliding window")
            plt.scatter(target_middle_df["A value"],target_middle_df["M value"],s=20,
            facecolor = "#FF0000",edgecolors="black",linewidths=0.2,alpha = 0.9,
            label = "Proteins in sliding window used\nfor model fitting")        
        
            plt.legend()
            plt.xlabel("log$_2$ ratio",size = 12)
            plt.ylabel("log$_2$ intensity" , size = 12)
            plt.savefig(maplot_dir + r"/window_"+"%s_vaplot.png"%(str(left)),dpi=300)
            plt.savefig(maplot_dir + r"/window_"+"%s_vaplot.pdf"%(str(left)),dpi=300)
            plt.close('all')  
            """
            pp_list = MS_plotting_position(1,gene_number_in_window+1,gene_number_in_window,0.375)
        
            theoretical_quantiles = [scipy.stats.norm.ppf(j) for j in pp_list] 
        
            theoretical_quantiles_used_to_fit = [theoretical_quantiles[int(i)] for i in middle_protein_index]
        
        #print (theoretical_quantiles_used_to_fit)
        #r2,p_function = linear_regression(theoretical_quantiles_used_to_fit,target_middle_df["M value"])
            p = np.polyfit(theoretical_quantiles_used_to_fit,target_middle_df["M value"],1)
    
            regr = LinearRegression()
            model_x = np.asarray(theoretical_quantiles_used_to_fit).reshape(-1,1)
            model_y = target_middle_df["M value"].values.reshape(-1,1)
            regr.fit(model_x,model_y)
            r2 = regr.score(model_x, model_y)
            r2_list.append(r2)
        
            all_window_variance_list.append(p[0]**2)
        
            all_window_intercept_list.append(p[1])
            
            """
            if p[1] > 0:
                p_function = "y= %s x + %s, r-square = %s" %("%.2f"%p[0],"%.2f"%p[1], "%.2f"%r2)
            elif p[1] < 0:
                p_function = "y= %s x - %s, r-square = %s" %("%.2f"%p[0],"%.2f"%p[1], "%.2f"%r2)
            else:
                p_function = "y= %s x, r-square = %s" %("%.2f"%p[0],"%.2f"%r2)
            plt.figure(figsize=(5,5))
            plt.scatter([theoretical_quantiles[int(i)] for i in edge_protein_index],target_edge_df["M value"],
                facecolor = "grey",edgecolors= "black",alpha = 0.7,linewidths=0.2,label = "The other proteins in sliding window")
            plt.scatter(theoretical_quantiles_used_to_fit,target_middle_df["M value"],
            facecolor = "#FF0000",edgecolors="black",linewidths=0.2,alpha = 0.8,label = "Proteins in sliding window used\nfor model fitting")
        
            start = int(len(theoretical_quantiles) * 0.1)
            end = int(len(theoretical_quantiles) * 0.9)
            plt.plot(theoretical_quantiles[start : end+1], 
                               regr.predict(np.array(theoretical_quantiles)[start:end+1].reshape(end+1-start,1)), 
                               color='orange',label ="Fitted linear model",linewidth=2)    
            plt.legend()
            plt.text(2.0,0.1,r"$R^2$=%s"%(str('%.3f'%r2)),fontsize =15)
        
            plt.title(p_function,size=15)
            plt.xlabel("Theoretical quantile of N(0,1)",fontsize=10)
            plt.ylabel("log$_2$-ratio",fontsize=10)
            plt.savefig(linear_dir + r"/linear_regression_window_" + 
                        "%s_linear_regression.png"%(str(left)),dpi=300)
            plt.savefig(linear_dir + r"/linear_regression_window_" + 
                        "%s_linear_regression.pdf"%(str(left)),dpi=300)
            plt.close('all')
            """
            left += step
        else:
            
            left += step
    mean_a_scope_df = pd.DataFrame([all_window_mean_value_list,all_window_variance_list,all_window_intercept_list],
                                   index = ["mean_A_value","variance","intercept"]).T
    mean_a_scope_df.to_csv(sample_outdir + r"/mean_A_estimated_variance.xls",sep ="\t")
    
    
    plt.figure(figsize=(10,5))
    plt.subplot(121)
    plt.boxplot(r2_list)
    plt.ylabel("R^2")
    plt.title("Boxplot for R^2 \nof linear regression",size=15)
        
    plt.subplot(122)
    plt.scatter(sort_normed_intensity_df["A value"],sort_normed_intensity_df["M value"],s=8,alpha=0.2)
    plt.plot([5,25],[0,0],color="black")
    plt.xlim(5,25)
    plt.xlabel("Log$_2$ intensity",size=15)
    plt.ylabel("Log$_2$ ratio",size=15)    
    
    plt.plot(all_window_mean_value_list,all_window_intercept_list,color="red")
    
    plt.savefig(sample_outdir + os.path.split(sample_outdir[:-1])[-1]+"_distribution_r2_intercept_of_linear_regression.png",dpi=300)
    plt.close("all")
        
    return all_window_mean_value_list,all_window_variance_list,all_window_intercept_list

    
              
def func_nl_lsq(x,a,b,c):
    return np.exp(a + (b*x))+ c
def func_nl_lsp_no_constant_term(x,a,b):
    return np.exp(a + (b*x))


################################fit non_linear model###########################################
from . import natural_cubic_model

def fit_natural_model(sample_outdir,log2_normalized_intensity_df,all_window_mean_value_list,all_window_variance_list,
                  used_for_model_fitting=[0.05,0.95]):
    
    nonlinear_dir = sample_outdir + "/nonlinear_model"
    mkdir(nonlinear_dir)
    plt.figure(figsize=(5,4))
    
    a_value = np.asarray(log2_normalized_intensity_df["A value"].sort_values(ascending=True))
    
    start = int(len(all_window_mean_value_list)*used_for_model_fitting[0])
    end = int(len(all_window_mean_value_list)*used_for_model_fitting[1])
    
    x_model = np.asarray(all_window_mean_value_list[start:end])
    y_model = np.asarray(all_window_variance_list[start:end])
    
    nonlinear_model_4 = natural_cubic_model.get_natural_cubic_spline_model(x_model,y_model,minval=min(x_model),maxval=max(x_model),n_knots=4)
    
    R_square_4 = getIndexes(nonlinear_model_4.predict(x_model),y_model)

    nonlinear_model_3 = natural_cubic_model.get_natural_cubic_spline_model(x_model,y_model,minval=min(x_model),maxval=max(x_model),n_knots=3)

    R_square_3 = getIndexes(nonlinear_model_3.predict(x_model),y_model)

    if R_square_3 > R_square_4:
        nonlinear_model = nonlinear_model_3
        R_square = R_square_3
    else:
        nonlinear_model = nonlinear_model_4
        R_square = R_square_4
    
    plt.scatter(all_window_mean_value_list[start:end],all_window_variance_list[start:end],facecolors='none', edgecolors="green",label = "Used for model fitting")
    plt.scatter(all_window_mean_value_list[0:start]+all_window_mean_value_list[end:],all_window_variance_list[0:start]+all_window_variance_list[end:],facecolors='none', edgecolors="red",label = "Not used for model fitting")

    plt.plot(a_value[100:-100],nonlinear_model.predict(a_value[100:-100]),c="darkorange",linewidth=3)
    plt.ylim(0,max(nonlinear_model.predict(a_value[100:-100]))+0.5)    
    plt.xlabel("log$_2$-intensity",size=15)
    plt.ylabel("$Variance\ \sigma^2$",size=15)
    
    #plt.text(12,0.1,"$R^2$=%s"%(str("%.3f"%(R_square))),fontsize = 15)
    plt.legend()
    
    plt.title(sample_outdir.split("/")[-2]+"$R^2$=%s"%(str("%.3f"%(R_square))))
    plt.savefig(nonlinear_dir + "/"+ sample_outdir.split("/")[-2] +"_nonlinear_nature_cubic_model.png",bbox_inches='tight')
    plt.close("all")
    
    return nonlinear_model,R_square   

def func_nl_lsq(x,a,b,c):
    return np.exp(a + (b*x))+ c
def func_nl_lsp_no_constant_term(x,a,b):
    return np.exp(a + (b*x))

##########################fit exp model#################################################
def fit_exp_model(sample_outdir,log2_normalized_intensity_df,all_window_mean_value_list,all_window_variance_list,used_for_model_fitting=[0.05,0.95]):
    
    nonlinear_dir = sample_outdir + "/nonlinear_model"
    mkdir(nonlinear_dir)    
    
    a_value = np.asarray(log2_normalized_intensity_df["A value"].sort_values(ascending=True))
    
    start = int(len(all_window_mean_value_list)*used_for_model_fitting[0])
    end = int(len(all_window_mean_value_list)*used_for_model_fitting[1])
    x_model = np.asarray(all_window_mean_value_list[start:end])
    y_model = np.asarray(all_window_variance_list[start:end])
    popt, pcov = curve_fit(func_nl_lsq, x_model, y_model,bounds=([-50,-50,-1], [50, 50, 1]))
 #   popt, pcov = curve_fit(func_nl_lsq, x_model, y_model)
    if popt[2] > 0:
        R_square = getIndexes(func_nl_lsq(x_model,*popt),y_model)
        plt.figure(figsize=(5,5))
        plt.scatter(all_window_mean_value_list[0:start] + all_window_mean_value_list[end:],
                    all_window_variance_list[0:start]+all_window_variance_list[end:],
                    facecolors="none",s = 50,edgecolors = "red",label = "Not used for model fitting" )
        
        plt.scatter(x_model, y_model,s=50,facecolors='none', edgecolors="green",label = "Used for model fitting" )
        plt.plot(a_value, func_nl_lsq(a_value, *popt), '--',linewidth = 3,color = "darkorange",
             label = "Fitted exponential funcation:\n$\sigma^2=exp(%s-%s*S)+%s$"%(str("%.3f"%(popt[0])),
                                                                   str("%.3f"%(-popt[1])),str("%.3f"%(popt[2]))))
    else:
        popt, pcov = curve_fit(func_nl_lsp_no_constant_term, x_model, y_model,bounds=([-50,-50], [50, 50]))
        R_square = getIndexes(func_nl_lsp_no_constant_term(x_model,*popt),y_model)
        plt.figure(figsize=(5,5))
        plt.scatter(all_window_mean_value_list[0:start]+all_window_mean_value_list[end:],
                    all_window_variance_list[0:start]+all_window_variance_list[end:],
                    facecolors="none",s = 50,edgecolors = "red",label = "Not used for model fitting" )
        
        plt.scatter(x_model, y_model,s=50,facecolors='none', edgecolors="green",label = "Used for model fitting" )
        plt.plot(a_value, func_nl_lsp_no_constant_term(a_value, *popt), '--',linewidth = 3,color = "darkorange",
             label = "Fitted exponential funcation:\n$\sigma^2=exp(%s-%s*S)$"%(str("%.3f"%(popt[0])),str("%.3f"%(-popt[1]))))
        
    plt.xlabel("log$_2$-intensity",size=15)
    plt.ylabel("$Variance\ \sigma^2$",size=15)
  #  plt.title("%s $R^2$=%s"%(batch,str("%.3f"%(R_square))),fontsize = 15)
    plt.ylim(0,1)
    plt.legend()
    plt.savefig(nonlinear_dir + r"/"+sample_outdir.split("/")[-2]+"_nonlinear_exp_model.pdf",dpi=300,bbox_inches='tight')
    plt.savefig(nonlinear_dir + r"/"+sample_outdir.split("/")[-2]+"_nonlinear_exp_model.png",dpi=300,bbox_inches='tight')
    plt.close("all")
    return popt,R_square




def theoretical_quantiles_scatterplot(zmap_result_df,result_dir = "",batch=""):
    z_score_columns=[]
    z_transformed_columns=[]
    for column in zmap_result_df.columns:
        if "z statistic" in column:
            z_transformed_columns.append(column)
    
    plt.figure(figsize=(8,10))
    plt.title(batch,size=20)
    for z_transformed_column in z_transformed_columns:
        z_transformed = zmap_result_df[z_transformed_column].sort_values(ascending=True)
        #z_transformed_quantiles_percent = [stats.percentileofscore(z_transformed,x) for x in z_transformed]
        
        z_transformed_quantiles_percent= MS_plotting_position(1,len(z_transformed)+1,len(z_transformed),0.375)
        
        theoretical_quantiles = [scipy.stats.norm.ppf(j) for j in z_transformed_quantiles_percent]
        plt.scatter(theoretical_quantiles,z_transformed,label=z_transformed_column)
    plt.xlabel("Theoretical quantiles of N(0,1)",size=15)
    plt.ylabel("Z-transformed of $log_2$-intensity",size=15)
    plt.plot(theoretical_quantiles,theoretical_quantiles,'k--',label="y = x")
    plt.legend()
    plt.savefig(result_dir+"/"+"%s_theoretical_quantiles.pdf"%(batch),dpi=300)
    plt.savefig(result_dir+"/"+"%s_theoretical_quantiles.png"%(batch),dpi=300)
    plt.close("all")

