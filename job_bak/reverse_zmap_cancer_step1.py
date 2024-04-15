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
from . import reverse_zmap_cancer_qc
from . import reverse_zmap_cancer_samplecluster
from . import reverse_zmap_cancer_combine
from . import reverse_zmap_cancer_pca

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

       
def normalization(raw_intensity_df,sample_outdir,l=1.5):
    
    ### log2 转化时为避免0值，丰度�?
    #raw_intensity_df = read_protein_intensity_file(file_for_map) 
    sample = raw_intensity_df.columns[0]
    ref = raw_intensity_df.columns[1]
    log2_intensity_df = np.log2(raw_intensity_df+1)
    plt.figure(figsize=(11,11))
    plt.subplot(221)
    plt.scatter(log2_intensity_df[sample],log2_intensity_df[ref],s=8,alpha=0.1)
    plt.xlim(5,25)
    plt.ylim(5,25)
    plt.xlabel(sample,size=15)
    plt.ylabel(ref,size=15)
    plt.plot([5,25],[5,25])
    R_square = getIndexes(log2_intensity_df[sample].values,log2_intensity_df[ref].values)
    plt.title("Raw intensity\nR^2=%s"%("%.2f"%R_square),size=20)
    
    plt.subplot(223)
    plt.scatter((log2_intensity_df[sample]+log2_intensity_df[ref])/2,
                log2_intensity_df[sample]-log2_intensity_df[ref],
                s=8,alpha=0.2)
    plt.xlabel("Log$_2$ intensity",size=15)
    plt.ylabel("Log$_2$ ratio",size=15)    
    plt.title("M-A plot of raw intensity",size=20)
        
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
    
    log2_normalized_intensity_df = np.log2(normalized_intensity_df+0.5)
    plt.subplot(222)
    plt.scatter(log2_normalized_intensity_df[sample],
                log2_normalized_intensity_df[ref],s=8,alpha=0.1)
    plt.xlim(5,25)
    plt.ylim(5,25)
    plt.xlabel(sample,size=15)
    plt.ylabel(ref,size=15)
    plt.plot([5,25],[5,25])
    R_square = getIndexes(log2_normalized_intensity_df[sample].values,
                          log2_normalized_intensity_df[ref].values)    
    
    plt.title("Normalized intensity\nR^2=%s"%("%.2f"%R_square),size=20)
     
    
    plt.subplot(224)
    plt.scatter((log2_normalized_intensity_df[sample]+log2_normalized_intensity_df[ref])/2,
                log2_normalized_intensity_df[sample]-log2_normalized_intensity_df[ref],
                s=8,alpha=0.2)
    
    plt.xlabel("Log$_2$ intensity",size=15)
    plt.ylabel("Log$_2$ ratio",size=15)    
    plt.title("M-A plot of normalized intensity",size=20)
    plt.savefig(sample_outdir+"intensity_MA_scatter.png",dpi=300)   
    log2_normalized_intensity_df["M value"] = log2_normalized_intensity_df[sample]-log2_normalized_intensity_df[ref]
    log2_normalized_intensity_df["A value"] = (log2_normalized_intensity_df[sample]+log2_normalized_intensity_df[ref])/2
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
    
    plt.savefig(sample_outdir +r"distribution_r2_intercept_of_linear_regression.png",dpi=300)
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
    if R_square > 0.5:
        plt.scatter(all_window_mean_value_list[start:end],all_window_variance_list[start:end],facecolors='none', edgecolors="green",label = "Used for model fitting")
        plt.scatter(all_window_mean_value_list[0:start]+all_window_mean_value_list[end:],all_window_variance_list[0:start]+all_window_variance_list[end:],facecolors='none', edgecolors="red",label = "Not used for model fitting")

        plt.plot(a_value[100:-100],nonlinear_model.predict(a_value[100:-100]),c="darkorange",linewidth=3)
        plt.ylim(0,max(nonlinear_model.predict(a_value[100:-100]))+0.5)    
        plt.xlabel("log$_2$-intensity",size=15)
        plt.ylabel("$Variance\ \sigma^2$",size=15)
        
        #plt.text(12,0.1,"$R^2$=%s"%(str("%.3f"%(R_square))),fontsize = 15)
        plt.legend()
        
        plt.title(sample_outdir.split("/")[-2]+"$R^2$=%s"%(str("%.3f"%(R_square))))
        plt.savefig(nonlinear_dir + r"/nonlinear_model.png",bbox_inches='tight')
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
    plt.savefig(nonlinear_dir + r"/"+"nonlinear_exp_model.pdf",dpi=300,bbox_inches='tight')
    plt.savefig(nonlinear_dir + r"/"+"nonlinear_exp_model.png",dpi=300,bbox_inches='tight')
    plt.close("all")
    return popt,R_square


# if __name__=="__main__":
#     parser=argparse.ArgumentParser(description="")
#     parser.add_argument("--filein",help="protein intensity file",required=True)
#     parser.add_argument("--outdir",help="outdir for result",required=True)
#     parser.add_argument("--window_size",help="window size for linear regression",required=True)
#     parser.add_argument("--step_size",help="step for sliding window",required=True)
#     parser.add_argument("--percent",help="middle percent proteins used to linear fitting",required=True)
#     parser.add_argument("--method",help="expotential_function or natural_cubic_spline",required = True)
#     argv=vars(parser.parse_args())
#     filein=argv['filein']
#     outdir=argv['outdir']
#     window_size=int(argv['window_size'])
#     step_size=int(argv['step_size'])
#     percent = int(argv['percent'])/100
#     method = argv["method"]
    
#    # filein = r"H:\project\protein\web_server\fig_for_reversezmap_workflow\reverse_zmap_input_final_1.xlsx"
#    # outdir = r"H:\project\protein\web_server\fig_for_reversezmap_workflow\test"
#     mkdir(outdir) 
#     raw_intensity_data = pd.ExcelFile(filein)
#     sheet_names = raw_intensity_data.sheet_names
    
#     protein_intensity_df = raw_intensity_data.parse(sheet_names[0],index_col=0)
    
#     batch_df = raw_intensity_data.parse(sheet_names[1],index_col=0)
#     final_protein_intensity_df = check_df_names(protein_intensity_df,batch_df)
#     reverse_zmap_qc.data_qc(final_protein_intensity_df,fig_dir = outdir)
    
#     patient_column_list = []
#     for column in batch_df.columns:
#         if column.startswith("ref"):
#             pass
#         else:
#             patient_column_list.append(column)
#     paired_df_list = []
#     for index in batch_df.index:
#         ref_sample = batch_df.loc[index,"ref_sample"]
#         patient_samples = batch_df.loc[index,patient_column_list]
#         for sample in patient_samples:
#             if isinstance(sample,str):
#                 paired_df_list.append(final_protein_intensity_df[[sample,ref_sample]])
#             else:
#                 pass
#     sample_results_list = []
#     for paired_df in paired_df_list[20:]:
#         print (paired_df.columns)
#         column_list = paired_df.columns
#         new_paired_df = paired_df.loc[(paired_df[column_list[0]]!=0) & (paired_df[column_list[1]]!=0)]
#         #new_paired_df = paired_df.loc[paired_df[""]]
#         log2_normalized_intensity_df,sample_outdir = normalization(new_paired_df) 
#        # """
#         all_window_mean_value_list,all_window_variance_list,all_window_intercept_list = \
#         MA_linear_regression(log2_normalized_intensity_df,wsize=window_size,step=step_size,\
#                       middle_percent = percent)  
#         if method == "natural_cubic_spline":
#             nonlinear_model,R_square = fit_natural_model(log2_normalized_intensity_df,all_window_mean_value_list,\
#                                 all_window_variance_list,used_for_model_fitting=[0.05,0.95])
#  #   r_square_list.append(R_square)   
    
#             if R_square > 0.5:
#                 log2_normalized_intensity_df["estimated variance"] =  \
#                 nonlinear_model.predict(log2_normalized_intensity_df["A value"].values)  
#                 for index in log2_normalized_intensity_df.index:
#                     if log2_normalized_intensity_df.loc[index,"A value"]< all_window_mean_value_list[0]:
#                         log2_normalized_intensity_df.loc[index,"estimated variance"] = all_window_variance_list[0]
#                     elif log2_normalized_intensity_df.loc[index,"A value"] > all_window_mean_value_list[-1]:
#                         log2_normalized_intensity_df.loc[index,"estimated variance"] = all_window_variance_list[-1]
#                     else:
#                         pass
#             else:
#                 log2_normalized_intensity_df["estimated variance"] = \
#                 [np.mean(all_window_variance_list)] * len(log2_normalized_intensity_df)
                
#         else:
#             popt,exp_R_square = fit_exp_model(log2_normalized_intensity_df,all_window_mean_value_list,all_window_variance_list,used_for_model_fitting=[0.05,0.95])
            
#             if len(popt) == 3:
#                 exp_function = func_nl_lsq
#             else:
#                 exp_function = func_nl_lsp_no_constant_term
            
#             if exp_R_square > 0.5:
                
#                 log2_normalized_intensity_df["estimated variance"] = exp_function(log2_normalized_intensity_df["A value"].values,*popt) 
                
#                 for index in log2_normalized_intensity_df.index:
                    
#                     if log2_normalized_intensity_df.loc[index,"A value"]< all_window_mean_value_list[0]:
#                         log2_normalized_intensity_df.loc[index,"estimated variance"] = exp_function(all_window_variance_list[0],*popt)
                        
#                     elif log2_normalized_intensity_df.loc[index,"A value"] > all_window_mean_value_list[-1]:
#                         log2_normalized_intensity_df.loc[index,"estimated variance"] = exp_function(all_window_variance_list[-1],*popt)
                        
#                     else:
#                         pass
#             else:
#                 log2_normalized_intensity_df["estimated variance"] = \
#                 [np.mean(all_window_variance_list)] * len(log2_normalized_intensity_df)
                
#     ### adjust M-value
#         x_model = np.asarray(all_window_mean_value_list)
#         y_model = np.asarray(all_window_intercept_list)
#         intercept_model = natural_cubic_model.get_natural_cubic_spline_model(x_model,y_model,minval=min(x_model),maxval=max(x_model),n_knots=4)
#         intercept_rsquare = getIndexes(intercept_model.predict(x_model),y_model)
#     #intercept_R_square_list.append(intercept_rsquare)
#         from scipy.stats import norm
#         pvalue_list = []
#         z_statiatic_list = []
#         adjust_mvalue =[]
#         fig = plt.figure(figsize=(10,5))
#         ax1= fig.add_subplot(121)
#         ax2 = fig.add_subplot(122)
    
#         ax1.scatter(log2_normalized_intensity_df["A value"],log2_normalized_intensity_df["M value"],c="blue",alpha=0.2,s=8)
#         ax1.plot([5,25],[0,0],c="grey")
#         ax1.scatter(all_window_mean_value_list,all_window_intercept_list,s=15)
#         ax1.plot(x_model,intercept_model.predict(x_model),c="red")
#     #ax1.plot([5,25],[0,0],c="grey")
#         ax2.scatter(log2_normalized_intensity_df["A value"].values,
#                 [log2_normalized_intensity_df.loc[index,"M value"]\
#                  -intercept_model.predict(log2_normalized_intensity_df.loc[index,"A value"])[0]\
#                  for index in log2_normalized_intensity_df.index],c="blue",alpha=0.2,s=8)
#         ax2.plot([5,25],[0,0],c="grey")
#         for index in log2_normalized_intensity_df.index:
    
#         #mvalue_abs = abs(log2_normalized_intensity_df.loc[index,"M value"])
#             mvalue = log2_normalized_intensity_df.loc[index,"M value"] \
#             - intercept_model.predict(log2_normalized_intensity_df.loc[index,"A value"])[0]
#             adjust_mvalue.append(mvalue)
#             mvalue_abs = abs(mvalue)
#         #ax2.scatter(log2_normalized_intensity_df.loc[index,"A value"],mvalue,c="blue",alpha=0.2,s=8)    
        
#             standard_variance = log2_normalized_intensity_df.loc[index,"estimated variance"]**0.5
    
#             pvalue_list.append(2*norm.sf(mvalue_abs,loc=0,scale=standard_variance))
#             z_statiatic_list.append(mvalue/standard_variance)
#         log2_normalized_intensity_df["adjust mvalue"] = adjust_mvalue
#         log2_normalized_intensity_df["pvalue"] = pvalue_list
#         log2_normalized_intensity_df["z statistic"] = z_statiatic_list 
#         log2_normalized_intensity_df.to_csv(sample_outdir+"map_output_results.xls",sep="\t")
#         sample_results_list.append(sample_outdir+"map_output_results.xls")
#         plt.savefig(sample_outdir+"mascatter.png",dpi=300)
#         plt.close("all")
#     all_sample_df  = reverse_zmap_combine.combine_all_sample_results(sample_results_list,outdir)
#     z_sta_table = reverse_zmap_combine.get_z_sta_df(all_sample_df,outdir).copy(deep=True)
#     new_z_sta_table = z_sta_table.copy()
#     new_z_sta_table.columns = [column.split(" ")[0] for column in new_z_sta_table.columns]
    
#     batch_color_dict = reverse_zmap_samplecluster.get_batch_color_dict(batch_df)
#     color_list_df = reverse_zmap_samplecluster.get_color_list(new_z_sta_table.columns,batch_color_dict)
#     pearsonr_cor = reverse_zmap_samplecluster.correlation(new_z_sta_table,method="pearsonr")
#     reverse_zmap_samplecluster.cluster_map(new_z_sta_table,pearsonr_cor,color_list_df,vmin=-0.5,vmax=1,fig=outdir+"/clustermap")
#     pca,pc_df = reverse_zmap_pca.pca_analysis(pearsonr_cor)            
#     reverse_zmap_pca.pc_scatter(pc_df,pca,batch_color_dict,outdir)        
#     final_z_sta_pvalue_df = reverse_zmap_combine.calculate_chi2_pvalue(z_sta_table,outdir)

# #plt.figure(figsize=(10,5))
# #plt.hist(r_square_list,bins=20)
# #plt.xlabel("R square of natural cubic spline",size=15)
# #plt.ylabel("Sample numbers",size=15)
# #plt.title("Estimated variance",size=20)
# #plt.savefig(outdir+"r_square.png",dpi=300)

# #plt.figure(figsize=(10,5))
# #plt.hist(intercept_R_square_list,bins=20)
# #plt.xlabel("R square of natural cubic spline",size=15)
# #plt.ylabel("Sample numbers",size=15)
# #plt.title("Intercept of linear regression",size=20)
# #plt.savefig(outdir+"intercept_r_square.png",dpi=300)





    



    
    
    

    



















