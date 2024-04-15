# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 13:54:18 2020

@author: guixiuqi@gmail.com
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import scipy
import scipy.stats as stats
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import argparse
import sys
import re
from . import reverse_zmap_qc
#import reverse_zmap_samplecluster
from . import reverse_zmap_combine
#import reverse_zmap_pca
from . import reverse_zmap_r2_distribution_of_linear_fitting
from . import reverse_zmap_r2_distribution_of_nonlinear_fitting
from . import reverse_zmap_r2_distribution_of_nonlinear_fitting_for_intercept_u
#import matplotlib
#matplotlib.use('Agg') 
#import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
#matplotlib.rcParams['font.family'] ='Arial'
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

       
def normalization(raw_intensity_df,outdir,l=1.5):
    
    ### log2(intensity+1)
    #raw_intensity_df = read_protein_intensity_file(file_for_map) 
    sample = raw_intensity_df.columns[0]
    ref = raw_intensity_df.columns[1]
    sample_outdir = outdir + "/" + sample + "/"
    mkdir(sample_outdir)
    log2_intensity_df = np.log2(raw_intensity_df+1)
    plt.figure(figsize=(11,11))
    plt.subplot(221)
    plt.scatter(log2_intensity_df[sample],log2_intensity_df[ref],s=8,alpha=0.1)
#    plt.xlim(0,30)
#    plt.ylim(0,30)
    plt.xlabel(sample,size=15)
    plt.ylabel(ref,size=15)
    min_value = min(min(log2_intensity_df[sample]),min(log2_intensity_df[ref]))
    max_value = max(max(log2_intensity_df[sample]),max(log2_intensity_df[ref]))
    plt.plot([min_value-1,max_value+1],[min_value-1,max_value+1],c="black",linestyle="--")
    plt.xlim([min_value-1,max_value+1])
    plt.ylim([min_value-1,max_value+1])
    R_square = getIndexes(log2_intensity_df[sample].values,log2_intensity_df[ref].values)
    plt.title("Raw intensity\n$R^2$=%s"%("%.2f"%R_square),size=20)
    
    plt.subplot(223)
    plt.scatter((log2_intensity_df[sample]+log2_intensity_df[ref])/2,
                log2_intensity_df[sample]-log2_intensity_df[ref],
                s=8,alpha=0.2)
    plt.xlim([min_value-1,max_value+1])
 #   plt.ylim([min_value-1,max_value+1])
    plt.plot([min_value-1,max_value+1],[0,0],color="black",linestyle="--")
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
    plt.plot([min_value-1,max_value+1],[min_value-1,max_value+1],c="black",linestyle="--")
   # plt.xlim(0,30)
   # plt.ylim(0,30)
    plt.xlim([min_value-1,max_value+1])
    plt.ylim([min_value-1,max_value+1])
    plt.xlabel(sample,size=15)
    plt.ylabel(ref,size=15)
   # plt.plot([0,30],[0,30])
    R_square = getIndexes(log2_normalized_intensity_df[sample].values,
                          log2_normalized_intensity_df[ref].values)    
    
    plt.title("Normalized intensity\n$R^2$=%s"%("%.2f"%R_square),size=20)
     
    
    plt.subplot(224)
    plt.scatter((log2_normalized_intensity_df[sample]+log2_normalized_intensity_df[ref])/2,
                log2_normalized_intensity_df[sample]-log2_normalized_intensity_df[ref],
                s=8,alpha=0.2)
    plt.xlim([min_value-1,max_value+1])
#    plt.ylim([min_value-1,max_value+1])
    plt.plot([min_value-1,max_value+1],[0,0],color="black",linestyle="--")
    
    plt.xlabel("Log$_2$ intensity",size=15)
    plt.ylabel("Log$_2$ ratio",size=15)    
    plt.title("M-A plot of normalized intensity",size=20)
    plt.savefig(sample_outdir+"intensity_MA_scatter.png",dpi=300)   
    log2_normalized_intensity_df["M value"] = log2_normalized_intensity_df[sample]-log2_normalized_intensity_df[ref]
    log2_normalized_intensity_df["A value"] = log2_normalized_intensity_df[ref]
    return log2_normalized_intensity_df,sample,sample_outdir


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

def MA_linear_regression(log2_normalized_intensity_df,sample_outdir,wsize=400,step=100,middle_percent = 0.5):
    sort_normed_intensity_df = log2_normalized_intensity_df.sort_values(by="A value",ascending = True)
    sliding_region_left = 0
    sliding_region_right = len(sort_normed_intensity_df)
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
    mean_a_scope_df.to_csv(sample_outdir + r"/mean_A_estimated_variance.txt",sep ="\t")
    
    
    plt.figure(figsize=(10,5))
    plt.subplot(121)
    plt.boxplot(r2_list)
    plt.ylabel("$R^2$",size=15)
    plt.title("Boxplot for $R^2$ \nof linear regression",size=15)
    out=open(sample_outdir+"linear_model_r2.txt","w")
    for i in r2_list:
        out.write(str(i)+"\n")
    out.close()
        
    plt.subplot(122)
    plt.scatter(sort_normed_intensity_df["A value"],sort_normed_intensity_df["M value"],s=8,alpha=0.2)
    min_value = min(sort_normed_intensity_df["A value"])
    max_value = max(sort_normed_intensity_df["A value"])
    plt.plot([min_value-1,max_value+1],[0,0],color="black",linestyle="--")
    plt.xlim([min_value-1,max_value+1])
    plt.xlabel("Log$_2$ intensity",size=15)
    plt.ylabel("Log$_2$ ratio",size=15)    
    
    plt.scatter(all_window_mean_value_list,all_window_intercept_list,edgecolors="tomato",facecolors="none")
    
    plt.savefig(sample_outdir +r"distribution_r2_intercept_of_linear_regression.png",dpi=300)
    plt.close("all")
        
    return all_window_mean_value_list,all_window_variance_list,all_window_intercept_list

    
              
def func_nl_lsq(x,a,b,c):
    return np.exp(a + (b*x))+ c
def func_nl_lsp_no_constant_term(x,a,b):
    return np.exp(a + (b*x))


################################fit non_linear model###########################################
from . import natural_cubic_model
def fit_natural_model(log2_normalized_intensity_df,all_window_mean_value_list,all_window_variance_list,sample_outdir,
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
#        print(R_square_3,R_square_4,3)
    else:
        nonlinear_model = nonlinear_model_4
        R_square = R_square_4
#        print(R_square_3,R_square_4,4)
    
    plt.scatter(all_window_mean_value_list[start:end],all_window_variance_list[start:end],facecolors='none', edgecolors="green",label = "Used for model fitting")
    plt.scatter(all_window_mean_value_list[0:start]+all_window_mean_value_list[end:],all_window_variance_list[0:start]+all_window_variance_list[end:],facecolors='none', edgecolors="red",label = "Not used for model fitting")

    plt.plot(a_value[400:-100],nonlinear_model.predict(a_value[400:-100]),"--",c="darkorange",linewidth=3)
    plt.ylim(0,max(nonlinear_model.predict(a_value[400:-100]))+0.5)    
    plt.xlabel("log$_2$-intensity",size=15)
    plt.ylabel("$Variance\ \sigma^2$",size=15)
    
    #plt.text(12,0.1,"$R^2$=%s"%(str("%.3f"%(R_square))),fontsize = 15)
    plt.legend()
    
    plt.title(sample_outdir.split("/")[-2]+" $R^2$=%s"%(str("%.3f"%(R_square))))
    plt.savefig(nonlinear_dir + r"/nonlinear_model.png",dpi=300,bbox_inches='tight')
    plt.savefig(nonlinear_dir + r"/nonlinear_model.pdf",dpi=300,bbox_inches='tight')
    plt.close("all")
    
    return nonlinear_model,R_square   

def func_nl_lsq(x,a,b,c):
    return np.exp(a + (b*x))+ c
def func_nl_lsp_no_constant_term(x,a,b):
    return np.exp(a + (b*x))

##########################fit exp model#################################################
def fit_exp_model(log2_normalized_intensity_df,all_window_mean_value_list,all_window_variance_list,sample_outdir,used_for_model_fitting=[0.05,0.95]):
    
    nonlinear_dir = sample_outdir + "/nonlinear_model"
    mkdir(nonlinear_dir)    
    
    a_value = np.asarray(log2_normalized_intensity_df["A value"].sort_values(ascending=True))
    
    start = int(len(all_window_mean_value_list)*used_for_model_fitting[0])
    end = int(len(all_window_mean_value_list)*used_for_model_fitting[1])
    x_model = np.asarray(all_window_mean_value_list[start:end])
    y_model = np.asarray(all_window_variance_list[start:end])
    popt, pcov = curve_fit(func_nl_lsq, x_model, y_model,bounds=([0,-10,-1], [10, 10, 1]))
 #   popt, pcov = curve_fit(func_nl_lsq, x_model, y_model)
    if popt[2] > 0:
        R_square = getIndexes(func_nl_lsq(x_model,*popt),y_model)
        plt.figure(figsize=(5,5))
        plt.scatter(all_window_mean_value_list[0:start] + all_window_mean_value_list[end:],
                    all_window_variance_list[0:start]+all_window_variance_list[end:],
                    facecolors="none",s = 50,edgecolors = "red",label = "Not used for model fitting" )
        
        plt.scatter(x_model, y_model,s=50,facecolors='none', edgecolors="green",label = "Used for model fitting" )
        plt.plot(a_value[400:-100], func_nl_lsq(a_value[400:-100], *popt), '--',linewidth = 3,color = "darkorange",
             label = "Fitted exponential funcation:\n$\sigma^2=exp(%s-%s*S)+%s$"%(str("%.3f"%(popt[0])),
                                                                   str("%.3f"%(-popt[1])),str("%.3f"%(popt[2]))))
        plt.ylim(0,max(func_nl_lsq(a_value[400:-100], *popt))+0.5)
    else:
        popt, pcov = curve_fit(func_nl_lsp_no_constant_term, x_model, y_model,bounds=([-50,-50], [50, 50]))
        R_square = getIndexes(func_nl_lsp_no_constant_term(x_model,*popt),y_model)
        plt.figure(figsize=(5,5))
        plt.scatter(all_window_mean_value_list[0:start]+all_window_mean_value_list[end:],
                    all_window_variance_list[0:start]+all_window_variance_list[end:],
                    facecolors="none",s = 50,edgecolors = "red",label = "Not used for model fitting" )
        
        plt.scatter(x_model, y_model,s=50,facecolors='none', edgecolors="green",label = "Used for model fitting" )
        plt.plot(a_value[400:-100], func_nl_lsp_no_constant_term(a_value[400:-100], *popt), '--',linewidth = 3,color = "darkorange",
             label = "Fitted exponential funcation:\n$\sigma^2=exp(%s-%s*S)$"%(str("%.3f"%(popt[0])),str("%.3f"%(-popt[1]))))
        plt.ylim(0,max(func_nl_lsp_no_constant_term(a_value[400:-100], *popt))+0.5)
        
    plt.xlabel("log$_2$-intensity",size=15)
    plt.ylabel("$Variance\ \sigma^2$",size=15)
   # plt.title("%s $R^2$=%s"%(batch,str("%.3f"%(R_square))),fontsize = 15)
    plt.title(sample_outdir.split("/")[-2]+" $R^2$=%s"%(str("%.3f"%(R_square))))
#    plt.ylim(0,1)
    plt.legend()
    plt.savefig(nonlinear_dir + r"/"+"nonlinear_exp_model.pdf",dpi=300,bbox_inches='tight')
    plt.savefig(nonlinear_dir + r"/"+"nonlinear_exp_model.png",dpi=300,bbox_inches='tight')
    plt.close("all")
    return popt,R_square


def check_df_names(intensity_df,batch_df,intensity_data):
    not_uniq_gene_list = []
    nan_sample = 0
    nan_gene = 0
    not_uniq_sample_list = []  
    final_gene_list = []
    final_sample_list = []

    ### check gene symbol
    for index in intensity_df.index:
        if index in final_gene_list:
            not_uniq_gene_list.append(index)
        elif isinstance(index,float):
  #          print(index)
            nan_gene +=1
        else:
            final_gene_list.append(index)

    ### check sample id
    with open(intensity_data) as f:
        firstline = f.readline()
    sample_columns = firstline.strip().split("\t")[1:]
   # print(sample_columns) 
    for column in sample_columns:
        if column in final_sample_list:
            not_uniq_sample_list.append(column)
        else:
            final_sample_list.append(column)

    if len(not_uniq_gene_list) >0 :
        print ("%s gene symbols were not unique:%s,Only the first one will be retained."%(str(len(not_uniq_gene_list)),";".join(list(set(not_uniq_gene_list)))))
    if len(not_uniq_sample_list) > 0:
        print ("%s sample ids were not unique:%s,please check your inputs."%(str(len(not_uniq_sample_list)),";".join(list(set(not_uniq_sample_list)))))
        sys.exit(0)
    
    no_intensity_sample = []
    for sample in batch_df.index:
        if sample in intensity_df.columns:
            pass
        else:
            no_intensity_sample.append(sample)

    if len(no_intensity_sample) >0:

        print ("No protein profiles was found for the sample %s,please check your inputs"%(";".join(no_intensity_sample)))
        sys.exit(0)
    else:
        new_intensity_df = intensity_df.loc[final_gene_list,final_sample_list]

        return new_intensity_df.loc[~new_intensity_df.index.duplicated(keep='first')].copy()

"""
if __name__=="__main__":
    parser=argparse.ArgumentParser(description="")
    parser.add_argument("--intensity_data",help="protein intensity file",required=True)
    parser.add_argument("--sample_information",help="sample information",required=True)
    parser.add_argument("--outdir",help="outdir for result",required=True)
    parser.add_argument("--window_size",help="window size for linear regression",required=True)
    parser.add_argument("--step_size",help="step for sliding window",required=True)
    parser.add_argument("--percent",help="middle percent proteins used to linear fitting,30<=percent<=50 ",required=True)
    parser.add_argument("--method",help="exponential_function or natural_cubic_spline",required = True)
    argv=vars(parser.parse_args())
    intensity_data=argv['intensity_data']
    sample_information = argv['sample_information']
    outdir=argv['outdir']
    window_size=int(argv['window_size'])
    step_size=int(argv['step_size'])
    percent = int(argv['percent'])/100
    method = argv["method"]
    mkdir(outdir) 
    import math
"""
def reverse_zmap_(intensity_data,sample_information,outdir,window_size,step_size,percent,method):    
    print ("start processing.")
    print("read data.")
    protein_intensity_df = pd.read_csv(intensity_data,index_col=0,sep="\t")    
    batch_df = pd.read_csv(sample_information,sep='\t',index_col=0)
    ### Check the input data
    print("check the input data")
    final_protein_intensity_df = check_df_names(protein_intensity_df,batch_df,intensity_data)


    ### make a dict,key is MS_run,value is internal ref in corresponding MS_run
    ms_internal_dict = {}
    for index in batch_df.index:
        sample_type = batch_df.loc[index,'internal_ref']
        ms_run_ = batch_df.loc[index,'MS_run']
        if sample_type=="Yes":
        
            ms_internal_dict[ms_run_]=index


    print("summarize the input data")
    reverse_zmap_qc.data_qc(final_protein_intensity_df,fig_dir = outdir)
    

    paired_df_list = []
    for index in batch_df.index:
        ms_run = batch_df.loc[index,"MS_run"]

        ref_sample = ms_internal_dict[ms_run]

        if batch_df.loc[index,'internal_ref']!="Yes":
            paired_df_list.append(final_protein_intensity_df[[index,ref_sample]])
           
    sample_results_list = []
    nonlinear_r2=[]
    sample_list=[]
    intercept_R_square_list = []
    print("start pairwise comparison to reference sample.")
    for paired_df in paired_df_list:
        column_list = paired_df.columns
        
        new_paired_df = paired_df.dropna(how="any",axis=0)
        print("%s vs %s"%(column_list[0],column_list[1]))
        print("start normalizing.")
        log2_normalized_intensity_df,sample,sample_outdir = normalization(new_paired_df,outdir) 
        sample_list.append(sample)
        print(sample+":linear regression with %s"%(int(percent*100))+"% proteins of each sliding window.")        

        all_window_mean_value_list,all_window_variance_list,all_window_intercept_list = \
        MA_linear_regression(log2_normalized_intensity_df,sample_outdir,wsize=window_size,step=step_size,\
                      middle_percent = percent)  

        final_method = method

        if final_method == "natural_cubic_spline":
            print(sample+":start mean variance curve(MVC) fitting using natural cubic spline model.")
            nonlinear_model,R_square = fit_natural_model(log2_normalized_intensity_df,all_window_mean_value_list,\
                                all_window_variance_list,sample_outdir,used_for_model_fitting=[0,1])
            nonlinear_r2.append(R_square)
            if R_square > 0.5:
                log2_normalized_intensity_df["estimated variance"] =  \
                nonlinear_model.predict(log2_normalized_intensity_df["A value"].values)  
                for index in log2_normalized_intensity_df.index:
                    if log2_normalized_intensity_df.loc[index,"A value"]< all_window_mean_value_list[0]:
                        log2_normalized_intensity_df.loc[index,"estimated variance"] = all_window_variance_list[0]
                    elif log2_normalized_intensity_df.loc[index,"A value"] > all_window_mean_value_list[-1]:
                        log2_normalized_intensity_df.loc[index,"estimated variance"] = all_window_variance_list[-1]
                    else:
                        pass
            else:
                log2_normalized_intensity_df["estimated variance"] = \
                [np.mean(all_window_variance_list)] * len(log2_normalized_intensity_df)
                
        else:
            print(sample+":start mean variance curve(MVC) fitting based on exponential decay function.")
            popt,exp_R_square = fit_exp_model(log2_normalized_intensity_df,all_window_mean_value_list,all_window_variance_list,sample_outdir,used_for_model_fitting=[0,1])
            nonlinear_r2.append(exp_R_square)
            if len(popt) == 3:
                exp_function = func_nl_lsq
            else:
                exp_function = func_nl_lsp_no_constant_term
            
            if exp_R_square > 0.5:
                
                log2_normalized_intensity_df["estimated variance"] = exp_function(log2_normalized_intensity_df["A value"].values,*popt) 
                
                for index in log2_normalized_intensity_df.index:
                    
                    if log2_normalized_intensity_df.loc[index,"A value"]< all_window_mean_value_list[0]:
                        log2_normalized_intensity_df.loc[index,"estimated variance"] = exp_function(all_window_variance_list[0],*popt)
                        
                    elif log2_normalized_intensity_df.loc[index,"A value"] > all_window_mean_value_list[-1]:
                        log2_normalized_intensity_df.loc[index,"estimated variance"] = exp_function(all_window_variance_list[-1],*popt)
                        
                    else:
                        pass
            else:
                log2_normalized_intensity_df["estimated variance"] = \
                [np.mean(all_window_variance_list)] * len(log2_normalized_intensity_df)
                
        x_model = np.asarray(all_window_mean_value_list)
        y_model = np.asarray(all_window_intercept_list)
        print(sample+": fitting M-A curve to account for normalization biases based on natural cubic spline.")
        intercept_model = natural_cubic_model.get_natural_cubic_spline_model(x_model,y_model,minval=min(x_model),maxval=max(x_model),n_knots=4)
        intercept_rsquare = getIndexes(intercept_model.predict(x_model),y_model)
        intercept_R_square_list.append(intercept_rsquare)
        from scipy.stats import norm
        pvalue_list = []
        z_statiatic_list = []
        adjust_mvalue =[]
        fig = plt.figure(figsize=(10,5))
        ax1= fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
    
        ax1.scatter(log2_normalized_intensity_df["A value"],log2_normalized_intensity_df["M value"],alpha=0.2,s=8)
        min_value = min(log2_normalized_intensity_df["A value"])
        max_value = max(log2_normalized_intensity_df["A value"])

        ax1.plot([min_value-1,max_value+1],[0,0],c="black",linestyle="--",linewidth=0.5)
        ax1.scatter(all_window_mean_value_list,all_window_intercept_list,s=15,facecolors='none', edgecolors="tomato")
        ax1.plot(x_model,intercept_model.predict(x_model),"--",c="orange")
        ax1.set_xlabel("Log$_2$ intensity",size=15)
        ax1.set_ylabel("Log$_2$ ratio",size=15)  
        ax2.scatter(log2_normalized_intensity_df["A value"].values,
                [log2_normalized_intensity_df.loc[index,"M value"]\
                 -intercept_model.predict(log2_normalized_intensity_df.loc[index,"A value"])[0]\
                 for index in log2_normalized_intensity_df.index],alpha=0.2,s=8)
        ax2.plot([min_value-1,max_value+1],[0,0],c="black",linestyle='--',linewidth=0.5)
        ax2.set_xlabel(" Log$_2$ intensity",size=15)
        ax2.set_ylabel("Adjusted Log$_2$ ratio",size=15)
        for index in log2_normalized_intensity_df.index:
    
            mvalue = log2_normalized_intensity_df.loc[index,"M value"] \
            - intercept_model.predict(log2_normalized_intensity_df.loc[index,"A value"])[0]
            adjust_mvalue.append(mvalue)
            mvalue_abs = abs(mvalue)
        
            standard_variance = log2_normalized_intensity_df.loc[index,"estimated variance"]**0.5
    
            pvalue_list.append(2*norm.sf(mvalue_abs,loc=0,scale=standard_variance))
            z_statiatic_list.append(mvalue/standard_variance)
        print(sample+":calculate z-statistic.")
        log2_normalized_intensity_df["adjust mvalue"] = adjust_mvalue
        log2_normalized_intensity_df["pvalue"] = pvalue_list
        log2_normalized_intensity_df["z statistic"] = z_statiatic_list 
        log2_normalized_intensity_df.to_csv(sample_outdir+"map_output_results.txt",sep="\t")
        sample_results_list.append(sample_outdir+"map_output_results.txt")
        plt.savefig(sample_outdir+"mascatter.png",dpi=300)
        plt.close("all")

    a_df = pd.DataFrame([sample_list,nonlinear_r2],index=["sample","r^2_of_nonlinear_model"]).T
    a_df.to_csv(outdir+"/r2_of_nonlinear_model_fitting_for_estimated_variance.csv",sep="\t")

    u_r2_df = pd.DataFrame([sample_list,intercept_R_square_list],index=["sample","r^2_of_u_natural_cubic_model"]).T
    u_r2_df.to_csv(outdir+"/r2_of_natural_cubic_model_fitting_for_intercept_u.csv",sep="\t")

     
    print("summarize and output.")
    all_sample_df  = reverse_zmap_combine.combine_all_sample_results(sample_results_list,outdir,batch_df)
    z_sta_table = reverse_zmap_combine.get_z_sta_df(all_sample_df,outdir).copy(deep=True)
    new_z_sta_table = z_sta_table.copy()
    new_z_sta_table.columns = [column.split(" ")[0] for column in new_z_sta_table.columns]
    new_z_sta_table.to_csv(outdir+"/z_statistic_table.txt",sep="\t")
    
    final_z_sta_pvalue_df = reverse_zmap_combine.calculate_chi2_pvalue(new_z_sta_table,outdir)
   
    reverse_zmap_r2_distribution_of_linear_fitting.r2_linear_boxplot(outdir)
    reverse_zmap_r2_distribution_of_nonlinear_fitting.barplot_r2_estimated_variance(outdir)
    reverse_zmap_r2_distribution_of_nonlinear_fitting_for_intercept_u.barplot_for_u_fitting(outdir)
    print("All finished.")
    

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="")
    parser.add_argument("--intensity_data",help="protein intensity file",required=True)
    parser.add_argument("--sample_information",help="sample information",required=True)
    parser.add_argument("--outdir",help="outdir for result",required=True)
    parser.add_argument("--window_size",help="window size for linear regression",required=True)
    parser.add_argument("--step_size",help="step for sliding window",required=True)
    parser.add_argument("--percent",help="middle percent proteins used to linear fitting,30<=percent<=50 ",required=True)
    parser.add_argument("--method",help="exponential_function or natural_cubic_spline",required = True)
    argv=vars(parser.parse_args())
    intensity_data=argv['intensity_data']
    sample_information = argv['sample_information']
    outdir=argv['outdir']
    window_size=int(argv['window_size'])
    step_size=int(argv['step_size'])
    percent = int(argv['percent'])/100
    method = argv["method"]
    mkdir(outdir)
    import math
    reverse_zmap_(intensity_data,sample_information,outdir,window_size,step_size,percent,method)



    



    
    
    

    



















