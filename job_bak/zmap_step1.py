# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 18:08:21 2019

@author: guixiuqi
E-mail:  guixiuqi@126.com
"""

import scipy
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import os
import math
import scipy.stats as stats
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import glob
import argparse
import pandas as pd
from . import natural_cubic_model

####################log2 normalization##################

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

######################## mkdir ################################

def mkdir(path_dir):
    isexists = os.path.exists(path_dir)
    if not isexists:
        os.makedirs(path_dir)
    else:
        pass
		
########### calculate mean log2-intensity among samples and observed variance#####
def calculate_a_variance_value(log2_normalized_intensity_df):
    #log2_normalized_intensity_df = pd.read_csv(new_table_txt,sep="\t", index_col = 'Gene Name')
    sample_num = len(log2_normalized_intensity_df.columns)
    A_value_list = log2_normalized_intensity_df.mean(axis=1).values
    variance_list = log2_normalized_intensity_df.var(axis=1).values
    log2_normalized_intensity_df["A-value"] = A_value_list
    log2_normalized_intensity_df["Observed Variance"] = variance_list
    return log2_normalized_intensity_df,sample_num

################## sort log2-intensity table by A-value###########################
def get_all_sliding_region(table):
    protein_num = len(table)
    sorted_table = table.sort_values(by = "A-value")
 
    return sorted_table


################################## Michael and Schucany ##########################
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

#####################linear regression for sliding window#########################
def window_vaplot(sorted_table,num_of_sample,window_size = 400,sliding_size = 100,result_dir="",top_gene_used_to_fit_linear_model = 0.3,batch=""):
    
    sliding_region_left = 0
    sliding_region_right = len(sorted_table)
    
  #  vaplot_dir = result_dir + "/VAplot"
  #  linear_regression = result_dir + "/linear_model"
    mkdir(result_dir)
    #mkdir(vaplot_dir)
    #mkdir(linear_regression)
    window_mean_A_list = []
    window_slope_list = []
    r2_list = []
    left = sliding_region_left
    while left < sliding_region_right:
        if left + window_size > sliding_region_right:
            non_target_window_df = pd.concat([sorted_table.iloc[0:left,:],sorted_table.iloc[sliding_region_right:,:]],axis=0)
        else:
            non_target_window_df = pd.concat([sorted_table.iloc[0:left,:],sorted_table.iloc[left+window_size:,:]],axis=0)
        
        target_window_df = sorted_table.iloc[left:min(left+window_size,sliding_region_right),:]
        
        target_window_df_variance = target_window_df["Observed Variance"]
        
        gene_number_in_window = len(target_window_df_variance)
        
        if gene_number_in_window > 100:
            
            top_percent_num = int(gene_number_in_window * top_gene_used_to_fit_linear_model)
        
            sorted_target_window_df_variance = target_window_df_variance.sort_values(ascending=True)
        
            gene_used_to_fit_curve = sorted_target_window_df_variance.index[0:top_percent_num]
        
            other_genes = sorted_target_window_df_variance.index[top_percent_num:]
        
            target_window_bottom_gene_df = target_window_df.loc[gene_used_to_fit_curve]
        
            target_window_up_gene_df = target_window_df.loc[other_genes]
        
            window_mean_A = np.mean(target_window_df["A-value"])
        
            window_mean_A_list.append(window_mean_A)
        
            degree_of_freedom = num_of_sample - 1
            
            pp_list = MS_plotting_position(1,gene_number_in_window+1,gene_number_in_window,0.375)

            theoretical_quantiles = [scipy.stats.chi2.ppf(j,degree_of_freedom)/degree_of_freedom for j in pp_list] 
        
            p = np.polyfit(theoretical_quantiles[0:top_percent_num],sorted_target_window_df_variance[0:top_percent_num],1)
            regr = LinearRegression()
            model_x = np.asarray(theoretical_quantiles[0:top_percent_num]).reshape(len(theoretical_quantiles[0:top_percent_num]),1)
            model_y = np.asarray(sorted_target_window_df_variance[0:top_percent_num]).reshape(len(sorted_target_window_df_variance[0:top_percent_num]),1)
            regr.fit(model_x,model_y)
            r2 = regr.score(model_x, model_y)
            r2_list.append(r2)
            window_slope_list.append(p[0])
            if p[1] > 0:
                p_function = "y= %s x + %s, r-square = %s" %(str(p[0]), str(p[1]), str(r2))
            elif p[1] < 0:
                p_function = "y= %s x - %s, r-square = %s" %(str(p[0]), str(-p[1]), str(r2))
            else:
                p_function = "y= %s x, r-square = %s" %(str(p[0]), str(r2))
        left += sliding_size
    mean_a_scope_df = pd.DataFrame([window_mean_A_list,window_slope_list],index = ["mean_A_value","estimated_variance"]).T
    mean_a_scope_df.to_csv(result_dir + "/"+"%s_mean_A_estimated_variance.xls"%(batch),sep ="\t")
    plt.figure(figsize= (5,6))
    plt.boxplot(r2_list)
    plt.ylabel("$R^2$")
    plt.title("Boxplot for %s $R^2$ \nof linear regression"%(batch),size=15)
    
    
    plt.savefig(result_dir + "/"+ batch + "_"+result_dir.split("/")[-2]+"_boxplot_for_r2_of_linear_regression.png",dpi=300,bbox_inches="tight")
    
    plt.savefig(result_dir + "/"+ batch + "_"+result_dir.split("/")[-2]+"_boxplot_for_r2_of_linear_regression.pdf",dpi=300,bbox_inches="tight")
    plt.close("all")
    return window_mean_A_list,window_slope_list


#############calculate R^2###################
def getIndexes(y_predict, y_data):
    
    n = y_data.size
    
    SSE = ((y_data-y_predict)**2).sum()
    
    MSE = SSE/n
    
    RMSE = np.sqrt(MSE)
    
    u = y_data.mean()
    SST = ((y_data -u)**2).sum()
    SSR = SST -SSE
    R_square = SSR/SST
    
    return SSE,MSE,RMSE,R_square


###########expotential_function##########
def func_nl_lsq(x,a,b,c):
    return np.exp(a + (b*x))+ c
def func_nl_lsp_no_constant_term(x,a,b):
    return np.exp(a + (b*x))


################################fit non_linear model###########################################
def fit_exp_model(table,mean_a_list,slope_list,used_for_model_fitting=[0.05,0.95],result_dir="",batch=""):
    
    #nonlinear_dir = result_dir+"/nonlinear_model"
    #mkdir(nonlinear_dir)
    
    a_value = np.asarray(table["A-value"].sort_values(ascending=True))
    
    start = int(len(mean_a_list)*used_for_model_fitting[0])
    end = int(len(mean_a_list)*used_for_model_fitting[1])
    x_model = np.asarray(mean_a_list[start:end])
    y_model = np.asarray(slope_list[start:end])
    #popt, pcov = curve_fit(func_nl_lsq, x_model, y_model,bounds=([-1,-5,-1], [2, 0, 0.1]))
    popt, pcov = curve_fit(func_nl_lsq, x_model, y_model,bounds=([-50,-50,-1], [50, 50, 1]))
    if popt[2] > 0:
        SSE,MSE,RMSE,R_square = getIndexes(func_nl_lsq(x_model,*popt),y_model)
        plt.figure(figsize=(5,5))
        plt.scatter(mean_a_list[0:start]+mean_a_list[end:],slope_list[0:start]+slope_list[end:],
                    facecolors="none",s = 50,edgecolors = "red",label = "Not used for model fitting" )
        plt.scatter(x_model, y_model,s=50,facecolors='none', edgecolors="green",label = "Used for model fitting" )
        plt.plot(a_value, func_nl_lsq(a_value, *popt), '--',linewidth = 3,color = "darkorange",
             label = "Fitted exponential funcation:\n$\sigma^2=exp(%s-%s*S)+%s$"%(str("%.3f"%(popt[0])),
                                                                   str("%.3f"%(-popt[1])),str("%.3f"%(popt[2]))))
    else:
        popt, pcov = curve_fit(func_nl_lsq, x_model, y_model,bounds=([-50,-50], [50, 50]))
        SSE,MSE,RMSE,R_square = getIndexes(func_nl_lsp_no_constant_term(x_model,*popt),y_model)
        plt.figure(figsize=(5,5))
        plt.scatter(mean_a_list[0:start]+mean_a_list[end:],slope_list[0:start]+slope_list[end:],
                    facecolors="none",s = 50,edgecolors = "red",label = "Not used for model fitting" )
        plt.scatter(x_model, y_model,s=50,facecolors='none', edgecolors="green",label = "Used for model fitting" )
        plt.plot(a_value, func_nl_lsp_no_constant_term(a_value, *popt), '--',linewidth = 3,color = "darkorange",
             label = "Fitted exponential funcation:\n$\sigma^2=exp(%s-%s*S)$"%(str("%.3f"%(popt[0])),str("%.3f"%(-popt[1]))))
        
    plt.xlabel("$Average\ log_2-intensity\ S$",size=15)
    plt.ylabel("$Estimated\ variance\ \sigma^2$",size=15)
    plt.title("%s $R^2$=%s"%(batch,str("%.3f"%(R_square))),fontsize = 15)
    plt.legend()
    plt.savefig(result_dir + r"/"+"%s_nonlinear_exp_model.pdf"%(batch),dpi=300,bbox_inches='tight')
    plt.savefig(result_dir + r"/"+"%s_nonlinear_exp_model.png"%(batch),dpi=300,bbox_inches='tight')
    plt.close("all")
    return popt,R_square

#########################################fit natural cubic spline model##########################
def fit_natural_cubic_spline_model(table,mean_a_list,slope_list,used_for_model_fitting=[0.05,0.95],result_dir="",batch=""):
   # nonlinear_dir = result_dir+"/nonlinear_model"
   # mkdir(nonlinear_dir)
    a_value = np.asarray(table["A-value"].sort_values(ascending=True))
    start = int(len(mean_a_list)*used_for_model_fitting[0])
    end = int(len(mean_a_list)*used_for_model_fitting[1])
    x_model = np.asarray(mean_a_list[start:end])
    y_model = np.asarray(slope_list[start:end])
    nonlinear_model = natural_cubic_model.get_natural_cubic_spline_model(x_model,y_model,minval=min(x_model),maxval=max(x_model),n_knots=4)

    R_square = getIndexes(nonlinear_model.predict(x_model),y_model)[-1]
 #   print (R_square)
#    print (type(R_square))
    plt.figure(figsize=(5,5))
    plt.scatter(mean_a_list[0:start]+mean_a_list[end:],slope_list[0:start]+slope_list[end:],
                    facecolors="none",s = 50,edgecolors = "red",label = "Not used for model fitting" )
    plt.scatter(x_model, y_model,s=50,facecolors='none', edgecolors="green",label = "Used for model fitting" )

    plt.plot(a_value,nonlinear_model.predict(a_value),c="darkorange",linewidth=3,label = "Fitted natural cubic spline")
        
    plt.xlabel("Average log$_2$-intensity",size=15)
    plt.ylabel("$Estimated\ variance\ \sigma^2$",size=15)
    plt.title("%s $R^2$=%s"%(batch,str("%.3f"%(R_square))),fontsize = 15)
    plt.legend()
    plt.savefig(result_dir + r"/"+"%s_nonlinear_natural_cubic_spline_model.pdf"%(batch),dpi=300,bbox_inches='tight')
    plt.savefig(result_dir + r"/"+"%s_nonlinear_natural_cubic_spline_model.png"%(batch),dpi=300,bbox_inches='tight')
    plt.close("all")
    return nonlinear_model,R_square


    
#########################################calculate chi^2 pvalue based on expotential_function model###################################
def calculate_exp_chi2(sorted_table,popt,result_dir="",num_of_sample=4,batch=""):
    
    sample_columns = sorted_table.columns[0:num_of_sample]
    for sample_column in sample_columns:
        sorted_table["z-score of " + sample_column] = (sorted_table[sample_column]-sorted_table["A-value"])/sorted_table["Observed Variance"]**0.5
    if len(popt)==3:
        exp_function = func_nl_lsq
    else:
        exp_function = func_nl_lsp_no_constant_term
    sorted_table["Estimated Variance"] = [exp_function(mean_value,*popt) for mean_value in sorted_table["A-value"].values]
    for sample_column in sample_columns:
        sorted_table["z-transformed of " + sample_column] = (sorted_table[sample_column]-sorted_table["A-value"])/sorted_table["Estimated Variance"]**0.5
    
    sorted_table["Chi^2-statistic"] = sorted_table["Observed Variance"]*(num_of_sample-1)/sorted_table["Estimated Variance"]
    sorted_table["pvalue"] = [scipy.stats.chi2.sf(x,num_of_sample-1) for x in sorted_table["Chi^2-statistic"].values]
    sorted_by_pvalue_table = sorted_table.sort_values(by = "pvalue")
    sorted_by_pvalue_table["BH-corrected P-value"] =[i if i<1 else 1 for i in sorted_by_pvalue_table["pvalue"]/sorted_by_pvalue_table["pvalue"].rank(axis=0)*len(sorted_by_pvalue_table)]
    sorted_by_pvalue_table.to_csv(result_dir+"/"+"%s_zmap_result_output.xls"%(batch),sep="\t")
    return sorted_by_pvalue_table

##########################################calculate chi^2 pvalue based on natural cubic spline model##################################
def calculate_natural_chi2(sorted_table,natural_cubic_spline_model,result_dir="",num_of_sample=4,batch=""):

    sample_columns = sorted_table.columns[0:num_of_sample]

    for sample_column in sample_columns:
        sorted_table["z-score of " + sample_column] = (sorted_table[sample_column]-sorted_table["A-value"])/sorted_table["Observed Variance"]**0.5

    sorted_table["Estimated Variance"] = natural_cubic_spline_model.predict(sorted_table["A-value"].values)

    for sample_column in sample_columns:
        sorted_table["z-transformed of " + sample_column] = (sorted_table[sample_column]-sorted_table["A-value"])/sorted_table["Estimated Variance"]**0.5

    sorted_table["Chi^2-statistic"] = sorted_table["Observed Variance"]*(num_of_sample-1)/sorted_table["Estimated Variance"]
    sorted_table["pvalue"] = [scipy.stats.chi2.sf(x,num_of_sample-1) for x in sorted_table["Chi^2-statistic"].values]
    sorted_by_pvalue_table = sorted_table.sort_values(by = "pvalue")
    sorted_by_pvalue_table["BH-corrected P-value"] =[i if i<1 else 1 for i in sorted_by_pvalue_table["pvalue"]/sorted_by_pvalue_table["pvalue"].rank(axis=0)*len(sorted_by_pvalue_table)]
    sorted_by_pvalue_table.to_csv(result_dir+"/"+"%s_zmap_result_output.xls"%(batch),sep="\t")
    return sorted_by_pvalue_table
	

	

###################################scatterplot of theoretical quantiles###############################
def theoretical_quantiles_scatterplot(zmap_result_df,result_dir = "",batch=""):
    z_score_columns=[]
    z_transformed_columns=[]
    for column in zmap_result_df.columns:
        if "z-score" in column:
            z_score_columns.append(column)
        elif "z-trans" in column:
            z_transformed_columns.append(column)
    
    plt.figure(figsize=(8,10))
    plt.subplot(211) 
    plt.title(batch,size=20)
    for z_score_column in z_score_columns:
        z_score = zmap_result_df[z_score_column].sort_values(ascending=True)
        #z_score_quantiles_percent = [stats.percentileofscore(z_score,x) for x in z_score]
        z_score_quantiles_percent = MS_plotting_position(1,len(z_score)+1,len(z_score),0.375)
        theoretical_quantiles = [scipy.stats.norm.ppf(j) for j in z_score_quantiles_percent]
        plt.scatter(theoretical_quantiles,z_score,label=z_score_column)
    plt.xlabel("Theoretical quantiles of N(0,1)",size=15)
    plt.ylabel("Z-score of $log_2$-intensity",size=15)
    plt.plot(theoretical_quantiles,theoretical_quantiles,'k--',label="y = x")
    plt.title(batch,size=20)
    plt.legend()
    plt.subplot(212)
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

def combine_result(zmap_output_results,sheet_names,sample_num=4,out_dir=""):
    results_df_list=[]
    batch_id_list = []
    i=0
    for result in zmap_output_results:
        sample_id_list = []
        result_df = pd.read_csv(result,sep="\t",index_col = 0)
        batch_id = sheet_names[i]
        i+=1
        batch_id_list.append(batch_id)
        sample_id_list = [i for i in result_df.columns[0:sample_num]]
        result_df.columns = [np.array([batch_id]*len(result_df.columns)),result_df.columns]
        results_df_list.append(result_df)
    all_df = pd.concat(results_df_list,axis=1,sort=False)

    summary_df = pd.DataFrame(index = all_df.index)
    for sample_id in  sample_id_list:
        summary_df["mean " + sample_id] = all_df.loc(axis=1)[:,"z-transformed of " + sample_id].mean(axis=1)
    summary_df["Best pvalue"] = all_df.loc(axis=1)[:,"pvalue"].min(axis=1)

    second_pvalue_list=[]
    for index in all_df.index:
        second_pvalue_list.append(all_df.loc(axis=1)[:,"pvalue"].loc[index].sort_values()[1])

    summary_df["Second best pvalue"] = second_pvalue_list
    summary_df["Average Chi^2-statistic"] = all_df.loc(axis=1)[:,"Chi^2-statistic"].mean(axis=1)

    N = sample_num

    summary_df["Number of detected"] = all_df.loc(axis=1)[:,"pvalue"].count(axis=1)
    summary_df["Average pvalue"]= [scipy.stats.gamma.sf(summary_df["Average Chi^2-statistic"].loc[index],
                       summary_df["Number of detected"].loc[index]*(N-1)/2,
                       scale = 2/summary_df["Number of detected"].loc[index]) for index in summary_df.index]
    all_sort_df = summary_df.sort_values(by = "Average pvalue")
    all_sort_df["BH-corrected average pvalue"] =[i if i<1 else 1 for i in \
               all_sort_df["Average pvalue"]/all_sort_df["Average pvalue"].rank(axis=0)*len(all_sort_df)]

    all_sort_df.columns = [np.array(["-"]*len(all_sort_df.columns)),all_sort_df.columns]
    all_batch_final_results = pd.concat([all_sort_df,all_df],axis=1,sort = False)

    all_batch_final_results.to_csv(out_dir+"/"+"all_batch_final_results.xls",sep="\t")

def pvalue_scatterplot(table,result_dir="",batch=""):
	
	#table = pd.read_csv(table,sep="\t",index_col=0)
    sort_table = table.sort_values(by = "Observed Variance",ascending=True)
    plt.figure(figsize=(6,6))
    plt.scatter(sort_table["A-value"],sort_table["Observed Variance"],s=20,c = -np.log10(sort_table["pvalue"].values),edgecolors="black",cmap="Reds",vmin=0, vmax=7,linewidth=0.2,alpha=0.8)
    plt.xlabel("Mean of $log_2-intensities$",size = 12)
    plt.ylabel("Variance of $log_2-intensities$" , size = 12)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=10)
    cbar.solids.set_edgecolor("face")
    cbar.set_label('$-log_{10}{(P-value)}$',size=10)
    plt.title(batch,size=20)
    #all_hemoglobin_protein = ["HBE1","HBM","HBB","HBD","HBQ1","HBA1","AHSP","HBG1","HBG2","HBZ"]
    #for gene in all_hemoglobin_protein:
    #    plt.annotate(s=gene,xy=(table.loc[gene]["A-value"],table.loc[gene]["Observed Variance"]),
    #                 xytext = (table.loc[gene]["A-value"]+1,table.loc[gene]["Observed Variance"]+0.5),
    #                 arrowprops={'arrowstyle':'fancy'})
    plt.savefig(result_dir+"/"+"%s_pvalue.pdf"%(batch),dpi=300,bbox_inches="tight")
    plt.savefig(result_dir+"/"+"%s_pvalue.png"%(batch),dpi=300,bbox_inches="tight")
    plt.close("all")

    
    
    
    
    
    
