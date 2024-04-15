from __future__ import with_statement
import signal
from contextlib import contextmanager

class TimeoutException(Exception): pass

@contextmanager
def time_limit(seconds):
    def signal_handler(signum, frame):
        raise TimeoutException
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)



from app import celery
import os
from flask import current_app
import zipfile
from .reverse_zmap_step1_for_linux_for_xlsx_1_5_2023_11_15 import *
import numpy as np
import json,random

import smtplib,mimetypes

from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email import encoders
from email.utils import formataddr
from jobs.zMap import send_async_email
from celery.task.control import revoke
def zip_ya(start_dir):
    start_dir = start_dir #compress directory path
    file_news = start_dir + '.zip'

    z = zipfile.ZipFile(file_news, 'w', zipfile.ZIP_DEFLATED)
    for dir_path, dir_names, file_names in os.walk(start_dir):
        f_path = dir_path.replace(start_dir, '')  # from current directory start to copy
        f_path = f_path and f_path + os.sep or '' # 
        for filename in file_names:
            z.write(os.path.join(dir_path, filename), f_path + filename)
    z.close()
    return file_news


@celery.task(name='async_reverse_zMap',bind=True)
def async_reverse_zMap(self,intensity_file,sample_info_file,window_size,step_size,percent,method,email):
    try:
        with time_limit(60*60*6):
            reverse_zMap_status = ["start processing..."]
            self.update_state(state='PROGRESS',meta={'status':reverse_zMap_status})
            files_dir = os.path.dirname(intensity_file)
            os.chdir(files_dir)
            
            outdir = os.path.join(files_dir,"result")
            os.mkdir(outdir)
            reverse_zMap_status.append("read data...")
            self.update_state(state='PROGRESS',meta={'status':reverse_zMap_status})
            protein_intensity_df = pd.read_csv(intensity_file,index_col=0,sep="\t")    
            batch_df = pd.read_csv(sample_info_file,sep='\t',index_col=0)

            final_protein_intensity_df = check_df_names(protein_intensity_df,batch_df,intensity_file)

            ms_internal_dict = {}
            for index in batch_df.index:
                sample_type = batch_df.loc[index,'internal_ref']
                ms_run_ = batch_df.loc[index,'MS_run']
                if sample_type=="Yes":
                
                    ms_internal_dict[ms_run_]=index

            reverse_zMap_status.append("summarize the input data")
            self.update_state(state='PROGRESS',meta={'status':reverse_zMap_status})
            reverse_zmap_qc.data_qc(final_protein_intensity_df,fig_dir = os.path.join(files_dir,outdir))

    
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
            reverse_zMap_status.append("start pairwise comparison to reference sample.")
        
            for paired_df in paired_df_list:
                column_list = paired_df.columns
                reverse_zMap_status.append([column_list[0] + " vs reference",column_list[0] + " : normalize..."])
                self.update_state(state='PROGRESS',meta={'status': reverse_zMap_status})

                new_paired_df = paired_df.dropna(how="any",axis=0)
                log2_normalized_intensity_df,sample,sample_outdir = normalization(new_paired_df,outdir)

                reverse_zMap_status.append(sample + ':linear regression with ' + str(percent * 100) + '% proteins of each sliding window...')
                self.update_state(state='PROGRESS',meta={'status': reverse_zMap_status})

                sample_list.append(sample)
                all_window_mean_value_list,all_window_variance_list,all_window_intercept_list = \
                MA_linear_regression(log2_normalized_intensity_df,sample_outdir,wsize=window_size,step=step_size,\
                              middle_percent = percent)  
                

                final_method = method

                if final_method == "natural_cubic_spline":
                    nonlinear_model,R_square = fit_natural_model(log2_normalized_intensity_df,all_window_mean_value_list,\
                                        all_window_variance_list,sample_outdir,used_for_model_fitting=[0,1])
                    reverse_zMap_status.append(sample + ' : fit mean variance curve(MVC) with natural cubic spline...')
                    reverse_zMap_status.append(sample + ' : regression coefficient (' + str(R_square) + ')')
                    self.update_state(state='PROGRESS',meta={'status': reverse_zMap_status})
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
                    popt,exp_R_square = fit_exp_model(log2_normalized_intensity_df,all_window_mean_value_list,all_window_variance_list,sample_outdir,used_for_model_fitting=[0,1])
                    reverse_zMap_status.append(sample + ' : fit mean variance curve(MVC) with exponential function...')
                    reverse_zMap_status.append(sample + ' : regression coefficient (' + str(exp_R_square) + ')')
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
                reverse_zMap_status.append(sample+": fitting M-A curve to account for normalization biases based on natural cubic spline.")
                self.update_state(state='PROGRESS',meta={'status': reverse_zMap_status})

                x_model = np.asarray(all_window_mean_value_list)
                y_model = np.asarray(all_window_intercept_list)
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

                reverse_zMap_status.append(sample + ' : calculate p value and z-transformed statistics...')
                self.update_state(state='PROGRESS',meta={'status': reverse_zMap_status})

            
                for index in log2_normalized_intensity_df.index:
            
                    mvalue = log2_normalized_intensity_df.loc[index,"M value"] \
                    - intercept_model.predict(log2_normalized_intensity_df.loc[index,"A value"])[0]
                    adjust_mvalue.append(mvalue)
                    mvalue_abs = abs(mvalue)
                
                    standard_variance = log2_normalized_intensity_df.loc[index,"estimated variance"]**0.5
            
                    pvalue_list.append(2*norm.sf(mvalue_abs,loc=0,scale=standard_variance))
                    z_statiatic_list.append(mvalue/standard_variance)
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

             
            reverse_zMap_status.append("summarize and output...")
            self.update_state(state='PROGRESS',meta={'status':reverse_zMap_status})
            all_sample_df  = reverse_zmap_combine.combine_all_sample_results(sample_results_list,outdir,batch_df)

            z_sta_table = reverse_zmap_combine.get_z_sta_df(all_sample_df,outdir).copy(deep=True)
            new_z_sta_table = z_sta_table.copy()
            new_z_sta_table.columns = [column.split(" ")[0] for column in new_z_sta_table.columns]
            new_z_sta_table.to_csv(outdir+"/z_statistic_table.txt",sep="\t")
            
            final_z_sta_pvalue_df = reverse_zmap_combine.calculate_chi2_pvalue(new_z_sta_table,outdir)
            reverse_zmap_r2_distribution_of_linear_fitting.r2_linear_boxplot(outdir)
            reverse_zmap_r2_distribution_of_nonlinear_fitting.barplot_r2_estimated_variance(outdir)
            reverse_zmap_r2_distribution_of_nonlinear_fitting_for_intercept_u.barplot_for_u_fitting(outdir)

            file_zip = zip_ya(files_dir)
            if email != "":
                send_async_email.apply_async(args=[file_zip,email])
            reverse_zMap_status.append("reverse_zMap finished")
            self.update_state(state='PROGRESS',meta={'status': reverse_zMap_status})
            return {'status': reverse_zMap_status}
    except TimeoutException:
        reverse_zMap_status.append("Error : Excution time out of limit")
        # for one_task in all_one_reverse_zMap_task:
        #     revoke(one_task.task_id,terminate=True)
        return {'status':reverse_zMap_status}
    except Exception as e:
        reverse_zMap_status.append("Error : " + str(e) + "please check your file and feel free to contact us if you have any question.")
        return {'status':reverse_zMap_status}




