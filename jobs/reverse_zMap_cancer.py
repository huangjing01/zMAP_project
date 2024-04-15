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
from .reverse_zmap_cancer_step1 import *
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

@celery.task(name="reverse_zMap_cancer",bind=True)
def one_reverse_zMap_cancer(self,file,window_size,step_size,percent,method):
    paired_df = pd.read_csv(file,sep="\t",index_col=0)
    os.remove(file)
    column_list = paired_df.columns
    sample = column_list[0]
    one_reverse_zMap_cancer_info = [sample + " vs reference"]
    self.update_state(state='PROGRESS',meta={'status': one_reverse_zMap_cancer_info})
    new_paired_df = paired_df.loc[(paired_df[column_list[0]]!=0) & (paired_df[column_list[1]]!=0)]
    sample_outdir = os.path.dirname(file) + '/' + column_list[0] + '/'
    mkdir(sample_outdir)
    one_reverse_zMap_cancer_info.append(sample + " : normalize...")
    self.update_state(state='PROGRESS',meta={'status': one_reverse_zMap_cancer_info})
    log2_normalized_intensity_df = normalization(new_paired_df,sample_outdir) 

    one_reverse_zMap_cancer_info.append(sample + ' : model with ' + str(percent * 100) + '% proteins of each sliding window...')
    self.update_state(state='PROGRESS',meta={'status': one_reverse_zMap_cancer_info})
    all_window_mean_value_list,all_window_variance_list,all_window_intercept_list = MA_linear_regression(sample_outdir,log2_normalized_intensity_df,wsize=window_size,step=step_size,middle_percent = percent)


    if method == "best":
        nonlinear_model,R_square = fit_natural_model(sample_outdir,log2_normalized_intensity_df,all_window_mean_value_list,all_window_variance_list,used_for_model_fitting=[0.05,0.95])
        one_reverse_zMap_cancer_info.append(sample + ' : fit mean variance curve(MVC) with natural cubic spline...')
        one_reverse_zMap_cancer_info.append(sample + ' : regression coefficient (' + str(R_square) + ')')
        self.update_state(state='PROGRESS',meta={'status': one_reverse_zMap_cancer_info})
        popt,exp_R_square = fit_exp_model(sample_outdir,log2_normalized_intensity_df,all_window_mean_value_list,all_window_variance_list,used_for_model_fitting=[0.05,0.95])
        one_reverse_zMap_cancer_info.append(sample + ' : fit mean variance curve(MVC) with exponential function...')
        one_reverse_zMap_cancer_info.append(sample + ' : regression coefficient (' + str(exp_R_square) + ')')
        self.update_state(state='PROGRESS',meta={'status': one_reverse_zMap_cancer_info})
        if R_square > exp_R_square:
            final_method="natural_cubic_spline"
            one_reverse_zMap_cancer_info.append(sample + ' : best fitted curve (natural cubic spline)')
            self.update_state(state='PROGRESS',meta={'status': one_reverse_zMap_cancer_info})
        else:
            final_method="exponential_function"
            one_reverse_zMap_cancer_info.append(sample + ' : best fitted curve (exponential function)')
            self.update_state(state='PROGRESS',meta={'status': one_reverse_zMap_cancer_info})

    elif method == "natural_cubic_spline":
        nonlinear_model,R_square = fit_natural_model(sample_outdir,log2_normalized_intensity_df,all_window_mean_value_list,all_window_variance_list,used_for_model_fitting=[0.05,0.95])
        one_reverse_zMap_cancer_info.append(sample + ' : fit mean variance curve(MVC) with natural cubic spline...')
        one_reverse_zMap_cancer_info.append(sample + ' : MVC regression coefficient (' + str(R_square) + ')')
        self.update_state(state='PROGRESS',meta={'status': one_reverse_zMap_cancer_info})
        final_method = method
    
    else:
        popt,exp_R_square = fit_exp_model(sample_outdir,log2_normalized_intensity_df,all_window_mean_value_list,all_window_variance_list,used_for_model_fitting=[0.05,0.95])
        one_reverse_zMap_cancer_info.append(sample + ' : fit mean variance curve(MVC) with exponential function...')
        one_reverse_zMap_cancer_info.append(sample + ' : MVC regression coefficient (' + str(R_square) + ')')
        self.update_state(state='PROGRESS',meta={'status': one_reverse_zMap_cancer_info})
        final_method = method


    one_reverse_zMap_cancer_info.append(sample + ' : calculate estimated variance...')
    self.update_state(state='PROGRESS',meta={'status': one_reverse_zMap_cancer_info})

    if final_method == "natural_cubic_spline":
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
    one_reverse_zMap_cancer_info.append(sample + ' : calculate normalization factor bias...')
    self.update_state(state='PROGRESS',meta={'status': one_reverse_zMap_cancer_info})
    x_model = np.asarray(all_window_mean_value_list)
    y_model = np.asarray(all_window_intercept_list)
    intercept_model = natural_cubic_model.get_natural_cubic_spline_model(x_model,y_model,minval=min(x_model),maxval=max(x_model),n_knots=4)
    intercept_rsquare = getIndexes(intercept_model.predict(x_model),y_model)
#intercept_R_square_list.append(intercept_rsquare)
    from scipy.stats import norm
    pvalue_list = []
    z_statiatic_list = []
    adjust_mvalue =[]
    fig = plt.figure(figsize=(10,5))
    ax1= fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.set_title("MA plot with normalization factors")
    ax2.set_title("MA plot after removing normalization factors")
    ax1.scatter(log2_normalized_intensity_df["A value"],log2_normalized_intensity_df["M value"],c="blue",alpha=0.2,s=8)
    ax1.plot([5,25],[0,0],c="grey")
    ax1.scatter(all_window_mean_value_list,all_window_intercept_list,s=15)
    ax1.plot(x_model,intercept_model.predict(x_model),c="red")
#ax1.plot([5,25],[0,0],c="grey")
    ax2.scatter(log2_normalized_intensity_df["A value"].values,
            [log2_normalized_intensity_df.loc[index,"M value"]\
             -intercept_model.predict(log2_normalized_intensity_df.loc[index,"A value"])[0]\
             for index in log2_normalized_intensity_df.index],c="blue",alpha=0.2,s=8)
    ax2.plot([5,25],[0,0],c="grey")

    one_reverse_zMap_cancer_info.append(sample + ' : calculate p value and z-transformed statistics...')
    self.update_state(state='PROGRESS',meta={'status': one_reverse_zMap_cancer_info})

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
    log2_normalized_intensity_df.to_csv(sample_outdir+"map_output_results.xls",sep="\t")
    plt.savefig(sample_outdir+"mascatter.png",dpi=300)
    plt.close("all")
    return {'status':one_reverse_zMap_cancer_info}


@celery.task(name='async_reverse_zMap_cancer',bind=True)
def async_reverse_zMap_cancer(self,filepath,window_size,step_size,percent,method,email):
    try:
        with time_limit(60*60*6):
            reverse_zMap_cancer_status = ["start processing..."]
            self.update_state(state='PROGRESS',meta={'status':reverse_zMap_cancer_status})
            files_dir = os.path.dirname(filepath)
            os.chdir(files_dir)
            outdir="result"
            os.mkdir(outdir)
            os.chdir(outdir)
            reverse_zMap_cancer_status.append("read data...")
            self.update_state(state='PROGRESS',meta={'status':reverse_zMap_cancer_status})
            raw_intensity_data = pd.ExcelFile(filepath)
            sheet_names = raw_intensity_data.sheet_names
            protein_intensity_df = raw_intensity_data.parse(sheet_names[0],index_col=0)
            batch_df = raw_intensity_data.parse(sheet_names[1],index_col=0)
            reverse_zMap_cancer_status.append("quality controlling...")
            self.update_state(state='PROGRESS',meta={'status':reverse_zMap_cancer_status})
            reverse_zmap_cancer_qc.data_qc(protein_intensity_df,fig_dir = os.getcwd())


            reverse_zMap_cancer_status.append("start one vs one comparison...")
            self.update_state(state='PROGRESS',meta={'status':reverse_zMap_cancer_status})
            all_one_reverse_zMap_cancer_task = []
            sample_results_list = []
            for index in batch_df.index:
                ref_sample = batch_df.loc[index,"ref_sample"]
                patient_samples = batch_df.loc[index,[i for i in batch_df.columns.tolist() if i != "ref_sample"]]
                for sample in patient_samples:
                    protein_intensity_df[[sample,ref_sample]].to_csv(sample+".xls",sep="\t")
                    one_file = os.path.join(os.getcwd(),sample+".xls")
                    all_one_reverse_zMap_cancer_task.append(one_reverse_zMap_cancer.apply_async(args=[one_file,window_size,step_size,percent,method]))
                    sample_results_list.append(os.path.join(files_dir,outdir,sample,"map_output_results.xls"))
                    print(os.path.join(files_dir,outdir,sample,"map_output_results.xls"))

            while len(all_one_reverse_zMap_cancer_task) != 0:
                one_reverse_zMap_cancer_task = all_one_reverse_zMap_cancer_task[random.randint(0,len(all_one_reverse_zMap_cancer_task)-1)]
                if one_reverse_zMap_cancer_task.state == 'SUCCESS':
                    for one_reverse_zMap_cancer_info in one_reverse_zMap_cancer_task.info.get('status',''):
                        if one_reverse_zMap_cancer_info not in reverse_zMap_cancer_status:
                            reverse_zMap_cancer_status.append(one_reverse_zMap_cancer_info)
                            self.update_state(state='PROGRESS',meta={'status':reverse_zMap_cancer_status})
                    all_one_reverse_zMap_cancer_task.remove(one_reverse_zMap_cancer_task)
                elif one_reverse_zMap_cancer_task.state == 'FAILURE':
                    all_one_reverse_zMap_cancer_task.remove(one_reverse_zMap_cancer_task)
                    revoke(one_reverse_zMap_cancer_task.task_id,terminate=True)
                    raise Exception("Unexpected error.")
                elif one_reverse_zMap_cancer_task.state != 'PENDING':
                    for one_reverse_zMap_cancer_info in one_reverse_zMap_cancer_task.info.get('status',''):
                        if one_reverse_zMap_cancer_info not in reverse_zMap_cancer_status:
                            reverse_zMap_cancer_status.append(one_reverse_zMap_cancer_info)
                            self.update_state(state='PROGRESS',meta={'status':reverse_zMap_cancer_status})
                else:
                    pass
            
            reverse_zMap_cancer_status.append("combine all sample results...")
            self.update_state(state='PROGRESS',meta={'status':reverse_zMap_cancer_status})
            all_sample_df  = reverse_zmap_cancer_combine.combine_all_sample_results(sample_results_list,os.path.join(files_dir,outdir))

            reverse_zMap_cancer_status.append("plot correlations between samples...")
            self.update_state(state='PROGRESS',meta={'status':reverse_zMap_cancer_status})
            z_sta_table = reverse_zmap_cancer_combine.get_z_sta_df(all_sample_df,os.path.join(files_dir,outdir)).copy(deep=True)
            new_z_sta_table = z_sta_table.copy()
            new_z_sta_table.columns = [column.split(" ")[0] for column in new_z_sta_table.columns]
            batch_color_dict = reverse_zmap_cancer_samplecluster.get_batch_color_dict(batch_df)
            color_list_df = reverse_zmap_cancer_samplecluster.get_color_list(new_z_sta_table.columns,batch_color_dict)
            pearsonr_cor = reverse_zmap_cancer_samplecluster.correlation(new_z_sta_table,method="pearsonr")
            reverse_zmap_cancer_samplecluster.cluster_map(new_z_sta_table,pearsonr_cor,color_list_df,vmin=-0.5,vmax=1,fig=os.path.join(files_dir,outdir,"clustermap"))

            reverse_zMap_cancer_status.append("draw pca plot of samples...")
            self.update_state(state='PROGRESS',meta={'status':reverse_zMap_cancer_status})
            pca,pc_df = reverse_zmap_cancer_pca.pca_analysis(pearsonr_cor)
            reverse_zMap_cancer_status.append('explained variance ratio (first three components): %s' % str(pca.explained_variance_ratio_))
            self.update_state(state='PROGRESS',meta={'status':reverse_zMap_cancer_status})
            reverse_zmap_cancer_pca.pc_scatter(pc_df,pca,batch_color_dict,os.path.join(files_dir,outdir))        
            final_z_sta_pvalue_df = reverse_zmap_cancer_combine.calculate_chi2_pvalue(z_sta_table,os.path.join(files_dir,outdir))

            
            file_zip = zip_ya(files_dir)
            if email != "":
                send_async_email.apply_async(args=[file_zip,email])
            reverse_zMap_cancer_status.append("reverse_zMAP.cancer finished")
            return {'status': reverse_zMap_cancer_status}
    except TimeoutException:
        reverse_zMap_cancer_status.append("Error : Excution time out of limit")
        for one_task in all_one_reverse_zMap_cancer_task:
            revoke(one_task.task_id,terminate=True)
        return {'status':reverse_zMap_cancer_status}
    except Exception as e:
        reverse_zMap_cancer_status.append("Error : " + str(e) + "please check your file and feel free to contact us if you have any question.")
        return {'status':reverse_zMap_cancer_status}