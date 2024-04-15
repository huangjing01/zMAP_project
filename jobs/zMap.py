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
from .zmap_step1_quantity_anaysis_2023_9_14 import *
import numpy as np
import json,random

import smtplib,mimetypes

from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email import encoders
from email.utils import formataddr
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






@celery.task(name='send_async_email')
def send_async_email(filepath,email):

    message = MIMEMultipart()
    sender = "bioinfo1519@163.com"
    message['From'] = formataddr(['RSG',"bioinfo1519@163.com"])
    message['To'] = email
    message['Subject'] = "[zMAP result]"
    message['Passwd'] = '###'
    os.chdir(filepath[:-4])
    if os.path.isfile("zMAP_input.info.log"):
        content = """<p>hello zMAP result, the accessible url is """ + """<a>http://bioinfo.sibs.ac.cn/shaolab/zMAP/zMAP/result/""" + os.path.split(filepath[:-4])[-1] + """</a></p>"""
    elif os.path.isfile("reverse_zMAP_input.info.log"):
        content = """<p>hello reverse_zMAP result, the accessible url is """ + """<a>http://bioinfo.sibs.ac.cn/shaolab/zMAP/reverse_zMAP/result/""" + os.path.split(filepath[:-4])[-1] + """</a></p>"""
    elif os.path.isfile("reverse_zMAP.cancer_input.info.log"):
        content = """<p>hello reverse_zMAP.cancer result, the accessible url is """ + """<a>http://bioinfo.sibs.ac.cn/shaolab/zMAP/reverse-zMAP.CancerPro/result/""" + os.path.split(filepath[:-4])[-1] + """</a></p>"""
    else:
        pass
    msg = MIMEText(content,'html','utf-8')
    message.attach(msg)
    if os.path.getsize(filepath)/float(1024*1024) < 50:    
        ctype, encoding = mimetypes.guess_type(filepath)
        if ctype is None or encoding is not None:
            ctype = 'application/octet-stream'
        maintype, subtype = ctype.split('/', 1)
        file_msg = MIMEBase(maintype, subtype)


        file = open(filepath,'rb')
        file_msg.set_payload(file.read())
        file.close()
        encoders.encode_base64(file_msg)
        file_msg.add_header('Content-Disposition','attachment',filename=os.path.split(filepath)[-1])
        message.attach(file_msg)
    else:
        pass
    try:
        s = smtplib.SMTP_SSL("smtp.163.com", 465)
        s.login(sender, message['Passwd'])
        s.sendmail(sender, message['To'], message.as_string())
    except s.SMTPException as e:
        print(e)
    finally:
        s.quit()




@celery.task(name='async_zMap',bind=True)
def async_zMap(self,intensity_file,sample_info_file,window_size,step_size,percent,method,email):
    try:
        with time_limit(60*60*6):
            
            zMap_status = ["start zMAP processing..."]
            self.update_state(state='PROGRESS',meta={'status': zMap_status})
            files_dir = os.path.dirname(intensity_file)
            os.chdir(files_dir)

            
            zMap_status.append("read data...")
            self.update_state(state='PROGRESS',meta={'status': zMap_status})
            outdir = 'result'
            mkdir(outdir)
            raw_intensity_data = pd.read_csv(intensity_file,sep="\t",index_col=0)
            sample_info_df = pd.read_csv(sample_info_file,sep="\t",index_col=0)
            ms_run_list = list(set(list(sample_info_df['MS_run'].values)))
    
            all_rep_zmap_result_xls_list =[]
            
            for ms_run in ms_run_list:
                batch = ms_run
                output_dir = outdir +"/"+batch
                mkdir(output_dir)
                all_rep_zmap_result_xls_list.append(output_dir+"/"+"%s_zmap_result_output.txt"%(batch))
                # """
                zMap_status.append(batch + ' : read data...')
                self.update_state(state='PROGRESS',meta={'status': zMap_status})
                sample_ = list(sample_info_df.loc[sample_info_df['MS_run']==batch].index)
                raw_intensity_df = raw_intensity_data[sample_]
                raw_intensity_df = raw_intensity_df.dropna(how="any")

                zMap_status.append(batch + ' : normalize...')
                self.update_state(state='PROGRESS',meta={'status': zMap_status})

                raw_intensity_df = remove_bad_gene(raw_intensity_df,cutoff=10)
                normalized_df = normalization(raw_intensity_df,output_dir,l=1.5,batch=batch)

                zMap_status.append(batch + ' : model with ' + str(percent * 100) + '% proteins of each sliding window...')
                self.update_state(state='PROGRESS',meta={'status': zMap_status})

                table_df,sample_num = calculate_a_variance_value(normalized_df)

                sorted_table = get_all_sliding_region(table_df)

            
                mean_a_list,slope_list = window_vaplot(sorted_table,sample_num,window_size = window_size,
                                                   sliding_size = step_size,output_dir = output_dir,
                                                   top_gene_used_to_fit_linear_model = percent,batch=batch)

                if method == "exponential_function":
                    popt, exp_R_square = fit_exp_model(sorted_table,mean_a_list,slope_list,
                                 used_for_model_fitting=[0.05,0.95],result_dir=output_dir,batch=batch)
                    zMap_status.append(batch + ' : fit mean variance curve(MVC) with exponential function...')
                    zMap_status.append(batch + ' : MVC regression coeffeicient (' + str(exp_R_square) + ')')
                    self.update_state(state='PROGRESS',meta={'status': zMap_status})

                    zmap_result_df = calculate_exp_chi2(sorted_table,popt,result_dir = output_dir,num_of_sample=sample_num,batch=batch)

                else:# method == "natural_cubic_spline":
                    natural_cubic_spline_model,ncs_R_square = fit_natural_cubic_spline_model(sorted_table,mean_a_list,slope_list,used_for_model_fitting=[0.05,0.95],result_dir=output_dir,batch=batch)
                    zMap_status.append(batch + ' : fit mean variance curve(MVC) with natural cubic spline...')
                    zMap_status.append(batch + ' : MVC regression coeffeicient (' + str(ncs_R_square) + ')')
                    self.update_state(state='PROGRESS',meta={'status': zMap_status})

                    zmap_result_df = calculate_natural_chi2(sorted_table,natural_cubic_spline_model,result_dir = output_dir,num_of_sample = sample_num,batch=batch)

                zMap_status.append(batch + ' : calculate p value and z-transformed statistics...')
                self.update_state(state='PROGRESS',meta={'status': zMap_status})
                    
                zMap_status.append(batch + ' : VA-plot of all proteins with p-values...')
                self.update_state(state='PROGRESS',meta={'status': zMap_status})
                pvalue_scatterplot(zmap_result_df,result_dir=output_dir,batch=batch)
                
                zMap_status.append(batch + ' : plot modeling performance...')
                self.update_state(state='PROGRESS',meta={'status': zMap_status})
                theoretical_quantiles_scatterplot(zmap_result_df,result_dir = output_dir,batch=batch)
                
             
            combine_result(all_rep_zmap_result_xls_list,sample_info_df,ms_run_list,out_dir=outdir)
            zMap_status.append("combine completed")
            self.update_state(state='PROGRESS',meta={'status': zMap_status})
           
            file_zip = zip_ya(files_dir)
            if email != "":
                send_async_email.apply_async(args=[file_zip,email])
            zMap_status.append("zMAP finished...")
            self.update_state(state='PROGRESS',meta={'status': zMap_status})

            return {'status': zMap_status}
    except TimeoutException:
        zMap_status.append("Error : Excution time out of limit. please check your file and feel free to contact us if you have any question.")
        return {'status':zMap_status}
    except Exception as e:
        zMap_status.append("Error : " + str(e) + ". please check your file and feel free to contact us if you have any question.")
        return {'status':zMap_status}



