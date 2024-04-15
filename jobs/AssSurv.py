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
from celery.task.control import revoke
import zipfile
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


import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt

import argparse
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12

from lifelines import CoxPHFitter

@celery.task(name='cox_regression',bind=True)
def cox_regression(self,z_statistic,survival_f,outdir):
    try:
        with time_limit(60*60*6):
            z_statistic_f = pd.read_csv(z_statistic,sep="\t",index_col=0)
            survival_df = pd.read_csv(survival_f,sep="\t",index_col=0)
            overlap_sample = list(set(z_statistic_f.columns)&set(survival_df.index))
            n=len(overlap_sample)

            status = ["Only %s common samples in input files were included in cox proportional hazards regression analysis."%str(n)]
            self.update_state(state='PROGRESS',meta={'status': status})

            z_statistic_f = z_statistic_f[overlap_sample]
            survival_df = survival_df.loc[overlap_sample]
            status.append("Starting cox proportional hazards regression analysis.This may take a few minutes.")
            self.update_state(state='PROGRESS',meta={'status': status})
            subgroup_df = survival_df

            half_n = int(len(z_statistic_f.columns)/2)

            z_statistic_f = z_statistic_f.dropna(thresh=half_n)

            protein_list = []
            results_list = []

            for index in z_statistic_f.index:
                gene_df = z_statistic_f.loc[index].to_frame()
                gene_df["survival_time"] = subgroup_df.loc[gene_df.index,"survival_time"]
                gene_df["death_or_not"] = subgroup_df.loc[gene_df.index,"death_or_not"]
                gene_df = gene_df.dropna(how="any")
            
                cph = CoxPHFitter()
                cph.fit(gene_df, duration_col='survival_time', event_col='death_or_not')
                results = cph.summary
                results_list.append(results)

            status.append("combine results and make multiple hypothesis correction...")
            self.update_state(state='PROGRESS',meta={'status': status})
            all_results_df = pd.concat(results_list)
            all_results_df["BH-corrected pvalue"] =all_results_df["p"]/all_results_df["p"].rank(axis=0) * all_results_df["p"].count()
            
            sub_all_results_df = all_results_df[['coef','coef lower 95%','coef upper 95%','p','BH-corrected pvalue']]
            sub_all_results_df.columns = ['log(HR)','log(HR) lower 95%','log(HR) upper 95%','pvalue','BH-corrected pvalue']
            protein_type = {}
            for index in sub_all_results_df.index:
                if sub_all_results_df.loc[index,'BH-corrected pvalue']<0.05:
                    if sub_all_results_df.loc[index,'log(HR) upper 95%'] > 0:
                        protein_type[index] = "unfavorable"
                    else:
                        protein_type[index] = "favorable"
                else:
                    protein_type[index] = "not prognostic"
            new_df = pd.DataFrame.from_dict(protein_type,orient="index")
            new_df.columns=["prognostic association"]
        #    print(new_df)
            all_df = pd.concat([sub_all_results_df,new_df],axis=1)

            all_df.to_csv(outdir+"/results_cox_regression.txt",sep = "\t")

            zip_ya(outdir)    
            status.append("All finished.")
            self.update_state(state='PROGRESS',meta={'status': status})

            return {'status': status}
    except TimeoutException:
        status.append("Error : Excution time out of limit. please check your file and feel free to contact us if you have any question.")
        return {'status':status}
    except Exception as e:
        status.append("Error : " + str(e) + ". please check your file and feel free to contact us if you have any question.")
        return {'status':status}

