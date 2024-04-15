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
from flask import render_template, request, jsonify,current_app,redirect, url_for,send_file,send_from_directory
from . import Surv
from werkzeug.utils import secure_filename

import os,string,random
import numpy as np
import pandas as pd

import pandas as pd
import numpy as np
import matplotlib
import os
matplotlib.use('Agg') 
import matplotlib.pyplot as plt

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
#matplotlib.rcParams['font.family'] ='Arial'
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.statistics import multivariate_logrank_test
import argparse
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec


def survival_analysis(survival_df,outdir):
    survival_df = survival_df.dropna(how="any")
    kmf = KaplanMeierFitter()
    fig = plt.figure(figsize=(5,5),dpi=300)

    color_dict = {1:"tab:green",2:"tab:blue",3:"tab:red",4:"tab:orange",5:'tab:purple',6:'tab:brown',7:'tab:pink',8:'tab:gray',9:'tab:olive',10:'tab:cyan'}

    ax = fig.gca()

    for name, grouped_df in survival_df.groupby("group"):
        #print (group_df)
        kmf.fit(grouped_df["survival_time"], grouped_df["death_or_not"], label="group"+str(name)+"(%s)"%(len(grouped_df)),
                alpha =0.1)

        kmf.plot(ax=ax,show_censors=True,ci_show=False,color=color_dict[name],
               censor_styles={'ms': 6},linewidth=3)
    from lifelines.statistics import multivariate_logrank_test

    results = multivariate_logrank_test(survival_df['survival_time'],
                                   survival_df['group'], 
                                   survival_df['death_or_not'])
    pvalue = results.p_value
    plt.rcParams['xtick.labelsize']=12
    plt.rcParams['ytick.labelsize']=12

    plt.ylim([0,1.1])
    plt.text(5,0.5,"P-value:" +'{:.2e}'.format(pvalue),size=12)
    ax.set_xlabel("Survival time",size=15)
    ax.set_ylabel("Survival ratio",size=15)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    plt.savefig(outdir+"/survival_analysis.pdf",dpi=300,bbox_inches="tight")
def mkdir(path_dir):
    isexists = os.path.exists(path_dir)
    if not isexists:
        os.makedirs(path_dir)
    else:
        pass

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


#from jobs.pca import async_pca
def processID(length=8,chars=string.ascii_letters+string.digits):
    return ''.join([random.choice(chars) for i in range(length)])


@Surv.route('/', methods=['GET', 'POST'])
def index():
    return render_template('Surv_index.html')


@Surv.route('/result', methods=['POST'])
def result():
    if request.files['input-1']:
        survival_df = pd.read_csv(request.files['input-1'],sep="\t",index_col=0)

        proc_id = processID(8)
        os.chdir(current_app.config['UPLOAD_FOLDER'])
        while(os.path.exists(proc_id)):
            proc_id = processID(8)
        proc_id = proc_id + "_Surv"
        outdir = os.path.join(os.getcwd(),proc_id)
        mkdir(outdir)
        log_info = ["Plots a figure of the Kaplan–Meier estimate model"]

        survival_analysis(survival_df,outdir)
        log_info.append("All finished.")

        zip_ya(outdir)

        return render_template('Surv_result.html',log_info=log_info,filedir_id=proc_id+".zip")
    else:
        flash('please specify valid input files.')
        return redirect(url_for('Surv.index'))

@Surv.route('/result/<filedir>', methods=['GET'])
def fetch_filedir(filedir):
    if filedir[-4:] == '.zip':
        return send_from_directory(os.path.abspath(current_app.config['UPLOAD_FOLDER']),filedir,as_attachment=True)
