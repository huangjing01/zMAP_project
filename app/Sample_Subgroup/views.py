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
from . import Sample_Subgroup
from werkzeug.utils import secure_filename

import os,string,random
import numpy as np
import pandas as pd




import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import matplotlib
#matplotlib.use('Agg') 
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
#matplotlib.rcParams['font.family'] ='Arial'
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
import os
import math
from scipy import stats
import seaborn as sns
import argparse
import matplotlib.gridspec as gridspec

def mkdir(path_dir):
    isexists = os.path.exists(path_dir)
    if not isexists:
        os.makedirs(path_dir)
    else:
        pass



def samplelist_for_cluster(sample_info_df,sample_type):
    
    sample_type=sample_type.strip().split(',')
    sample_list = []
    for index in sample_info_df.index:
        sample_condition = sample_info_df.loc[index,"Sample_condition"]
        if sample_condition in sample_type:
            sample_list.append(index)
        else:
            pass
     
    return sample_list


def tumor_mean_variance(tumor_z_df,outdir):
    
    mean_ = tumor_z_df.mean(axis=1)
    
    variance_ = tumor_z_df.var(axis=1)
    x=mean_
    y=variance_
    
    sns.scatterplot(x=x, y=y, s=5, color=".15")
    sns.histplot(x=x, y=y, bins=50, pthresh=.1, cmap="mako")
    sns.kdeplot(x=x, y=y, levels=5, color="w", linewidths=1)
    plt.xlabel("Average of z-statistic")
    plt.ylabel("Variance of z-statistic")
    plt.savefig(outdir+"/mean_variance_scatterplot.png",dpi=300)
    plt.savefig(outdir+"/mean_variance_scatterplot.pdf",dpi=300)
    

def top_tumor_variance_table(z_df,tumor_z_df,outdir,top=3000):
    
    mean_ = tumor_z_df.mean(axis=1)
    
    variance_ = tumor_z_df.var(axis=1).to_frame()
    variance_sort = variance_.sort_values(by=0,ascending=False)
    
    top_index = variance_sort.index[0:top]
        
#    top_variance_all_z_df = tumor_z_df.loc[top_index]

    top_variance_z_df_dropna = z_df.loc[top_index].dropna(how="any")

    top_variance_tumor_z_df_dropna = top_variance_z_df_dropna[tumor_z_df.columns]
    
   

    file_ = outdir+"/top_variance_tumor_z_df_dropna_%s_%s.txt"%(top,len(top_variance_tumor_z_df_dropna))

    top_variance_tumor_z_df_dropna.to_csv(file_,sep="\t")
        
    return file_

def prepare_ConsensusClusterPlus(file_,outdir): 
    out = open(outdir+"/ConsensusClusterPlus.r","w")   
    out.write("setwd('%s')\n"%(outdir))
    out.write("library(ConsensusClusterPlus)\n")
    out.write("dc = read.csv('%s',sep='\t',row.names = 1,check.names=FALSE)\n"%(file_))
    out.write("dc = as.matrix(dc)\n")
    out.write("rcc = ConsensusClusterPlus(dc,maxK=5,reps=1000,pItem=0.8,pFeature=1,title='euclidean_km', \
                           distance='euclidean',clusterAlg='km',plot='pdf',seed=1262118322)\n")
    out.write("cluster <- rcc[[3]]$consensusClass\n")
    out.write("write.csv(cluster,file='%s', quote = FALSE)\n"%(outdir+'/cluster_3.csv'))
    out.write("cluster <- rcc[[4]]$consensusClass\n")
    out.write("write.csv(cluster,file='%s', quote = FALSE)\n"%(outdir+'/cluster_4.csv'))
    out.write("cluster <- rcc[[5]]$consensusClass\n")
    out.write("write.csv(cluster,file='%s', quote = FALSE)\n"%(outdir+'/cluster_5.csv'))
    out.close()
    os.system("Rscript %s"%(outdir+"/ConsensusClusterPlus.r"))

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


@Sample_Subgroup.route('/', methods=['GET', 'POST'])
def index():
    return render_template('Sample_Subgroup_index.html')

@Sample_Subgroup.route('/result', methods=['POST'])
def result():
    if request.files['input-1'] and request.files['input-2'] and request.form.get('condition_name') and request.form.get('protein_number'):

        z_df = pd.read_csv(request.files['input-1'],sep="\t",index_col=0)

        sample_info_df = pd.read_csv(request.files['input-2'],sep="\t",index_col=0)

        sample_type = request.form.get('condition_name')

        top_n = int(request.form.get('protein_number'))

        sample_list = samplelist_for_cluster(sample_info_df,sample_type)



        overlap_sample = list(set(sample_list)&set(z_df.columns))

        proc_id = processID(8)
        os.chdir(current_app.config['UPLOAD_FOLDER'])
        while(os.path.exists(proc_id)):
            proc_id = processID(8)
        proc_id = proc_id + "_Subgrouping"
        outdir = os.path.join(os.getcwd(),proc_id)
        mkdir(outdir)
        log_info = ["%s common samples in input files were included in subsequant analysis."%(len(overlap_sample))]

        tumor_z_df = z_df[overlap_sample]

        tumor_mean_variance(tumor_z_df,outdir)
        log_info.append("Generating mean-variance scatterplot.")

        top_variance_z_file = top_tumor_variance_table(z_df,tumor_z_df,outdir,top=top_n)
        log_info.append("Perform clustering on %s %s samples."%(len(overlap_sample),sample_type))

        prepare_ConsensusClusterPlus(top_variance_z_file,outdir)
        log_info.append("All finished.")

        zip_ya(outdir)

        return render_template('Sample_Subgroup_result.html',log_info=log_info,filedir_id=proc_id+".zip")
    else:
        flash('please specify valid input files.')
        return redirect(url_for('Sample_Subgroup.index'))

@Sample_Subgroup.route('/result/<filedir>', methods=['GET'])
def fetch_filedir(filedir):
    if filedir[-4:] == '.zip':
        return send_from_directory(os.path.abspath(current_app.config['UPLOAD_FOLDER']),filedir,as_attachment=True)


