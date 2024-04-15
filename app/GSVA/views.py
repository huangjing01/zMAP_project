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
from flask import render_template, request, jsonify,current_app,redirect, url_for,send_file,send_from_directory,flash
from . import GSVA
from werkzeug.utils import secure_filename

import os,string,random
import numpy as np
import pandas as pd
#from jobs.pca import async_pca

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import argparse
import os
import random
import matplotlib.pyplot as plt

def get_batch_color_dict(batch_df,col_=""):
    
    num = len(set(batch_df[col_].values))
#    print(num)
    batch_color_dict = {}
    color_list = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])for i in range(num)]
                  
    i=0
    list_ = list(set(batch_df[col_].values))
    list_ = sorted(list_)
   
    for index in list_:
        
        batch_color_dict[index] = color_list[i]
        i+=1
    return batch_color_dict


def mkdir(path_dir):
    isexists = os.path.exists(path_dir)
    if not isexists:
        os.makedirs(path_dir)
    else:
        pass

def prepare_gsva_kegg(z_statistic_matrix,sample_info,gmt,outdir,type_="kegg",top_n=50,adj_P_Val=0.05):
    os.chdir(outdir)
    r_dir = outdir+'/gsva_%s.r'%type_
    out = open(r_dir,"w")
    out.write("library(GSVA)\n")
    out.write("library(limma)\n")
    out.write("library(GSEABase,quietly=TRUE)\n")
    out.write('group_info <- read.csv("%s",sep="\\t",row.names = 1,stringsAsFactors = FALSE)\n'%sample_info)
  
    out.write('data <- read.csv("%s",sep="\\t",row.names = 1)\n'%z_statistic_matrix)

    out.write('data <- na.omit(as.matrix(data))\n')
    out.write('keggset <- getGmt("%s")\n'%gmt)
    out.write('keggs_kcdf_none <- gsva(data,gset.idx.list=keggset,kcdf="Gaussian", parallel.sz=1,min.sz=2)\n')
    out.write('write.csv(keggs_kcdf_none,file="%s")\n'%("%s/%s_gsva_kcdf_Gaussian.csv"%(outdir,type_)))

    out.write('design <- factor(c(group_info$Sample_condition))\n')
    out.write('design1 <- model.matrix(~0+design)\n')

    out.write('colnames(design1) <- levels(design)\n')

    out.write('fit <- lmFit(keggs_kcdf_none, design1)\n')
    
    condition_df = pd.read_csv(sample_info,sep="\t",index_col=0)
    condition_list = list(set(list(condition_df['Sample_condition'].values)))
    str_ = ""
    n_condition = len(condition_list)

    for i in range(0,n_condition-1):
        for j in range(i+1,n_condition):
            str_ = str_  +condition_list[i]+"-"+condition_list[j]+","
    
    

    out.write('contrast.matrix <- makeContrasts('+str_+' levels = design1)\n')

    out.write('fit2 <- contrasts.fit(fit, contrast.matrix)\n')
    out.write('fit2 <- eBayes(fit2)\n')
    out.write("top <- topTable(fit2, number= %s)\n"%(top_n))

    out.write('write.csv(fit2,file="%s")\n'%("%s/%s_gsva_differential_pathway_activity.csv"%(outdir,type_)))

    out.write('write.csv(top,file="%s")\n'%("%s/%s_gsva_differential_pathway_activity_top_%s.csv"%(outdir,type_,str(top_n))))
    out.close()

    os.system("Rscript %s"%r_dir)

    gsva_score_df = pd.read_csv("%s"%("%s/%s_gsva_kcdf_Gaussian.csv"%(outdir,type_)),sep=",",index_col=0)

    gsva_score_df.columns = [i.strip('"') for i in gsva_score_df.columns]
    gsva_score_df.index = [i.strip('"') for i in gsva_score_df.index]
    
    diff_pathway_df = pd.read_csv("%s"%("%s/%s_gsva_differential_pathway_activity_top_%s.csv"%(outdir,type_,str(top_n))),sep=",",index_col=0)
    diff_pathway_df.columns = [i.strip('"') for i in diff_pathway_df.columns]
    diff_pathway_df.index = [i.strip('"') for i in diff_pathway_df.index]
    sub_diff_pathway_df = diff_pathway_df.loc[diff_pathway_df['adj.P.Val'] < adj_P_Val]

    list_ = list(condition_df.index)
    list_ = sorted(list_)

    diff_score_df = gsva_score_df.loc[sub_diff_pathway_df.index,list_]

    condition_color_dict = get_batch_color_dict(condition_df,col_="Sample_condition")
    sample_color_list = [condition_color_dict[condition_df.loc[i,'Sample_condition']] for i in list_]

    sns.clustermap(diff_score_df,col_cluster=False,cmap="bwr",col_colors=sample_color_list,xticklabels=True,yticklabels=True,row_cluster=True,
                   cbar_kws = {"shrink": .5,'label':"GSVA score"})

    plt.savefig("%s"%("%s/%s_gsva_differential_pathway_activity_top_%s_adj.P.Val_%s.pdf"%(outdir,type_,str(top_n),adj_P_Val)),dpi=300,bbox_inches="tight")

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


def processID(length=8,chars=string.ascii_letters+string.digits):
    return ''.join([random.choice(chars) for i in range(length)])


@GSVA.route('/', methods=['GET', 'POST'])
def index():
    return render_template('GSVA_index.html')


@GSVA.route('/result', methods=['POST'])
def result():
    if request.files['z_statistic_matrix'] and request.files['sample_info'] and request.form.get('cluster_number_for_top_proteins') and request.form.get('fdr'):
        z_statistic_matrix = request.files['z_statistic_matrix']
        sample_info = request.files['sample_info']
        top_n = int(request.form.get('cluster_number_for_top_proteins'))
        adj_P_Val = float(request.form.get('fdr'))

        proc_id = processID(8)
        os.chdir(current_app.config['UPLOAD_FOLDER'])
        while(os.path.exists(proc_id)):
            proc_id = processID(8)
        proc_id = proc_id + "_GSVA"
        outdir = os.path.join(os.getcwd(),proc_id)
        mkdir(outdir)

        z_statistic_matrix_filename = secure_filename(z_statistic_matrix.filename)
        z_statistic_matrix.save(os.path.join(outdir,z_statistic_matrix_filename))
        sample_info_filename = secure_filename(sample_info.filename)
        sample_info.save(os.path.join(outdir,sample_info_filename))
        
        log_info = ["start gene sets variation analysis."]
        kegg_gmt = "/home/shao/zMap/app/templates/input_file/reverse_zMAP/gmt/enrichr.KEGG_2016.gmt"
        prepare_gsva_kegg(z_statistic_matrix_filename,sample_info_filename,kegg_gmt,outdir,type_="kegg",top_n=50,adj_P_Val=adj_P_Val)

        go_gmt = '/home/shao/zMap/app/templates/input_file/reverse_zMAP/gmt/human_GO_Biological_Process_2015.gmt'
        prepare_gsva_kegg(z_statistic_matrix_filename,sample_info_filename,go_gmt,outdir,type_="go",top_n=50,adj_P_Val=adj_P_Val)
        zip_ya(outdir)
        log_info.append("All finished.")


        return render_template('GSVA_result.html',log_info=log_info,filedir_id=proc_id+".zip")
    else:
        flash('please specify valid input files.')
        return redirect(url_for('GSVA.index'))



@GSVA.route('/result/<filedir>', methods=['GET'])
def fetch_filedir(filedir):
    if filedir[-4:] == '.zip':
        return send_from_directory(os.path.abspath(current_app.config['UPLOAD_FOLDER']),filedir,as_attachment=True)



