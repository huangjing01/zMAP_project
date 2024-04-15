from __future__ import with_statement
import signal
from contextlib import contextmanager
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
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
import os,zipfile,json,random
from flask import current_app
import numpy as np

from celery.task.control import revoke
# from .zmap_step2_downstream_analysis import *

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



def pathway_info():
    kegg_pathway_file = "/home/shao/zMap/jobs/gmt_files/enrichr.KEGG_2016.gmt"
    go_bp_pathway_file = "/home/shao/zMap/jobs/gmt_files/human_GO_Biological_Process_2015.gmt"
    go_cc_pathway_file = "/home/shao/zMap/jobs/gmt_files/human_GO_Cellular_Component_2015.gmt"
    go_mf_pathway_file = "/home/shao/zMap/jobs/gmt_files/human_GO_Molecular_Function_2015.gmt"

    pathway_dict={}
    for gmt_file in [kegg_pathway_file,go_bp_pathway_file,go_cc_pathway_file,go_mf_pathway_file]:
        f = open(gmt_file,"r")
        gene_list = f.readlines()
        for line in gene_list:
            line_elements = line.strip().split("\t")
            term = line_elements[0]
            for gene in line_elements[1:]:
                if gene == "":
                    pass
                else:
                    if term in pathway_dict.keys():
                        pathway_dict[term].append(gene)
                    else:
                        pathway_dict[term]=[gene]
    return pathway_dict


def enrich(pathway_dict,gene_list):
    bg_gene = []
    [bg_gene.extend(i) for i in pathway_dict.values()]
    enrich_result=[]
    for key in pathway_dict.keys():
        target_pathway_gene_list = pathway_dict[key]
        a = len(set(gene_list) & set(target_pathway_gene_list)) ###interest gene numbers in target path  
        b = len(gene_list) - a   ### insterest gene numbers not in target path
        c = len(target_pathway_gene_list) - a
        d = len(bg_gene) - (a+b) - c
        # print(a)
        if a > 0:
            pvalue = fisher_exact([[a,b],[c,d]],alternative="greater")[1]
            genes=",".join(list(set(gene_list)&set(target_pathway_gene_list)))
            enrich_result.append([key,str(a)+"/"+str(len(target_pathway_gene_list)),genes,pvalue])
        else:
            pass
    enrich_result_df = pd.DataFrame(enrich_result,columns=["pathway","overlap","gene list","pvalue"])
    sorted_df = enrich_result_df.sort_values(by=["pvalue"])
    padjust_value = [len(enrich_result_df)*i if len(enrich_result_df)*i<=1 else 1  for i in sorted_df["pvalue"]]
    sorted_df["padjust_by_bonferroni"]=padjust_value
    sorted_df["padjust_by_BH"] =[i if i<1 else 1 for i in sorted_df["pvalue"]/sorted_df["pvalue"].rank(axis=0)*len(sorted_df)]
    sorted_df = sorted_df.sort_values(by=["padjust_by_BH"])
   # sorted_df = sorted_df.sort_values(by="padjust_by_BH")
    sorted_df.index = np.arange(len(sorted_df))
    # print(sorted_df)
    return sorted_df

def enrich_plot(enrichment_result,fig_file):
    sig_enrichment_df = enrichment_result.iloc[0:20,:].loc[enrichment_result["padjust_by_BH"]<=0.05]
    sig_enrichment_df = sig_enrichment_df.sort_values(by=["padjust_by_BH"],ascending=False)
    if len(sig_enrichment_df) >0:
        plt.figure(figsize=(5,len(sig_enrichment_df)/2 ))
        plt.barh([i for i in range(len(sig_enrichment_df))],-np.log(sig_enrichment_df["padjust_by_BH"].values))
        
        plt.xlabel('$-log_{10}{ (padjust\ by\ BH)}$',size=20)
        plt.yticks([i for i in range(len(sig_enrichment_df))],sig_enrichment_df["pathway"] + " ["+sig_enrichment_df["overlap"]+"]",size=15)
        plt.title("enrichment in GO, KEGG",size=20)
        plt.savefig(fig_file,dpi=300,bbox_inches="tight")
        plt.close("all")
    else:
        pass


@celery.task(name='async_enrichment',bind=True)
def async_enrichment(self,filepath,species,id_type):
    if os.path.split(filepath)[-1] == 'single_gene_list.txt':
        self.update_state(state='PROGRESS',meta={'status':'start enrichment processing'})
        gene_list = []
        with open(filepath,'r') as f:
            for line in f:
                gene_list.append(line.rstrip())
        outdir = os.path.dirname(filepath) + "/result"
        os.mkdir(outdir)
        pathway_dict = pathway_info()
        enrichment_result = enrich(pathway_dict,gene_list)
        enrich_file = outdir + "/enrichment_result.png"
        enrich_plot(enrichment_result,enrich_file)
        enrichment_result.to_csv(outdir+"/enrichment_result.txt",sep = "\t")
        zip_ya(os.path.dirname(filepath))
        return {'status': 'enrichment finished'}                   
