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
import os,zipfile,json,random
from flask import current_app
import numpy as np

from celery.task.control import revoke
from .zmap_step2_downstream_analysis import *

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



def pca_scatter_2(condition_df,pc_df,pca,condition_dict,info_dict,outdir):
    pc1_ratio = "%.2f"% (pca.explained_variance_ratio_[0]*100)
    pc2_ratio = "%.2f"% (pca.explained_variance_ratio_[1]*100)
    plt.figure(figsize=(5,5))
    for index in pc_df.index:
        plt.scatter(pc_df.loc[index,"PC1"],pc_df.loc[index,"PC2"],edgecolors = condition_dict[condition_df.loc[index,condition_df.columns[0]]],
                    linewidths=2,facecolors="white",marker=info_dict[condition_df.loc[index,condition_df.columns[1]]],s=70)

    plt.xlabel("PC1(%s"%(pc1_ratio)+"%)",size=15)
    plt.ylabel("PC2(%s"%(pc2_ratio)+"%)",size=15)
    legend_elements = []

    for key in info_dict.keys():
        legend_elements.append(Line2D([0], [0],color="black", 
                                      markeredgewidth = 2,
                                      markerfacecolor="white",
                                      marker = info_dict[key],label=key,
                                      linestyle='None',
                                      markersize=8))
    for key in condition_dict.keys():
        legend_elements.append(Line2D([0], [0], color = condition_dict[key], lw=3, label=key))

    plt.legend(handles=legend_elements, loc="best")
    plt.savefig(outdir + "/sample_cluster_pca.png",dpi=300,bbox_inches="tight")


def pca_scatter_1(condition_df,pc_df,pca,condition_dict,outdir):
    pc1_ratio = "%.2f"% (pca.explained_variance_ratio_[0]*100)
    pc2_ratio = "%.2f"% (pca.explained_variance_ratio_[1]*100)
    plt.figure(figsize=(5,5))
    # print(pc_df)
    for index in pc_df.index:
        plt.scatter(pc_df.loc[index,"PC1"],pc_df.loc[index,"PC2"],edgecolors = condition_dict[condition_df.loc[index,condition_df.columns[0]]],
                    linewidths=2,facecolors="white",s=70)
    plt.xlabel("PC1(%s"%(pc1_ratio)+"%)",size=15)
    plt.ylabel("PC2(%s"%(pc2_ratio)+"%)",size=15)
    legend_elements = []

    for key in condition_dict.keys():
        legend_elements.append(Line2D([0], [0], color = condition_dict[key], lw=3, label=key))
    plt.legend(handles=legend_elements, loc="best")
    plt.savefig(outdir + "/sample_cluster_pca.png",dpi=300,bbox_inches="tight")

@celery.task(name='async_pca',bind=True)
def async_pca(self,filepath,task_type):
    if task_type == 'upload':
        self.update_state(state='PROGRESS',meta={'status':'start pca processing'})
        raw_intensity_data = pd.ExcelFile(filepath)
        sheet_names = raw_intensity_data.sheet_names
        intensity_df = raw_intensity_data.parse(sheet_names[0],index_col=0)
        condition_df = raw_intensity_data.parse(sheet_names[1],index_col=0)
        if len(set(intensity_df.index))!=len(intensity_df.index):
            raise Exception("Gene IDs were not uniq")
        elif len([i for i in list(intensity_df.index) if isinstance(i,float)]) > 0:
            raise Exception("Nans are included in the gene ID")
        elif len(set(intensity_df.columns))!=len(intensity_df.columns):
            raise Exception("Sample IDs were not uniq")

        elif len([i for i in list(intensity_df.index) if str(i).startswith('Unnamed')]) > 0:
            raise Exception("Nans are included in the sample ID")

        elif set(intensity_df.columns) != set(condition_df.index):
            raise Exception("Sample list not exactly the same")
        elif condition_df.isna().any().any():
            raise Exception("Nans are included in condition info.")
        else:
            self.update_state(state='PROGRESS',meta={'status':'calculate sample correlation'})
            os.chdir(os.path.dirname(filepath))
            os.mkdir("result")
            res_dir = os.path.join(os.path.dirname(filepath),"result")
            sample_pcc = pcc(intensity_df,res_dir)
            self.update_state(state='PROGRESS',meta={'status':'start pca'})
            sample_pca, pca_res = pca_analysis(sample_pcc,res_dir)
            self.update_state(state='PROGRESS',meta={'status':'explained variance ratio (first three components): %s' % str(sample_pca.explained_variance_ratio_)})
            if len(condition_df.columns) > 1:

                condition_label = condition_df.columns
                condition_group = list(set(condition_df[condition_label[0]]))
                condition_group_num = len(condition_group)
                condition_group_color_list = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(condition_group_num)]
                condition_dict = dict(zip(condition_group,condition_group_color_list))

                info_group = list(set(condition_df[condition_label[1]]))
                info_group_num = len(info_group)
                markers = ['.', ',', 'o', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p', 'P', '*', 'h', 'H', '+', 'x', 'X', 'D','d', '|', '_']
                info_dict = dict(zip(info_group,markers[0:info_group_num]))
                pca_scatter_2(condition_df,pca_res,sample_pca,condition_dict,info_dict,res_dir)
                zip_ya(os.path.dirname(filepath))
            else:
                condition_label = condition_df.columns
                condition_group = list(set(condition_df[condition_label[0]]))
                condition_group_num = len(condition_group)
                condition_group_color_list = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(condition_group_num)]
                condition_dict = dict(zip(condition_group,condition_group_color_list))
                pca_scatter_1(condition_df,pca_res,sample_pca,condition_dict,res_dir)
                zip_ya(os.path.dirname(filepath))
            return {'status': 'pca finished'}                   
