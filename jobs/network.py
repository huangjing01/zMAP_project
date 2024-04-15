from app import celery
import os,zipfile,json
from flask import current_app

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
    
# @celery.task(name='async_network',bind=True)
# def async_network(self,filepath,module_detect):
#     zMap_status = ["start network processing..."]
#     self.update_state(state='PROGRESS',meta={'status': zMap_status})
#     if module_detect == 'static':
#         raw_intensity_data = pd.read_excel(filepath,index_col=0)
#         if len(set(raw_intensity_data.index))!=len(raw_intensity_data.index):
#             raise Exception("Gene IDs were not uniq")
#         elif len([i for i in list(raw_intensity_data.index) if isinstance(i,float)]) > 0:
#             raise Exception("Nans are included in the gene ID")
#         elif len(set(raw_intensity_data.columns))!=len(raw_intensity_data.columns):
#             raise Exception("Sample IDs were not uniq")
#         else:
#             zMap_status.append("calculate sample correlation")
#             self.update_state(state='PROGRESS',meta={'status':'calculate sample correlation'})
#             os.chdir(os.path.dirname(filepath))
#             os.mkdir("result")
#             res_dir = os.path.join(os.path.dirname(filepath),"result")
#             sample_pcc = pcc(intensity_df,res_dir)
#             self.update_state(state='PROGRESS',meta={'status':'start cluster'})
#             condition_color_list = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(len(set(condition_df[condition_df.columns[0]])))]
#             rep_color_list = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(len(set(condition_df[condition_df.columns[1]])))]
#             condition_color_dict = dict(zip(list(set(condition_df[condition_df.columns[0]])),condition_color_list))
#             rep_color_dict = dict(zip(list(set(condition_df[condition_df.columns[1]])),rep_color_list))
            
#             condition_color = map(condition_color_dict.get,list(condition_df[condition_df.columns[0]]))
#             rep_color = map(rep_color_dict.get,list(condition_df[condition_df.columns[1]]))
#             df = pd.DataFrame([condition_color,rep_color],index=["condition","batch"],columns=condition_df.index)
#             df = df.T
#             sample_order = hierarchical_clustering(sample_pcc,df,res_dir)
#             sub_protein_clusters,cluster_list = gene_cluster(intensity_df,sample_order,8,50,-6,6,res_dir)
#             #condition order
#             zip_ya(os.path.dirname(filepath))
#             return {'status': 'cluster finished'}                   




@celery.task(name='async_network',bind=True)
def async_network(self,filepath,module_detect):
    zMap_status = ["start network processing..."]
    self.update_state(state='PROGRESS',meta={'status': zMap_status})
    os.system("Rscript ~/zMap/jobs/net_construct.R "+filepath+ " static")

    if module_detect == 'static':
        os.system("Rscript ~/zMap/jobs/mod_detect.R "+ os.path.split(filepath)[0]+"/Protein_partial_correlation.txt"+" static")
    else:
        os.system("Rscript ~/zMap/jobs/mod_detect.R "+ os.path.split(filepath)[0]+"/Protein_partial_correlation.txt"+" dynamic")
    zip_ya(os.path.dirname(filepath))
    zMap_status.append("Network finished")
    return {'status': zMap_status} 
