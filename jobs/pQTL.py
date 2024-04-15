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
import os
matplotlib.use('Agg') 
import matplotlib.pyplot as plt

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
#matplotlib.rcParams['font.family'] ='Arial'
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
import argparse
import os


def prepare_pQTL(SNP_file_name,expression_file_name,covariates_file_name,outdir):
    
    os.chdir(outdir)
    out=open('pQTL.r',"w")
    out.write("library(MatrixEQTL)\n")


    out.write('useModel = modelLINEAR;\n\n')    
    out.write("output_file_name = tempfile();\n\n")

    # Only associations significant at this level will be saved
    out.write("pvOutputThreshold = 0.1;\n\n")

    # Error covariance matrix
    # Set to numeric() for identity.
    out.write("errorCovariance = numeric();\n\n\n")
    # errorCovariance = read.table("Sample_Data/errorCovariance.txt");


    ## Load genotype data

    out.write('snps = SlicedData$new();\nsnps$fileDelimiter = "\\t";\nsnps$fileOmitCharacters = "NA";\nsnps$fileSkipRows = 1;\nsnps$fileSkipColumns = 1;\nsnps$fileSliceSize = 2000;\nsnps$LoadFile("%s");\n\n\n'%SNP_file_name)

    ## Load gene expression data

    out.write('gene = SlicedData$new();\n \
    gene$fileDelimiter = "\\t";      # the TAB character\n \
    gene$fileOmitCharacters = "NA"; # denote missing values;\n \
    gene$fileSkipRows = 1;          # one row of column labels\n \
    gene$fileSkipColumns = 1;       # one column of row labels\n \
    gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows\ngene$LoadFile("%s");\n\n\n'%expression_file_name)

    ## Load covariates

    out.write('cvrt = SlicedData$new();\n \
cvrt$fileDelimiter = "\\t";      # the TAB character\n \
cvrt$fileOmitCharacters = "NA"; # denote missing values;\n \
cvrt$fileSkipRows = 1;          # one row of column labels\n \
cvrt$fileSkipColumns = 1;       # one column of row labels\n \
if(length("%s")>0) {\n \
cvrt$LoadFile("%s");\n \
}\n\n\n'%(covariates_file_name,covariates_file_name))

    ## Run the analysis

    out.write('me = Matrix_eQTL_engine(\n \
    snps = snps,\n \
  gene = gene,\n \
  cvrt = cvrt,\n \
  output_file_name = output_file_name,\n \
  pvOutputThreshold = pvOutputThreshold,\n \
  useModel = useModel,\n \
  errorCovariance = errorCovariance,\n \
  verbose = TRUE,\n \
  pvalue.hist = TRUE,\n \
  min.pv.by.genesnp = FALSE,\n \
  noFDRsaveMemory = FALSE);\n\n\n')

    out.write('unlink(output_file_name);\n\n\n')

    ## Results:

    out.write("#cat('Analysis done in: ', me$time.in.sec, ' seconds', '\\n');\n\n\n")
    out.write("#cat('Detected eQTLs:', '\\n');\n\n\n")
    out.write("#show(me$all$eqtls)\n\n\n")

## Plot the histogram of all p-values

    out.write("plot(me)\n\n\n")
 
    out.write("results <- me$all$eqtls\n\n\n")

    out.write('write.table (results, file ="results.txt", sep ="\\t",row.names =FALSE, col.names =TRUE,quote=FALSE)\n')
    out.close()
    os.system("Rscript pQTL.r")
    pqtl_results = outdir+"/results.txt"
    return pqtl_results

def mkdir(path_dir):
    isexists = os.path.exists(path_dir)
    if not isexists:
        os.makedirs(path_dir)
    else:
        pass

def add_bh_pvalue(z_df,pqtl_results,fdr=0.05):

    new_pqtl_df = pd.read_csv(pqtl_results,sep="\t")
    remain_pqtl_df = new_pqtl_df
    p_num = len(z_df)

    all_snps_gene_list = set(remain_pqtl_df["snps"])

    df_list = []
    for snp_gene in all_snps_gene_list:
        snp_gene_df = remain_pqtl_df.loc[remain_pqtl_df["snps"]==snp_gene]
        sort_snp_gene_df = snp_gene_df.sort_values(by="pvalue",ascending = True)
        pvalue_list = []
        i=0
        for index in sort_snp_gene_df.index:
            pvalue = sort_snp_gene_df.loc[index,"pvalue"]
            if i==0:    
            #pvalue = sort_snp_gene_df.loc[index,"pvalue"]
                pvalue_list.append(min(1,pvalue*p_num))
                i+=1
            else:
                bh_pvalue = pvalue*p_num/(i+1)
          #  print(bh_pvalue)
                i+=1
                if bh_pvalue < pvalue_list[-1]:
                    pvalue_list.append(pvalue_list[-1])
                else:
                    pvalue_list.append(min(1,bh_pvalue))
        sort_snp_gene_df["FDR"] = pvalue_list
    
        df_list.append(sort_snp_gene_df)
    all_pqtl_bh_df = pd.concat(df_list)
    all_pqtl_bh_df = all_pqtl_bh_df.loc[all_pqtl_bh_df['FDR']<fdr]
    all_pqtl_bh_df.to_csv("results_association_z_statistic_with_mutation.txt",sep="\t",index=False)
    return all_pqtl_bh_df

def generate_figure(gene_tss_location,chr_length,pQTL_out,fdr=0.05):

    chromo_size_df = pd.read_csv(chr_length,sep="\t",index_col=0)
    chrom_dict = {}
    axis=0
    for i in chromo_size_df.index:
        chrom_dict[i]=axis
        axis += chromo_size_df.loc[i,"length"]
    
    gene_tss_df = pd.read_csv(gene_tss_location,sep="\t",index_col=0)
    autosome_gene_tss_df = gene_tss_df.loc[[ i in chrom_dict.keys() for i in gene_tss_df["chromosome"]]]

    pqtl_df = pQTL_out
    pair_index_list = []
    snps_location_list = []
    gene_location_list = []
    color_list = []
    pvalue_list = []
    snps_list = []
    gene_list = []
    pqtl_df = pqtl_df.loc[pqtl_df["FDR"]<fdr]
    #pqtl_df = pqtl_df.loc[pqtl_df["pvalue"]<0.05]
    for index in pqtl_df.index:
        snps = pqtl_df.loc[index,"snps"]
        gene = pqtl_df.loc[index,"gene"]
        if snps in autosome_gene_tss_df.index and gene in autosome_gene_tss_df.index:
            snps_chrom = gene_tss_df.loc[snps,"chromosome"]
            gene_chrom = gene_tss_df.loc[gene,"chromosome"]
            snps_location = chrom_dict[snps_chrom] + gene_tss_df.loc[snps,"tss"]
            gene_location = chrom_dict[gene_chrom] + gene_tss_df.loc[gene,"tss"]
            if pqtl_df.loc[index,"beta"] >0 :
                color="red"
            else:
                color="blue"
            color_list.append(color)
            pvalue_list.append(pqtl_df.loc[index,"FDR"])
            pair_index_list.append(index)
            snps_location_list.append(snps_location)
            gene_location_list.append(gene_location)
            snps_list.append(snps)
            gene_list.append(gene)
        else:
            pass
    final_pair_df = pd.DataFrame([snps_list,gene_list,snps_location_list,gene_location_list,color_list,pvalue_list,-np.log10(pvalue_list)],
                             index=["mutation","protein","mutation_site","protein_loation","color","pvalue","-log10pvalue"]).T

    plt.figure(figsize=(10,10),dpi=300)
    plt.subplot(211)
    color_list = ["orange","white"]
    m=0
    for i in chrom_dict.keys():
        plt.bar(chrom_dict[i],10**7.5,width = chromo_size_df.loc[i,"length"],bottom=-10**7.5,align="edge",color=color_list[(m+1)%2])
        plt.bar(0,chromo_size_df.loc[i,"length"],bottom=chrom_dict[i], width=10**7.3,align="edge",color=color_list[(m+1)%2])
        if (m+1)%2 == 0:
            plt.text(chrom_dict[i],0,str(i).strip("chr"),va="top",ha="left",size=10)
            plt.text(0,chrom_dict[i],str(i).strip("chr"),va="bottom",ha="right",size=10)
        m+=1
    plt.scatter(final_pair_df["mutation_site"],final_pair_df["protein_loation"],
            color = final_pair_df["color"],s = [int(i*2) for i in final_pair_df["-log10pvalue"].values])
    plt.xticks([])    
    plt.box(False)
    plt.yticks([])

#    print (chrom_dict.keys())    
    total_length = chrom_dict[list(chrom_dict.keys())[-1]]+chromo_size_df.loc[list(chrom_dict.keys())[-1],"length"]

    plt.plot([0,total_length],[total_length,total_length],'k--')
    plt.plot([total_length,total_length],[0,total_length],'k--')

    plt.scatter(0,total_length+10**8,c="red",s=50)
    plt.text(10**7.5,total_length+10**8,"Beta>0",va="center",size=20)

    plt.scatter(10**8.77,total_length+10**8,c="blue",s=50)
    plt.text(10**8.8,total_length+10**8,"Beta<0",va="center",size=20)

    plt.text(10**9.1,total_length+10**8,"FDR<%s"%fdr,va="center",size=20)
#    plt.savefig("pQTL_scatter_fdr_%s.pdf"%fdr,dpi=300)
    ax=plt.subplot(212)
    import random
    def randomcolor():
        colorArr = ['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
        color = ""
        for i in range(6):
            color += colorArr[random.randint(0,14)]
        return "#"+color

    import collections
    pair_dict = collections.Counter(snps_list)
#    print(pair_dict)
    number_df = pd.DataFrame.from_dict(pair_dict,orient='index')
    number_df.columns=["associated_protein_number"]
    number_df = number_df.sort_values(by="associated_protein_number",ascending=False)

 #   plt.figure(figsize=(10,3),dpi=300)
 #   ax=plt.subplot()
    for key in pair_dict.keys():
        chrom = gene_tss_df.loc[key,"chromosome"]
        chrom_location = chrom_dict[chrom]
        key_location = chrom_location + gene_tss_df.loc[key,"tss"]
        femu = 10**7
        plt.bar(key_location/femu,pair_dict[key],color="black",edgecolor="black",width=0.85)
    
        if pair_dict[key] > 25:
            plt.text(key_location/femu,pair_dict[key],key,size=10,c=randomcolor())
    
    color_list = ["orange","white"]

    m=0    
    for i in chrom_dict.keys():
        plt.bar(chrom_dict[i]/femu,2,width = chromo_size_df.loc[i,"length"]/femu,align="edge",color=color_list[(m+1)%2])
        if (m+1)%2 == 0:
            plt.text(chrom_dict[i]/femu,0,str(i).strip("chr"),va="top",ha="left",size=10)
        m+=1
        
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    plt.xlim(0,total_length/femu)
    plt.xticks([])  
    plt.ylabel("Count",size=10)
    plt.xlabel("pQTL location",size=10)
    number_df.to_csv("mutation_associated_protein_count.txt",sep="\t")
    plt.savefig("pQTL_results_FDR_%s.pdf"%fdr,dpi=300)

def delete_intermediate_files():
    os.system("rm results.txt")
    os.system("rm z_df.txt")
    os.system("rm mutation_df.txt")
    os.system("rm covariates.txt")

@celery.task(name='async_pqtl',bind=True)
def async_pqtl(self,mutation_f,z_statistic_matrix,covariates_f,gene_tss_location,chr_length,fdr,outdir):
    try:
        with time_limit(60*60*6):
            mutation_df = pd.read_csv(mutation_f,sep="\t",index_col=0)
            z_df = pd.read_csv(z_statistic_matrix,sep='\t',index_col=0)
            covariates_df = pd.read_csv(covariates_f,sep="\t",index_col=0)

            overlap_sample = list(set(mutation_df.columns)&set(z_df.columns)&set(covariates_df.columns))

            status = [("Only %s common samples in three input files were included in the analysis."%(str(len(overlap_sample))))]
            self.update_state(state='PROGRESS',meta={'status': status})

            half_ = int(len(overlap_sample)/2)
        #    print('%s Proteins detected in at least half of the samples are included in the mutation-association analysis.')
            

            mutation_df = mutation_df[overlap_sample]
            z_df = z_df[overlap_sample]
            z_df = z_df.dropna(thresh=half_)
            status.append('%s Proteins detected in at least half of the samples are included in the mutation-association analysis.'%(str(len(z_df))))
            self.update_state(state='PROGRESS',meta={'status': status})
            
            covariates_df = covariates_df[overlap_sample]

            os.chdir(outdir)

            mutation_df.to_csv("mutation_df.txt",sep="\t")
            z_df.to_csv("z_df.txt",sep="\t")
            covariates_df.to_csv("covariates.txt",sep="\t")
            
            SNP_file_name = "mutation_df.txt"
            expression_file_name = "z_df.txt"
            covariates_file_name = "covariates.txt"
            status.append("Starting pQTL analysis.")
            self.update_state(state='PROGRESS',meta={'status': status})
            pqtl_results = prepare_pQTL(SNP_file_name,expression_file_name,covariates_file_name,outdir)
            status.append("pQTL analysis was done.")
            self.update_state(state='PROGRESS',meta={'status': status})
            os.chdir(outdir)
            all_pqtl_bh_df = add_bh_pvalue(z_df,pqtl_results,fdr=fdr)
            status.append("Generating pQTL results chart.")
            self.update_state(state='PROGRESS',meta={'status': status})
            generate_figure(gene_tss_location,chr_length,all_pqtl_bh_df,fdr=fdr)
            delete_intermediate_files()

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

