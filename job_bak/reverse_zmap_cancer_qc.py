# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 10:49:24 2020

@author: guixiuqi
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import six
def render_mpl_table(data, col_width=6.0, row_height=0.625, font_size=14,
                     header_color='#5B9BD5', row_colors=['#f1f1f2', 'w'], edge_color='black',
                     bbox=[0, 0, 1, 1], header_columns=0,
                     ax=None, **kwargs):
    if ax is None:
        size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
        fig, ax = plt.subplots(figsize=size)
        ax.axis('off')

    mpl_table = ax.table(cellText=data.values, cellLoc = 'center',bbox=bbox, colLabels=data.columns, **kwargs)

    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(font_size)

    for k, cell in six.iteritems(mpl_table._cells):
        cell.set_edgecolor(edge_color)
        if k[0] == 0 or k[1] < header_columns:
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
        else:
            cell.set_facecolor(row_colors[k[0]%len(row_colors) ])
    ax.axis('off')
    return ax

def data_qc(final_intensity_df,fig_dir=""):
    sample_num = len(final_intensity_df.columns)
    total_gene_number = len(final_intensity_df) 
    replace_intensity_df = final_intensity_df.replace(0,np.nan)
    
    detected_sample_number_each_gene = replace_intensity_df.count(axis=1)
    detected_gene_number_each_sample = replace_intensity_df.count(axis=0)
    
    boxprops = dict(linestyle = '-', linewidth = 2, color = 'darkorange')
    medianprops = dict(linestyle = '-', linewidth = 1.5, color = 'firebrick')
    flierprops = dict(marker = 'o', markerfacecolor = 'grey',
                  linestyle = 'none')
    whiskerprops = dict(linestyle = '-', linewidth = 2, color = 'grey')

    capprops = dict(linestyle='-', linewidth=2, color='grey')
    plt.figure(figsize=(3,5))
    ax1 = plt.subplot(111)
    ax1.boxplot(detected_gene_number_each_sample,
                vert = True,widths=0.7,
            boxprops = boxprops,medianprops=medianprops,whiskerprops=whiskerprops,
            capprops=capprops,flierprops=flierprops)
    ax1.set_ylabel("Detected gene numbers",size=15)
    ax1.set_xlabel("Sample",size=15)
    ax1.set_title("Number of gene detected by \nat least one sample:%s"%(str(total_gene_number)),size=15)
    plt.savefig(fig_dir+"/"+"qc_boxplot.png",dpi=300,bbox_inches="tight")
    plt.figure(figsize=(10,4))
    ax2=plt.subplot(111)
    
    half_detected_num = np.sum(detected_sample_number_each_gene > 0.5*sample_num)
    three_quarters_detected_number = np.sum(detected_sample_number_each_gene > 0.75*sample_num)
    all_detected_number = np.sum(detected_sample_number_each_gene > (sample_num-1))
    
    col_labels = ['Description','Counts']
    table_values = [["Proteins quantified in at least\nhalf samples at gene level",half_detected_num],
                    ["Proteins quantified in at least three\nquarters samples at gene sample",
                     three_quarters_detected_number],
                    ["Proteins quantified in all \nsamples at gene level",all_detected_number]]
    table_df = pd.DataFrame(table_values,columns=col_labels)  
    ax = render_mpl_table(table_df, header_columns=0, col_width=6.0,ax=ax2,font_size = 16)  
    plt.savefig(fig_dir+"/"+"qc_table.png",dpi=300,bbox_inches="tight")
    plt.close("all")