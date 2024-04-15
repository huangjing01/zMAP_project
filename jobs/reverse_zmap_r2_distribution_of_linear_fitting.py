import glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
#import matplotlib
#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
#matplotlib.rcParams['font.family'] ='Arial'

def r2_linear_boxplot(outdir):
	files = glob.glob( outdir + "/*/linear_model_r2.txt")
	r2_list=[]
	for file in files:
	    file_df = pd.read_csv(file,header=None)
	    r2 = list(file_df[0].values)
	    r2_list += r2
        
	plt.figure(figsize=(2,6),dpi=300)

	boxprops = dict(linestyle = '-', linewidth = 2, color = 'darkorange')
	medianprops = dict(linestyle = '-', linewidth = 1.5, color = 'firebrick')
	flierprops = dict(marker = 'o', markerfacecolor = 'grey',linestyle = 'none')
	whiskerprops = dict(linestyle = '-', linewidth = 2, color = 'grey')
	capprops = dict(linestyle='-', linewidth=2, color='grey')

 
	plt.boxplot(r2_list,
                vert = True,widths=0.7,showfliers=False,
            boxprops = boxprops,medianprops=medianprops,whiskerprops=whiskerprops,
            capprops=capprops,flierprops=flierprops)

	plt.xticks([1],["Sample"])

	plt.yticks([0.985,0.99,0.995,1])

	plt.ylim([0.985,1])

	plt.ylabel("R$^2$ of linear model",size=10)
	plt.savefig(outdir+"/r2_linear_model_boxplot.pdf",dpi=300,bbox_inches="tight")

#r2_linear_boxplot(outdir)
