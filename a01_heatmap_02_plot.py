import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import json
from scipy.stats import zscore
import numpy as np

#----------------------------------------------------------------
f_gl='./out/a01_heatmap_01_pp/gene.json'
f_ada='./out/a00_p7_04_anno/p7.h5ad'
f_meta='./raw/meta_info.csv'
fd_out='./out/a01_heatmap_02_plot'

l_cell=['IHC', 'OHC', 'Pillar', 'Deiter']

#--------------------------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

#color map
df_meta=pd.read_csv(f_meta, index_col=0)
dic_cmap=df_meta.to_dict()['color']

#---------------------------------------------------------
def plot_top(df, f_out, l_cell=l_cell, sz=(6,6), dic_cmap=dic_cmap):
	cmap=[dic_cmap[i] for i in l_cell]
	#clean
	df=df.loc[:, ['anno']]
	df['anno']=df.anno.cat.codes
	df=df.loc[:, ['anno']].T
	#plot
	fig, ax=plt.subplots(figsize=sz)
	ax=sns.heatmap(df, cmap=cmap, cbar=False)	
	#adjust
	ax.xaxis.label.set_visible(False)
	ax.yaxis.label.set_visible(False)
	plt.xticks([])
	plt.yticks([])
	#save
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()	
	return

def plot_hm(df, f_out, size=(10,15), vmax=None, vmin=-0.3, y=14):
	df=df.drop_duplicates()
	#2. heatmap
	fig, ax=plt.subplots(figsize=size)
	ax=sns.heatmap(df, cmap='Purples', vmax=vmax, vmin=vmin, cbar=True)
	#3. adjust
	ax.xaxis.label.set_visible(False)
	ax.yaxis.label.set_visible(False)
	plt.xticks([])
	plt.yticks(np.arange(0.5, df.shape[0]+0.5, 1), df.index.tolist(), fontsize=y, rotation=0, weight='medium')
	#4. save
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return
	
	
############################################################################
#load
ada=sc.read(f_ada)
with open(f_gl, 'r') as f:
	dic_gene=json.load(f)

#get gene list
l_gene=[]
l=l_cell+['Others']
for cell in l:
	l_gene.extend(dic_gene[cell])

#make df
ada=ada[:, l_gene]
df=pd.DataFrame(ada.X, index=ada.obs.index, columns=ada.var.index)
df['anno']=ada.obs['anno']

df=df.loc[df['anno'].isin(l_cell), :]
df['anno']=pd.Categorical(df['anno'], categories=l_cell, ordered=True)
df=df.sort_values('anno')

##plot top
#f_out=f'{fd_out}/top.png'
#plot_top(df, f_out)

#plot
df=df.drop('anno', axis=1).T

f_out=f'{fd_out}/hm_bar.png'
plot_hm(df, f_out)





