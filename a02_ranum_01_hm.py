import pandas as pd
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#-----------------------------------------------
f_in='./out/a02_ranum_00_pp/ada.h5ad'
f_gl='./raw/genelist/gl_ssnhl_2020-10-07.txt'
f_meta='./raw/meta_info.csv'
fd_out='./out/a02_ranum_01_hm'
l_cell=['IHC', 'OHC', 'Deiter']
#cmap_top=['#4287f5', '#f5a142', '#4bf542']

#-----------------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
ada=sc.read(f_in)

#color map
df_meta=pd.read_csv(f_meta, index_col=0)
dic_cmap=df_meta.to_dict()['color']

#----------------------------------------------
def plt_top(df, f_out):
	cmap=[dic_cmap[i] for i in l_cell]
	fig, ax=plt.subplots(figsize=(6,5))
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


def sort_df(df):
	df=df.loc[~df.index.duplicated()]
	#groupby
	dfi=df.copy().T
	dfi['cell']=dfi.index.str.split('_').str[0]
	dfi=dfi.groupby('cell').mean()
	#sort
	dfi=dfi.T
	dfi=dfi.sort_values('Deiter', ascending=False)
	#sort original
	l_order=dfi.index
	df=df.reindex(l_order)
	return df


def plt_hm(df, f_out, sz=(8,10), hide_gene=False, vmax=8, vmin=-0.1, y=10, cbar=False):
	#plot
	fig, ax=plt.subplots(figsize=sz)
	ax=sns.heatmap(df, cmap='Purples', vmax=vmax, vmin=vmin, cbar=cbar)
	#adjust
	ax.xaxis.label.set_visible(False)
	ax.yaxis.label.set_visible(False)	
	plt.xticks([])
	plt.yticks(np.arange(0.5, df.shape[0]+0.5, 1), df.index.tolist(), fontsize=y, rotation=0, weight='medium')
	if hide_gene:
		plt.yticks([])
	#save
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return


###############################################################
#gene list
l_gene=Path(f_gl).read_text().strip().split('\n')
l_gene=[i.capitalize() for i in l_gene]
l_gene=[i for i in l_gene if i in ada.var.index]


#make df
ada=ada[:, l_gene]
df=ada.to_df()

#sort df
df['cell']=df.index.str.split('_').str[0]
df['cell']=pd.Categorical(df['cell'], categories=l_cell, ordered=True)
df=df.sort_values('cell')

#split
df_top=df.loc[:, ['cell']].copy()
df=df.drop('cell', axis=1)
df=df.T


#top bar
df_top['cell']=df_top['cell'].cat.codes
df_top=df_top.T

f_out=f'{fd_out}/top.png'
plt_top(df_top, f_out)


##------------------------------------------------
##filter low exp
#df=df.loc[df.sum(axis=1)>0.5, :]

##sort df
#df=sort_df(df)

##plot
#f_out=f'{fd_out}/hm.png'
#plt_hm(df, f_out)





