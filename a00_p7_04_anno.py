import pandas as pd
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt

#------------------------------------------
f_ada='./out/a00_p7_03_cluster/clus.h5ad'
fd_out='./out/a00_p7_04_anno'

dic_cmap={'IHC': 'orange', 
		  'OHC': 'Purple', 
		  'Pillar': 'Blue',
		  'Deiter': 'Green',
		  'Others': 'grey'}

#------------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

#------------------------------------------
def anno_ada(ada, dic_anno):
	#0. get number of clusters
	n=ada.obs['leiden'].astype('int').max()
	#1.pp
	l_anno=list(range(n+1))   #hold cell annotation
	l_cell=[]                 #hold cell order
	l_temp=[]                 #hold all cluster number to check duplicate
	#2. replace cell name
	for cell, l_clus in dic_anno.items():
		#1. append
		if len(l_clus)==0:
			continue
		l_cell.append(cell)
		l_temp.extend(l_clus)
		#2. replace by cell name
		l_anno=[(i, cell)[i in l_clus] for i in l_anno]
	#3. chekc duplicate clus
	for i in l_temp:
		if l_temp.count(i)>1:
			print(f'duplicate cluster number: {i}\n\n')
			return ada
	#4. add Unknown
	l_anno=[(i, 'Others')[isinstance(i, int)] for i in l_anno]
	#5. annotate adata
	ada.obs['anno']=ada.obs['leiden'].astype('int')
	ada.obs['anno']=ada.obs['anno'].replace(list(range(n+1)), l_anno)
	ada.obs['anno']=pd.Categorical(ada.obs['anno'], categories=l_cell+['Others'], ordered=True)
	return ada

def plot_anno(ada, f_out, title=None, dic_cmap=dic_cmap, sz=(8,6), col='anno'):
	#1. get cmap
	idx_cell=ada.obs[col].cat.categories
	cmap=[dic_cmap[i] for i in idx_cell]
	#2. plot
	ax=sc.pl.umap(ada, color=[col], s=8, alpha=1, show=False, legend_loc='right margin', legend_fontsize=10, legend_fontweight='medium', frameon=False, palette=cmap)
	#3. adjust
	fig = plt.gcf()
	fig.set_size_inches(sz)
	ax.set_title(title, fontsize=20, pad=10, weight='semibold')
	ax.xaxis.labelpad = 10
	ax.yaxis.labelpad = 10
	plt.legend(loc=(1.01, 0.1), frameon=False, prop={'size': 15, 'weight': 'semibold'})
	#4. save
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return


####################################################################
ada=sc.read(f_ada)

dic_anno={'IHC': [9], 
		  'OHC': [6],
		  'Pillar': [10],
		  'Deiter': [8]}

#anno
ada=anno_ada(ada, dic_anno)
ada.write(f'{fd_out}/p7.h5ad')

plot_anno(ada, f'{fd_out}/p7.png', title='')




















