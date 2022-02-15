import pandas as pd
import scanpy as sc
from pathlib import Path
from scipy.stats import zscore
import json

#---------------------------------------------------------
r=0.2

f_rss='./out/a01_heatmap_00_rss/rss.csv'
f_ada='./out/a00_p7_04_anno/p7.h5ad'
fd_out='./out/a01_heatmap_01_pp'

l_cell=['IHC', 'OHC', 'Pillar', 'Dieter']

#--------------------------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
ada=sc.read(f_ada)
df_rss=pd.read_csv(f_rss, index_col=0).fillna(0)

###################################################################
#make dic
l_all=df_rss.index.tolist()
dic_gene={}

for cell in l_cell:
	dfi=df_rss.loc[df_rss[cell]>=r, :]	
	l_gene=dfi.index.tolist()
	dic_gene[cell]=l_gene
	l_all=[i for i in l_all if i not in l_gene]

dic_gene['Others']=l_all

#save
f_out=f'{fd_out}/gene.json'
with open(f_out, 'w') as f:
	json.dump(dic_gene, f)
