import pandas as pd
import scanpy as sc
from pathlib import Path
from pyscenic.rss import regulon_specificity_scores

#-------------------------------------------
prj='/mnt/c/Users/gus5/Desktop/a01_proj/a08_ssnhl'
f_gl=f'{prj}/a00_raw/genelist/gl_ssnhl_2020-10-07.txt'
f_ada='./out/a00_p7_04_anno/p7.h5ad'
fd_out='./out/a01_heatmap_00_rss'

#------------------------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
l_gene=Path(f_gl).read_text().strip().split('\n')
l_gene=[i.capitalize() for i in l_gene]

##############################################################
#load
ada=sc.read(f_ada)
l_gl=[i for i in l_gene if i in ada.var.index]
l_gl=list(set(l_gl))

#make df
df=pd.DataFrame(ada.X, columns=ada.var.index, index=ada.obs.index)
df=df.loc[:, l_gl]

#rss
l_anno=ada.obs['anno']


df_rss=regulon_specificity_scores(df, l_anno).T

df_rss.to_csv(f'{fd_out}/rss.csv')
