import pandas as pd
from pathlib import Path
import scanpy as sc
import anndata as ad

#-----------------------------------------------------
prj='/mnt/c/Users/gus5/Desktop/a01_proj/a08_ssnhl'
f_in=f'{prj}/a00_raw/count/ranum/GSE114157_p15_Expression_Matrix.csv'
fd_out='./out/a02_ranum_00_pp'


#----------------------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)


#####################################################################
#load
df=pd.read_csv(f_in, index_col=0)
df=df.sort_index()
df=df.loc[~df.index.duplicated(), :]

#sort
df=df.T
df=df.sort_index()

#convert
ada=ad.AnnData(df)

#norm
sc.pp.normalize_total(ada, exclude_highly_expressed=True)
sc.pp.log1p(ada)

ada.write(f'{fd_out}/ada.h5ad')

