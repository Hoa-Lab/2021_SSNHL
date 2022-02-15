import pandas as pd
import scanpy as sc
import numpy as np
import anndata as ad
from pathlib import Path

#--------------------------------------------
prj='/mnt/c/Users/gus5/Desktop/a01_proj/a08_ssnhl'
f_in=f'{prj}/a00_raw/count/Kolla/GSE137299_P7_Normalized_Counts.txt'
fd_out='./out/a00_p7_00_load'

#-------------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

######################################################################
#pp
df=pd.read_csv(f_in, index_col=0, sep='\t')
df=df.T

#2. convert to ada
ada=ad.AnnData(df)
ada.obs['sample']='P7'
ada.write(f'{fd_out}/p7.h5ad')
	
