import scanpy as sc
import pandas as pd
from pathlib import Path

#---------------------------------------
f_ada='./out/a00_p7_00_load/p7.h5ad'
fd_out='./out/a00_p7_01_embed'

#----------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

####################################################
#load
ada=sc.read(f_ada)

#embed
sc.tl.pca(ada, svd_solver='arpack')
sc.pp.neighbors(ada)
sc.tl.umap(ada, n_components=2, random_state=42)

ada.write(f'{fd_out}/embed.h5ad')
