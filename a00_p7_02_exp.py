import scanpy as sc
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from multiprocessing import Pool

#----------------------------------------------
f_ada='./out/a00_p7_01_embed/embed.h5ad'
fd_out='./out/a00_p7_02_exp'

l_gene=['Anxa5', 'Fabp7', 'S100a6', 'Pla2g7', 'Prox1', 'Tuba1b', 'Cd44', 'Slc17a8', 'Slc26a5']

#----------------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
ada=sc.read(f_ada)

#-----------------------------------------------
def plot_exp(gene, fd_out=fd_out, ada=ada, vmin=-0.8, vmax=None, raw=False, sz=(6,5)):
	try:
		#plot
		ax=sc.pl.umap(ada, color=[gene], show=False, s=10, alpha=1, frameon=False, cmap='BuPu', use_raw=raw, vmin=vmin, vmax=vmax)
		#2. adjust
		fig = plt.gcf()
		fig.set_size_inches(sz)
		ax.set_title(gene, fontsize=16, pad=10, weight='semibold')
		#3. save
		plt.tight_layout()
		plt.savefig(f'{fd_out}/{gene}.png', dpi=300)
		plt.close()
	except Exception as e:
		plt.close()
		print(e)
	return


#########################################################
with Pool(6) as pool:
	pool.map(plot_exp, l_gene)

