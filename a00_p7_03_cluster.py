import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt
from pathlib import Path

#------------------------------------------
res=0.5

f_ada='./out/a00_p7_01_embed/embed.h5ad'
fd_out='./out/a00_p7_03_cluster'

#-----------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)


####################################################################
#cluster
ada=sc.read(f_ada)
sc.tl.leiden(ada, resolution=res, random_state=42) 
ada.write(f'{fd_out}/clus.h5ad')

#plot
ada=sc.read(f'{fd_out}/clus.h5ad')
title=None

ax=sc.pl.umap(ada, color=['leiden'], legend_loc='on data', s=12, alpha=0.8, palette='tab20', show=False, frameon=False)
fig = plt.gcf()
fig.set_size_inches((6,6))

#4. adjust
ax.set_title(title, fontsize=24, pad=10, weight='semibold')

#5. save plot
plt.tight_layout()
plt.savefig(f'{fd_out}/cluster.png', dpi=300)
plt.close()
