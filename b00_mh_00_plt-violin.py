import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

#----------------------------------------------------
prj='/mnt/c/Users/gus5/Desktop/a01_proj/a08_ssnhl'
f_sc=f'{prj}/a00_raw/count/sc_P30/sc_P30.h5ad'
f_sn=f'{prj}/a00_raw/count/sNuc_P30/sNuc_P30.h5ad'
f_p7='./out/a00_p7_04_anno/p7.h5ad'
f_ranum='./out/a02_ranum_00_pp/ada.h5ad'
fd_out='./out/b00_mh_00_plt-violin'

#----------------------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

#----------------------------------------------------
def get_df(ada, l_cell, l_gene, col='anno'):
	ada=ada[ada.obs[col].isin(l_cell), :]
	ada=ada[:, l_gene]
	#make df
	df=ada.to_df()
	df=df.merge(ada.obs, left_index=True, right_index=True)
	#transform df
	l_col=['cnt', 'Cell']
	df1=df.loc[:, [l_gene[0], 'anno']]
	df2=df.loc[:, [l_gene[1], 'anno']]
	df1.columns=l_col
	df2.columns=l_col
	df1['Gene']=l_gene[0]
	df2['Gene']=l_gene[1]
	df=pd.concat([df1, df2])
	#make categories
	df['Cell']=pd.Categorical(df['Cell'], categories=l_cell, ordered=True)
	df['Gene']=pd.Categorical(df['Gene'], categories=l_gene, ordered=True)
	return df


def plt_viol(df, f_out, sz=(9,6), title=None, scale='width'):
	#plot
	sns.set()
	fig, ax=plt.subplots(figsize=sz)
	ax=sns.violinplot(x='Cell', y='cnt', hue='Gene', data=df, cut=0, scale=scale)
	#adjust
	ax.set_title(title, weight='semibold', fontsize='22', pad=15)
	plt.xlabel('')
	plt.ylabel('Expression level\n (normalized counts)', fontsize=20, labelpad=15, weight='semibold')
	plt.xticks(fontsize=20, rotation=0, weight='semibold')
	plt.legend(loc=9, frameon=True, ncol=2, prop={'size': 15, 'weight': 'semibold'})
	#save
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return


#################################################################
#-------------------sc-------------------------
sample='sc'
title='P30 Stria Vascularis (scRNA-Seq)'
ada=sc.read(f_sc)

l_cell=['Marginal', 'Intermediate', 'Basal']
l_gl1=['Gpx1', 'Gpx3']
l_gl2=[ 'Nr3c1', 'Nr3c2']

#make df
df1=get_df(ada, l_cell, l_gl1)
df2=get_df(ada, l_cell, l_gl2)

#plot
f_out=f'{fd_out}/{sample}_1.png'
plt_viol(df1, f_out, title=title)

f_out=f'{fd_out}/{sample}_2.png'
plt_viol(df2, f_out, title=title)


#-------------------sn-------------------------
sample='sn'
title="P30 Stria Vascularis (snRNA-Seq)"
ada=sc.read(f_sn)

l_cell=['Marginal', 'Intermediate', 'Basal']
l_gl1=['Gpx1', 'Gpx3']
l_gl2=[ 'Nr3c1', 'Nr3c2']

#make df
df1=get_df(ada, l_cell, l_gl1)
df2=get_df(ada, l_cell, l_gl2)

#plot
f_out=f'{fd_out}/{sample}_1.png'
plt_viol(df1, f_out, title=title)

f_out=f'{fd_out}/{sample}_2.png'
plt_viol(df2, f_out, title=title)


#----------------p7------------------------------
sample='p7'
title="P7 Organ of Corti (scRNA-Seq)"
#sz=(12, 6)
ada=sc.read(f_p7)

l_cell=['IHC', 'OHC', 'Pillar', 'Deiter']
l_gl1=['Gpx1', 'Gpx3']
l_gl2=[ 'Nr3c1', 'Nr3c2']

#make df
df1=get_df(ada, l_cell, l_gl1)
df2=get_df(ada, l_cell, l_gl2)

#plot
f_out=f'{fd_out}/{sample}_1.png'
plt_viol(df1, f_out, title=title)#, sz=sz)

f_out=f'{fd_out}/{sample}_2.png'
plt_viol(df2, f_out, title=title)#, sz=sz)


#--------------Ranum--------------------------
sample='ranum'
title='P15 Organ of Corti (scRNA-Seq)'
ada=sc.read(f_ranum)
ada.obs['anno']=ada.obs.index.str.split('_').str[0]

l_cell=['IHC', 'OHC', 'Deiter']
l_gl1=['Gpx1', 'Gpx3']
l_gl2=[ 'Nr3c1', 'Nr3c2']

#make df
df1=get_df(ada, l_cell, l_gl1)
df2=get_df(ada, l_cell, l_gl2)

#plot
f_out=f'{fd_out}/{sample}_1.png'
plt_viol(df1, f_out, title=title)

f_out=f'{fd_out}/{sample}_2.png'
plt_viol(df2, f_out, title=title)






