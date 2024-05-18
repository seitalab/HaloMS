import numpy as np
import glob
import os
import pandas as pd
import scanpy as sc
import json
 
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt

import plotly.graph_objs as go
from plotly.offline import plot, init_notebook_mode


emb_name = '../analysis/data/project_af2_mono__emb.csv'
full_res = '../analysis/data/adsc_af2_score_full_lists.csv'
Kegg_list_path = '../analysis/data/adsc_kegg_lists.csv'
adata_output = '../analysis/data/adata_adsc.h5ad'
bait_list = "../analysis/data/bait_list.csv"



adata = sc.read_csv(emb_name)

qqlist = pd.read_csv(bait_list)['AccB'].to_list()
rlist = pd.read_csv(full_res)
rlist_full = pd.DataFrame(columns=['ID'])



for i in qqlist:
    rlist1_A = rlist[rlist['A_ID'] == i][['B_ID','model_1_multimer_score(iptm+ptm)', 'model_1_multimer_score(iptm)', 'model_1_multimer_score(ptm)']]
    rlist1_B = rlist[rlist['B_ID'] == i][['A_ID','model_1_multimer_score(iptm+ptm)', 'model_1_multimer_score(iptm)', 'model_1_multimer_score(ptm)']]
    rlist1_A = rlist1_A.rename(columns={'B_ID': 'ID'})
    rlist1_B = rlist1_B.rename(columns={'A_ID': 'ID'})
    rlist1AB = pd.concat([rlist1_A, rlist1_B])
    rlist1AB = rlist1AB.rename(columns={'model_1_multimer_score(iptm+ptm)' : i + '_model_1_multimer_score(iptm+ptm)'})
    rlist1AB = rlist1AB.rename(columns={'model_1_multimer_score(iptm)' : i + '_model_1_multimer_score(iptm)'})
    rlist1AB = rlist1AB.rename(columns={'model_1_multimer_score(ptm)': i + '_model_1_multimer_score(ptm)'})
    
    rlist_full = pd.merge(rlist_full, rlist1AB, on='ID', how='outer')
    

count_per_column = rlist_full.count()
total_distribution = count_per_column.value_counts()


df = rlist.copy()
df_a = df[['A_ID', 'preA_ID']]
df_a = df_a.rename(columns={'A_ID':'ID', 'preA_ID':'gene_name'})
df_b = df[['B_ID', 'preB_ID']]
df_b = df_b.rename(columns={'B_ID':'ID', 'preB_ID':'gene_name'})
gene_name_list = pd.concat([df_a, df_b]).drop_duplicates().reset_index(drop=True)

sample_id = adata.obs.index.tolist()
sample_id_list = pd.DataFrame(sample_id, columns=['ID'])
result_df = pd.merge(gene_name_list, sample_id_list, on='ID', how='left')
result_df = pd.merge(result_df, rlist_full, on='ID', how='left')


obs_df = adata.obs
obs_df['ID'] = obs_df.index.tolist()
merged_df = pd.merge(obs_df, result_df, on='ID', how='left')
merged_df.set_index('ID', inplace=True)
merged_df.fillna('None', inplace=True)

for name in merged_df.columns.tolist():
    adata.obs[name] = merged_df[name].values
    
adata = adata[adata.obs.index.isin(gene_name_list['ID']), :]


#add kegg info
obs_df = adata.obs['gene_name'].copy()
kegg_list = pd.read_csv(Kegg_list_path)
kegg_list = kegg_list.rename(columns={'Unnamed: 0':'gene_name'})

merged_df = pd.merge(obs_df, kegg_list, on='gene_name', how='left')
merged_df.fillna(0, inplace=True)

adata = adata.copy()
for name in merged_df.columns.tolist():
    adata.obs[name] = merged_df[name].values

adata.write(adata_output)