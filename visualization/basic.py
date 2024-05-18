import numpy as np
import glob
import os
import pandas as pd
import scanpy as sc
import json
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
import seaborn as sns
import ipywidgets as widgets
from plotly.subplots import make_subplots
import plotly.graph_objs as go
from IPython.display import HTML, display
import ipywidgets as widgets
from plotly.offline import plot, iplot, init_notebook_mode
import plotly.io as pio


def basic_input(input_file):
    adata = sc.read(input_file)
    
    return adata
    
    
def basic_computing(adata):
    seed=0
    sc.tl.pca(adata, random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=15)
    sc.tl.umap(adata, random_state=seed)
    sc.tl.louvain(adata, random_state=seed) 
    
    return adata


def easy_heatmap(adata, output_df, target, target_name, project_name, fig_num):
    hsa_cols = [col for col in adata.obs.columns if col.startswith('hsa')] #kegg pathway
    hsa_df = adata.obs[hsa_cols].copy()
    hsa_df['gene_name'] = adata.obs['gene_name'].copy()
    m_hsa_df = pd.merge(output_df['gene_name'], hsa_df, on='gene_name', how='left')
    m_hsa_df = m_hsa_df.set_index('gene_name').iloc[::-1]

    m_hsa_df.to_csv("../figures/{}/heatmap_results/{}_sorted_kegg_pathway_names.csv".format(project_name, target))

    fig, ax = plt.subplots(figsize=(120,150))
    cmap = (['#440154FF','#FDE725FF'])
    ax.set_yticks(np.arange(len(m_hsa_df.index)))
    ax.set_yticklabels(m_hsa_df.index, rotation=90, size=60)
    ax.set_xticks(np.arange(len(m_hsa_df.columns)))
    ax.set_xticklabels(m_hsa_df.columns, rotation=90, size=80, ha="center")
    ax.text(-0.05, 1.01, 'Figure S6-{}-B {} {}'.format(fig_num, 'Adipocyte', target_name), fontsize=120, va='bottom', ha='left', transform=ax.transAxes)

    sns.heatmap(m_hsa_df, cmap=cmap, ax=ax, cbar=False)
    plt.subplots_adjust(left=0.15, right=0.85, bottom=0.15, top=0.85)
    plt.savefig("../figures/{}/heatmap_results/{}.pdf".format(project_name, target), bbox_inches="tight")

def easy_dendrogram(adata, target, target_name, project_name, fig_num):
    a_euc = linkage(adata.X, metric = 'euclidean', method = 'single')
    fig = plt.figure(figsize=(15, 20))  # figを取得
    dendrogram(a_euc, orientation='right',
               labels=list(adata.obs['gene_name'].to_list()), color_threshold=0.01) 
    dendro_labels = [int(x) for x in dendrogram(a_euc, no_plot=True)['ivl']]
    sorted_gene_names = adata.obs['gene_name'][dendro_labels]
    sorted_df = pd.DataFrame({'gene_name': sorted_gene_names, 'accession_id': adata.obs.index[dendro_labels],
                              'gene_name_bait': target_name, 'accession_id_bait': target})

    output_df = pd.merge(sorted_df, adata.obs[['gene_name','{}_model_1_multimer_score(iptm+ptm)'.format(target),'{}_model_1_multimer_score(iptm)'.format(target),                                     '{}_model_1_multimer_score(ptm)'.format(target)]], on='gene_name', how='left')
    output_df.to_csv("../figures/{}/dendrogram_results/{}_sorted_gene_names.csv".format(project_name, target), index=False)

    plt.tick_params(bottom=False, labelbottom=False)
    plt.title('Figure S6-{}-C {} {}'.format(fig_num, 'Adipocyte', target_name), fontsize=16, ha='left', va='top', x=-0.1, y=1.02)
    
    plt.subplots_adjust(top=0.95)  # Add this line to adjust the top margin.
    plt.savefig("../figures/{}/dendrogram_results/{}.pdf".format(project_name, target))
    plt.show()


    easy_heatmap(adata, output_df, target, target_name, project_name, fig_num)
    

def on_click(trace, points, state):
    gene_name = adata_target.obs['gene_name'][point]
    info = f"{gene_name}'s info"
    display(Javascript(f"alert('{info}');"))

        
def adata_target_view(adata, target, project_name, fig_num):
    adata_target = adata[adata.obs['{}_model_1_multimer_score(iptm+ptm)'.format(target)] != 'None']
    adata_target = adata_target[adata_target.obs['{}_model_1_multimer_score(iptm+ptm)'.format(target)] != 'error']
    obs_copy = adata_target.obs.copy()
    obs_copy['{}_model_1_multimer_score(iptm+ptm)'.format(target)] = obs_copy['{}_model_1_multimer_score(iptm+ptm)'.format(target)].astype('float64')
    obs_copy['{}_model_1_multimer_score(iptm)'.format(target)] = obs_copy['{}_model_1_multimer_score(iptm)'.format(target)].astype('float64')
    obs_copy['{}_model_1_multimer_score(ptm)'.format(target)] = obs_copy['{}_model_1_multimer_score(ptm)'.format(target)].astype('float64')
    adata_target.obs = obs_copy
    
    #target only
    seed=0
    sc.tl.pca(adata_target, random_state=seed)
    if 50 > len(adata_target):
        sc.pp.neighbors(adata_target, n_neighbors=3)
    else:
        sc.pp.neighbors(adata_target, n_neighbors=4)

    sc.tl.umap(adata_target, random_state=seed)
    sc.tl.louvain(adata_target, random_state=seed) 
    print(adata[adata.obs['ID'] == target].obs['gene_name'][0])
    gene_name =  adata[adata.obs['ID'] == target].obs['gene_name'][0]
    easy_dendrogram(adata_target, target, gene_name, project_name, fig_num)
    print(len(adata_target))

    init_notebook_mode()
    fig = go.Figure(data=go.Scatter(x=adata_target.obsm['X_umap'][:, 0], y=adata_target.obsm['X_umap'][:, 1],mode='markers+text', text=adata_target.obs.index))
    color_array = adata_target.obs['{}_model_1_multimer_score(iptm+ptm)'.format(target)]
    fig = go.Figure(data=go.Scatter(
            x=adata_target.obsm['X_umap'][:, 0],
            y=adata_target.obsm['X_umap'][:, 1],
            mode='markers+text',
            text=['<a href="https://www.uniprot.org/uniprot/{}" target="_blank">{}</a>'.format(i, j) for i, j in zip(adata_target.obs.index, adata_target.obs['gene_name'])],
            marker=dict(
                size=20,
                color=color_array, 
                colorbar=dict(title='score(iptm+ptm)'),
                # colorbar=dict(title='{}_model_1_multimer_score(iptm+ptm)'.format(target)),
                colorscale='ylgn',
                showscale=True
            ),
            textfont=dict(color='#447ADB')
        ),layout=go.Layout(
            title='Figure S6-{}-A {} {}'.format(fig_num, 'Adipocyte', gene_name),
            plot_bgcolor="white",
            xaxis=dict(linecolor='black',mirror=True, title='UMAP 1', tickfont=dict(color='white')),
            yaxis=dict(linecolor='black',mirror=True, title='UMAP 2', tickfont=dict(color='white')),
            # yaxis=dict(showticklabels=False),
            width=1000,
            height=1000,
        
        ))
    


    fig.update_traces(mode='markers+text', textfont=dict(size=7))
    fig.data[0].on_click(on_click)
    plot(fig, filename='../figures/{}/umap_results/{}_only_umap.html'.format(project_name, target))
    iplot(fig, filename='../figures/{}/umap_results/{}_only_umap.html'.format(project_name, target))
    
    pio.write_image(fig, '../figures/{}/umap_results/{}_only_umap.pdf'.format(project_name, target))

    