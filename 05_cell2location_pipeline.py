import numpy as np
import pandas as pd
import cell2location
import scanpy as sc
import subprocess
import matplotlib.pyplot as plt
import anndata
import scvi
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--sample",type=str, default='1-1')
args = parser.parse_args()

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns
from scipy.sparse import csr_matrix
from cell2location.utils.filtering import filter_genes

sample = args.sample
adata_scrna_raw = sc.read_10x_h5("/fs/ess/PCON0022/guoqi/NC/Revision/Output/scRNA_6.h5")
adata_vis = sc.read_visium("/fs/ess/PCON0022/Yuzhou/NSF_2022/Spatial/" + sample)
adata_vis.obs['sample'] = list(adata_vis.uns['spatial'].keys())[0]

celltype_anno = pd.read_csv("/fs/ess/PCON0022/guoqi/NC/Revision/Output/scRNA_6_annotation.csv", index_col=0)

overlap_barcode = np.intersect1d(adata_scrna_raw.obs.index.tolist(), celltype_anno.index.tolist())
celltype_anno = celltype_anno.loc[overlap_barcode, :]
adata_scran_raw = adata_scrna_raw[overlap_barcode, :]
adata_scrna_raw.obs['celltype'] = celltype_anno['celltype'].copy()

celltype_key = 'celltype'
adata_scrna_raw = adata_scrna_raw[~adata_scrna_raw.obs[celltype_key].isin(np.array(adata_scrna_raw.obs[celltype_key].value_counts()[adata_scrna_raw.obs[celltype_key].value_counts() <=1].index))]
# remove cells and genes with 0 counts everywhere
adata_scrna_raw.var_names_make_unique()
sc.pp.filter_genes(adata_scrna_raw,min_cells=1)
sc.pp.filter_cells(adata_scrna_raw,min_genes=1)
adata_scrna_raw.obs[celltype_key] = pd.Categorical(adata_scrna_raw.obs[celltype_key])
adata_scrna_raw = adata_scrna_raw[~adata_scrna_raw.obs[celltype_key].isna(), :]
from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_scrna_raw, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_scrna_raw = adata_scrna_raw[:, selected].copy()
cell2location.models.RegressionModel.setup_anndata(adata=adata_scrna_raw,
                        # 10X reaction / sample / batch
                        # cell type, covariate used for constructing signatures
                        labels_key=celltype_key,
                       )
from cell2location.models import RegressionModel
mod = RegressionModel(adata_scrna_raw)
mod.view_anndata_setup()
mod.train(max_epochs=250, use_gpu=True)
adata_snrna_raw = mod.export_posterior(
    adata_scrna_raw, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)
# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_snrna_raw.varm.keys():
    inf_aver = adata_scrna_raw.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' for i in adata_snrna_raw.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_scrna_raw.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_snrna_raw.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_scrna_raw.uns['mod']['factor_names']
adata_vis.var_names_make_unique()
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")
# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)
mod.view_anndata_setup()
mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True)

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training'])
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)
print(adata_vis)
output_file_path = "/fs/ess/PCON0022/guoqi/NC-snrna/Spatial/Nature_MTG_deconvolution/cell2location/output"
cell2loc_results = adata_vis.obsm['q05_cell_abundance_w_sf'].copy()
cell2loc_results.columns = [c.split('q05cell_abundance_w_sf_')[1] for c in cell2loc_results.columns]
cell2loc_results = cell2loc_results.loc[:,np.unique(cell2loc_results.columns)]
cell2loc_results = cell2loc_results.div(cell2loc_results.sum(axis=1), axis=0)
cell2loc_results.to_csv(output_file_path + sample + '_Cell2location_results.csv')

