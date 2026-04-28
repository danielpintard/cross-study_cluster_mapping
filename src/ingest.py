# Description: takes in directories containing NS-Forest results and respective datasets and produces matrix, barcodes, features and cell type gene
# marker information to be used in producing SingleCellExperiment a object for FR-Match. 

# [ ] TODO: refactor to accept sample sheet from args and encode associated args with dataset runs in sheet

import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import os
import argparse
import warnings
import gc

warnings.filterwarnings('ignore')
seed = 42 

# considering submitting a sample sheet

parser = argparse.ArgumentParser(description="Ingest scRNA-seq data")
###########
parser.add_argument("--data_path", required=True, type=str, help="root/parent path to access all datasets")
parser.add_argument("--data_ids", type=str, required=True, help="Array of strings used to ID the datasets. Should include unique strings for identifying respective results folder")
# NOTE: I think this should allow for multiple inputs at once ^
parser.add_argument("--markers_path", type=str, required=True, help="root/parent path to access all datasets marker information")
parser.add_argument("--include_fscores", action='store_true', help="Indicate if an additional file containing NS-Forest F-Beta scores should be written")
parser.add_argument("--include_dendrogram", action='store_true', help="Indicate if an additional file containing a vector with the order of cell types from hiearchical clustering. When flag used, sc.tl.dendrogram is ran.")
parser.add_argument("--cluster_header", type=str, required=True, help = "Column name of adata.obs that contains cell type labels of interest")
parser.add_argument("--tmpdir", type=str, required=True, help = "Temporary space for holding intermediate files. On Biowulf, set $TMPDIR to lscratch space.")
parser.add_argument("--cxg", action="store_true", help="Indicate whether or not data is sourced from CellxGene. Omit if data not from CellxGene. This is to deal with how CellxGene organizes their adata.var")
parser.add_argument("--var_col", type=str, default="", help="Column in adata.var where gene symbols are held")

args = parser.parse_args()

data_ids = args.data_ids # from arg parse
data_path = args.data_path # data_path from argparse
tmpdir = args.tmpdir # from argparse, it'll take TMPDIR to place data in scratch space
var_col = args.var_col # coming from argparse, colname - if left empty, it is var_names
cxg = args.cxg # cxg from argparse
cluster_header = args.cluster_header #from argparse
markers_path = args.markers_path
include_fscores = args.include_fscores
include_dendrogram = args.include_dendrogram

def adata_checks(adata, cxg = cxg,):
    # ensure .var index is gene symbols
    if cxg:
        adata.var['ensembl_id'] = adata.var_names
        adata.var.index = adata.var['feature_name'].astype(str)
        adata.var_names_make_unique()
        adata.var.index.name = None
    elif var_col == "":
        var_col = adata.var_names
    else:
        var_col = var_col
        adata.var.index = adata.var[var_col].astype(str)

    ds_X = adata.X[:100].copy()
    ds_arr = ds_X.toarray() if sp.issparse(ds_X) else ds_X
    is_integer = np.allclose(ds_arr, np.round(ds_arr), rtol=1e-5, atol=1e-5)
    max_val = ds_arr.max()

    # NOTE: logic needs testing 
    # check .X to see if its transformed already
    if 'log1p' in adata.uns:
        print("Metadata Check: Found 'log1p' in adata.uns. Data is already transformed.")
    elif (not is_integer) and max_val < 30: # if matrix is not raw counts andmax value less than 30, then it is transformed
        print(f"Heuristic Check: 'log1p' metadata missing, but data appears transformed (contains floats, max={max_val:.2f}). Skipping.")
    else:
        print(f"Heuristic Check: Data appears not to be transformed (contains ints)")
        print('Transforming data...')
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    
    return adata

def make_sce_obj_files(adata, cluster_header, out_dir, data_id, mkr_info_dir, fscores_included = include_fscores, include_dendrogram = include_dendrogram):
    
    frmatch_files_dir =  os.path.join(out_dir, f'{data_id}_FRMatch_files/')
    # write cluster labels for each barcode to file
    adata.obs[cluster_header].to_csv(os.path.join(frmatch_files_dir, f'{data_id}_clusters.csv', index = False))
    
    # write cxg matrix to file - writing in chunks to combat memory spikes
    chunk_sz = 2500
    total_cells = adata.shape[0]
    with open(os.path.join(frmatch_files_dir, f'{data_id}_matrix.csv', 'w')) as f:
        f.write("," + ",".join(adata.var_names) + "\n")
        for i in range(0, total_cells, chunk_sz):
            end = min(i + chunk_sz, total_cells)   
            chunk_X = adata.X[i:end] 
            dense_chunk = chunk_X.toarray() if sp.issparse(chunk_X) else chunk_X
            chunk_df = pd.DataFrame(dense_chunk, index=adata.obs_name[i:end])
            chunk_df.to_csv(f, header=False, index=True)
    
    # write marker genes per cluster to csv
    markers = pd.read_csv(os.path.join(mkr_info_dir, f'{data_id}_results', 'tables', 'combined_markers_eval_results.csv')) 
    markers_series = pd.Series(markers['markers'].values, index = markers['clusterName'])
    markers_series.to_csv(os.path.join(frmatch_files_dir, f'{data_id}_markers.csv'))
    if include_fscores: 
       markers_fscores = pd.Series(markers['f_score'].values, index = markers['clusterName'])
       markers_fscores.to_csv(os.path.join(frmatch_files_dir, f'{data_id}_fscores.csv'))
    if include_dendrogram:
        if "X_pca" not in adata.obsm:
            print("No `X_pca` is .obsm, calculating...")
            sc.pp.pca(adata, n_comps=30, random_state=seed)
        else:
            print("`X_pca` space present...")
        
        print("running sc.tl.dendrogram...")
        sc.tl.dendrogram(adata, groupby=cluster_header, use_rep="X_pca", use_raw=False)
        dendro_file_path = os.path.join(frmatch_files_dir, f'{data_id}_dendrogram.csv')
        dendro_key = f"dendrogram_{cluster_header}"
        order = adata.uns[dendro_key]['categories_order']
        with open(dendro_file_path, 'w') as f:
            # NOTE: KEEP IN MIND that dendrogram order file being written with no header
            for label in order:
                f.write(f"{label}\n")
    
for dataset in data_ids:
    # for loop that loops thru data ids and reads in data and associated args with that data for ingest
        # yeah I think doing this with a sample sheet will be easier and we can probably bypass needing argparse
    pass 