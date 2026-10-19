import sys
import numpy as np
import polars as pl
import scanpy as sc
from pathlib import Path
sys.path.append(f'{Path.home()}/brisc-benchmarking')
from utils_local import MemoryTimer, system_info, remove_self_neighbors

import warnings
# For sc.tl.rank_genes_groups
warnings.filterwarnings('ignore', message='DataFrame is highly fragmented')

DATA_NAME = sys.argv[1]
DATA_PATH = sys.argv[2]
OUTPUT_PATH_TIME = sys.argv[3]
OUTPUT_PATH_EMBEDDING = sys.argv[4]
OUTPUT_PATH_PCS = sys.argv[5]
OUTPUT_PATH_NEIGHBORS = sys.argv[6]
OUTPUT_PATH_HVGS = sys.argv[7]
OUTPUT_PATH_MARKERS = sys.argv[8]
REP = int(sys.argv[-1])

if __name__ == '__main__':

    system_info()
    print('--- Params ---')
    print('scanpy basic')
    print(f'{DATA_PATH=}')

    timers = MemoryTimer(
        silent=False, csv_path=OUTPUT_PATH_TIME,
        csv_columns={'library': 'scanpy', 'test': 'basic',
                     'dataset': DATA_NAME, 'rep': REP})

    with timers('Load data'):
        data = sc.read_h5ad(DATA_PATH)

    with timers('Quality control'):
        data.var['mt'] = data.var_names.str.upper().str.startswith('MT-')
        data.var['malat1'] = data.var_names.str.upper() == 'MALAT1'
        sc.pp.calculate_qc_metrics(
            data, qc_vars=['mt', 'malat1'],
            percent_top=None, log1p=False, inplace=True)
        keep = ((data.obs['n_genes_by_counts'].values >= 100) &
                (data.obs['pct_counts_mt'].values <= 5) &
                (data.obs['pct_counts_malat1'].values > 0))
        data = data[keep].copy()

    with timers('Normalization'):
        sc.pp.normalize_total(data)
        sc.pp.log1p(data)

    with timers('Feature selection'):
        sc.pp.highly_variable_genes(
            data, n_top_genes=2000, batch_key='donor')

    with timers('PCA'):
        sc.tl.pca(data)

    with timers('Nearest neighbors'):
        sc.pp.neighbors(data)

    with timers('Clustering'):
        for res in [0.25, 0.5, 1.0, 1.5, 2.0]:
            sc.tl.leiden(
                data,
                resolution=res,
                flavor='igraph', # Future default
                n_iterations=2, # Future default
                key_added=f'leiden_res_{res:4.2f}')

    with timers('Embedding'):
        sc.tl.umap(data)

    with timers('Find markers'):
        sc.tl.rank_genes_groups(data, groupby='cell_type')

    # Save markers
    pl.from_pandas(sc.get.rank_genes_groups_df(data, group=None))\
        .rename({'group': 'cell_type', 'names': 'gene', 'scores': 'score',
                 'logfoldchanges': 'logFC', 'pvals': 'p_value',
                 'pvals_adj': 'p_value_adj'})\
        .with_columns(rank=pl.int_range(1, pl.len() + 1).over('cell_type'))\
        .select('cell_type', 'gene', 'rank', 'score', 'logFC', 'p_value',
                'p_value_adj')\
        .write_parquet(OUTPUT_PATH_MARKERS)

    # Save HVGs
    pl.DataFrame({
        'gene': data.var_names[data.var['highly_variable'].values].tolist()
    }).write_csv(OUTPUT_PATH_HVGS)

    # Save PCs
    pcs = data.obsm['X_pca']
    pc_df = pl.DataFrame({
        f'PC_{i+1}': pcs[:, i] for i in range(pcs.shape[1])
    })
    pc_df.write_parquet(OUTPUT_PATH_PCS)

    # Save neighbors
    dist = data.obsp['distances']
    neighbors = remove_self_neighbors(
        dist.indices.reshape(data.n_obs, -1)).astype(np.uint32)
    neighbors_df = pl.DataFrame({
        f'neighbor_{i+1}': neighbors[:, i]
        for i in range(neighbors.shape[1])
    })
    neighbors_df.write_parquet(OUTPUT_PATH_NEIGHBORS)

    # Save embeddings
    ct_col = 'cell_type_fine' if DATA_NAME == 'PanSci' else 'cell_type'
    embedding_df = pl.DataFrame({
        'cell_id': data.obs_names.tolist(),
        'embed_1': data.obsm['X_umap'][:, 0],
        'embed_2': data.obsm['X_umap'][:, 1],
        'cell_type': data.obs[ct_col].tolist(),
        'cell_type_broad': data.obs['cell_type_broad'].tolist(),
        'cluster_res_0.25': data.obs['leiden_res_0.25'].tolist(),
        'cluster_res_0.5': data.obs['leiden_res_0.50'].tolist(),
        'cluster_res_1.0': data.obs['leiden_res_1.00'].tolist(),
        'cluster_res_1.5': data.obs['leiden_res_1.50'].tolist(),
        'cluster_res_2.0': data.obs['leiden_res_2.00'].tolist(),
    })
    embedding_df.write_csv(OUTPUT_PATH_EMBEDDING)

    timers.shutdown()
    print('--- Completed successfully ---')

    print('\n--- Session Info ---')
    sc.logging.print_header()
