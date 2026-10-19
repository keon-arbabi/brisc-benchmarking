import os
import gc
import sys
import h5py
import logging
import numpy as np
import polars as pl
import anndata as ad
import scanpy as sc
from pathlib import Path
sys.path.append(f'{Path.home()}/brisc-benchmarking')
from utils_local import MemoryTimer, system_info, remove_self_neighbors

os.environ.setdefault('CUPY_CACHE_DIR', '/tmp/cupy_cache')
os.environ.setdefault('CUDA_PATH', os.path.join(
    os.path.dirname(os.__file__), 'site-packages',
    'nvidia', 'cuda_runtime'))

from dask_cuda import LocalCUDACluster
from dask.distributed import Client

DATA_NAME = sys.argv[1]
DATA_PATH = sys.argv[2]
OUTPUT_PATH_TIME = sys.argv[3]
OUTPUT_PATH_EMBEDDING = sys.argv[4]
OUTPUT_PATH_PCS = sys.argv[5]
OUTPUT_PATH_NEIGHBORS = sys.argv[6]
OUTPUT_PATH_HVGS = sys.argv[7]
OUTPUT_PATH_MARKERS = sys.argv[8]
REP = int(sys.argv[-1])

CHUNK_SIZE = 20_000

if __name__ == '__main__':

    system_info()
    print('--- Params ---')
    print('rapids basic (managed)')
    print(f'{DATA_PATH=}')

    cluster = LocalCUDACluster(
        CUDA_VISIBLE_DEVICES='0,1,2,3',
        threads_per_worker=1,
        protocol='tcp',
        rmm_managed_memory=True,
        rmm_allocator_external_lib_list='cupy',
        silence_logs=logging.ERROR,
    )
    client = Client(cluster)
    print(client)

    import rapids_singlecell as rsc
    from anndata.experimental import read_elem_lazy as read_dask

    timers = MemoryTimer(
        silent=False, csv_path=OUTPUT_PATH_TIME,
        csv_columns={'library': 'rapids', 'test': 'basic',
                     'dataset': DATA_NAME, 'rep': REP})

    with timers('Load data'):
        f = h5py.File(DATA_PATH, 'r')
        shape = tuple(f['X'].attrs['shape'])
        data = ad.AnnData(
            X=read_dask(f['X'], (CHUNK_SIZE, shape[1])),
            obs=ad.io.read_elem(f['obs']),
            var=ad.io.read_elem(f['var']),
        )
        rsc.get.anndata_to_GPU(data)
        data.X = data.X.persist()
        data.X.compute_chunk_sizes()

    with timers('Quality control'):
        data.var['mt'] = data.var_names.str.upper().str.startswith('MT-')
        data.var['malat1'] = data.var_names.str.upper() == 'MALAT1'
        rsc.pp.calculate_qc_metrics(
            data, qc_vars=['mt', 'malat1'], log1p=False)
        keep = ((data.obs['n_genes_by_counts'].values >= 100) &
                (data.obs['pct_counts_mt'].values <= 5) &
                (data.obs['pct_counts_malat1'].values > 0))
        data = data[keep].copy()
        data.X = data.X.persist()
        data.X.compute_chunk_sizes()

    gc.collect()

    with timers('Normalization'):
        rsc.pp.normalize_total(data, target_sum=1e4)
        rsc.pp.log1p(data)
        data.X = data.X.persist()
        data.X.compute_chunk_sizes()

    with timers('Feature selection'):
        rsc.pp.highly_variable_genes(
            data, n_top_genes=2000, flavor="cell_ranger")
        rsc.pp.filter_highly_variable(data)

    with timers('PCA'):
        n_rows = data.shape[0]
        n_cols = data.shape[1]
        n_workers = len(client.scheduler_info()['workers'])
        rows_per_worker = (n_rows + n_workers - 1) // n_workers
        data.X = data.X.rechunk((rows_per_worker, n_cols)).persist()
        data.X.compute_chunk_sizes()

        rsc.pp.scale(data, zero_center=False)
        data.X = data.X.persist()
        data.X.compute_chunk_sizes()

        rsc.pp.pca(data, n_comps=50, mask_var=None)
        data.obsm['X_pca'] = data.obsm['X_pca'].persist()
        data.obsm['X_pca'].compute_chunk_sizes()
        data.obsm['X_pca'] = data.obsm['X_pca'].compute()

    with timers('Nearest neighbors'):
        rsc.pp.neighbors(
            data, n_neighbors=15, n_pcs=50, algorithm='mg_ivfflat')

    with timers('Clustering'):
        rsc.tl.leiden(data, resolution=[0.25, 0.5, 1.0, 1.5, 2.0])

    with timers('Embedding'):
        rsc.tl.umap(data)

    with timers('Find markers'):
        rsc.tl.rank_genes_groups(
            data, groupby='cell_type', method='logreg', use_raw=False)

    # Save markers
    pl.from_pandas(sc.get.rank_genes_groups_df(data, group=None))\
        .rename({'group': 'cell_type', 'names': 'gene', 'scores': 'score'})\
        .with_columns(rank=pl.int_range(1, pl.len() + 1).over('cell_type'))\
        .select('cell_type', 'gene', 'rank', 'score')\
        .write_parquet(OUTPUT_PATH_MARKERS)

    # Save HVGs
    pl.DataFrame({'gene': data.var_names.tolist()}).write_csv(OUTPUT_PATH_HVGS)

    pcs = data.obsm['X_pca'].get()
    pc_df = pl.DataFrame({
        f'PC_{i+1}': pcs[:, i] for i in range(pcs.shape[1])
    })
    pc_df.write_parquet(OUTPUT_PATH_PCS)

    dist = data.obsp['distances']
    neighbors = remove_self_neighbors(
        dist.indices.reshape(data.n_obs, -1)).astype(np.uint32)
    neighbors_df = pl.DataFrame({
        f'neighbor_{i+1}': neighbors[:, i]
        for i in range(neighbors.shape[1])
    })
    neighbors_df.write_parquet(OUTPUT_PATH_NEIGHBORS)

    ct_col = 'cell_type_fine' if DATA_NAME == 'PanSci' else 'cell_type'
    embedding_df = pl.DataFrame({
        'cell_id': data.obs_names.tolist(),
        'embed_1': data.obsm['X_umap'][:, 0],
        'embed_2': data.obsm['X_umap'][:, 1],
        'cell_type': data.obs[ct_col].tolist(),
        'cell_type_broad': data.obs['cell_type_broad'].tolist(),
        'cluster_res_0.25': data.obs['leiden_0.25'].tolist(),
        'cluster_res_0.5': data.obs['leiden_0.5'].tolist(),
        'cluster_res_1.0': data.obs['leiden_1.0'].tolist(),
        'cluster_res_1.5': data.obs['leiden_1.5'].tolist(),
        'cluster_res_2.0': data.obs['leiden_2.0'].tolist(),
    })
    embedding_df.write_csv(OUTPUT_PATH_EMBEDDING)

    timers.shutdown()
    print('--- Completed successfully ---')

    print('\n--- Session Info ---')
    sc.logging.print_header()
    sys.stdout.flush()
    os._exit(0)
