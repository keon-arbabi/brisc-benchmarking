import os
import sys
import polars as pl
from pathlib import Path
sys.path.append(f'{Path.home()}/brisc-benchmarking')
from brisc import SingleCell
from utils_local import MemoryTimer, system_info, remove_self_neighbors

GPU = os.uname().nodename.startswith('trig')
DATA_NAME = sys.argv[1]
DATA_PATH = sys.argv[2]
NUM_THREADS = int(sys.argv[3])
OUTPUT_PATH_TIME = sys.argv[4]
OUTPUT_PATH_EMBEDDING = sys.argv[5]
OUTPUT_PATH_PCS = sys.argv[6]
OUTPUT_PATH_NEIGHBORS = sys.argv[7]
OUTPUT_PATH_HVGS = sys.argv[8]
OUTPUT_PATH_MARKERS = sys.argv[9]
REP = int(sys.argv[-1])

if __name__ == '__main__':

    system_info()
    print('--- Params ---')
    print('brisc basic')
    print(f'{DATA_PATH=}')
    print(f'{NUM_THREADS=}')

    timers = MemoryTimer(
        silent=False, csv_path=OUTPUT_PATH_TIME,
        csv_columns={'library': 'brisc', 'test': 'basic',
                     'dataset': DATA_NAME, 'num_threads': NUM_THREADS,
                     'rep': REP})

    with timers('Load data'):
        data = SingleCell(DATA_PATH, num_threads=NUM_THREADS)

    with timers('Quality control'):
        data = data.qc(allow_float=True)

    # Rapids has degraded performance with hvg batch, match no batch
    with timers('Feature selection'):
        data = data.hvg() if GPU else data.hvg(batch_column='donor')

    with timers('Normalization'):
        data = data.normalize()

    with timers('PCA'):
        data = data.pca()

    with timers('Nearest neighbors'):
        data = data.neighbors().shared_neighbors()

    with timers('Clustering'):
        data = data.cluster(resolution=[0.25, 0.5, 1.0, 1.5, 2.0])

    with timers('Embedding (PaCMAP)'):
        data = data.pacmap()

    if not GPU:
        with timers('Embedding (LocalMAP)', exclude=True):
            data = data.localmap()

    if not GPU and NUM_THREADS != 1:
        with timers('Embedding (UMAP)', exclude=True):
            data = data.umap()

        with timers('Embedding (UMAP hogwild)', exclude=True):
            data = data.umap(
                hogwild=True, embedding_key='umap_hogwild')

    with timers('Find markers'):
        markers = data.find_markers('cell_type', all_genes=True)

    data = data.drop_X()
    data = data.filter_obs('passed_QC')

    # Save markers
    markers\
        .sort('cell_type', 'fold_change', descending=[False, True])\
        .with_columns(rank=pl.int_range(1, pl.len() + 1).over('cell_type'))\
        .select('cell_type', 'gene', 'rank', 'marker', 'detection_rate',
                'fold_change')\
        .write_parquet(OUTPUT_PATH_MARKERS)

    # Save HVGs
    pl.DataFrame({
        'gene': data.var_names.filter(data.var['highly_variable'])
    }).write_csv(OUTPUT_PATH_HVGS)

    # Save PCs
    pcs = data.obsm['pca']
    pc_df = pl.DataFrame({
        f'PC_{i+1}': pcs[:, i] for i in range(pcs.shape[1])
    })
    pc_df.write_parquet(OUTPUT_PATH_PCS)

    # Save neighbors
    neighbors = remove_self_neighbors(data.obsm['neighbors'])
    neighbors_df = pl.DataFrame({
        f'neighbor_{i+1}': neighbors[:, i]
        for i in range(neighbors.shape[1])
    })
    neighbors_df.write_parquet(OUTPUT_PATH_NEIGHBORS)

    # Save embeddings
    ct_col = 'cell_type_fine' if DATA_NAME == 'PanSci' else 'cell_type'
    embed_cols = {
        'cell_id': data.obs['cell_id'],
        'cell_type': data.obs[ct_col],
        'cell_type_broad': data.obs['cell_type_broad'],
        'cluster_res_0.25': data.obs['cluster_0'],
        'cluster_res_0.5': data.obs['cluster_1'],
        'cluster_res_1.0': data.obs['cluster_2'],
        'cluster_res_1.5': data.obs['cluster_3'],
        'cluster_res_2.0': data.obs['cluster_4'],
    }
    for key in ('pacmap', 'localmap', 'umap', 'umap_hogwild'):
        if key in data.obsm:
            embed_cols[f'{key}_1'] = data.obsm[key][:, 0]
            embed_cols[f'{key}_2'] = data.obsm[key][:, 1]
    pl.DataFrame(embed_cols).write_csv(OUTPUT_PATH_EMBEDDING)

    timers.shutdown()
    print('--- Completed successfully ---')
