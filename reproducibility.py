import os
import numpy as np
import polars as pl
from pathlib import Path
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree

OUTPUT_DIR = f'{Path.home()}/brisc-benchmarking/output'
SUMMARY_PATH = f'{OUTPUT_DIR}/reproducibility_summary.csv'
DATASETS = ['SEAAD', 'Parse', 'PanSci']
REPS = [1, 2, 3]
PAIRS = [(1, 2), (1, 3), (2, 3)]
RESOLUTION = '1.0'
EMBED_KNN_K = 15  # neighbourhood size for embedding-space overlap

CONFIGS = {
    'brisc_mt': ('basic_brisc_{ds}_-1_', ('localmap_1', 'localmap_2')),
    'brisc_st': ('basic_brisc_{ds}_1_', ('localmap_1', 'localmap_2')),
    'scanpy': ('basic_scanpy_{ds}_', ('embed_1', 'embed_2')),
    'seurat': ('basic_seurat_{ds}_', ('embed_1', 'embed_2')),
    'rapids': ('basic_rapids_{ds}_gpu_', ('embed_1', 'embed_2')),
}

def load(prefix, name, reader):
    paths = {r: f'{OUTPUT_DIR}/rep{r}/{prefix}{name}' for r in REPS}
    missing = [p for p in paths.values() if not os.path.exists(p)]
    if missing:
        print(f'  missing {missing[0]}; skipping', flush=True)
        return None
    return {r: reader(p) for r, p in paths.items()}

def encode(arr):
    return np.unique(arr, return_inverse=True)[1].astype(np.int32)

def adjusted_rand(a, b):
    na, nb = int(a.max()) + 1, int(b.max()) + 1
    n = a.shape[0]
    c = np.bincount(a.astype(np.int64) * nb + b.astype(np.int64),
                    minlength=na * nb).reshape(na, nb).astype(np.float64)
    sum_c = (c * (c - 1)).sum() / 2
    sum_a = (c.sum(1) * (c.sum(1) - 1)).sum() / 2
    sum_b = (c.sum(0) * (c.sum(0) - 1)).sum() / 2
    expected = sum_a * sum_b / (n * (n - 1) / 2)
    max_idx = (sum_a + sum_b) / 2
    return float((sum_c - expected) / (max_idx - expected)) \
        if max_idx != expected else 0.0

def cluster_flip_rate(a, b):
    # fraction of cells whose cluster changes after one-to-one matching of
    # clusters across runs (cells of unmatched clusters count as changed)
    na, nb = int(a.max()) + 1, int(b.max()) + 1
    c = np.bincount(a.astype(np.int64) * nb + b.astype(np.int64),
                    minlength=na * nb).reshape(na, nb)
    rows, cols = linear_sum_assignment(c, maximize=True)
    return float(1 - c[rows, cols].sum() / a.shape[0])

def max_pc_rel_rmse(A, B):
    # worst PC's RMS difference relative to that PC's spread, sign-aligned
    A, B = A.astype(np.float64), B.astype(np.float64)
    sign = np.sign((A * B).sum(axis=0))
    sign[sign == 0] = 1
    D = A - B * sign
    return float((np.sqrt((D ** 2).mean(axis=0)) / A.std(axis=0)).max())

def neighbor_overlap(A, B):
    # mean per-cell fraction of shared neighbours
    A = np.sort(A, axis=1)
    B = np.sort(B, axis=1)
    k = B.shape[1]
    new_value = np.empty_like(B, dtype=bool)
    new_value[:, 0] = True
    new_value[:, 1:] = B[:, 1:] != B[:, :-1]
    hits = np.zeros(A.shape[0], dtype=np.int32)
    for j in range(k):
        hits += np.any(A == B[:, j:j + 1], axis=1) & new_value[:, j]
    return float((hits / k).mean())

def embedding_knn(coords):
    # kNN within the 2D embedding itself (self excluded); overlap of on-map
    # neighbourhoods is invariant to the arbitrary global cluster arrangement
    return cKDTree(coords).query(coords, k=EMBED_KNN_K + 1,
                                 workers=-1)[1][:, 1:]

def score(config, ds):
    print(f'[{ds}/{config}]', flush=True)
    prefix, embed_cols = CONFIGS[config]
    prefix = prefix.format(ds=ds)
    cluster_col = f'cluster_res_{RESOLUTION}'
    emb = load(prefix, 'embedding.csv',
               lambda p: pl.read_csv(p, columns=['cell_id', cluster_col,
                                                 *embed_cols]))
    if emb is None:
        return []
    ids = {r: emb[r]['cell_id'] for r in REPS}
    clusters = {r: emb[r][cluster_col] for r in REPS}
    coords = {r: emb[r].select(embed_cols).to_numpy() for r in REPS}
    del emb

    rows = []
    def add(step, metric, pair, value, identical=None):
        rows.append(dict(config=config, dataset=ds, step=step, metric=metric,
                         pair=pair, value=value, identical=identical))

    for r in REPS:
        add('Clustering', 'n_clusters', f'rep{r}',
            float(clusters[r].n_unique()))

    knn = {}
    def emb_knn(r):
        if r not in knn:  # only ever computed when coordinates differ
            knn[r] = embedding_knn(coords[r])
        return knn[r]

    for a, b in PAIRS:
        pair = f'{a}v{b}'
        same = ids[a].equals(ids[b])
        add('Cells', 'cell_id_jaccard', pair, float(same), same)

        # identical runs short-circuit to the exact perfect value
        same = clusters[a].equals(clusters[b])
        if same:
            ari, flips = 1.0, 0.0
        else:
            ca = encode(clusters[a].to_numpy())
            cb = encode(clusters[b].to_numpy())
            ari, flips = adjusted_rand(ca, cb), cluster_flip_rate(ca, cb)
        add('Clustering', f'ari_res_{RESOLUTION}', pair, ari, same)
        add('Clustering', 'cluster_flip_rate', pair, flips, same)

        same = np.array_equal(coords[a], coords[b])
        add('Embedding', 'embedding_neighbor_overlap', pair,
            1.0 if same else neighbor_overlap(emb_knn(a), emb_knn(b)), same)
    del coords, knn

    hvgs = load(prefix, 'hvgs.csv', lambda p: set(pl.read_csv(p)['gene']))
    if hvgs is not None:
        for a, b in PAIRS:
            same = hvgs[a] == hvgs[b]
            add('Feature selection', 'hvg_jaccard', f'{a}v{b}',
                len(hvgs[a] & hvgs[b]) / len(hvgs[a] | hvgs[b]), same)
        del hvgs

    pcs = load(prefix, 'pcs.parquet', lambda p: pl.read_parquet(p).to_numpy())
    if pcs is not None:
        for a, b in PAIRS:
            same = np.array_equal(pcs[a], pcs[b])
            add('PCA', 'max_pc_rel_rmse', f'{a}v{b}',
                0.0 if same else max_pc_rel_rmse(pcs[a], pcs[b]), same)
        del pcs

    nn = load(prefix, 'neighbors.parquet',
              lambda p: pl.read_parquet(p).to_numpy())
    if nn is not None:
        for a, b in PAIRS:
            same = np.array_equal(nn[a], nn[b])
            add('Nearest neighbors', 'mean_neighbor_overlap', f'{a}v{b}',
                1.0 if same else neighbor_overlap(nn[a], nn[b]), same)
        del nn
    return rows

rows = []
for ds in DATASETS:
    for config in CONFIGS:
        rows += score(config, ds)
pl.DataFrame(rows)\
    .sort('dataset', 'config', 'step', 'metric', 'pair')\
    .write_csv(SUMMARY_PATH)
print(f'Saved {SUMMARY_PATH}', flush=True)
