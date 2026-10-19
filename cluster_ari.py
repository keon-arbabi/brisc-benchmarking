import os
import numpy as np
import polars as pl
from pathlib import Path

OUTPUT_DIR = f'{Path.home()}/brisc-benchmarking/output'
DATASETS = ['SEAAD', 'Parse', 'PanSci']
LIBRARIES = {
    'brisc': 'basic_brisc_{dataset}_-1_embedding.csv',
    'scanpy': 'basic_scanpy_{dataset}_embedding.csv',
    'seurat': 'basic_seurat_{dataset}_embedding.csv',
    'rapids': 'basic_rapids_{dataset}_gpu_embedding.csv',
}
LABEL_COLS = ['cell_type', 'cell_type_broad']
RESOLUTIONS = ['0.25', '0.5', '1.0', '1.5', '2.0']
REP = 1
MIN_CLUSTER_SIZE = 100

def log(msg):
    print(msg, flush=True)

def encode(arr):
    _, inv = np.unique(arr, return_inverse=True)
    return inv.astype(np.int32)

def fast_ari(a, b):
    na, nb = int(a.max()) + 1, int(b.max()) + 1
    n = a.shape[0]
    flat = a.astype(np.int64) * nb + b.astype(np.int64)
    c = np.bincount(flat, minlength=na * nb)\
        .reshape(na, nb).astype(np.float64)
    ai, bj = c.sum(1), c.sum(0)
    sum_c = (c * (c - 1)).sum() / 2
    sum_a = (ai * (ai - 1)).sum() / 2
    sum_b = (bj * (bj - 1)).sum() / 2
    n_pairs = n * (n - 1) / 2
    expected = sum_a * sum_b / n_pairs
    max_idx = (sum_a + sum_b) / 2
    return float((sum_c - expected) / (max_idx - expected)) \
        if max_idx != expected else 0.0

def balanced_ari(it, ip):
    # Balanced ARI of Maan et al. 2024 (Nat Biotechnol, doi 10/gtb3zf;
    # github.com/hsmaan/balanced-clustering): each row of the contingency
    # table is rescaled so every class contributes equally while the total
    # count is preserved, then the standard ARI is computed from it.
    n = it.shape[0]
    na, nb = int(it.max()) + 1, int(ip.max()) + 1
    c = np.bincount(it.astype(np.int64) * nb + ip.astype(np.int64),
                    minlength=na * nb).reshape(na, nb).astype(np.float64)
    c = (c * (round(c.sum(1).mean()) / c.sum(1))[:, None]).astype(np.int64)
    sq = int((c ** 2).sum())
    tp = sq - n
    fp = int((c * c.sum(0)).sum()) - sq
    fn = int((c * c.sum(1)[:, None]).sum()) - sq
    tn = n ** 2 - fp - fn - sq
    return float(2 * (tp * tn - fn * fp) /
                 ((tp + fn) * (fn + tn) + (tp + fp) * (fp + tn)))

all_summary_rows = []
rep_dir = os.path.join(OUTPUT_DIR, f'rep{REP}')

for dataset_name in DATASETS:
    for lib_name, fname in LIBRARIES.items():
        path = os.path.join(
            rep_dir, fname.format(dataset=dataset_name))
        if not os.path.exists(path):
            log(f'[{dataset_name}/{lib_name}] Embedding not found: {path}')
            continue

        cols = LABEL_COLS + [f'cluster_res_{r}' for r in RESOLUTIONS]
        log(f'[{dataset_name}/{lib_name}] Loading embedding...')
        df = pl.read_csv(path, columns=cols)

        for label_col in LABEL_COLS:
            it = encode(df[label_col].to_numpy())
            n_classes = int(it.max()) + 1
            for res in RESOLUTIONS:
                ip = encode(df[f'cluster_res_{res}'].to_numpy())
                n_clusters = int(ip.max()) + 1
                cluster_sizes = np.bincount(ip, minlength=n_clusters)
                n_singletons = int((cluster_sizes == 1).sum())
                n_clusters_big = int(
                    (cluster_sizes >= MIN_CLUSTER_SIZE).sum())
                ari = fast_ari(it, ip)
                bari = balanced_ari(it, ip)

                all_summary_rows.append({
                    'method': lib_name,
                    'dataset': dataset_name,
                    'label': label_col,
                    'resolution': float(res),
                    'n_classes': n_classes,
                    'n_clusters': n_clusters,
                    'n_singletons': n_singletons,
                    f'n_clusters_ge_{MIN_CLUSTER_SIZE}': n_clusters_big,
                    'ari': ari,
                    'balanced_ari': bari})

                log(f'[{dataset_name}/{lib_name}] {label_col} '
                    f'res={res}: ARI={ari:.4f}, '
                    f'balanced ARI={bari:.4f} '
                    f'(n_classes={n_classes}, n_clusters={n_clusters}, '
                    f'n_singletons={n_singletons}, '
                    f'n_clusters>={MIN_CLUSTER_SIZE}: {n_clusters_big})')

if all_summary_rows:
    out_path = f'{OUTPUT_DIR}/cluster_ari_summary.csv'
    pl.DataFrame(all_summary_rows).write_csv(out_path)
    log(f'Saved summary to {out_path}')
