import os
import time
import faiss
import numpy as np
import polars as pl
from pathlib import Path

OUTPUT_DIR = f'{Path.home()}/brisc-benchmarking/output'
DATASETS = ['SEAAD', 'Parse', 'PanSci']
LIBRARIES = {
    'brisc': {
        'pcs': 'basic_brisc_{dataset}_-1_pcs.parquet',
        'neighbors': 'basic_brisc_{dataset}_-1_neighbors.parquet',
    },
    'scanpy': {
        'pcs': 'basic_scanpy_{dataset}_pcs.parquet',
        'neighbors': 'basic_scanpy_{dataset}_neighbors.parquet',
    },
    'seurat': {
        'pcs': 'basic_seurat_{dataset}_pcs.parquet',
        'neighbors': 'basic_seurat_{dataset}_neighbors.parquet',
    },
    'rapids': {
        'pcs': 'basic_rapids_{dataset}_gpu_pcs.parquet',
        'neighbors': 'basic_rapids_{dataset}_gpu_neighbors.parquet',
    },
}

REP = 1  # score only the first replicate's outputs
K = 14   # test the top-14 nearest (min neighbors across libraries: rapids=14)
BATCH_SIZE = 500_000

TAIL_FRAC = 0.25

def log(msg):
    print(msg, flush=True)

def tail_mean(rc, frac=TAIL_FRAC):
    n = max(1, int(frac * rc.size))
    return float(np.sort(rc)[:n].mean())

def exact_knn(PCs, k=20, batch_size=100_000):
    # Exact k nearest neighbors (GPU brute force), self excluded.
    n, d = PCs.shape
    assert PCs.dtype == np.float32 and PCs.flags['C_CONTIGUOUS']
    if hasattr(faiss, 'omp_set_num_threads'):
        faiss.omp_set_num_threads(os.cpu_count())
    res = faiss.StandardGpuResources()
    index = faiss.index_cpu_to_gpu(res, 0, faiss.IndexFlatL2(d))
    index.add(PCs)

    all_I = np.empty((n, k + 1), dtype=np.int64)
    t0 = time.time()
    for s in range(0, n, batch_size):
        e = min(s + batch_size, n)
        _, all_I[s:e] = index.search(PCs[s:e], k + 1)
        log(f'  {e:,}/{n:,} ({time.time() - t0:.0f}s)')

    # Drop self from each row; flag rows where self was not found.
    self_mask = all_I == np.arange(n)[:, None]
    valid = self_mask.any(axis=1)
    n_missing = int((~valid).sum())
    if n_missing:
        log(f'  WARN: {n_missing:,}/{n:,} queries missing self in top-{k+1} '
            f'(excluded from recall)')
    keep = np.ones((n, k + 1), dtype=bool)
    keep[np.arange(n), self_mask.argmax(axis=1)] = False
    return all_I[keep].reshape(n, k).astype(np.uint32), valid

def recall_at_k(gt, approx):
    # Fraction of exact neighbors recovered per cell (approx ties deduplicated).
    n, k = gt.shape
    assert approx.shape == (n, k)
    gt_s = np.sort(gt, axis=1)
    ap_s = np.sort(approx, axis=1)
    is_first = np.empty_like(ap_s, dtype=bool)
    is_first[:, 0] = True
    is_first[:, 1:] = ap_s[:, 1:] != ap_s[:, :-1]
    n_dup_rows = int((~is_first).any(axis=1).sum())
    if n_dup_rows:
        log(f'  WARN: {n_dup_rows:,}/{n:,} rows have duplicate approx '
            f'indices (deduplicated for recall)')
    hits = np.zeros(n, dtype=np.int32)
    for j in range(k):
        hits += np.any(gt_s == ap_s[:, j:j + 1], axis=1) & is_first[:, j]
    return hits.astype(np.float32) / k

all_summary_rows = []

for dataset_name in DATASETS:
    dataset_cells = {}
    for lib_name, paths in LIBRARIES.items():
        tag = f'[{dataset_name}/{lib_name}]'
        rep_dir = os.path.join(OUTPUT_DIR, f'rep{REP}')
        pc_path = os.path.join(
            rep_dir, paths['pcs'].format(dataset=dataset_name))
        nn_path = os.path.join(
            rep_dir, paths['neighbors'].format(dataset=dataset_name))

        if not os.path.exists(pc_path):
            log(f'{tag} PCs not found: {pc_path}')
            continue
        if not os.path.exists(nn_path):
            log(f'{tag} Neighbors not found: {nn_path}')
            continue

        log(f'{tag} Loading approximate neighbors...')
        approx_neighbors = pl.read_parquet(nn_path).to_numpy().astype(np.uint32)
        if approx_neighbors.shape[1] < K:
            log(f'{tag} only {approx_neighbors.shape[1]} neighbors (< K={K}); '
                f'skipping')
            continue
        # neighbors are distance-sorted, so the first K are the nearest K
        approx_neighbors = approx_neighbors[:, :K]

        gt_cache = os.path.join(
            OUTPUT_DIR, f'knn_exact_{lib_name}_{dataset_name}.npz')
        if os.path.exists(gt_cache) and \
                os.path.getmtime(gt_cache) > os.path.getmtime(pc_path):
            log(f'{tag} Loading cached exact KNN...')
            cached = np.load(gt_cache)
            gt_neighbors, valid = cached['gt_neighbors'], cached['valid']
        else:
            log(f'{tag} Loading PCs...')
            PCs = np.ascontiguousarray(
                pl.read_parquet(pc_path).to_numpy(), dtype=np.float32)
            log(f'{tag} Computing exact KNN '
                f'({PCs.shape[0]:,} cells x {PCs.shape[1]} PCs)...')
            gt_neighbors, valid = exact_knn(PCs, k=K, batch_size=BATCH_SIZE)
            np.savez(gt_cache, gt_neighbors=gt_neighbors, valid=valid)

        gt_k = gt_neighbors[:, :K]
        n_invalid = int((~valid).sum())
        if n_invalid:
            log(f'{tag} Excluding {n_invalid:,}/{valid.shape[0]:,} '
                f'rows with missing self from recall')
            gt_k = gt_k[valid]
            approx_neighbors = approx_neighbors[valid]

        log(f'{tag} Computing recall@{K}...')
        rc = recall_at_k(gt_k, approx_neighbors)

        # full-length per-cell recall (NaN for cells excluded above), so all
        # libraries share one positionally-aligned per-dataset table
        rc_full = np.full(valid.shape[0], np.nan, dtype=np.float32)
        rc_full[valid] = rc
        dataset_cells[lib_name] = rc_full

        rc_path = os.path.join(
            OUTPUT_DIR, f'knn_recall_{lib_name}_{dataset_name}.csv')
        pl.DataFrame({'recall': rc}).write_csv(rc_path)
        log(f'{tag} Saved per-cell recall to {rc_path}')

        all_summary_rows.append(dict(
            method=lib_name,
            dataset=dataset_name,
            k=K,
            mean=float(rc.mean()),
            min=float(rc.min()),
            p5=float(np.quantile(rc, 0.05)),
            p25=float(np.quantile(rc, 0.25)),
            cvar25=tail_mean(rc)))

        log(f'{tag} Recall@{K}: mean={rc.mean():.4f}, min={rc.min():.4f}')

    # Cache every cell's recall for this dataset, one column per library
    if dataset_cells:
        lens = {k: len(v) for k, v in dataset_cells.items()}
        if len(set(lens.values())) == 1:
            cells_path = os.path.join(
                OUTPUT_DIR, f'knn_recall_cells_{dataset_name}.parquet')
            pl.DataFrame(dataset_cells).write_parquet(cells_path)
            log(f'[{dataset_name}] Saved per-cell recall for '
                f'{len(dataset_cells)} libraries to {cells_path}')
        else:
            log(f'[{dataset_name}] cell counts differ across libraries '
                f'{lens}; skipping consolidated per-cell recall file')

if all_summary_rows:
    pl.DataFrame(all_summary_rows)\
        .write_csv(f'{OUTPUT_DIR}/knn_recall_summary.csv')
    log(f'Saved summary to {OUTPUT_DIR}/knn_recall_summary.csv')
