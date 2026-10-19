import gc
import sys
import scanpy as sc
import polars as pl
from pathlib import Path
sys.path.append(f'{Path.home()}/brisc-benchmarking')
from utils_local import MemoryTimer, system_info, transfer_accuracy

DATA_NAME = sys.argv[1]
DATA_PATH = sys.argv[2]
OUTPUT_PATH_TIME = sys.argv[3]
OUTPUT_PATH_ACC = sys.argv[4]
OUTPUT_PATH_LABELS = sys.argv[5]
REP = int(sys.argv[-1])

if __name__ == '__main__':

    system_info()
    print('--- Params ---')
    print('scanpy transfer')
    print(f'{DATA_PATH=}')

    timers = MemoryTimer(
        silent=False, csv_path=OUTPUT_PATH_TIME,
        csv_columns={'library': 'scanpy', 'test': 'transfer',
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

    with timers('Split data'):
        data_ref = data[data.obs['is_ref'] == 1].copy()
        data_query = data[data.obs['is_ref'] == 0].copy()

    del data; gc.collect()

    with timers('Normalization'):
        sc.pp.normalize_total(data_ref)
        sc.pp.log1p(data_ref)
        sc.pp.normalize_total(data_query)
        sc.pp.log1p(data_query)

    with timers('Feature selection'):
        sc.pp.highly_variable_genes(
            data_ref, n_top_genes=2000, batch_key='donor')
        hvg_genes = data_ref.var_names[data_ref.var['highly_variable']]
        data_ref = data_ref[:, hvg_genes].copy()
        data_query = data_query[:, hvg_genes].copy()

    with timers('PCA'):
        sc.pp.pca(data_ref)

    cell_type_col = 'cell_type_fine' if DATA_NAME == 'PanSci' else 'cell_type'
    data_query.obs['cell_type_orig'] = data_query.obs[cell_type_col]

    with timers('Transfer labels'):
        sc.pp.neighbors(data_ref)
        sc.tl.ingest(
            data_query,
            data_ref,
            obs=cell_type_col,
            embedding_method='pca')

    print('--- Transfer Accuracy ---')
    transfer_accuracy(
        data_query.obs, 'cell_type_orig', cell_type_col)\
    .write_csv(OUTPUT_PATH_ACC)

    pl.DataFrame({
        'cell_id': data_query.obs_names.tolist(),
        'cell_type': data_query.obs['cell_type_orig'].tolist(),
        'cell_type_transferred': data_query.obs[cell_type_col].tolist(),
    }).write_csv(OUTPUT_PATH_LABELS)

    timers.shutdown()
    print('--- Completed successfully ---')

    print('\n--- Session Info ---')
    sc.logging.print_header()
