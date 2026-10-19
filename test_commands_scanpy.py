import gc
import sys
import scanpy as sc
import anndata as ad
from pathlib import Path
sys.path.append(f'{Path.home()}/brisc-benchmarking')
from utils_local import MemoryTimer, system_info

DATA_NAME = sys.argv[1]
DATA_PATH = sys.argv[2]
OUTPUT_PATH_TIME = sys.argv[3]
REP = int(sys.argv[-1])

if __name__ == '__main__':

    system_info()
    print('--- Params ---')
    print('scanpy commands')
    print(f'{DATA_PATH=}')

    timers = MemoryTimer(
        silent=False, csv_path=OUTPUT_PATH_TIME, summary_unit='ms',
        csv_columns={'library': 'scanpy', 'test': 'commands',
                     'dataset': DATA_NAME, 'rep': REP})

    # Setup
    data = sc.read_h5ad(DATA_PATH)
    data.var['mt'] = data.var_names.str.upper().str.startswith('MT-')
    data.var['malat1'] = data.var_names.str.upper() == 'MALAT1'
    sc.pp.calculate_qc_metrics(
        data, qc_vars=['mt', 'malat1'],
        percent_top=None, log1p=False, inplace=True)
    keep = ((data.obs['n_genes_by_counts'].values >= 100) &
            (data.obs['pct_counts_mt'].values <= 5) &
            (data.obs['pct_counts_malat1'].values > 0))
    data = data[keep].copy()
    sc.pp.normalize_total(data)
    sc.pp.log1p(data)
    sc.pp.highly_variable_genes(
        data, n_top_genes=2000, batch_key='donor')

    gene_idx = 0
    cell_type_select = data.obs['cell_type'].iloc[0]

    with timers('Get expression by gene'):
        data.X[:, gene_idx].toarray().ravel()

    with timers('Subset to one cell type'):
        data[data.obs['cell_type'] == cell_type_select].copy()

    with timers('Subset to highly variable genes'):
        data[:, data.var['highly_variable']].copy()

    with timers('Subsample to 10,000 cells'):
        sc.pp.sample(data, n=10_000, copy=True)

    with timers('Split by cell type'):
        data_split = [data[data.obs['cell_type_broad'] == ct].copy()
                      for ct in data.obs['cell_type_broad'].unique()]

    del data; gc.collect()

    with timers('Concatenate cell types'):
        data = ad.concat(data_split)

    timers.shutdown()
    print('--- Completed successfully ---')
