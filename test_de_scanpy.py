import os
import gc
import sys
import decoupler as dc
import numpy as np
import polars as pl
import scanpy as sc
from formulaic import Formula
from pydeseq2.ds import DeseqStats
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pathlib import Path
sys.path.append(f'{Path.home()}/brisc-benchmarking')
from utils_local import MemoryTimer, system_info

DATA_NAME = sys.argv[1]
DATA_PATH = sys.argv[2]
OUTPUT_PATH_TIME = sys.argv[3]
OUTPUT_PATH_DE = sys.argv[4]
REP = int(sys.argv[-1])

if __name__ == '__main__':

    system_info()
    print('--- Params ---')
    print('scanpy de')
    print(f'{DATA_PATH=}')

    if DATA_NAME == 'SEAAD':
        sample_cols = ['cond', 'apoe4_dosage', 'sex', 'age_at_death', 'pmi']
        contrast = ['cond', 'AD', 'Control']
        design = ('~ cond + apoe4_dosage + sex + age_at_death + pmi +'
                  'log2_psbulk_cells + log2_psbulk_counts')
    elif DATA_NAME == 'Parse':
        sample_cols = ['cond', 'donor']
        contrast = ['cond', 'IFN-gamma', 'PBS']
        design = '~ cond + donor + log2_psbulk_cells + log2_psbulk_counts'
    elif DATA_NAME == 'PanSci':
        sample_cols = ['cond', 'sex']
        contrast = ['cond', 'Aged', 'Young']
        design = '~ cond + sex + log2_psbulk_cells + log2_psbulk_counts'

    timers = MemoryTimer(
        silent=False, csv_path=OUTPUT_PATH_TIME,
        csv_columns={'library': 'scanpy', 'test': 'de',
                     'dataset': DATA_NAME, 'rep': REP})

    with timers('Load data'):
        data_sc = sc.read_h5ad(DATA_PATH)

    with timers('Quality control'):
        data_sc.var['mt'] = data_sc.var_names.str.upper().str.startswith('MT-')
        data_sc.var['malat1'] = data_sc.var_names.str.upper() == 'MALAT1'
        sc.pp.calculate_qc_metrics(
            data_sc, qc_vars=['mt', 'malat1'],
            percent_top=None, log1p=False, inplace=True)
        keep = ((data_sc.obs['n_genes_by_counts'].values >= 100) &
                (data_sc.obs['pct_counts_mt'].values <= 5) &
                (data_sc.obs['pct_counts_malat1'].values > 0) &
                (data_sc.obs['cond'].notna().values))
        data_sc = data_sc[keep].copy()

    with timers('Pseudobulk'):
        data_pb = sc.get.aggregate(
            data_sc, by=['sample', 'cell_type'], func='sum')
        data_pb.X = data_pb.layers.pop('sum')
        sample_meta = data_sc.obs[['sample'] + sample_cols]\
            .drop_duplicates('sample').set_index('sample')
        for col in sample_cols:
            data_pb.obs[col] = data_pb.obs['sample'].map(sample_meta[col])
        data_pb.obs['psbulk_cells'] = data_pb.obs['n_obs_aggregated']
        data_pb.obs['psbulk_counts'] = np.asarray(
            data_pb.X.sum(axis=1)).ravel()
        data_pb.obs['log2_psbulk_cells'] = np.log2(
            data_pb.obs['psbulk_cells'])
        data_pb.obs['log2_psbulk_counts'] = np.log2(
            data_pb.obs['psbulk_counts'])

    del data_sc; gc.collect()

    # z-score so exp(Xb) in IRLS doesn't overflow
    for col in ('age_at_death', 'pmi',
                'log2_psbulk_cells', 'log2_psbulk_counts'):
        if col in data_pb.obs:
            s = data_pb.obs[col]
            data_pb.obs[col] = (s - s.mean()) / s.std()

    with timers('Quality control'):
        dc.pp.filter_samples(data_pb, min_cells=10, min_counts=1000)

    # drop donors lacking both cond levels per cell type (post-sample-filter)
    if 'donor' in data_pb.obs:
        paired = (data_pb.obs.groupby(['cell_type', 'donor'], observed=True)
                  ['cond'].transform('nunique') == 2)
        data_pb = data_pb[paired.values].copy()

    with timers('Differential expression'):
        inference = DefaultInference(n_cpus=16)
        cell_types = data_pb.obs['cell_type'].unique()

        de = {}
        for ct in cell_types:
            data_pb_ct = data_pb[data_pb.obs['cell_type'] == ct].copy()

            # skip saturated or rank-deficient designs
            mm = Formula(design).get_model_matrix(data_pb_ct.obs)
            if mm.shape[0] <= mm.shape[1] \
                    or np.linalg.matrix_rank(mm.values) < mm.shape[1]:
                print(f'Skipping {ct}: design not full rank '
                      f'({mm.shape[0]} samples, {mm.shape[1]} cols)')
                continue

            dc.pp.filter_by_expr(
                data_pb_ct, group='cond', min_count=10)
            dds = DeseqDataSet(
                adata=data_pb_ct, design=design,
                refit_cooks=True, inference=inference, quiet=True)
            dds.deseq2()
            stat_res = DeseqStats(
                dds, contrast=contrast,
                inference=inference, quiet=True)
            stat_res.summary()
            de[ct] = stat_res.results_df

    de_df = pl.concat([
        pl.DataFrame({
            'cell_type': ct,
            'gene': df.index.tolist(),
            'logFC': df['log2FoldChange'].values,
            'p_value': df['pvalue'].values,
            'p_value_adj': df['padj'].values,
        })
        for ct, df in de.items()])
    de_df.write_csv(OUTPUT_PATH_DE)

    timers.shutdown()
    print('--- Completed successfully ---')

    print('\n--- Session Info ---')
    sc.logging.print_header()
