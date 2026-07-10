import gc
import sys
from pathlib import Path
import polars as pl
sys.path.append(f'{Path.home()}/sc-benchmarking')
from brisc import SingleCell
from utils_local import MemoryTimer, system_info

DATA_NAME = sys.argv[1]
DATA_PATH = sys.argv[2]
NUM_THREADS = int(sys.argv[3])
OUTPUT_PATH_TIME = sys.argv[4]
OUTPUT_PATH_DE = sys.argv[5]
REP = int(sys.argv[-1])

if __name__ == '__main__':

    system_info()
    print('--- Params ---')
    print('brisc de')
    print(f'{DATA_PATH=}')
    print(f'{NUM_THREADS=}')

    if DATA_NAME == 'SEAAD':
        formula = '~ 0 + cond + apoe4_dosage + sex + age_at_death + pmi + ' \
            'log2(num_cells) + log2(library_size)'
        contrasts = {'AD_vs_Control': '`condAD` - `condControl`'}
    elif DATA_NAME == 'Parse':
        formula = '~ 0 + cond + donor + log2(num_cells) + log2(library_size)'
        contrasts = {'IFN-gamma_vs_PBS': '`condIFN-gamma` - `condPBS`'}
    elif DATA_NAME == 'PanSci':
        formula = '~ 0 + cond + sex + log2(num_cells) + log2(library_size)'
        contrasts = {'Aged_vs_Young': '`condAged` - `condYoung`'}

    timers = MemoryTimer(
        silent=False, csv_path=OUTPUT_PATH_TIME,
        csv_columns={'library': 'brisc', 'test': 'de',
                     'dataset': DATA_NAME, 'num_threads': NUM_THREADS,
                     'rep': REP})

    with timers('Load data'):
        data_sc = SingleCell(DATA_PATH, num_threads=NUM_THREADS)

    with timers('Quality control'):
        data_sc = data_sc.qc(
            custom_filter=pl.col('cond').is_not_null(),
            allow_float=True)

    with timers('Pseudobulk'):
        data_pb = data_sc.pseudobulk('sample', 'cell_type')

    del data_sc; gc.collect()

    with timers('Quality control'):
        data_pb = data_pb.qc('cond', verbose=False)

    with timers('Differential expression'):
        de = data_pb\
            .library_size()\
            .de(formula,
                contrasts=contrasts,
                verbose=False)

    de_df = de.table\
        .select('cell_type', 'gene', 'logFC', 'p', 'FDR')\
        .rename({'p': 'p_value', 'FDR': 'p_value_adj'})
    de_df.write_csv(OUTPUT_PATH_DE)

    timers.shutdown()
    print('--- Completed successfully ---')
