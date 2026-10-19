import gc
import sys
from pathlib import Path
import numpy as np
import polars as pl
sys.path.append(f'{Path.home()}/brisc-benchmarking')
from brisc import SingleCell, concat_obs
from utils_local import MemoryTimer, system_info

DATA_NAME = sys.argv[1]
DATA_PATH = sys.argv[2]
NUM_THREADS = int(sys.argv[3])
OUTPUT_PATH_TIME = sys.argv[4]
REP = int(sys.argv[-1])

if __name__ == '__main__':

    system_info()
    print('--- Params ---')
    print('brisc commands')
    print(f'{DATA_PATH=}')
    print(f'{NUM_THREADS=}')

    timers = MemoryTimer(
        silent=False, csv_path=OUTPUT_PATH_TIME, summary_unit='ms',
        csv_columns={'library': 'brisc', 'test': 'commands',
                    'dataset': DATA_NAME, 'num_threads': NUM_THREADS,
                    'rep': REP})

    # Setup
    data = SingleCell(DATA_PATH, num_threads=NUM_THREADS)\
        .qc(allow_float=True)\
        .hvg(batch_column='donor')\
        .normalize()

    gene_name = data.var_names[0]
    cell_type_select = data.obs['cell_type'][0]

    with timers('Get expression by gene'):
        data.gene(gene_name)

    with timers('Subset to one cell type'):
        data.filter_obs(pl.col('cell_type').eq(cell_type_select))

    with timers('Subset to highly variable genes'):
        data.filter_var(pl.col('highly_variable'))

    with timers('Subsample to 10,000 cells'):
        data.subsample_obs(n=10_000)

    with timers('Split by cell type'):
        data_split = data.split_by_obs('cell_type_broad')

    del data; gc.collect()

    with timers('Concatenate cell types'):
        data = concat_obs(data_split.values())

    timers.shutdown()
    print('--- Completed successfully ---')
