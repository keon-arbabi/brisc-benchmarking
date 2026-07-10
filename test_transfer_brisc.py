import gc
import sys
from pathlib import Path
sys.path.append(f'{Path.home()}/sc-benchmarking')
from brisc import SingleCell
from utils_local import MemoryTimer, system_info, transfer_accuracy

DATA_NAME = sys.argv[1]
DATA_PATH = sys.argv[2]
NUM_THREADS = int(sys.argv[3])
OUTPUT_PATH_TIME = sys.argv[4]
OUTPUT_PATH_ACC = sys.argv[5]
OUTPUT_PATH_LABELS = sys.argv[6]
REP = int(sys.argv[-1])

if __name__ == '__main__':

    system_info()
    print('--- Params ---')
    print('brisc transfer')
    print(f'{DATA_PATH=}')
    print(f'{NUM_THREADS=}')

    timers = MemoryTimer(
        silent=False, csv_path=OUTPUT_PATH_TIME,
        csv_columns={'library': 'brisc', 'test': 'transfer',
                     'dataset': DATA_NAME, 'num_threads': NUM_THREADS,
                     'rep': REP})

    with timers('Load data'):
        data = SingleCell(DATA_PATH, num_threads=NUM_THREADS)

    with timers('Quality control'):
        data = data.qc(allow_float=True)

    with timers('Split data'):
        data_query, data_ref = data.split_by_obs('is_ref').values()

    del data; gc.collect()

    with timers('Feature selection'):
        data_ref, data_query = data_ref.hvg(data_query, batch_column='donor')

    with timers('Normalization'):
        data_ref = data_ref.normalize()
        data_query = data_query.normalize()

    with timers('PCA'):
        data_ref, data_query = data_ref.pca(data_query)

    cell_type_col = 'cell_type_fine' if DATA_NAME == 'PanSci' else 'cell_type'

    with timers('Transfer labels'):
        data_ref, data_query = data_ref.harmonize(data_query)
        data_query = data_query.label_transfer_from(
            data_ref,
            cell_type_col,
            cell_type_column='cell_type_transferred')

    print('--- Transfer Accuracy ---')
    transfer_accuracy(
        data_query.obs, cell_type_col, 'cell_type_transferred')\
    .write_csv(OUTPUT_PATH_ACC)

    data_query.obs\
        .select(['cell_id', cell_type_col, 'cell_type_transferred'])\
        .rename({cell_type_col: 'cell_type'})\
        .write_csv(OUTPUT_PATH_LABELS)

    timers.shutdown()
    print('--- Completed successfully ---')
