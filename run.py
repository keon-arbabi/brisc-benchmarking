#!/usr/bin/env python
import argparse
import os
import sys
import subprocess
from pathlib import Path
_HOME_DIR = Path.home()
sys.path.append(f'{_HOME_DIR}/sc-benchmarking')
from utils_local import run_slurm

DATA_DIR = f'{_HOME_DIR}/single-cell'
WORK_DIR = f'{_HOME_DIR}/sc-benchmarking'

OUTPUT_DIR = os.path.join(WORK_DIR, 'output')
LOG_DIR = os.path.join(WORK_DIR, 'logs')
FIGURES_DIR = os.path.join(WORK_DIR, 'figures')

RSCRIPT = '/home/wainberg/bin/Rscript-4.5.3'
PYTHON = '/home/wainberg/bin/python3.14'
PYTHON_RAPIDS = '/home/karbabi/miniforge3/envs/rapids_singlecell/bin/python'

DATASET_NAMES = [
    'SEAAD',
    'Parse',
    'PanSci',
]
DATASETS = {
    name: os.path.join(DATA_DIR, name, f'{name}_raw.h5ad')
    for name in DATASET_NAMES
}
THREADS = [-1, 1]
N_REPS = 3
# (file, tool, task, thread_params, gpu)
SCRIPTS = [
    ('test_basic_brisc.py', 'brisc', 'basic', THREADS, False),
    ('test_basic_brisc.py', 'brisc', 'basic', [-1], True),
    ('test_de_brisc.py', 'brisc', 'de', THREADS, False),
    ('test_transfer_brisc.py', 'brisc', 'transfer', THREADS, False),
    ('test_commands_brisc.py', 'brisc', 'commands', THREADS, False),

    ('test_basic_rapids.py', 'rapids', 'basic', None, True),

    ('test_basic_scanpy.py', 'scanpy', 'basic', None, False),
    ('test_de_scanpy.py', 'scanpy', 'de', None, False),
    ('test_transfer_scanpy.py', 'scanpy', 'transfer', None, False),
    ('test_commands_scanpy.py', 'scanpy', 'commands', None, False),

    ('test_basic_seurat.R', 'seurat', 'basic', None, False),
    ('test_de_seurat.R', 'seurat', 'de', None, False),
    ('test_transfer_seurat.R', 'seurat', 'transfer', None, False),
    ('test_commands_seurat.R', 'seurat', 'commands', None, False),
]
SHORT = {
    'basic': 'ba', 'de': 'de', 'transfer': 'tr', 'commands': 'cm',
    'brisc': 'br', 'scanpy': 'sc', 'seurat': 'sr', 'rapids': 'rp',
    'SEAAD': 'SE', 'Parse': 'PA', 'PanSci': 'PS',
}
OUTPUTS = {
    'basic': ['embedding', 'pcs', 'neighbors'],
    'de': ['de'],
    'transfer': ['accuracy', 'labels'],
    'commands': [],
}

_OUTPUT_EXT = {'pcs': 'parquet', 'neighbors': 'parquet'}

def out(job_name, suffix, rep):
    ext = _OUTPUT_EXT.get(suffix, 'csv')
    return os.path.join(
        OUTPUT_DIR, f'rep{rep}', f'{job_name}_{suffix}.{ext}')

def _job_name(head, dataset, threads, gpu, gpu_tag):
    parts = [head, dataset]
    if threads is not None:
        parts.append(str(threads))
    if gpu:
        parts.append(gpu_tag)
    return '_'.join(parts)

def _slurm_job_pending_or_running(job_name):
    try:
        result = subprocess.run(
            ['squeue', '--me', f'--name={job_name}', '--noheader'],
            capture_output=True, text=True, timeout=10)
        return bool(result.stdout.strip())
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False

def run_jobs(scripts):
    for script_file, tool, task, thread_params, gpu in scripts:
        if script_file.endswith('.R'):
            interpreter = RSCRIPT
        elif tool == 'rapids':
            interpreter = f'PYTHONPATH= {PYTHON_RAPIDS}'
        else:
            interpreter = PYTHON
        script_path = os.path.join(WORK_DIR, script_file)

        for d_name, d_path in DATASETS.items():
            for threads in (thread_params or [None]):
                for rep in range(1, N_REPS + 1):
                    job_name = _job_name(
                        f'{task}_{tool}', d_name, threads, gpu, 'gpu')

                    log = os.path.join(
                        LOG_DIR, f'rep{rep}', f'{job_name}.log')
                    if os.path.exists(log):
                        with open(log) as f:
                            if 'Completed successfully' in f.read():
                                print(f'Skipping completed job: '
                                      f'{job_name} rep{rep}')
                                continue

                    slurm_name = _job_name(
                        SHORT[task] + SHORT[tool], SHORT[d_name],
                        threads, gpu, 'g') + f'_r{rep}'

                    if _slurm_job_pending_or_running(slurm_name):
                        print(f'Skipping in-progress job: '
                              f'{job_name} rep{rep}')
                        continue

                    cmd = [interpreter, script_path, d_name, d_path]
                    if threads is not None:
                        cmd.append(str(threads))
                    cmd.append(out(job_name, 'timer', rep))
                    cmd.extend(out(job_name, s, rep) for s in OUTPUTS[task])
                    cmd.append(str(rep))

                    run_slurm(
                        ' '.join(cmd),
                        account='def-wainberg',
                        job_name=slurm_name, log_file=log, hours=24)

def filter_scripts_and_datasets(filters):
    tools = {t for _, t, _, _, _ in SCRIPTS}
    tasks = {t for _, _, t, _, _ in SCRIPTS}
    datasets = set(DATASET_NAMES)

    sel_tools, sel_tasks, sel_datasets = set(), set(), set()
    for a in filters:
        if a in tools:
            sel_tools.add(a)
        elif a in tasks:
            sel_tasks.add(a)
        elif a in datasets:
            sel_datasets.add(a)
        else:
            raise SystemExit(
                f'Unknown argument: {a!r}. '
                f'Tools: {sorted(tools)}, tasks: {sorted(tasks)}, '
                f'datasets: {sorted(datasets)}')

    scripts = SCRIPTS
    if sel_tools:
        scripts = [s for s in scripts if s[1] in sel_tools]
    if sel_tasks:
        scripts = [s for s in scripts if s[2] in sel_tasks]

    selected_datasets = (
        {n: DATASETS[n] for n in sel_datasets} if sel_datasets else DATASETS)
    return scripts, selected_datasets

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filters', nargs='*',
        help='Any combination of tool/task/dataset names to filter jobs')
    parser.add_argument('--log-dir', default=LOG_DIR)
    parser.add_argument('--output-dir', default=OUTPUT_DIR)
    parser.add_argument('--reps', type=int, default=N_REPS)
    args = parser.parse_args()

    LOG_DIR = args.log_dir
    OUTPUT_DIR = args.output_dir
    N_REPS = args.reps
    for rep in range(1, N_REPS + 1):
        os.makedirs(os.path.join(OUTPUT_DIR, f'rep{rep}'), exist_ok=True)
        os.makedirs(os.path.join(LOG_DIR, f'rep{rep}'), exist_ok=True)
    os.makedirs(FIGURES_DIR, exist_ok=True)

    is_gpu = os.environ.get('CLUSTER') == 'trillium-gpu'
    scripts, DATASETS = filter_scripts_and_datasets(args.filters)
    scripts = [s for s in scripts if s[4] == is_gpu]
    run_jobs(scripts)
