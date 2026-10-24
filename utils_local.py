import gc
import io
import os
import sys
import time
import socket
import tempfile
import subprocess
import numpy as np
import polars as pl
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from timeit import default_timer
from contextlib import contextmanager
from typing import ContextManager
from tabulate import tabulate

def plural(string: str, count: int) -> str:
    return string if abs(count) == 1 else f'{string}s'

def run(cmd, *, log_file=None, unbuffered=False, pipefail=True,
        num_threads=None, **kwargs):
    run_kwargs = dict(check=True, shell=True, executable='/bin/bash')
    run_kwargs.update(**kwargs)
    thread_str = (
        f'export MKL_NUM_THREADS={num_threads}; '
        f'export OMP_NUM_THREADS={num_threads}; '
        f'export OPENBLAS_NUM_THREADS={num_threads}; '
        f'export NUMEXPR_MAX_THREADS={num_threads}; '
        if num_threads is not None else '')
    return subprocess.run(
        f'{thread_str}'
        f'set -eu{"o pipefail" if pipefail else ""}; '
        f'{"stdbuf -i0 -o0 -e0 " if unbuffered else ""}{cmd}'
        f'{f" 2>&1 | tee {log_file}" if log_file is not None else ""}',
        **run_kwargs)

_MONITOR_MEM_SH_PATH = os.path.join(
    os.path.dirname(__file__), 'monitor_mem.sh')

POLLING_INTERVAL = '0.05'

class MemoryTimer:
    _TIME_CONVERSIONS = {
        's': 1, 'ms': 1000, 'us': 1e6, 'µs': 1e6,
        'ns': 1e9, 'm': 1/60, 'h': 1/3600, 'd': 1/86400
    }

    def __init__(self, silent=True, csv_path=None, csv_columns=None,
                 unit='s', summary_unit=None):
        import atexit
        import signal

        self.timings = {}
        self.silent = silent
        self._csv_path = csv_path
        self._csv_columns = csv_columns or {}
        self._unit = unit
        self._summary_unit = summary_unit or unit
        self._shutdown_done = False

        def _sigterm_handler(signum, frame):
            signal.signal(signal.SIGTERM, signal.SIG_IGN)
            signal.signal(signal.SIGABRT, signal.SIG_IGN)
            raise SystemExit(1)

        signal.signal(signal.SIGTERM, _sigterm_handler)
        signal.signal(signal.SIGABRT, _sigterm_handler)
        atexit.register(self.shutdown)

    def __call__(self, message: str, exclude=False) -> ContextManager[None]:
        pid = os.getpid()
        @contextmanager
        def timer():
            if not self.silent:
                print(f'{message}...', flush=True)
            monitor, tmpf_path = self._start_monitor(pid)
            start = default_timer()
            aborted = False
            try:
                yield
            except BaseException:
                aborted = True
                raise
            finally:
                duration = default_timer() - start
                memory_gb, percent_mem = self._stop_monitor(monitor, tmpf_path)
                if not self.silent:
                    status = 'aborted after' if aborted else 'took'
                    time_str = self._format_time(duration)
                    print(f'{message} {status} {time_str} using '
                          f'{memory_gb} GiB\n', flush=True)
                self._record(message, duration, memory_gb, percent_mem,
                             aborted, exclude)
                gc.collect()

        return timer()

    def _start_monitor(self, pid):
        # Write monitor output to a temp file instead of a pipe.
        # The 64 KB pipe buffer fills during long operations (the
        # main thread is busy and never reads), blocking the monitor
        # and missing late peaks.
        tmpf = tempfile.NamedTemporaryFile(
            prefix='memmon_', suffix='.csv', delete=False)
        tmpf.close()
        stdout_file = open(tmpf.name, 'w')
        monitor = subprocess.Popen(
            [_MONITOR_MEM_SH_PATH, '-p', str(pid), '-i', POLLING_INTERVAL],
            stdout=stdout_file
        )
        stdout_file.close()
        # Sync: wait for first sample (file must have content)
        deadline = time.monotonic() + 5
        while time.monotonic() < deadline:
            if os.path.getsize(tmpf.name) > 0:
                break
            time.sleep(0.005)
        return monitor, tmpf.name

    def _stop_monitor(self, monitor, tmpf_path):
        try:
            monitor.terminate()
            monitor.wait(timeout=5)
        except ProcessLookupError:
            pass
        except subprocess.TimeoutExpired:
            monitor.kill()
            monitor.wait()
        try:
            with open(tmpf_path) as f:
                output = f.read().strip()
        except OSError:
            output = ''
        finally:
            try:
                os.unlink(tmpf_path)
            except OSError:
                pass

        if not output:
            data = np.array([[0.0, 0.0]])
        else:
            try:
                data = np.loadtxt(io.StringIO(output), delimiter=',')
            except ValueError:
                data = np.array([[0.0, 0.0]])

        if data.size == 0:
            data = np.array([[0.0, 0.0]])
        elif data.ndim == 1:
            data = data[None, :]

        max_mem = np.max(data, axis=0)
        memory_gb = np.round(max_mem[0] / 1024 / 1024, 2)
        percent_mem = np.round(max_mem[1], 2)
        return memory_gb, percent_mem

    def _record(self, message, duration, memory_gb, percent_mem, aborted, exclude=False):
        if message in self.timings:
            self.timings[message]['duration'] += duration
            self.timings[message]['memory'] = max(
                self.timings[message]['memory'], memory_gb)
            self.timings[message]['%mem'] = max(
                self.timings[message]['%mem'], percent_mem)
            self.timings[message]['aborted'] = (
                self.timings[message]['aborted'] or aborted)
        else:
            self.timings[message] = {
                'duration': duration,
                'memory': memory_gb,
                '%mem': percent_mem,
                'aborted': aborted,
                'exclude': exclude,
            }
        self._write_csv()

    def _write_csv(self):
        if not self._csv_path or not self.timings:
            return
        df = self.to_dataframe(sort=False, unit=self._unit)
        if self._csv_columns:
            df = df.with_columns(
                [pl.lit(v).alias(k)
                 for k, v in self._csv_columns.items()])
        df.write_csv(self._csv_path)

    def shutdown(self):
        # Write timing summary and CSV. No-op if already called.
        if self._shutdown_done:
            return
        self._shutdown_done = True
        if not self.timings:
            return
        self.print_summary(unit=self._summary_unit)
        self._write_csv()

    def print_summary(self, sort=False, unit=None):
        print('\n--- Timing Summary ---')
        items = (sorted(self.timings.items(),
                        key=lambda x: x[1]['duration'], reverse=True)
                 if sort else list(self.timings.items()))

        total_time = sum(info['duration'] for _, info in items
                         if not info.get('exclude'))
        table_data = []
        duration_header = f'Duration ({unit})' if unit else 'Duration'
        headers = [
            'Operation', 'Status', duration_header,
            '% of Total', 'Memory (GiB)', '% of Avail']

        for message, info in items:
            duration = info['duration']
            memory = info['memory']
            pct = (duration / total_time * 100) if total_time else 0
            status = 'aborted' if info['aborted'] else 'completed'
            time_str = self._format_time(duration, unit)
            table_data.append([
                message,
                status,
                time_str,
                'N/A' if info.get('exclude') else f'{pct:.2f}%',
                f'{memory}',
                'N/A' if info.get('exclude') else f'{info["%mem"]:.2f}%'
            ])

        print(tabulate(table_data, headers=headers, tablefmt='simple'))
        print(f'\nTotal time: {self._format_time(total_time, unit)}')

    def _format_time(self, duration, unit=None):
        if unit:
            if unit not in self._TIME_CONVERSIONS:
                raise ValueError(f'Unsupported unit: {unit}')
            return f'{duration * self._TIME_CONVERSIONS[unit]:.2f}{unit}'

        units = [
            (86400, 'd'), (3600, 'h'), (60, 'm'), (1, 's'),
            (0.001, 'ms'), (1e-6, 'µs'), (1e-9, 'ns')
        ]

        parts = []
        for threshold, suffix in units:
            if duration >= threshold or (not parts and suffix == 'ns'):
                if threshold >= 1:
                    value = int(duration // threshold)
                    duration %= threshold
                else:
                    value = int((duration / threshold) % 1000)
                if value > 0 or (not parts and suffix == 'ns'):
                    parts.append(f'{value}{suffix}')
                if len(parts) == 2:
                    break

        return ' '.join(parts) if parts else 'less than 1ns'

    def to_dataframe(self, sort=True, unit=None):
        if not self.timings:
            return pl.DataFrame({
                'operation': [], 'duration': [], 'duration_unit': [],
                'aborted': [], 'percentage': [], 'memory': [],
                'memory_unit': [], 'percent_mem': []
            })

        items = (sorted(self.timings.items(),
                        key=lambda x: x[1]['duration'], reverse=True)
                 if sort else list(self.timings.items()))

        total = sum(info['duration'] for info in self.timings.values()
                    if not info.get('exclude'))

        conv = 1.0
        if unit:
            if unit not in self._TIME_CONVERSIONS:
                raise ValueError(f'Unsupported unit: {unit}')
            conv = self._TIME_CONVERSIONS[unit]

        rows = [
            {
                'operation': msg,
                'duration': info['duration'] * conv if unit else info['duration'],
                'duration_unit': unit or 's',
                'aborted': info['aborted'],
                'percentage': None if info.get('exclude') else ((info['duration'] / total * 100) if total else 0),
                'memory': info['memory'],
                'memory_unit': 'GiB',
                'percent_mem': None if info.get('exclude') else info['%mem'],
            }
            for msg, info in items
        ]

        return pl.DataFrame(rows)

def print_df(df, num_rows=-1, num_columns=-1):
    with pl.Config(tbl_rows=num_rows, tbl_cols=num_columns):
        print(df)

def system_info():
    hostname = socket.gethostname()
    user = os.environ.get('USER', 'N/A')
    cpu_cores = (os.environ.get('SLURM_CPUS_PER_TASK') or
                 os.environ.get('SLURM_CPUS_ON_NODE'))
    if not cpu_cores:
        try:
            cpu_cores = os.cpu_count()
        except NotImplementedError:
            cpu_cores = 'N/A'

    mem_gb = 'N/A'
    try:
        mem_mb = int(os.environ['SLURM_MEM_PER_NODE'])
        mem_gb = f'{mem_mb / 1024:.1f} GB'
    except (KeyError, ValueError):
        try:
            with open('/proc/meminfo', 'r') as f:
                for line in f:
                    if line.startswith('MemTotal:'):
                        mem_kb = int(line.split()[1])
                        mem_gb = f'{mem_kb / 1024 / 1024:.1f} GB'
                        break
        except FileNotFoundError:
            pass

    print(f'\n--- User Resource Allocation ---', flush=True)
    print(f'Node: {hostname}', flush=True)
    print(f'User: {user}', flush=True)
    print(f'CPU Cores Allocated: {cpu_cores}', flush=True)
    print(f'Memory Allocated: {mem_gb}', flush=True)
    print(f'Python Version: {sys.version}', flush=True)
    try:
        from importlib.metadata import version
        print(f'brisc Version: {version("brisc")}', flush=True)
    except ImportError:
        pass

    gpu_info = 'N/A'
    try:
        result = subprocess.run(
            ['nvidia-smi', '--query-gpu=name,memory.total',
             '--format=csv,noheader'],
            capture_output=True, text=True, timeout=5)
        if result.returncode == 0 and result.stdout.strip():
            gpu_info = result.stdout.strip()
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    print(f'GPU: {gpu_info}', flush=True)

    print(flush=True)

def run_slurm(cmd, *, job_name='job', log_file=None,
              CPUs=192 if os.environ.get('CLUSTER') == 'trillium' else
                   96 if os.environ.get('CLUSTER') == 'trillium-gpu' else 1,
              GPUs=4 if os.environ.get('CLUSTER') == 'trillium-gpu' else 0,
              days=None, hours=None, memory=None,
              partition='compute' if os.environ.get('CLUSTER') == 'trillium'
                         else 'compute_full_node' if os.environ.get('CLUSTER') == 'trillium-gpu'
                         else None,
              account=None, verbose=False):
    cluster = os.environ.get('CLUSTER')
    runtime = f'{days}-00:00:00' if days is not None else \
        f'{hours}:00:00' if hours is not None else '1-00:00:00' \
            if partition != 'debug' else '1:00:00'
    if verbose:
        memory_description = f' and {memory[:-1]} {memory[-1]}iB memory' \
            if memory is not None else ''
        partition_description = \
            f' on the {partition} partition' if partition is not None else ''
        print(f'Requesting {CPUs} {plural("CPU", CPUs)}, {GPUs} '
              f'{plural("GPU", GPUs)}{memory_description} for {runtime}'
              f'{partition_description}...')
    job_name = job_name.replace(' ', '_')
    from tempfile import NamedTemporaryFile
    try:
        with NamedTemporaryFile('w', dir=os.environ.get('SCRATCH', '.'),
                                suffix='.sh', delete=False) as temp_file:
            partition_settings = f'#SBATCH -p {partition}\n' \
                if partition is not None else ''
            account_settings = f'#SBATCH --account={account}\n' \
                if account is not None else ''
            cpu_settings = '' if cluster.startswith('trillium') else \
                f'#SBATCH -c {CPUs}\n'
            gpu_settings = f'#SBATCH --gpus-per-node=h100:{GPUs}\n' \
                if GPUs > 0 else ''
            memory_settings = f'#SBATCH --mem {memory}\n' \
                if memory is not None else ''
            print(
                f'#!/bin/bash\n'
                f'{partition_settings}'
                f'{account_settings}'
                f'#SBATCH -N 1\n'
                f'#SBATCH -n 1\n'
                f'{gpu_settings}'
                f'{cpu_settings}'
                f'{memory_settings}'
                f'#SBATCH -t {runtime}\n'
                f'#SBATCH -J {job_name}\n'
                f'{f"#SBATCH -o {log_file}" if log_file is not None else ""}\n'
                f'#SBATCH --signal=B:TERM@30\n'
                f'export PYTHONUNBUFFERED=1\n'
                f'export R_LIBS_USER=/home/wainberg/R/x86_64-pc-linux-gnu-library/4.4\n'
                f'export OMP_PLACES=cores\n'
                f'export OMP_PROC_BIND=spread\n'
                f'export NUMBA_NUM_THREADS=${{SLURM_CPUS_ON_NODE:-$(nproc)}}\n'
                f'export NUMBA_CACHE_DIR=${{SCRATCH:-/scratch/$USER}}/.cache/numba\n'
                f'export PYTHONFAULTHANDLER=1\n'
                f'export HDF5_USE_FILE_LOCKING=FALSE\n'
                f'set -m\n'
                f'{cmd} &\n'
                f'CHILD=$!\n'
                f'trap "kill -INT -$CHILD 2>/dev/null; '
                f'wait $CHILD 2>/dev/null; exit \\$?" TERM\n'
                f'wait $CHILD\n'
                f'exit $?\n',
                file=temp_file)
        sbatch = '.sbatch' if cluster.startswith('trillium') else 'sbatch'
        sbatch_message = run(f'{sbatch} {temp_file.name}',
                             stdout=subprocess.PIPE)\
            .stdout.decode().rstrip('\n')
        print(f'{sbatch_message} ("{job_name}")')
    finally:
        try:
            os.unlink(temp_file.name)
        except NameError:
            pass

def remove_self_neighbors(neighbors):
    # Drop each cell's self-neighbor (value == row index); keep all if none.
    n, m = neighbors.shape
    drop = neighbors == np.arange(n)[:, None]
    if not drop.any():
        return neighbors
    drop[~drop.any(axis=1), -1] = True
    return neighbors[~drop].reshape(n, m - 1)

def transfer_accuracy(obs, orig_col, trans_col):
    if isinstance(obs, pd.DataFrame):
        obs = pl.from_pandas(obs)
    groups = obs\
        .with_columns(pl.col(orig_col, trans_col).cast(pl.String))\
        .group_by(orig_col)\
        .agg(n_correct=pl.col(orig_col).eq(pl.col(trans_col)).sum(),
             n_total=pl.len())
    df = pl.concat([
            groups,
            groups.sum().with_columns(pl.lit('Total').alias(orig_col))])\
        .with_columns(
            percent_correct=pl.col('n_correct') / pl.col('n_total') * 100)\
        .sort(pl.col(orig_col).eq('Total'), orig_col)\
        .rename({orig_col: 'cell_type'})
    print_df(df)
    return df
