import os
import subprocess
from pathlib import Path
LOG_DIR = f'{Path.home()}/brisc-benchmarking/logs'

SHORT = {
    'basic': 'ba', 'de': 'de', 'transfer': 'tr', 'commands': 'cm',
    'brisc': 'br', 'scanpy': 'sc', 'seurat': 'sr', 'rapids': 'rp',
    'SEAAD': 'SE', 'Parse': 'PA', 'PanSci': 'PS',
}

def slurm_short_name(job_name):
    parts = job_name.split('_')
    if len(parts) < 3 or any(p not in SHORT for p in parts[:3]):
        return None
    short = [SHORT[parts[0]] + SHORT[parts[1]], SHORT[parts[2]]]
    rest = parts[3:]
    is_gpu = rest[-1:] == ['gpu']
    if is_gpu:
        rest = rest[:-1]
    if rest:
        short.append(rest[0])
    if is_gpu:
        short.append('g')
    return '_'.join(short)

try:
    result = subprocess.run(
        ['squeue', '--me', '--noheader', '-o', '%j'],
        capture_output=True, text=True, timeout=10)
    running_names = set(result.stdout.split())
except (FileNotFoundError, subprocess.TimeoutExpired):
    running_names = set()

running = []
failed = []
for rep_dir in sorted(os.listdir(LOG_DIR)):
    rep_path = os.path.join(LOG_DIR, rep_dir)
    if not (rep_dir.startswith('rep') and os.path.isdir(rep_path)):
        continue
    rep = rep_dir[3:]
    for log_file in os.listdir(rep_path):
        if not log_file.endswith('.log'):
            continue
        job_name = log_file[:-4]
        with open(os.path.join(rep_path, log_file)) as f:
            content = f.read()
        if 'Completed successfully' in content:
            continue
        short = slurm_short_name(job_name)
        if short:
            short = f'{short}_r{rep}'
        label = f'{job_name} rep{rep}'
        if short and short in running_names:
            running.append(label)
        elif content.strip():
            failed.append(label)

for label, jobs in [('Running', running), ('Failed', failed)]:
    if jobs:
        print(f'{label} jobs ({len(jobs)}):')
        for job in sorted(jobs):
            print(f'  {job}')

if not running and not failed:
    print('No failed or running jobs found.')
