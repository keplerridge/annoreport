#!/usr/bin/env python3
"""
benchmark_uhgg.py

Benchmarks annoreport on the UHGG human gut MAG subset (GFF-only mode).
Runs with and without UniProt enrichment, N_RUNS times each.
Reports mean and std dev for runtime and peak memory.
"""

import subprocess
import time
import os
import statistics
from datetime import datetime

ANNOREPORT = os.path.expanduser('~/groups/grp_ktr_research/annoreport/annotation_report.py') # path to annotation_report.py (default assumes same directory)
UHGG_DIR   = os.path.expanduser('/path/to/uhgg_gff_files') # directory containing UHGG .gff.gz files
OUTBASE    = os.path.expanduser('/path/to/benchmark_output') # directory where benchmark results will be written
N_RUNS     = 5

CONFIGURATIONS = [
    {
        'label':      'UHGG 1500 MAGs — GFF only, with UniProt',
        'no_uniprot': False,
    },
    {
        'label':      'UHGG 1500 MAGs — GFF only, no UniProt',
        'no_uniprot': True,
    },
]


def run_single(outdir, no_uniprot, run_num):
    if no_uniprot:
        cmd = [
            '/usr/bin/time', '-v',
            'python3', os.path.expanduser(ANNOREPORT),
            '--annotation_dir', UHGG_DIR,
            '--outdir', outdir,
            '--no_uniprot',
        ]
    else:
        cmd = [
            '/usr/bin/time', '-v',
            'python3', os.path.expanduser(ANNOREPORT),
            '--annotation_dir', UHGG_DIR,
            '--outdir', outdir,
        ]

    start   = time.time()
    result  = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - start

    peak_mem_kb = None
    for line in result.stderr.splitlines():
        if 'Maximum resident set size' in line:
            peak_mem_kb = int(line.split(':')[-1].strip())
            break

    peak_mem_mb = peak_mem_kb / 1024 if peak_mem_kb else None

    print(f"    Run {run_num}: {elapsed:.1f}s | {peak_mem_mb:.1f} MB")

    return elapsed, peak_mem_mb


def main():
    os.makedirs(OUTBASE, exist_ok=True)
    all_results = []

    print(f"\nannoreport UHGG Benchmark")
    print(f"{'=' * 60}")
    print(f"Dataset:           UHGG human gut MAGs (GFF-only)")
    print(f"Bins:              1500")
    print(f"Runs per config:   {N_RUNS}")
    print(f"Start time:        {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'=' * 60}\n")

    for i, config in enumerate(CONFIGURATIONS, 1):
        label      = config['label']
        no_uniprot = config['no_uniprot']

        print(f"[{i}/2] {label}")

        times = []
        mems  = []

        for run_num in range(1, N_RUNS + 1):
            outdir = os.path.join(OUTBASE, f"run_{i}_{run_num}")
            os.makedirs(outdir, exist_ok=True)

            elapsed, mem_mb = run_single(outdir, no_uniprot, run_num)
            times.append(elapsed)
            if mem_mb:
                mems.append(mem_mb)

        mean_time = statistics.mean(times)
        std_time  = statistics.stdev(times) if len(times) > 1 else 0
        mean_mem  = statistics.mean(mems) if mems else None
        std_mem   = statistics.stdev(mems) if len(mems) > 1 else 0

        result = {
            'label':     label,
            'mean_time': round(mean_time, 2),
            'std_time':  round(std_time, 2),
            'mean_mem':  round(mean_mem, 1) if mean_mem else None,
            'std_mem':   round(std_mem, 1),
            'runs':      list(zip([round(t, 2) for t in times],
                                  [round(m, 1) for m in mems]))
        }
        all_results.append(result)

        print(f"    Mean: {mean_time:.1f}s ± {std_time:.1f}s | "
              f"{mean_mem:.1f} ± {std_mem:.1f} MB\n")

    outfile = os.path.join(OUTBASE, 'benchmark_uhgg_summary.txt')
    with open(outfile, 'w') as f:
        f.write('annoreport UHGG Benchmark Results\n')
        f.write('=' * 60 + '\n')
        f.write(f'Date:              {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n')
        f.write(f'Dataset:           UHGG human gut MAGs (GFF-only mode)\n')
        f.write(f'Bins:              1500\n')
        f.write(f'Runs per config:   {N_RUNS}\n')
        f.write(f'top_n:             100\n')
        f.write('=' * 60 + '\n\n')

        for r in all_results:
            f.write(f"Configuration: {r['label']}\n")
            f.write(f"  Mean runtime:     {r['mean_time']}s ± {r['std_time']}s\n")
            f.write(f"  Mean peak memory: {r['mean_mem']} MB ± {r['std_mem']} MB\n")
            f.write(f"  Individual runs:  ")
            for run_num, (t, m) in enumerate(r['runs'], 1):
                f.write(f"run{run_num}: {t}s / {m}MB  ")
            f.write('\n\n')

    print(f"Benchmark complete.")
    print(f"Results written to: {outfile}")


if __name__ == '__main__':
    main()
