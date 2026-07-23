#!/usr/bin/env python3
"""
benchmark_annoreport.py

Runs annoreport across all input/mode combinations and records
runtime and peak memory for each. Reports mean and std dev across 5 runs.

Usage:
    1. Set BASE to the directory containing your annotation results
       (must contain subfolders: prokka, bakta, prokka_gff, bakta_gff)
    2. Set OUTBASE to the directory where the benchmark output will be written
    3. Run: python3 benchmark_annoreport.py
"""

import subprocess
import time
import os
import statistics
from datetime import datetime

ANNOREPORT = os.path.expanduser('annotation_report.py') # path to annotation_report.py (default assumes same directory)
BASE       = os.path.expanduser('/path/to/annotation_results') # directory containing prokka, bakta, prokka_gff, bakta_gff subfolders
OUTBASE    = os.path.expanduser('/path/to/benchmark_output') # directory where benchmark results will be written
N_RUNS     = 5

CONFIGURATIONS = [
    {
        'label':         'Prokka full output — with UniProt',
        'annotation_dir': os.path.join(BASE, 'prokka'),
        'no_uniprot':    False,
    },
    {
        'label':         'Prokka full output — no UniProt',
        'annotation_dir': os.path.join(BASE, 'prokka'),
        'no_uniprot':    True,
    },
    {
        'label':         'Prokka GFF only — with UniProt',
        'annotation_dir': os.path.join(BASE, 'prokka_gff'),
        'no_uniprot':    False,
    },
    {
        'label':         'Prokka GFF only — no UniProt',
        'annotation_dir': os.path.join(BASE, 'prokka_gff'),
        'no_uniprot':    True,
    },
    {
        'label':         'Bakta full output — with UniProt',
        'annotation_dir': os.path.join(BASE, 'bakta'),
        'no_uniprot':    False,
    },
    {
        'label':         'Bakta full output — no UniProt',
        'annotation_dir': os.path.join(BASE, 'bakta'),
        'no_uniprot':    True,
    },
    {
        'label':         'Bakta GFF only — with UniProt',
        'annotation_dir': os.path.join(BASE, 'bakta_gff'),
        'no_uniprot':    False,
    },
    {
        'label':         'Bakta GFF only — no UniProt',
        'annotation_dir': os.path.join(BASE, 'bakta_gff'),
        'no_uniprot':    True,
    },
]


def run_single(annotation_dir, outdir, no_uniprot, run_num):
    cmd = [
        '/usr/bin/time', '-v',
        'python3', os.path.expanduser(ANNOREPORT),
        '--annotation_dir', annotation_dir,
        '--outdir', outdir,
        '--no_uniprot' if no_uniprot else '--top_n', '100',
    ]

    # Clean up --top_n if no_uniprot since we don't need the extra arg
    if no_uniprot:
        cmd = [
            '/usr/bin/time', '-v',
            'python3', os.path.expanduser(ANNOREPORT),
            '--annotation_dir', annotation_dir,
            '--outdir', outdir,
            '--no_uniprot',
        ]
    else:
        cmd = [
            '/usr/bin/time', '-v',
            'python3', os.path.expanduser(ANNOREPORT),
            '--annotation_dir', annotation_dir,
            '--outdir', outdir,
        ]

    start = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - start

    # Parse peak memory from /usr/bin/time -v stderr
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

    print(f"\nannoreport Benchmark")
    print(f"{'=' * 60}")
    print(f"Runs per configuration: {N_RUNS}")
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'=' * 60}\n")

    for i, config in enumerate(CONFIGURATIONS, 1):
        label         = config['label']
        annotation_dir = config['annotation_dir']
        no_uniprot    = config['no_uniprot']

        print(f"[{i}/8] {label}")

        times   = []
        mems    = []

        for run_num in range(1, N_RUNS + 1):
            outdir = os.path.join(OUTBASE, f"run_{i}_{run_num}")
            os.makedirs(outdir, exist_ok=True)

            elapsed, mem_mb = run_single(annotation_dir, outdir, no_uniprot, run_num)
            times.append(elapsed)
            if mem_mb:
                mems.append(mem_mb)

        mean_time   = statistics.mean(times)
        std_time    = statistics.stdev(times) if len(times) > 1 else 0
        mean_mem    = statistics.mean(mems) if mems else None
        std_mem     = statistics.stdev(mems) if len(mems) > 1 else 0

        result = {
            'label':      label,
            'mean_time':  round(mean_time, 2),
            'std_time':   round(std_time, 2),
            'mean_mem':   round(mean_mem, 1) if mean_mem else None,
            'std_mem':    round(std_mem, 1),
            'runs':       list(zip([round(t, 2) for t in times],
                                   [round(m, 1) for m in mems]))
        }
        all_results.append(result)

        print(f"    Mean: {mean_time:.1f}s ± {std_time:.1f}s | "
              f"{mean_mem:.1f} ± {std_mem:.1f} MB\n")

    # Write results file
    outfile = os.path.join(OUTBASE, 'benchmark_summary.txt')
    with open(outfile, 'w') as f:
        f.write('annoreport Benchmark Results\n')
        f.write('=' * 60 + '\n')
        f.write(f'Date:              {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n')
        f.write(f'Runs per config:   {N_RUNS}\n')
        f.write(f'Bins:              206 (Antarctic soil MAGs)\n')
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

    print(f"\nBenchmark complete.")
    print(f"Results written to: {outfile}")


if __name__ == '__main__':
    main()
