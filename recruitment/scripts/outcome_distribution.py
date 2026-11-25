"""
Dataset-wide outcome distribution over the full recruitment classification.

Reads every per-particle recruitment file and tallies the outcome states,
named as the files name them (flag_meanings), so that the script serves
both the v2.0.0 and the v3.0.0 files. Writes a CSV and prints a summary suitable for the Validation table
in recruitment/README.md and for methodology.md section 7.

Usage
-----
    python outcome_distribution.py [--data DIR] [--out CSV]

Cohorts released on 29 February are included in the totals: they are part of
the classified dataset, and are excluded only from the season-day
climatologies in the figure repository.
"""

import argparse
import glob
import os

import numpy as np
import pandas as pd
import xarray as xr


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data', default='data',
                   help='directory of per-cohort recruitment files')
    p.add_argument('--out', default='outcome_distribution.csv')
    p.add_argument('--per-cohort', default=None,
                   help='optional CSV of per-cohort counts')
    args = p.parse_args()

    files = sorted(glob.glob(os.path.join(args.data, '*.nc')))
    if not files:
        raise SystemExit(f'No files found in {args.data}')
    print(f'{len(files)} cohorts', flush=True)
    with xr.open_dataset(files[0]) as d:
        LABELS = d.outcome.attrs['flag_meanings'].split()

    counts = np.zeros(len(LABELS), dtype=np.int64)
    rows = []

    for k, path in enumerate(files, 1):
        with xr.open_dataset(path) as d:
            if d.outcome.attrs['flag_meanings'].split() != LABELS:
                raise SystemExit(f'{path}: outcome labels differ from those '
                                 f'of {files[0]}')
            c = np.bincount(d.outcome.values.astype(int), minlength=len(LABELS))
        counts += c
        rows.append({'release_date': os.path.basename(path)[:-3].replace('_', '-'),
                     'n_particles': int(c.sum()),
                     **{label: int(n) for label, n in zip(LABELS, c)}})
        if k % 250 == 0:
            print(f'  {k}/{len(files)}', flush=True)

    total = counts.sum()

    df = pd.DataFrame({'outcome': LABELS,
                       'count': counts,
                       'fraction': counts / total})
    df.to_csv(args.out, index=False)

    if args.per_cohort:
        pd.DataFrame(rows).to_csv(args.per_cohort, index=False)
        print(f'per-cohort counts: {args.per_cohort}')

    print(f'\n{len(files)} cohorts, {total:,} particles\n')
    for label, n in zip(LABELS, counts):
        print(f'  {label:26s} {n:14,}  {n / total:6.2%}')

    censored = counts[LABELS.index('censored')]
    m6 = counts[[k for k, label in enumerate(LABELS) if 'M6' in label][0]]
    print(f'\n  M6 including censored    {(censored + m6) / total:15.2%}')
    print(f'\nWritten: {args.out}')


if __name__ == '__main__':
    main()
