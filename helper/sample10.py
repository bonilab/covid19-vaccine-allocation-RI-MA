#/usr/bin/env python3

import argparse
import pandas as pd
from pathlib import Path

def Export_Sample_10(csv_file, out_filename, index_col=None, random_state = 1):
    abspath = Path.resolve(Path(csv_file))
    if index_col is not None:
        f = pd.read_csv(abspath, index_col=index_col)
    else:
        f = pd.read_csv(abspath)
    f_10 = f.sample(n=10, axis=0, random_state=random_state)
    f_10.to_csv(Path(abspath.parent, out_filename))

if __name__ == "__main__":
    parser = argparse.ArgumentParser('Sample 10 rows from a csv file')
    parser.add_argument(
        '-i', '--input',
        default='input.csv',
        type=str,
        help='The csv where samples will be drawn from. Default: input.csv'
    )
    parser.add_argument(
        '-o', '--output',
        default='output.csv',
        type=str,
        help='The output csv file. Default: output.csv'
    )
    parser.add_argument(
        '-c', '--index_col',
        default=0,
        type=int,
        help='Number to indicate the index column. Default: 0'
    )
    parser.add_argument(
        '-r', '--random_state',
        default=1,
        type=int,
        help='The random state being used in selecting the 10 samples. Default: 1'
    )

    args = parser.parse_args()

    Export_Sample_10(args.input, args.output, args.index_col, args.random_state)

    print('--done--')