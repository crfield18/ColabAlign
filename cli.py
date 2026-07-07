'''Command-line argument parsing for the dendrogram generation script.'''

from argparse import ArgumentParser
from pathlib import Path


def script_args():
    '''Parse and return command-line arguments for the script.'''
    parser = ArgumentParser(
        description=(
            'Generate a dendrogram based on pairwise-structure alignments '
            'and select representative structures.'
        )
    )
    parser.add_argument(
        '-i', '--input',
        type=Path,
        required=True,
        nargs='+',
        metavar='PATH',
        help='Path(s) to input files.'
    )
    parser.add_argument(
        '-o', '--output',
        type=Path,
        required=True,
        metavar='PATH',
        help='Path to output directory. Will be created if it does not exist.'
    )
    parser.add_argument(
        '-c', '--cores',
        type=int,
        required=False,
        default=1,
        metavar='INTEGER',
        help='Number of CPU cores to use (Default = 1).'
    )
    parser.add_argument(
        '-t', '--threshold',
        type=float,
        required=False,
        nargs='+',
        default=[0.20],
        metavar='FLOAT',
        help=(
            'Set threshold(s) for clustering between 0.01 and 1.00 '
            '(Default = 0.20). Multiple values can be provided.'
        )
    )
    parser.add_argument(
        '-u', '--usalign',
        type=Path,
        required=False,
        default='USalign',
        metavar='PATH',
        help='Path to USalign executable. Not required if using conda install.'
    )
    parser.add_argument(
        '-b', '--beem',
        type=Path,
        required=False,
        default='BeEM',
        metavar='PATH',
        help=(
            'Path to BeEM executable. Not required if using conda install '
            'or working exclusively with .pdb files.'
        )
    )
    parser.add_argument(
        '-m', '--mode',
        type=str,
        required=False,
        default='all',
        choices=['all', 'first'],
        metavar='MODE',
        help='''Choose whether to align all chains from input structures or
        only the first chain (typically chain A). Valid arguments: %(choices)s'''
    )

    return parser.parse_args()
