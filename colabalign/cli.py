'''Command-line argument parsing for the colabalign command.'''

from argparse import ArgumentParser, ArgumentTypeError
from pathlib import Path

from . import __version__
from .panel import DEFAULT_LOW_COVERAGE_FRACTION, MAX_LOW_COVERAGE_FRACTION, parse_panel_size


def _panel_size(text):
    '''argparse type for -n: "N" or "MIN-MAX".'''
    try:
        return parse_panel_size(text)
    except ValueError as error:
        raise ArgumentTypeError(str(error)) from None


def _fraction(text):
    '''argparse type for -l: a fraction from 0 to MAX_LOW_COVERAGE_FRACTION.'''
    try:
        value = float(text)
    except ValueError:
        raise ArgumentTypeError(f'{text!r} is not a number') from None
    if not 0 <= value <= MAX_LOW_COVERAGE_FRACTION:
        raise ArgumentTypeError(f'must be between 0 and {MAX_LOW_COVERAGE_FRACTION}')
    return value


def build_parser() -> ArgumentParser:
    '''Build the argument parser (separate from parsing, for testing and docs).'''
    parser = ArgumentParser(
        prog='colabalign',
        description=(
            'Generate a dendrogram based on pairwise-structure alignments '
            'and select representative structures.'
        )
    )
    parser.add_argument(
        '-V', '--version',
        action='version',
        version=f'%(prog)s {__version__}',
    )
    parser.add_argument(
        '-i', '--input',
        type=Path,
        required=True,
        nargs='+',
        metavar='PATH',
        help='Path(s) to input .pdb/.cif files, directories of them, or .zip/.tar(.gz) archives.'
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
        default=2,
        metavar='INTEGER',
        help='Number of CPU cores to use (Default = 2).'
    )
    parser.add_argument(
        '-n', '--panel-size',
        type=_panel_size,
        required=False,
        nargs='+',
        default=None,
        metavar='SIZE',
        help=(
            'Desired panel size(s): the number of structures to select. A single '
            'number gives the achievable panel closest to it; a range such as 40-60 '
            'gives every achievable panel size within it. The clustering threshold '
            'is chosen for you. Panel sizes depend on the tree, so the result can '
            'differ from the number requested. Can be combined with -t.'
        )
    )
    parser.add_argument(
        '-t', '--threshold',
        type=float,
        required=False,
        nargs='+',
        default=None,
        metavar='FLOAT',
        help=(
            'Set threshold(s) for clustering between 0.01 and 1.00 '
            '(Default = 0.20 when -n is not given). Multiple values can be provided.'
        )
    )
    parser.add_argument(
        '-l', '--low-coverage-fraction',
        type=_fraction,
        required=False,
        default=DEFAULT_LOW_COVERAGE_FRACTION,
        metavar='FLOAT',
        help=(
            'Fraction of each panel given to low-coverage structures: those whose '
            'best TM-score to the panel is lowest, added one at a time, e.g. 0.05 '
            f'(Default = {DEFAULT_LOW_COVERAGE_FRACTION:g}, no additions; '
            f'maximum {MAX_LOW_COVERAGE_FRACTION:g}).'
        )
    )
    parser.add_argument(
        '--no-superpose',
        action='store_true',
        help=(
            'Do not write superposed (or copied singleton) structure files to '
            'results/clusters/. Representatives and panels are unchanged; this '
            'only saves time and disk space.'
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

    return parser


def script_args(argv=None):
    '''Parse command-line arguments. argv defaults to sys.argv[1:].'''
    return build_parser().parse_args(argv)
