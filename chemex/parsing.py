import argparse
import pkgutil
import sys
import re

import chemex.experiments
import chemex.version


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def arg_parse():
    description = (
        "ChemEx is an analysis program for chemical exchange detected by "
        "NMR. It is designed to take almost any kind of NMR data to aid the "
        "analysis, but the principle techniques are CPMG relaxation "
        "dispersion and Chemical Exchange Saturation Transfer."
    )

    parser = MyParser(
        description=description,
        prog='chemex',
        version='ChemEx version {}'.format(chemex.version.__version__)
    )

    subparsers = parser.add_subparsers(dest='commands', )

    # Parser Info
    parser_info = subparsers.add_parser(
        "info",
        help="Shows classes of experiments that can be fit",
        description="Enter a class of experiments.",
    )

    subparsers_info = parser_info.add_subparsers(dest='types')
    types = [name for _, name, ispkg in
             pkgutil.iter_modules(chemex.experiments.__path__) if ispkg]

    for type in types:

        type_help = __import__(
            '.'.join(['chemex', 'experiments', type, 'exp_help']),
            fromlist=['exp_help']
        )

        parser_info_exp = subparsers_info.add_parser(
            type,
            help=type_help.parse_line,
            description="Enter an experiment to obtain more info about it.",
        )

        subparsers_info_type = parser_info_exp.add_subparsers(
            dest='experiments',
        )

        path_experiments = __import__(
            '.'.join(['chemex', 'experiments', type]),
            fromlist=[type],
        ).__path__

        experiments = [name for _, name, ispkg in
                       pkgutil.iter_modules(path_experiments) if ispkg]

        for experiment in experiments:
            experiment_help = __import__(
                '.'.join(
                    ['chemex', 'experiments', type, experiment, 'exp_help']),
                fromlist=['exp_help']
            )

            subparsers_info_type.add_parser(
                '_'.join([experiment, type]),
                help=experiment_help.parse_line,
                add_help=False,
            )

    # Parser fit
    parser_fit = subparsers.add_parser(
        "fit",
        help="Starts a fit",
        prefix_chars='+-'
    )

    parser_fit.add_argument(
        '-e',
        dest='experiments',
        metavar='FILE',
        nargs='+',
        required=True,
        help='Input files containing experimental setup and data location'
    )

    parser_fit.add_argument(
        '-p',
        dest='parameters',
        metavar='FILE',
        required=True,
        help='Input file containing the fitting parameters'
    )

    parser_fit.add_argument(
        '-m',
        dest='method',
        metavar='FILE',
        help='Input file containing the fitting method'
    )

    parser_fit.add_argument(
        '-o',
        dest='out_dir',
        metavar='DIR',
        default='./output',
        help='Directory for output'
    )

    parser_fit.add_argument(
        '-i',
        '--info',
        action='store_true',
        help='List of experiments'
    )

    parser_fit.add_argument(
        '--noplot',
        action='store_true',
        help='No plots of the fits'
    )

    group_residue_selec = parser_fit.add_mutually_exclusive_group()

    group_residue_selec.add_argument(
        '+r',
        dest='res_incl',
        metavar='ID',
        nargs='+',
        help='residue(s) to include in the fit'
    )

    group_residue_selec.add_argument(
        '-r',
        dest='res_excl',
        metavar='ID',
        nargs='+',
        help='residue(s) to exclude from the fit'
    )

    group_simulation = parser_fit.add_mutually_exclusive_group()

    group_simulation.add_argument(
        '--mc',
        metavar='N',
        type=int,
        help='Run N Monte-Carlo simulation'
    )

    group_simulation.add_argument(
        '--bs',
        metavar='N',
        type=int,
        help='Run N Bootstrap simulation'
    )

    args = parser.parse_args()

    if args.commands == 'fit':
        if args.res_incl:
            args.res_incl = [res.lower() for res in args.res_incl]
        if args.res_excl:
            args.res_excl = [res.lower() for res in args.res_excl]

    return args


# Functions to parse Sparky-like assignment
# Code has been adapted from Sparky code

def parse_assignment(assignment):
    """
    Parse assignment of form g1a1-g2a2 to get ((g1, a1), (g2, a2))
    Or g1a1-a2 to get ((g1, a1), (g1, a2))
    A '?' component is translated to ('', '')
    """

    res = assignment.lower().split('-')
    assignment = []
    last_group = None
    for s in res:
        ga = split_group_atom(s)
        if ga == None:
            if last_group:
                ga = (last_group, s)
            else:
                return None
        assignment.append((parse_group_name(ga[0]) + ga[1:]))
        last_group = ga[0]

    return tuple(assignment)


def split_group_atom(group_atom):
    """
    Skip to first digit, then skip to first H|C|N to find start of atom name.
    """

    s = re.search('[0-9]', group_atom)
    if s:
        first_digit = s.start()
        s = re.search('[hHcCnNqQmM]', group_atom[first_digit:])
        if s:
            HCNQM_offset = s.start()
            d = first_digit + HCNQM_offset
            return (group_atom[:d], group_atom[d:])
    if group_atom == '?':
        return ('', '')
    return None


def parse_group_name(g):
    s = re.search('[0-9]+', g)
    if s:
        return int(s.group()), g[:s.start()]
    return g, None


