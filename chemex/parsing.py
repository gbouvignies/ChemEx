import importlib
import argparse
import pkgutil
import sys

from chemex import experiments
from chemex import version


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
        version='ChemEx version {}'.format(version.__version__)
    )

    subparsers = parser.add_subparsers(dest='commands', )

    parser_info = subparsers.add_parser(
        "info",
        help="Shows classes of experiments that can be fit",
        description="Enter a class of experiments.",
    )

    subparsers_info = parser_info.add_subparsers(dest='types')

    exp_types = [
        name
        for _, name, ispkg in pkgutil.iter_modules(experiments.__path__)
        if ispkg
    ]

    for exp_type in exp_types:

        package_name_type = '.'.join([experiments.__name__, exp_type])

        package_type = importlib.import_module(package_name_type)

        description_type = package_type.__doc__.split('\n')[0]

        parser_info_type = subparsers_info.add_parser(
            exp_type,
            help=description_type,
            description="Enter an experiment to obtain more info about it.",
        )

        subparsers_info_type = parser_info_type.add_subparsers(
            dest='experiments',
        )

        exps = [
            name
            for _, name, ispkg in pkgutil.iter_modules(package_type.__path__)
        ]

        for exp in exps:

            package_name_exp = '.'.join([package_name_type, exp])

            package_exp = importlib.import_module(package_name_exp)

            if hasattr(package_exp, 'Profile'):
                subparsers_info_type.add_parser(
                    exp,
                    help=package_exp.__doc__.split('\n')[0],
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
            args.res_incl = set([res.lower() for res in args.res_incl])
        if args.res_excl:
            args.res_excl = set([res.lower() for res in args.res_excl])

    return args
