"""
Created on Jul 1, 2011
Modified substantially Feb/March 2012 to support nicer help

@author: guillaume & Alex Hansen
"""

import argparse
import pkgutil
import sys

import chemex.experiments


def generate_experiment_subparser(prog, subparser, exp_pkg):
    """
    Generates a parsers with experimental arguments based on the directory tree

    Arguments
       prog       :   the program name
       subparser  :   the subparser you will be adding a parser to
       exp_pkg    :   the experimental class directory name

    Returns
       exp_parser :  the parser
    """

    class HelpAction(argparse.Action):
        """
        Special experiment-specific class for the parser arguments. Prints the
        nicely formatted help
        """

        def __call__(self, parser, namespace, values, option_string=None):

            try:
                exp = __import__(
                    '.'.join(['chemex', 'experiments', exp_pkg, values,
                              'exp_help']),
                    fromlist=['exp_help'],
                )
                dtp = __import__(
                    '.'.join(['chemex', 'experiments', exp_pkg, values,
                              'data_point']),
                    fromlist=['data_point'],
                )

            except ImportError:
                sys.stderr.write(
                    "\n ! {:s} is not a {:s} experiment\n\n"
                    .format(values, exp_pkg)
                )
                parser.print_help()
                exit(1)

            parse_line = exp.parse_line
            description = exp.description
            reference = exp.reference
            parameters = dtp.PAR_DICT

            format_experiment_help(parse_line, description, parameters)

            sys.exit(0)
            # ## end __call__

            # ## end HelpAction

    exp_class = __import__(
        ".".join(["chemex", "experiments", exp_pkg]),
        fromlist=[exp_pkg]
    )

    exp_help = __import__(
        ".".join(["chemex", "experiments", exp_pkg, "exp_help"]),
        fromlist=["exp_help"]
    )

    exp_parser = subparser.add_parser(
        exp_pkg,
        help=exp_help.parse_line,
        add_help=False,
        usage="{:s} -i {:s} [experiment]".format(prog, exp_pkg)
    )

    pkgs = [modname
            for _, modname, ispkg in pkgutil.iter_modules(exp_class.__path__)
            if ispkg]

    for pkg in sorted(pkgs):
        exp = __import__(
            ".".join(["chemex", "experiments", exp_pkg, pkg, "exp_help"]),
            fromlist=["exp_help"]
        )

        exp_parser.add_argument(pkg, action=HelpAction, help=exp.parse_line)

    return exp_parser


# ## end generate_experiment_subparser



# ## end format_experiment_help

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def arg_parse():
    description = (
    "ChemEx is an analysis program for chemical exchange detected by NMR. It "
    "is designed to take almost any kind of NMR data to aid the analysis, but "
    "the principle techniques are CPMG relaxation dispersion and Chemical "
    "Exchange Saturation Transfer.")

    parser = MyParser(description=description)

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
            fromlist=['exp_help'])

        parser_info_exp = subparsers_info.add_parser(
            type,
            help=type_help.parse_line,
            description="Enter an experiment to obtain more info about it.",
        )

        subparsers_info_type = parser_info_exp.add_subparsers(
            dest='experiments', )

        path_experiments = __import__(
            '.'.join(['chemex', 'experiments', type]),
            fromlist=[type]).__path__

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


# ######################### THE REAL PARSER IS HERE ##########################
# def arg_parse():
# with_info = set(['-i', '--info']) & set(sys.argv)
#
#     ######## PARSER TO GET HELP ########
#     if with_info:
#
#         args_to_parse = sys.argv[1:]
#
#         for flag in with_info:
#             args_to_parse.remove(flag)
#
#         parser = argparse.ArgumentParser(description="Classes of
# experiments", add_help=False,
#                                          usage="{:s} -i [experiment
# type]".format(sys.argv[0]))
#
#         subparser = parser.add_subparsers(help='Lists of fittable
# experiments')
#
#         exp_parsers = {}
#
#         exp_pkgs = [modname
#                     for _, modname, ispkg in pkgutil.iter_modules(
# chemex.experiments.__path__)
#                     if ispkg]
#
#         for exp_pkg in sorted(exp_pkgs):
#             exp_parsers[exp_pkg] = generate_experiment_subparser(
# parser.prog, subparser, exp_pkg)
#
#         if len(args_to_parse) == 0:
#             parser.print_help()
#
#         elif len(args_to_parse) == 1 and args_to_parse[0] in exp_parsers:
#             exp_parsers[args_to_parse[0]].print_help()
#
#         else:
#             parser.parse_args(args_to_parse)
#
#         exit(0)
#
#     ######## PARSER TO RUN THE CALCULATION ########
#     # Parses the command line
#
#     description = ("ChemEx is an analysis program for chemical exchange
# detected by NMR. It "
#                    "is designed to take almost any kind of NMR data to aid
#  the analysis, but "
#                    "the principle techniques are CPMG relaxation
# dispersion and Chemical "
#                    "Exchange Saturation Transfer.")
#
#     parser = argparse.ArgumentParser(description=description,
# prefix_chars='+-')
#
#     parser.add_argument('-e', dest='experiments', metavar='FILE',
# nargs='+', required=True,
#                         help='Input files containing experimental setup
# and data location')
#     parser.add_argument('-p', dest='parameters', metavar='FILE',
# required=True,
#                         help='Input file containing the fitting parameters')
#     parser.add_argument('-m', dest='method', metavar='FILE', help='Input
# file containing the fitting method')
#     parser.add_argument('-o', dest='out_dir', metavar='DIR',
# default='./output', help='Directory for output')
#     parser.add_argument('-i', '--info', action='store_true', help='List of
#  experiments')
#
#     group_residue_selec = parser.add_mutually_exclusive_group()
#     group_residue_selec.add_argument('+r', dest='res_incl', metavar='ID',
# nargs='+',
#                                      help='residue(s) to include in the fit')
#     group_residue_selec.add_argument('-r', dest='res_excl', metavar='ID',
# nargs='+',
#                                      help='residue(s) to exclude from the
# fit')
#
#     group_simulation = parser.add_mutually_exclusive_group()
#     group_simulation.add_argument('--mc', metavar='N', type=int, help='Run
#  N Monte-Carlo simulation')
#     group_simulation.add_argument('--bs', metavar='N', type=int, help='Run
#  N Bootstrap simulation')
#
#     parser.add_argument('--noplot', action='store_true', help='No plots of
#  the fits')
#
#     if len(sys.argv) == 1:
#         parser.print_help()
#         sys.exit(1)
#
#     args = parser.parse_args()
#
#     if args.res_incl:
#         args.res_incl = [res.lower() for res in args.res_incl]
#
#     if args.res_excl:
#         args.res_excl = [res.lower() for res in args.res_excl]
#
#     return args
