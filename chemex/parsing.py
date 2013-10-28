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
                exp = __import__('.'.join(['chemex', 'experiments',
                    exp_pkg, values, 'exp_help']), fromlist=['exp_help'])
                dtp = __import__('.'.join(['chemex', 'experiments',
                    exp_pkg, values, 'data_point']), fromlist=['data_point'])

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

            format_experiment_help(parse_line, description, reference, parameters)

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


def format_experiment_help(pline='unknown experiment',
                           desc='unknown experiment',
                           ref={'journal': '', 'year': 1900, 'volume': 0, 'pages': ''},
                           par={'spec': [], 'float': [], 'fix': []}):

    sys.stdout.write("""
    ---------------------------------------------
    {:s}
    ---------------------------------------------

    {:s}

    {:s} ({:d}) v.{:d}, p.{:s}

    Spectrometer parameters
    -----------------------
    """.format(pline, desc, ref['journal'], ref['year'],
                                        ref['volume'], ref['pages']))

    for p in par['exp']:
        sys.stdout.write('    {:s}\n'.format(p))

    sys.stdout.write("""
    Fitted parameters
    -----------------\n""")

    for p in par['fit']:
        sys.stdout.write('    {:s}\n'.format(p))

    if len(par['fix']):
        sys.stdout.write("""
    Fixed parameters
    ----------------\n""")

        for p in par['fix']:
            sys.stdout.write('    {:s}\n'.format(p))

    sys.stdout.write("\n")
# ## end format_experiment_help

########################## THE REAL PARSER IS HERE ##########################
def arg_parse():

    with_info = set(['-i', '--info']) & set(sys.argv)

    ######## PARSER TO GET HELP ########
    if with_info:

        args_to_parse = sys.argv[1:]

        for flag in with_info:
            args_to_parse.remove(flag)

        parser = argparse.ArgumentParser(description="Classes of experiments", add_help=False,
                                         usage="{:s} -i [experiment type]".format(sys.argv[0]))

        subparser = parser.add_subparsers(help='Lists of fittable experiments')

        exp_parsers = {}

        exp_pkgs = [modname
                    for _, modname, ispkg in pkgutil.iter_modules(chemex.experiments.__path__)
                    if ispkg]

        for exp_pkg in sorted(exp_pkgs):
            exp_parsers[exp_pkg] = generate_experiment_subparser(parser.prog, subparser, exp_pkg)

        if len(args_to_parse) == 0:
            parser.print_help()

        elif len(args_to_parse) == 1 and args_to_parse[0] in exp_parsers:
            exp_parsers[args_to_parse[0]].print_help()

        else:
            parser.parse_args(args_to_parse)

        exit(0)

    ######## PARSER TO RUN THE CALCULATION ########
    # Parses the command line
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
     description="""
     ChemEx is an analysis program for chemical exchange detected by NMR. It
     is designed to take almost any kind of NMR data to aid the analysis, but
     the principle techniques are CPMG relaxation dispersion and Chemical
     Exchange Saturation Transfer.
     """, prefix_chars='+-')

    parser.add_argument('-e', '--experiments', action='store', metavar='FILE',
                        nargs='+', help='Input files containing experimental setup and data')

    parser.add_argument('-p', '--parameters', action='store', metavar='FILE',
                        help='Input file containing the fitting parameters.')

    parser.add_argument('-m', '--method', action='store', metavar='FILE',
                        help='Input file containing the fitting method.')

    parser.add_argument('+r', '--res_incl', action='store', metavar='RESID',
                        nargs='+', help='residue(s) to include in the fit')

    parser.add_argument('-r', '--res_excl', action='store', metavar='RESID',
                        nargs='+', help='residue(s) to exclude from the fit')

    parser.add_argument('-o', '--out_dir', action='store', metavar='DIR',
                        help='Directory for output')

    parser.add_argument('-i', '--info', action='store_true', help='List of experiments')

    parser.add_argument('--noplot', action='store_true', help='No plots of the fits')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args()
