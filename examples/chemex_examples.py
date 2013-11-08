#!/usr/bin/env python
"""
Created on Mar 10, 2012

@author: Alex Hansen
         Guillaume Bouvignies (May 6, 2013)
"""

import os.path
import sys


def main():
    all_examples = set()
    for exp_type in os.listdir('./'):
        if os.path.isdir(exp_type):
            exp_list = os.listdir(exp_type)
            for exp in exp_list:
                all_examples.add(os.path.join(exp_type, exp))

    errors = 0
    error_list = set()
    for example in all_examples:
        exp_cfg = os.path.join(example, "experiments")
        par_cfg = os.path.join(example, "parameters")
        met_cfg = os.path.join(example, "methods")
        out_dir = os.path.join(example, "output")

        com_string = "chemex_fit.py -e {:s}/* -p {:s}/* -m {:s}/* -o {:s}".format(exp_cfg, par_cfg, met_cfg, out_dir)

        sys.argv = com_string.split()

        print ''.join(["#"] * 60)
        print "Testing {:s}".format(example).center(60)
        print ''.join(["#"] * 60)
        print
        print com_string

        try:
            os.system(com_string)
        except:
            print "Couldn't run:\n{:s}".format(sys.exc_info()[1])
            errors += 1
            error_list.add(example + " : " + sys.exc_info()[1])
            pass

    if errors:
        print "\n\nSome experiments had errors:\n\n{:s}".format("\n".join(error_list))
    else:
        print "\n\nAll systems a-ok!"


if __name__ == '__main__':
    main()
