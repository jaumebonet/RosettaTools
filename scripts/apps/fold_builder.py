# -*- coding: utf-8 -*-
# @Author: Jaume Bonet
# @Date:   2016-03-17 14:41:38
# @Last Modified by:   Jaume Bonet
# @Last Modified time: 2016-03-23 15:23:46

import os
import sys
import argparse

# Including rotools without adding them to the PYTHONPATH
scrdir  = os.path.dirname(os.path.realpath(__file__))
rotools = os.path.join(scrdir, "../..")
sys.path.append(rotools)

from rotools.forms import FormFabric


def get_options(*args, **kwds):

    parser = argparse.ArgumentParser(description="Plot Rosetta's fragments as an image")

    parser.add_argument('-in:form', dest='inform', type=str, action='store',
                        help='JSON Form definition')
    parser.add_argument('-in:id',   dest='inid', type=str, action='store',
                        help='Form identifier', default="form")

    options = parser.parse_args()

    return options


if __name__ == "__main__":
    options = get_options()

    fabric  = FormFabric()
    fabric.build(options.inid, options.inform)
    fabric.print_structures()
    fabric.dump()
