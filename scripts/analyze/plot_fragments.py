# -*- coding: utf-8 -*-
# @Author: Jaume Bonet
# @Date:   2016-03-15 15:06:45
# @Last Modified by:   Jaume Bonet
# @Last Modified time: 2016-03-16 18:49:18

import os
import sys
import argparse

scrdir  = os.path.dirname(os.path.realpath(__file__))
rotools = os.path.join(scrdir, "../..")
sys.path.append(rotools)
from rotools.fragments import FragSet


def get_options(*args, **kwds):

    parser = argparse.ArgumentParser(description="Plot Rosetta's fragments as an image")

    parser.add_argument('-in', dest='infile', type=str, action='store',
                        help='Fragment file name')
    parser.add_argument('-show', dest='show', action='store_true', default=False,
                        help='Show image with visualizer (def:false)')
    parser.add_argument('-out:prefix', dest='outpref', type=str, action='store',
                        help='Image output file prefix (def:input name_img)', default=None)
    parser.add_argument('-out:type', dest='outtype', type=str, action='store',
                        help='Image output file type (png/svg) (def:png)', default="png")

    options = parser.parse_args()
    if options.outpref is None:
        options.outpref = os.path.split(options.infile)[-1] + "_img"

    return options

if __name__ == "__main__":
    options = get_options()

    fset = FragSet.parse(options.infile)
    fset.plot(fileprefix= options.outpref, format=options.outtype, show=options.show )
