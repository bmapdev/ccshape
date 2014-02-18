#!/usr/local/epd/bin/python

__author__ = "Brandon Ayers"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "ayersb@ucla.edu"

import argparse
import curvematch
import numpy as np
from ccshape.corpus_callosum import CorpusCallosum
import os


def main():
    parser = argparse.ArgumentParser(description='Compute thickness of the Corpus Callosum using top and bottom curves.\n')
    parser.add_argument('-subjname', dest='subjname', help='subject name', required=True)
    parser.add_argument('-top', dest='top_ucf', help='top ucf curve', required=True)
    parser.add_argument('-bottom', dest='bottom_ucf', help='bottom ucf curve', required=True)
    parser.add_argument('-odir', dest='odir', help='output directory', required=True)
    args = parser.parse_args()
    corpus_callosum_thickness(args.subjname, args.top_ucf, args.bottom_ucf, args.odir)


def corpus_callosum_thickness(subjname, top_ucf, bottom_ucf, odir):

    callosal_curve = CorpusCallosum(os.path.join(odir, subjname), top_ucf, bottom_ucf)
    callosal_curve.compute_thickness()
    callosal_curve.output_thickness_ucf()
    callosal_curve.plot_thicknesses(False)

if __name__ == '__main__':
    main()
