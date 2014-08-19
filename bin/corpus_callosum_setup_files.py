#!/usr/local/epd/bin/python

__author__ = "Brandon Ayers"
__copyright__ = "Copyright 2013, Brandon Ayers, Shantanu H. Joshi \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "ayersb@ucla.edu"

"""
Quickly creates text files for running corpus_callosum_anaylze
from top and bottom callosal segementations formated as:
subjectid_TOP.ucf
subjectid_bot.ucf
"""

import argparse
import os


def main():
    parser = argparse.ArgumentParser(description="Quickly creates text files for running corpus_callosum_analyze \
                                                  from top and bottom callosal segementations formated \
                                                  as:\n subjectid_TOP.ucf \n subjectid_bot.ucf")
    parser.add_argument('source_directory', help='Directory containing properly formatted ucf files.\n')
    parser.add_argument('output_directory', help='Output directory for text files.')

    args = parser.parse_args()
    create_setup_files(args.source_directory, args.output_directory)


def create_setup_files(source_directory, output_directory):
    source_directory = os.path.abspath(source_directory)
    output_directory = os.path.abspath(output_directory)
    source_list = os.listdir(source_directory)
    top_curves = []
    bot_curves = []
    subject_ids = []
    for curve in source_list:
        if '_top' in curve.lower():
            top_curves.append(curve)
            subject_ids.append(curve[:curve.lower().find('_top')])
        elif '_bot' in curve.lower():
            bot_curves.append(curve)
    if len(top_curves) > len(bot_curves):
        raise ValueError("Error! There are more top curves than bottom curves!")
    elif len(top_curves) < len(bot_curves):
        raise ValueError("Error! There are more bottom curves than top curves!")
    top_curve_paths_file = open(os.path.join(output_directory, "top_curve_paths.txt"), 'w')
    bottom_curve_paths_file = open(os.path.join(output_directory, "bottom_curve_paths.txt"), "w")
    subject_ids_file = open(os.path.join(output_directory, "subjectIDs.txt"), 'w')

    for i in xrange(len(top_curves)):
        top_curve_paths_file.write(os.path.join(source_directory, top_curves[i])+'\n')
        bottom_curve_paths_file.write(os.path.join(source_directory, bot_curves[i])+'\n')
        subject_ids_file.write(subject_ids[i]+'\n')
    print "Setup Files Created Successfully.\n"

if __name__ == '__main__':
    main()