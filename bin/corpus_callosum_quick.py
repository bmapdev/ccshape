#!/usr/local/epd/bin/python

__author__ = "Brandon Ayers"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "ayersb@ucla.edu"

import argparse
import os


def main():
    parser = argparse.ArgumentParser(description="Quickly runs corpus_callosum_analyze \
                                                  from top and bottom callosal segementations formated \
                                                  as:\n subjectid_TOP.ucf \n subjectid_bot.ucf")
    parser.add_argument("source_directory", help="Directory containing properly formatted ucf files.\n")
    parser.add_argument("output_directory", help="Output directory for analyzed curves.\n")
    parser.add_argument("-test_level", dest="test_level", help="min/mid/max", default="mid", required=False)
    parser.add_argument("-template_id", dest="template_id", help="subjectID to use as template.",
                        default=False, required=False)
    args = parser.parse_args()
    corpus_callosum_quick(args.source_directory, args.output_directory, args.test_level, args.template_id)


def corpus_callosum_quick(source_directory, output_directory, test_level, template_id):
    """ Analyze properly named corpus callosum UCF segmentations by searching source directory for appropriate files, and
        automatically outputting setup information and analysis output to the output directory.

        File naming is assumed to be of the form:
        subjectid_top.ucf subjectid_bot.ucf
        If registration to a template curve is desired it should be specified by subject id.
            e.g. template_id = subject_j NOT template_id = subject_j_top.ucf
    """
    if test_level != "min" and test_level != "mid" and test_level != "max":
        raise ValueError("Error! test_lever must be 'min', 'mid' or 'max'")

    # Import models after function has been called to avoid unnecessary file imports on every call from the terminal.
    import corpus_callosum_analyze
    import corpus_callosum_setup_files
    source_directory = os.path.abspath(source_directory)
    output_directory = os.path.abspath(output_directory)
    paths = [os.path.join(output_directory, "setup_files"),
             os.path.join(output_directory, "output_files"),
             os.path.join(output_directory, "output_files", "elastic"),
             os.path.join(output_directory, "output_files", "linear"),
             os.path.join(output_directory, "output_files", "reg_elastic"),
             os.path.join(output_directory, "output_files", "reg_linear"),
             os.path.join(output_directory, "output_files", "uniform_reg_elastic"),
             os.path.join(output_directory, "output_files", "uniform_reg_linear")
             ]
    for p in paths:
        if 'reg' in p and test_level == 'min':
            break
        if not os.path.exists(p):
            os.mkdir(p)

    corpus_callosum_setup_files.create_setup_files(source_directory, os.path.join(output_directory, "setup_files"))
    subject_ids = os.path.join(output_directory, "setup_files", "subject_ids.txt")
    top_curves = os.path.join(output_directory, "setup_files", "top_curve_paths.txt")
    bot_curves = os.path.join(output_directory, "setup_files", "bottom_curve_paths.txt")

    if test_level == "mid":
        if not template_id:
            test_level = "min"
        else:
            corpus_callosum_analyze.corpus_callosum_analyze(subject_ids, top_curves, bot_curves,
                                                            os.path.join(output_directory, "output_files", "elastic"),
                                                            template_id=False, linear=False, linear_template_matching=False)
            corpus_callosum_analyze.corpus_callosum_analyze(subject_ids, top_curves, bot_curves,
                                                            os.path.join(output_directory, "output_files", "linear"),
                                                            template_id=False, linear=True, linear_template_matching=False)
            corpus_callosum_analyze.corpus_callosum_analyze(subject_ids, top_curves, bot_curves,
                                                            os.path.join(output_directory, "output_files", "reg_elastic"),
                                                            template_id, linear=False, linear_template_matching=False)

            corpus_callosum_analyze.corpus_callosum_analyze(subject_ids, top_curves, bot_curves,
                                                            os.path.join(output_directory, "output_files", "reg_linear"),
                                                            template_id, linear=True, linear_template_matching=False)
            corpus_callosum_analyze.corpus_callosum_analyze(subject_ids, top_curves, bot_curves,
                                                            os.path.join(output_directory, "output_files", "uniform_reg_elastic"),
                                                            template_id, linear=False, linear_template_matching=True)
            corpus_callosum_analyze.corpus_callosum_analyze(subject_ids, top_curves, bot_curves,
                                                            os.path.join(output_directory, "output_files", "uniform_reg_linear"),
                                                            template_id, linear=True, linear_template_matching=True)
    if test_level == "min":
            corpus_callosum_analyze.corpus_callosum_analyze(subject_ids, top_curves, bot_curves,
                                                            os.path.join(output_directory, "output_files", "elastic"),
                                                            template_id=False, linear=False, linear_template_matching=False)
            corpus_callosum_analyze.corpus_callosum_analyze(subject_ids, top_curves, bot_curves,
                                                            os.path.join(output_directory, "output_files", "linear"),
                                                            template_id=False, linear_template_matching=False)

if __name__ == "__main__":
    main()
