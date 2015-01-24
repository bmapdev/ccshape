#! /usr/local/epd/bin/python

__author__ = "Brandon Ayers"
__copyright__ = "Copyright 2013, Brandon Ayers, Shantanu H. Joshi \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "ayersb@ucla.edu"

import os
import sys
import subprocess


def main():
    print "\nStarting setup...\n"
    print "Testing rpy2..."
    try:
        import rpy2.robjects
        print "Success! R and rpy2 are installed.\n"
    except:
        print "rpy2 test failure! R or rpy2 are not installed correctly!\n"
        print "Stats portion will not function until this is fixed.\n"

    print "Testing corpus callosum thickness package...\n"
    test_dir_location = os.path.dirname(os.path.realpath(__file__))
    print test_dir_location
    os.chdir(test_dir_location)
    names = open("names.txt", "w")
    top_curves = open("topCurvePaths.txt", "w")
    bot_curves = open("botCurvePaths.txt", "w")
    curve_counter = 1
    for curve in os.listdir(os.getcwd()+"/dataSample/"):
        if 'top' in curve.lower():
            names.write('curve'+str(curve_counter)+'\n')
            curve_counter += 1
            top_curves.write(os.getcwd()+"/dataSample/"+curve+'\n')
        if 'bot' in curve.lower():
            bot_curves.write(os.getcwd()+"/dataSample/"+curve+'\n')
    names.close()
    top_curves.close()
    bot_curves.close()
    bindir = os.path.dirname(sys.executable)
    print "Test data will output to a specified directory"
    print "Blank input will output to the user's desktop"
    outdir = raw_input("Specify test output directory: ")
    if outdir == '':
        outdir = "~/Desktop"
    outdir = os.path.expanduser(outdir) + "/testOutput"
    print "Writing files to ", outdir

    # Test unregistered Matching
    c1 = subprocess.call([sys.executable, bindir + "/corpus_callosum_analyze.py", os.getcwd() + "/names.txt",
                    os.getcwd() + "/topCurvePaths.txt",
                    os.getcwd() + "/botCurvePaths.txt",
                    "-odir", outdir])

    # Test matching registered to Template
    subprocess.call([sys.executable, bindir + "/corpus_callosum_analyze.py", os.getcwd() + "/names.txt",
                    os.getcwd() + "/topCurvePaths.txt",
                    os.getcwd() + "/botCurvePaths.txt",
                    "-odir", outdir, "-templateID", "curve1"])

if __name__ == "__main__":
    main()