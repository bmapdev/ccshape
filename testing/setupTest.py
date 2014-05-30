#! /usr/local/epd/bin/python

__author__ = "Brandon Ayers"
__copyright__ = "Copyright 2013, Brandon Ayers, Shantanu H. Joshi \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "ayersb@ucla.edu"

import os
import sys

def main():
    print "'\nStarting setup...'\n"
    names = open("names.txt", "w")
    top_curves = open("topCurvePaths.txt", "w")
    bot_curves = open("botCurvePaths.txt", "w")
    for curve in os.listdir(os.getcwd()+"/dataSample/"):
        if 'top' in curve:
            names.write(curve.split('_top.ucf')[0]+'\n')
            top_curves.write(os.getcwd()+"/dataSample/"+curve+'\n')
        if 'bot' in curve:
            bot_curves.write(os.getcwd()+"/dataSample/"+curve+'\n')
    names.close()
    top_curves.close()
    bot_curves.close()
    print "Setup Completed! To run the test, run the following command,\n"
    bindir = os.path.dirname(sys.executable)
    print sys.executable, " ", bindir+ "/corpus_callosum_analyze.py", " ", os.getcwd() + \
          "/names.txt", os.getcwd() + "/topCurvePaths.txt", os.getcwd() + "/botCurvePaths.txt" +" "+ "-odir " +os.getcwd()+"/testOutput"

if __name__ == "__main__":
    main()