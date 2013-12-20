#!/usr/local/epd/bin/python


__author__ = "Brandon Ayers"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "ayersb@ucla.edu"

from math import pi
import os
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as LA
from curvematch.match import match_curve_pair
from curvematch.curve import Curve
from curvematch.qshape import QShape
from curvematch import geodesics



class CCThickness():
    def __init__(self, subject_name, curvefile_path_top, curvefile_path_bottom, resample_siz=100, geodesic_steps=10):
        self.settings = geodesics.Geodesic()
        self.settings.steps = geodesic_steps
        self.settings.closed = False

        self.subject_name = subject_name
        self.return_shapes = True
        self.resample_siz = resample_siz
        self.rotation = False

        self.curvefile_path_top = os.path.abspath(curvefile_path_top)
        self.curvefile_path_bottom = os.path.abspath(curvefile_path_bottom)
        self.gamma = []
        self.curve_top = []
        self.curve_bot = []
        self.curve_bot_gamma_adjusted = []
        self.naive_thickness = []
        self.thickness = []
        self.main_method()

    def main_method(self):
        self.get_matched_curve_pair_data()
        self.curve_bot_gamma_adjusted = self.get_curve_of_gamma(self.curve_bot)
        self.compute_thickness()

    def get_matched_curve_pair_data(self):
        curve_pair_data = match_curve_pair(self.curvefile_path_top, self.curvefile_path_bottom,
                                           self.settings, self.rotation,
                                           self.resample_siz, self.return_shapes)
        self.curve_top = curve_pair_data[1]
        self.curve_bot = curve_pair_data[2]
        self.gamma = curve_pair_data[0].gamma

        # WARNING! This assumes slice coordinate is in first pos of input array
        self.curve_top.coords = self.curve_top.coords[1:2+1]
        self.curve_top.dim -= 1
        self.curve_bot.coords = self.curve_bot.coords[1:2+1]
        self.curve_bot.dim -= 1

    def get_curve_of_gamma(self, curve1):
        ##Please check that this is correct!
        newcoords = np.zeros((curve1.dim, curve1.siz))
        for i in range(0, curve1.dim):
            newcoords[i, :] = np.interp(self.gamma, np.linspace(0, 2*pi, curve1.siz), curve1.coords[i, :])
        newcurve = Curve()
        newcurve.coords = newcoords
        newcurve.dim = curve1.dim
        newcurve.siz = curve1.siz
        return newcurve

    def compute_thickness(self, include_naive=True):
        #Just finds the euclidean distance between points
        if include_naive:
            self.naive_thickness = np.sqrt(np.sum((self.curve_top.coords - self.curve_bot.coords)**2, axis=0))
        self.thickness.append(np.sqrt(np.sum((self.curve_top.coords - self.curve_bot_gamma_adjusted.coords)**2, axis=0)))
        print "Naive: ", np.average(self.naive_thickness), '\n'
        print "Gamma: ", np.average(self.thickness), '\n'


def add_thickness_plot_given_curves(curve1, curve2):
        plt.plot(curve1.coords[0], curve1.coords[1])
        plt.plot(curve2.coords[0], curve2.coords[1])
        for i in xrange(0,len(curve1.coords[0])):
            plt.plot([curve1.coords[0][i], curve2.coords[0][i]],
                     [curve1.coords[1][i], curve2.coords[1][i]])


def plot_thicknesses(CCThickness_object_list, individual_plots=False, include_naive_plot=True):
    figure_count = 1
    for subject in CCThickness_object_list:
        plt.figure(figure_count)
        if include_naive_plot == True:
            plt.title(subject.subject_name)
            plt.subplot(211)
            plt.subplots_adjust(hspace = 0.4)
            plt.title(subject.subject_name)
            plt.xlabel("Naive Plot")
            add_thickness_plot_given_curves(subject.curve_top, subject.curve_bot)
            plt.subplot(212)
            plt.xlabel("Intelligent Plot")
            add_thickness_plot_given_curves(subject.curve_top, subject.curve_bot_gamma_adjusted)
        else:
            plt.title(subject.subject_name)
            add_thickness_plot_given_curves(subject.curve_top, subject.curve_bot_gamma_adjusted)
        plt.savefig("test"+str(figure_count)+".pdf")
        figure_count += 1

def analyze_thicknesses(ucf_curve_paths):
    if len(ucf_curve_paths) % 2 != 0:
        raise ValueError("Uneven Number of Segmentations! Can not process")
    subjects = []
    for pair_index in xrange(1, len(ucf_curve_paths), 2):
        current_subject = CCThickness(pair_index/2+1, ucf_curve_paths[pair_index],ucf_curve_paths[pair_index+1])
        subjects.append(current_subject)
    plot_thicknesses(subjects)


a = CCThickness("Subject 1", "002_S_0295_TOP.ucf", "002_S_0295_BOT.ucf")
b=  CCThickness("Subject 2", "002_S_0413_TOP.ucf", "002_S_0413_BOT.ucf")

plot_thicknesses([a,b])