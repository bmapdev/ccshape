#!/usr/local/epd/bin/python


__author__ = "Brandon Ayers"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "ayersb@ucla.edu"

import os
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as LA
from curvematch.match import match_curve_pair
from curvematch.curve import Curve
from curvematch import geodesics


class CCThickness():
    def __init__(self, curvefile_path_top, curvefile_path_bottom, resample_siz=150, geodesic_steps=10, rotation=False):
        self.settings = geodesics.Geodesic()
        self.settings.steps = geodesic_steps
        self.settings.closed = False
        self.return_shapes = True

        self.resample_siz = resample_siz
        self.rotation = False
        self.curvefile_path_top = os.path.abspath(curvefile_path_top)
        self.curvefile_path_bottom = os.path.abspath(curvefile_path_bottom)
        self.curve_top = []
        self.curve_bot = []
        self.matched_curve_pair_geodesic = []
        self.compute_thickness()

    def get_matched_curve_pair_data(self):
        curve_pair_data = match_curve_pair(self.curvefile_path_top, self.curvefile_path_bottom,
                                           self.settings, self.rotation,
                                           self.resample_siz, self.return_shapes)


        self.curve_top = curve_pair_data[1]
        self.curve_bot = curve_pair_data[2]
        self.matched_curve_pair_geodesic = curve_pair_data[0]

        # WARNING! This assumes slice coordinate is in first pos of input arrays
        self.curve_top.coords = self.curve_top.coords[1:2+1]
        self.curve_bot.coords = self.curve_bot.coords[1:2+1]

    def plot_thickness(self, curve1, curve2, labels=[]):
        plt.plot(curve1.coords[0], curve1.coords[1])
        plt.plot(curve2.coords[0], curve2.coords[1])
        for i in xrange(len(curve1.coords[1])):
            plt.plot([curve1.coords[0][i], curve2.coords[0][i]],
                     [curve1.coords[1][i], curve2.coords[1][i]])


    def get_curve_of_gamma(self, curve1, geodesic_with_gamma):
        gamma_coords = np.zeros(curve1.coords.shape)
        gamma = geodesic_with_gamma.gamma
        total_length_curve = curve1.length()
        for counter in xrange(len(gamma)):
            cur_curve[0] = Curve(curve1.coords[0][:counter+1])
            cur_curve[1] = Curve(curve1.coords[1][:counter+1])


    def compute_thickness(self):
        self.get_matched_curve_pair_data()
        #print self.curve_top.coords[0]
        #print self.curve_top.coords[1] - self.curve_bot.coords[1],'\n'
        #print self.curve_top.coords[2] - self.curve_bot.coords[2],'\n'
        plt.subplot(211)
        self.plot_thickness(self.curve_top,self.curve_bot)
        curve_bot_gamma = Curve()
        curve_bot_gamma.coords = (self.curve_bot.coords - self.matched_curve_pair_geodesic.gamma)
        plt.subplot(212)
        #print self.matched_curve_pair_geodesic.gamma
        self.plot_thickness(self.curve_top, curve_bot_gamma)
        plt.show()



a = CCThickness("002_S_0295_BOT.ucf","002_S_0295_TOP.ucf")