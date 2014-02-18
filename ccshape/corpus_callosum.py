#!/usr/local/epd/bin/python

__author__ = "Brandon Ayers"
__copyright__ = "Copyright 2013, Brandon Ayers, Shantanu H. Joshi \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "ayersb@ucla.edu"

import numpy as np
from curvematch import match
from curvematch import geodesics
from curvematch.curve import Curve
import os
from shapeio import curveio
import matplotlib.pyplot as plt


class CorpusCallosum:

    def __init__(self, subject_name, curvefile_path_top, curvefile_path_bottom, resample_siz=500,
                 geodesic_steps=7, elastic=True):

        self.settings = geodesics.Geodesic()
        self.settings.steps = geodesic_steps
        self.settings.closed = False
        self.elastic = elastic

        self.subject_name = subject_name
        self.return_shapes = True
        self.resample_siz = resample_siz
        self.rotation = False

        self.curvefile_path_top = os.path.abspath(curvefile_path_top)
        self.curvefile_path_bottom = os.path.abspath(curvefile_path_bottom)

        self.gamma = []
        self.curve_top = []
        self.curve_bot = []
        self.medial_curve = []
        self.plane_dim_data = []

        self.curve_bot_gamma_adjusted = []
        self.nonelastic_thickness = []
        self.elastic_thickness = []

        self.joined_elastic_coords = []
        self.joined_nonelastic_coords = []

    def match_top_and_bottom_curves(self):
        curve_pair_data = match.match_curve_pair(self.curvefile_path_top, self.curvefile_path_bottom,
                                           self.settings, self.rotation, self.resample_siz, self.return_shapes)
        self.gamma = curve_pair_data[0].gamma
        self.curve_top = curve_pair_data[1]
        self.curve_bot = curve_pair_data[2]
        diff_between_dims = np.abs(np.sum(self.curve_top.coords - self.curve_bot.coords, 1))
        if np.min(diff_between_dims) > 1:
            print ("Warning! These segmentations don't seem to be from the same plane.")
        least_var_dim = np.argmin(diff_between_dims)
        coords_top, coords_bot = [], []
        for i in xrange(self.curve_top.dim):
            if i == least_var_dim:
                self.plane_dim_data = [i, self.curve_top.coords[i], self.curve_bot.coords[i]]
            else:
                coords_top.append(self.curve_top.coords[i])
                coords_bot.append(self.curve_bot.coords[i])
        self.curve_top.dim -= 1
        self.curve_bot.dim -= 1
        self.curve_top.coords = np.array(coords_top)
        self.curve_bot.coords = np.array(coords_bot)

    def compute_thickness(self):
        self.match_top_and_bottom_curves()
        self.curve_bot_gamma_adjusted = self.curve_bot.return_reparameterized_by_gamma(self.gamma)
        #Compute nonelastic thickness by finding the euclidean distance between uniformly sampled points
        self.nonelastic_thickness = np.sqrt(np.sum((self.curve_top.coords - self.curve_bot.coords)**2, axis=0))
        #Compute elastic thickness by finding the euclidean distance between elastically matched points
        self.elastic_thickness = np.sqrt(np.sum((self.curve_top.coords - self.curve_bot_gamma_adjusted.coords)**2, axis=0))
        self.compute_medial_curve()

    def compute_medial_curve(self):
        medial_curve_coords = self.curve_top.coords - (0.5)*(self.curve_top.coords - self.curve_bot_gamma_adjusted.coords)
        self.medial_curve = Curve(medial_curve_coords)

    def output_thickness_ucf(self):
        self.join_top_and_bottom()
        joined_elastic_thickness = np.array(list(self.elastic_thickness) + list(self.elastic_thickness)[::-1])
        joined_nonelastic_thickness = np.array(list(self.nonelastic_thickness) + list(self.nonelastic_thickness)[::-1])
        curveio.WriteUCF(self.joined_elastic_coords, "thickness", joined_elastic_thickness, self.subject_name + "_elastic_thickness.ucf")
        curveio.WriteUCF(self.joined_nonelastic_coords, "thickness", joined_nonelastic_thickness, self.subject_name + "_nonelastic_thickness.ucf")

    def join_top_and_bottom(self):
        output_plane_data = np.array(list(self.plane_dim_data[1]) + list(self.plane_dim_data[2]))
        output_plane_dim = self.plane_dim_data[0]

        for coord_index in xrange(len(self.curve_top.coords)):
            self.joined_nonelastic_coords.append(np.array(list(self.curve_top.coords[coord_index]) +
                                              list(self.curve_bot.coords[coord_index][::-1])))
            self.joined_elastic_coords.append(np.array(list(self.curve_top.coords[coord_index]) +
                                              list(self.curve_bot_gamma_adjusted.coords[coord_index][::-1])))

        self.joined_nonelastic_coords.insert(output_plane_dim, output_plane_data)
        self.joined_nonelastic_coords = np.array(self.joined_nonelastic_coords).transpose()
        self.joined_elastic_coords.insert(output_plane_dim, output_plane_data)
        self.joined_elastic_coords = np.array(self.joined_elastic_coords).transpose()

    def plot_thicknesses(self, plot_title=False, include_naive_plot=True):
        set_fontsize = 10
        if include_naive_plot:
            test_code_only = "\n% Difference of Mean Thickness Values: " + str( round(100 * (np.average(self.nonelastic_thickness-self.elastic_thickness))/(.5*(np.average(np.average(self.nonelastic_thickness)+np.average(self.elastic_thickness)))),2))
            plt.subplot(211)
            if plot_title:
                plt.title(plot_title+test_code_only, fontsize=set_fontsize)
            else:
                plt.title(self.subject_name+test_code_only, fontsize=set_fontsize+1)
            plt.subplots_adjust(hspace=0.2)
            plt.tick_params(labelsize=set_fontsize, labelbottom=False, labelleft=False)
            plt.xlabel("Uniform Matching", fontsize=set_fontsize)
            self.add_thickness_plot_given_curves(self.curve_top, self.curve_bot)
            plt.plot(self.medial_curve.coords[0], self.medial_curve.coords[1])
            plt.subplot(212)
            plt.tick_params(labelsize=set_fontsize, labelbottom=False, labelleft=False)
            plt.xlabel("Elastic Matching ", fontsize=set_fontsize)
            self.add_thickness_plot_given_curves(self.curve_top, self.curve_bot_gamma_adjusted)
            plt.plot(self.medial_curve.coords[0], self.medial_curve.coords[1])
            plt.savefig(self.subject_name + ".pdf")
        else:
            plt.title(self.subject_name)
            plt.subplot(111)
            self.add_thickness_plot_given_curves(self.curve_top, self.curve_bot_gamma_adjusted)
            plt.plot(self.medial_curve.coords[0], self.medial_curve.coords[1])
            plt.savefig(self.subject_name + ".pdf")

    def add_thickness_plot_given_curves(self, curve1, curve2):
            plt.plot(curve1.coords[0], curve1.coords[1])
            plt.plot(curve2.coords[0], curve2.coords[1])
            for i in xrange(0, len(curve1.coords[0])):
                plt.plot([curve1.coords[0][i], curve2.coords[0][i]],
                         [curve1.coords[1][i], curve2.coords[1][i]])
