#!/usr/local/epd/bin/python

__author__ = "Brandon Ayers"
__copyright__ = "Copyright 2013, Brandon Ayers, Shantanu H. Joshi \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "ayersb@ucla.edu"

import numpy as np
from curvematch import match
from curvematch import geodesics
from curvematch.curve import Curve
from curvematch import plotting
import os
from shapeio import curveio
from shapeio import convert
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt


class CorpusCallosum:

    def __init__(self, subject_name, curvefile_path_top, curvefile_path_bottom, resample_siz=100,
                 geodesic_steps=5, linear=False, template_curve=False, linear_template_matching=False, outdir=''):

        self.settings = geodesics.Geodesic()
        self.settings.steps = geodesic_steps
        self.settings.closed = False
        self.linear = linear
        self.linear_template_matching = linear_template_matching
        if type(template_curve) is str:
            self.template_curve = os.path.abspath(template_curve)
        elif template_curve:
            self.template_curve = template_curve
        else:
            self.template_curve = False

        self.subject_name = subject_name
        self.return_shapes = True
        self.resample_siz = resample_siz
        self.rotation = False

        self.curvefile_path_top = os.path.abspath(curvefile_path_top)
        self.curvefile_path_bottom = os.path.abspath(curvefile_path_bottom)
        self.outdir = outdir

        self.gamma = []
        self.curve_top = Curve(file=self.curvefile_path_top)
        self.curve_bot = Curve(file=self.curvefile_path_bottom)
        self.medial_curve = []
        self.plane_dim_data = []

        self.curve_bot_elastic = []
        self.nonelastic_thickness = []
        self.elastic_thickness = []

        self.joined_elastic_curve = []
        self.joined_elastic_thickness = []
        self.joined_nonelastic_curve = []
        self.joined_nonelastic_thickness = []


    def match_top_and_bottom_curves(self):
        self.curve_top.resample_curve_uniform(self.resample_siz)
        self.curve_bot.resample_curve_uniform(self.resample_siz)
        diff_between_dims = np.abs(np.sum(self.curve_top.coords - self.curve_bot.coords, 1))
        least_var_dim = np.argmin(diff_between_dims)
        coords_top, coords_bot = [], []
        for i in xrange(self.curve_top.dim()):
            if i == least_var_dim:
                self.plane_dim_data = [i, self.curve_top.coords[i], self.curve_bot.coords[i]]
            else:
                coords_top.append(self.curve_top.coords[i])
                coords_bot.append(self.curve_bot.coords[i])
        self.curve_top.coords = np.array(coords_top)
        self.curve_bot.coords = np.array(coords_bot)
        self.curve_bot_elastic = match.elastic_curve_matching(self.curve_top, self.curve_bot, self.settings, linear=self.linear)
        self.gamma = self.curve_bot_elastic.geodesic.gamma



    def compute_thickness(self):
        self.match_top_and_bottom_curves()
        #Compute nonelastic thickness by finding the euclidean distance between uniformly sampled points
        self.nonelastic_thickness = np.sqrt(np.sum((self.curve_top.coords - self.curve_bot.coords)**2, axis=0))
        #Compute elastic thickness by finding the euclidean distance between elastically matched points
        self.elastic_thickness = np.sqrt(np.sum((self.curve_top.coords - self.curve_bot_elastic.coords)**2, axis=0))
        self.compute_medial_curve()
        if self.template_curve:
            if type(self.template_curve) == str:
                self.template_curve = Curve(self.template_curve)
            self.join_top_and_bottom()
            self.joined_nonelastic_curve = match.elastic_curve_matching(self.template_curve,
                                                                        self.joined_nonelastic_curve, self.settings, linear=self.linear_template_matching)
            self.joined_nonelastic_thickness = (Curve(np.array(self.joined_nonelastic_thickness)).return_reparameterized_by_gamma(
                self.joined_nonelastic_curve.geodesic.gamma)).coords
            self.joined_elastic_curve = match.elastic_curve_matching(self.template_curve,
                                                                     self.joined_elastic_curve, self.settings, linear=self.linear_template_matching)
            self.joined_elastic_thickness = (Curve(np.array(self.joined_elastic_thickness)).return_reparameterized_by_gamma(
                self.joined_elastic_curve.geodesic.gamma)).coords

    def compute_medial_curve(self):
        medial_curve_coords = self.curve_top.coords - (0.5)*(self.curve_top.coords - self.curve_bot_elastic.coords)
        self.medial_curve = Curve(medial_curve_coords)

    def output_thickness_ucf(self, convert_to_vtp=True):
        if self.joined_elastic_curve == [] or self.joined_nonelastic_curve == []:
            self.join_top_and_bottom()
        if self.linear:
            output_file_name = os.path.join(self.outdir, self.subject_name +
                                          "_uniform" * bool(self.linear_template_matching)*bool(self.template_curve) +
                                          "_reg" * bool(self.template_curve) + "_linear_thickness")
            curveio.WriteUCF(self.joined_nonelastic_curve.coords.T, "thickness", self.joined_nonelastic_thickness,
                             output_file_name + ".ucf")
        else:
            output_file_name = os.path.join(self.outdir, self.subject_name +
                                          "_uniform" * bool(self.linear_template_matching)*bool(self.template_curve) +
                                          "_reg" * bool(self.template_curve) + "_elastic_thickness")

            curveio.WriteUCF(self.joined_elastic_curve.coords.T, "thickness", self.joined_elastic_thickness,
                             output_file_name + ".ucf")
        if convert_to_vtp:
            convert.curve_format(output_file_name+'.ucf', output_file_name+'.vtp')




    def join_top_and_bottom(self):
        if self.joined_nonelastic_curve != [] or self.joined_elastic_curve != []:
            return
        output_plane_data = np.array(list(self.plane_dim_data[1]) + list(self.plane_dim_data[2]))
        output_plane_dim = self.plane_dim_data[0]
        joined_nonelastic_coords = []
        joined_elastic_coords = []
        for coord_index in xrange(len(self.curve_top.coords)):
            joined_nonelastic_coords.append(np.array(list(self.curve_top.coords[coord_index]) +
                                              list(self.curve_bot.coords[coord_index][::-1])))
            joined_elastic_coords.append(np.array(list(self.curve_top.coords[coord_index]) +
                                              list(self.curve_bot_elastic.coords[coord_index][::-1])))

        joined_nonelastic_coords.insert(output_plane_dim, output_plane_data)
        joined_nonelastic_coords = np.array(joined_nonelastic_coords).transpose()
        self.joined_nonelastic_curve = Curve(joined_nonelastic_coords)
        self.joined_nonelastic_thickness = np.array(list(self.nonelastic_thickness) +
                                                    list(self.nonelastic_thickness)[::-1])
        joined_elastic_coords.insert(output_plane_dim, output_plane_data)
        joined_elastic_coords = np.array(joined_elastic_coords).transpose()
        self.joined_elastic_curve = Curve(joined_elastic_coords)
        self.joined_elastic_thickness = np.array(list(self.elastic_thickness) + list(self.elastic_thickness)[::-1])


    def plot_thicknesses(self, plot_title=False, plot_linear=False, plot_both=False):
        set_fontsize = 10
        if bool(self.template_curve):
            if self.linear:
                plotting.plot_matching(self.subject_name + '_uniform_' * self.linear_template_matching *
                                       bool(self.template_curve) +
                                       '_template_matching_' * bool(self.template_curve) + 'linear',
                                       self.template_curve, self.joined_nonelastic_curve,
                                       outdir=self.outdir)
            else:
                plotting.plot_matching(self.subject_name + '_uniform_' * self.linear_template_matching *
                                       bool(self.template_curve)
                                       + '_template_matching_' * bool(self.template_curve) + 'elastic',
                                       self.template_curve, self.joined_elastic_curve,
                                       outdir=self.outdir)
        else:
            fig = plt.figure(1)
            if plot_title:
                plt.title(plot_title)
            else:
                plt.title(self.subject_name + '\nLinear Matching' * plot_linear +
                          '\nElastic Matching' * (not plot_linear))
            if plot_linear:
                self.add_thickness_plot_given_curves(self.curve_top, self.curve_bot)
            else:
                self.add_thickness_plot_given_curves(self.curve_top, self.curve_bot_elastic)
            plt.plot(self.medial_curve.coords[0], self.medial_curve.coords[1])
            plt.axis('equal')
            plt.savefig(os.path.join(self.outdir, self.subject_name + "_linear"*plot_linear +
                                     "_elastic"*(not plot_linear) + ".pdf"))
            plt.close()

    def plot_thickness_comparison(self):
        set_fontsize = 10
        ax = plt.subplot(111)
        plt.title(self.subject_name+" Elastic vs Linear thickness")
        plt.xlabel("Point number")
        plt.ylabel("Thickness value")
        ax.plot(self.nonelastic_thickness, label="Linear Thickness")
        ax.plot(self.elastic_thickness, label="Elastic Thickness")
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)
        plt.savefig(os.path.join(self.outdir, self.subject_name + "_thickness_comparison" + ".pdf"))
        plt.close()

    def save_run_data(self):
        np.savetxt(os.path.join(self.outdir, self.subject_name + "_gamma.csv"), self.gamma, delimiter=",")


    def add_thickness_plot_given_curves(self, curve1, curve2):
            plt.plot(curve1.coords[0], curve1.coords[1])
            plt.plot(curve2.coords[0], curve2.coords[1])
            for i in xrange(0, len(curve1.coords[0])):
                plt.plot([curve1.coords[0][i], curve2.coords[0][i]],
                         [curve1.coords[1][i], curve2.coords[1][i]])