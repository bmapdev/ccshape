#!/usr/local/epd/bin/python

__author__ = "Brandon Ayers"
__copyright__ = "Copyright 2013, Brandon Ayers, Shantanu H. Joshi \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "ayersb@ucla.edu"

"""
Compute and compare thickness of callosal curve segmenations naively or using elastic shape matching methods.

Allows many
"""

import os

import numpy as np
import matplotlib.pyplot as plt

from curvematch import match
from curvematch import geodesics
from curvematch.curve import Curve
from curvematch import plotting
from curvematch.utils import reparameterize_by_gamma
from shapeio import curveio
from shapeio import convert


class CorpusCallosum(object):
    """ Hold and process callosal curve segmentation information"""
    def __init__(self, subject_name, curvefile_path_top, curvefile_path_bottom, resample_siz=100,
                 geodesic_steps=5, linear=False, template_curve=False, linear_template_matching=False,
                 outdir='', alt_registration=False):
        self.settings = geodesics.Geodesic()
        self.settings.steps = geodesic_steps
        self.settings.closed = False

        self.linear = linear
        self.linear_template_matching = linear_template_matching
        self.alt_registration = alt_registration
        self.subject_name = subject_name
        self.return_shapes = True
        self.resample_siz = resample_siz
        self.rotation = False
        self.outdir = outdir
        if type(template_curve) is str:
            self.template_curve = os.path.abspath(template_curve)
        elif template_curve:
            self.template_curve = template_curve
        else:
            self.template_curve = False

        self.curvefile_path_top = os.path.abspath(curvefile_path_top)
        self.curvefile_path_bottom = os.path.abspath(curvefile_path_bottom)
        self.gamma = []
        self.curve_top = Curve(file=self.curvefile_path_top)
        self.curve_bot = Curve(file=self.curvefile_path_bottom)
        self.medial_curve = []
        self.least_var_dim = -1
        self._plane_dim_data = []

        self.curve_bot_elastic = []
        self.attributes_linear_matching = []
        self.attributes_elastic_matching = []

        self.joined_curve_elastic = []
        self.joined_curve_linear = []
        self.joined_attributes_elastic = []
        self.joined_attributes_linear = []
        self.template_geodesic_elastic = []
        self.template_gamma_linear = []
        self._compute()

    def plot(self, plot_title=False, plot_both=False):
        if bool(self.template_curve):
            if self.linear:
                plotting.plot_matching(self.subject_name + '_uniform' * self.linear_template_matching *
                                       bool(self.template_curve) +
                                       '_template_matching_' * bool(self.template_curve) + 'linear',
                                       self.template_curve, self.joined_curve_linear,
                                       outdir=self.outdir)
            else:
                plotting.plot_matching(self.subject_name + '_uniform' * self.linear_template_matching *
                                       bool(self.template_curve)
                                       + '_template_matching_' * bool(self.template_curve) + 'elastic',
                                       self.template_curve, self.joined_curve_elastic,
                                       outdir=self.outdir)
        else:
            plt.figure(1)
            if plot_title:
                plt.title(plot_title)
            else:
                plt.title(self.subject_name + '\nLinear Matching' * self.linear +
                          '\nElastic Matching' * (not self.linear))
            self._attribute_plot()
            #plt.plot(self.medial_curve.coords[0], self.medial_curve.coords[1])
            plt.axis('equal')
            plt.savefig(os.path.join(self.outdir, self.subject_name + "_linear"*self.linear +
                                     "_elastic"*(not self.linear) + ".pdf"))
            plt.close()

    def plot_comparison(self):
        if self.linear:
            return
        ax = plt.subplot(111)
        plt.title(self.subject_name+" Elastic vs Linear Thickness")
        plt.xlabel("Point Number")
        plt.ylabel("Thickness Value")
        ax.plot(self.attributes_linear_matching, label="Linear Thickness", color='r')
        ax.plot(self.attributes_elastic_matching, label="Elastic Thickness", color='g')
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)
        plt.savefig(os.path.join(self.outdir, self.subject_name + "_thickness_comparison" + ".pdf"))
        plt.close()

    def output_ucf(self, convert_to_vtp=True):
        if self.joined_curve_elastic == [] or self.joined_curve_linear == []:
            self._join_top_and_bottom()
        if self.linear:
            output_file_name = os.path.join(self.outdir, self.subject_name +
                                          "_uniform" * bool(self.linear_template_matching)*bool(self.template_curve) +
                                          "_reg" * bool(self.template_curve) + "_linear_thickness")
            curveio.WriteUCF(self.joined_curve_linear.coords.T, "thickness", self.joined_attributes_linear,
                             output_file_name + ".ucf")
        else:
            output_file_name = os.path.join(self.outdir, self.subject_name +
                                          "_uniform" * bool(self.linear_template_matching)*bool(self.template_curve) +
                                          "_reg" * bool(self.template_curve) + "_elastic_thickness")

            curveio.WriteUCF(self.joined_curve_elastic.coords.T, "thickness", self.joined_attributes_elastic,
                             output_file_name + ".ucf")
        if convert_to_vtp:
            convert.curve_format(output_file_name+'.ucf', output_file_name+'.vtp')

    def save_matching_data(self):
        if self.linear:
            return
        np.savetxt(os.path.join(self.outdir, self.subject_name + "_thickness_gamma.csv"),
                   self.gamma, fmt="%f", delimiter=',')
        if self.template_curve:
            np.savetxt(os.path.join(self.outdir, self.subject_name + "_registration_gamma.csv"),
                   self.template_gamma_elastic, fmt="%f", delimiter=',')
        plt.figure(1)
        plt.title(self.subject_name + '\n' + 'Template ' * bool(self.template_curve) +
                          "Gamma")
        plt.xlabel("Point Number")
        plt.ylabel("Gamma Value")
        if not self.template_curve:
            plt.plot(self.gamma, label="Gamma", color='r')
            plt.savefig(os.path.join(self.outdir, self.subject_name + '_template' * bool(self.template_curve) +
                    "_gamma" + ".pdf"))
        plt.close()

    def _compute(self):
        self._match_top_and_bottom_curves()

        # If a template curve was specified, match the joined output curves to it.
        if self.template_curve:
            self._join_top_and_bottom()
            if type(self.template_curve) == str:
                self.template_curve = Curve(self.template_curve)
            self.joined_curve_linear = \
                match.elastic_curve_matching(self.template_curve, self.joined_curve_linear,
                                             self.settings, linear=self.linear_template_matching)

            self.joined_curve_elastic = \
                match.elastic_curve_matching(self.template_curve, self.joined_curve_elastic,
                                             self.settings, linear=self.linear_template_matching)
            self.template_gamma_linear = self.joined_curve_linear.geodesic.gamma
            self.template_gamma_elastic = self.joined_curve_elastic.geodesic.gamma
        self._compute_attributes()
        self._compute_medial_curve()

    def _match_top_and_bottom_curves(self):
        """ Preform elastic shape matching between top and bottom portion of callosal curves."""

        # Resample curves to an equal parameterization for baseline comparison
        self.curve_top.resample_curve_uniform(self.resample_siz)
        self.curve_bot.resample_curve_uniform(self.resample_siz)

        # Input Callosal curve segmentations are three dimensional arrays but callosal thickness is only related to two
        # of these dimensions. Determine which dimension is least variant and store that information for latter, while
        # keeping the relevant dimensions as part of the curves.
        diff_between_dims = np.abs(np.sum(self.curve_top.coords - self.curve_bot.coords, 1))
        self.least_var_dim = np.argmin(diff_between_dims)
        coords_top, coords_bot = [], []
        for i in xrange(self.curve_top.dim()):
            if i == self.least_var_dim:
                # [Position data should be inserted in array, stored top curve data, bot curve data]
                self._plane_dim_data = [i, self.curve_top.coords[i], self.curve_bot.coords[i]]
            else:
                coords_top.append(self.curve_top.coords[i])
                coords_bot.append(self.curve_bot.coords[i])

        self.curve_top.coords = np.array(coords_top)
        self.curve_bot.coords = np.array(coords_bot)
        self.curve_bot_elastic = match.elastic_curve_matching(self.curve_top, self.curve_bot,
                                                              self.settings, linear=self.linear)
        self.gamma = self.curve_bot_elastic.geodesic.gamma

    def _join_top_and_bottom(self):
        output_plane_data = np.array(list(self._plane_dim_data[1]) + list(self._plane_dim_data[2]))
        output_plane_dim = self._plane_dim_data[0]
        joined_nonelastic_coords = []
        joined_elastic_coords = []
        for coord_index in xrange(len(self.curve_top.coords)):
            joined_nonelastic_coords.append(np.array(list(self.curve_top.coords[coord_index]) +
                                            list(self.curve_bot.coords[coord_index][::-1])))
            joined_elastic_coords.append(np.array(list(self.curve_top.coords[coord_index]) +
                                         list(self.curve_bot_elastic.coords[coord_index][::-1])))

        joined_nonelastic_coords.insert(output_plane_dim, output_plane_data)
        joined_nonelastic_coords = np.array(joined_nonelastic_coords).transpose()
        self.joined_curve_linear = Curve(joined_nonelastic_coords)
        self.joined_attributes_linear = np.array(list(self.attributes_linear_matching) +
                                                    list(self.attributes_linear_matching)[::-1])
        joined_elastic_coords.insert(output_plane_dim, output_plane_data)
        joined_elastic_coords = np.array(joined_elastic_coords).transpose()
        self.joined_curve_elastic = Curve(joined_elastic_coords)
        self.joined_attributes_elastic = np.array(list(self.attributes_elastic_matching) + list(self.attributes_elastic_matching)[::-1])

    def _compute_medial_curve(self):
        """ Compute a curve medial to the top and bottom curves. """
        medial_curve_coords = self.curve_top.coords - (0.5)*(self.curve_top.coords - self.curve_bot_elastic.coords)
        self.medial_curve = Curve(medial_curve_coords)

    def _compute_attributes(self):
        pass

    def _attribute_plot(self):
            pass


class CorpusCallosumThickness(CorpusCallosum):
    def __init__(self, subject_name, curvefile_path_top, curvefile_path_bottom, resample_siz=100,
                 geodesic_steps=5, linear=False, template_curve=False, linear_template_matching=False,
                 outdir='', alt_registration=False):
        super(CorpusCallosumThickness, self).__init__(subject_name, curvefile_path_top, curvefile_path_bottom, resample_siz,
                 geodesic_steps, linear, template_curve, linear_template_matching,
                 outdir, alt_registration)

    def _compute_attributes(self):

        # Compute nonelastic thickness by finding the euclidean distance between points
        self.attributes_linear_matching = \
            np.sqrt(np.sum((self.curve_top.coords - self.curve_bot.coords)**2, axis=0))
        self.attributes_elastic_matching = \
            np.sqrt(np.sum((self.curve_top.coords - self.curve_bot_elastic.coords)**2, axis=0))
        self._join_top_and_bottom()

        # Compute thickness for registered curve.
        if self.alt_registration and self.template_curve:
                # Change to compute the Euclidean distance on template match curve to get registered thickness
                nonelastic_reg_curve_top = \
                    self.joined_curve_linear.coords[:, :self.joined_curve_elastic.siz()//2]
                nonelastic_reg_curve_bot = \
                    self.joined_curve_linear.coords[:, self.joined_curve_elastic.siz()//2:][:, ::-1]
                elastic_reg_curve_top = \
                    self.joined_curve_linear.coords[:, :self.joined_curve_elastic.siz()//2]
                elastic_reg_curve_bot = \
                    self.joined_curve_linear.coords[:, self.joined_curve_elastic.siz()//2:][:, ::-1]
                nonelastic_thickness = np.sqrt(np.sum((nonelastic_reg_curve_top - nonelastic_reg_curve_bot)**2, axis=0))
                elastic_thickness = np.sqrt(np.sum((elastic_reg_curve_top - elastic_reg_curve_bot)**2, axis=0))
                self.joined_attributes_linear = np.array(list(nonelastic_thickness)+list(nonelastic_thickness)[::-1])
                self.joined_attributes_elastic = np.array(list(elastic_thickness)+list(elastic_thickness)[::-1])
        elif self.template_curve:
                # Obtain new thickness values by reparameterizing values for thickness computed using linear matching
                # with elastic matching output values.
                self.joined_attributes_linear = \
                    reparameterize_by_gamma(self.joined_attributes_linear, self.template_gamma_linear)
                self.joined_attributes_elastic = \
                    reparameterize_by_gamma(self.joined_attributes_linear, self.template_gamma_elastic)

    def _attribute_plot(self):
        curve1 = self.curve_top
        if self.linear:
            curve2 = self.curve_bot
        else:
            curve2 = self.curve_bot_elastic

        # Plot the top and bottom callosal curves.
        plt.plot(curve1.coords[0], curve1.coords[1])
        plt.plot(curve2.coords[0], curve2.coords[1])

        # Plot the distance between matched points on the curves
        for i in xrange(0, len(curve1.coords[0])):
            plt.plot([curve1.coords[0][i], curve2.coords[0][i]],
                     [curve1.coords[1][i], curve2.coords[1][i]])