#!/usr/local/epd/bin/python


__author__ = "Brandon Ayers"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "ayersb@ucla.edu"

import os
import match
from geodesics import Geodesic


class CCThickness():
    def __init__(self, curvefile_path_top, curvefile_path_bottom, resample_siz=100, geodesic_steps=100, rotation=False):
        self.setting = Geodesic()
        self.setting.steps = geodesic_steps
        self.resample_siz = resample_siz
        self.rotation = False
        self.curvefile_path_top = os.path.abspath(curvefile_path_top)
        self.curvefile_path_bottom = os.path.abspath(curvefile_path_bottom)


    def get_matched_curve_pair_geodesic(self):
        return match.match_curve_pair(self.curvefile_path_top, self.curvefile_path_bottom,
                                      self.setting, self.rotation, self.resample_siz)

    def compute_thickness(self):
        matched_curve_pair_geodesic = self.get_matched_curve_pair_geodesic()
        print matched_curve_pair_geodesic.gamma


