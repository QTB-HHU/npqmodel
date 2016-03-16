#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" A universal mathematical model of non-photochemical quenching

Copyright (C) 2015-2016  Anna Matuszyńska, Oliver Ebenhöh

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (license.txt).  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np

def pH(x):
    """ returns the pH of the solution for given hydrogen ions concentration"""
    return (-np.log(x*2.5e-4)/np.log(10))

def pHinv(ph):
    """ returns the hydrogenions concentration for given pH"""
    return (4e3*10**(-ph))

def load_colours():
    """ returns list of 8 colours that are unambiguous both to colorblinds and non-colorblinds
    source: http://jfly.iam.u-tokyo.ac.jp/color/#see
    """
    rgb = [(0, 0, 0), (230, 159, 0), (86, 180, 233), (0, 158, 115), (240, 228, 66),
           (0, 114, 178), (213, 94, 0), (204, 121, 167)]

    # Scale the RGB values to the [0, 1] range (acceptable by matplotlib)
    for i in range(len(rgb)):
        r, g, b = rgb[i]
        rgb[i] = (r / 255., g / 255., b / 255.)

    return rgb
