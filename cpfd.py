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
def cpfd(specie, PFD):
    """
    converts the given light intensity into internal activation rate
    calibrated for 3 light intensities for two species
    """
    if specie == 'Arabidopsis':
        light =  0.0005833 * PFD**2 + 0.2667 * PFD + 187.5
    elif specie == 'Pothos':
        light = 0.0004167 * PFD**2 + 0.3333 * PFD + 862.5

    return light