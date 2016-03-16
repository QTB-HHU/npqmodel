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
import matplotlib.pyplot as plt
import scipy.optimize

import parameters
import npqModel
import npqResults
import simulate
import misc

import dataAnalysis
import visualizeexperiment

db = dataAnalysis.DB()

# =============================================================== #
#              load parameters and create model object            #
# =============================================================== #
p = parameters.NPQmodelParameterSet()
m = npqModel.NPQModel(p)
s = simulate.Sim(m)

# =============================================================== #
# initial concentrations corresponding to a dark adapted organism #
# =============================================================== #
y0 = np.array([0., p.bH*misc.pHinv(p.pHstroma), 0, 1, 1, 0])


