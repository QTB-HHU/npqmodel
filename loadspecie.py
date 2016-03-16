#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" A universal mathematical model of non-photochemical quenching

File to import everything needed to start tests from console

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

import parameters
import npqModel
import simulate
import lightProtocol
import misc
import cpfd
import npqResults

def loadspecie(specie = 'Arabidopdis', organism='wt', mechanism=2):
    """
    :return: model object and starting vector
    """

    p = parameters.NPQmodelParameterSet()
    m = npqModel.NPQModel(p, mechanism)

    if organism == 'wt':
        # wild type parameters for Arabidopsis as default
        y0 = np.array([  0.,   6.32975752e-05,   p.APtot/2, 1., 1., 0.])

    elif organism == 'npq4':
        # npq4 mutant lacks the PsbS protein so no change in the y0
        m.par.kProtonationL = 0
        y0 = np.array([  0.,   6.32975752e-05,   p.APtot/2, 1., 1., 0.])

    elif organism == 'npq1':
        # npq1 mutant lacks the deepoxidase
        m.par.kDeepoxV = 0
        y0 = np.array([  0.,   6.32975752e-05,   p.APtot/2, 1., 1., 0.])

    elif organism == 'npq2':
        # npq2 mutant lacks the epoxidase
        m.par.kEpoxZ = 0
        y0 = np.array([  0.,   6.32975752e-05,   p.APtot/2, 1., 0., 0.])

    if specie == 'Pothos':
        # increase the gamma2 coefficient for Pothos
        m.par.gamma2 = 1.


    return m, y0


if __name__ == '__main__':

    m, y0 = loadspecie()