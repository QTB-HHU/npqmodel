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
from numpy import log, exp


class NPQmodelParameterSet:
    """ Parameters for the model NPQ for specie Arabidopsis thaliana """

    default = {
        # lumenal pool sizes [mmol/molChl]
        'PSIItot': 2.5,         # [mmol/molChl] total concentration of PSII
        'PQtot': 20.,           # [mmol/molChl]
        'APtot': 50.,           # [mmol/molChl] Bionumbers ~2.55mM (=81mmol/molChl)
        'xO2': 8.,              # external oxygen, kept constant, corresponds to 250 microM, corr. to 20%

        # artificial, normalized pool sizes
        'PsbStot': 1.,             # [relative] LHCs that get phosphorylated and protonated
        'Xtot': 1.,             # [relative] xanthophylls - both still relative

        # pH and protons
        'pHstroma': 7.8,
        'kLeak': 1000,          # [1/s] leakage rate -- inconsistency with Kathrine
        'bH': 100.,             # proton buffer: ratio total / free protons
        'HPR': 14./3.,

        # rate constants of the light reactions of zero dimension [ unit: s-1]
        'kH': 5.e9,
        'kF': 6.25e8,           # fluorescence 16ns
        'k2': 5.e9,             # original 5e9 (charge separation limiting step ~ 200ps) - made this faster for higher Fs fluorescence

        # quencher fitted parameters
        'gamma0': 0.1,          # slow quenching of Vx present despite lack of protonation
        'gamma1': 0.25,         # fast quenching present due to the protonation
        'gamma2': 0.6,          # slow quenching of Zx present despite lack of protonation
        'gamma3': 0.15,         # fastest possible quenching

        # ATP parameters
        'kActATPase': 0.01,     # paramter relating the rate constant of activation of the ATPase in the light
        'kDeactATPase': 0.002,  # paramter relating the deactivation of the ATPase at night
        'kATPsynth': 20.,       # taken from MATLAB
        'kATPcons': 10.,        # taken from MATLAB
        'Pi_mol': 0.01,
        'DeltaG0_ATP': 30.6,    # 30.6kJ/mol / RT

        # non-photochemical quenching PROTONATION
        'kDeprotonation': 0.0096,
        'kProtonationL': 0.0096,
        'kphSatLHC': 5.8,

        # non-photochemical quenching XANTOPHYLLS
        'kDeepoxV': 0.0024,
        'kEpoxZ': 0.00024,      # 6.e-4,        # converted to [1/s]
        'kphSat': 5.8,          # [-] half-saturation pH value for activity de-epoxidase, highest activity at ~pH 5.8
        'kHillX': 5.,     # [-] hill-coefficient for activity of de-epoxidase
        'kHillL': 3.,     # [-] hill-coefficient for activity of de-epoxidase
        'kZSat': 0.12,          # [-], half-saturation constant (relative conc. of Z) for quenching of Z

        # standard redox potentials (at pH=0) in V
        'E0_QA': -0.140,
        'E0_PQ': 0.354,
        'E0_PC': 0.380,

        # physical constants
        'F': 96.485,            # Faraday constant
        'R': 8.3e-3,            # universal gas constant
        'T': 298.,              # Temperature in K - for now assumed to be constant at 25 C

        # rate constants
        'kPQred': 250.,         # [1/(s*(mmol/molChl))]
        'kPTOX': .01,           # ~ 5 electrons / seconds. This gives a bit more (~20)
        'kCytb6f': 0.104        # a rough estimate of the transfer from PQ to cyt that is equal to ~ 10ms
                                # [1/s*(mmol/(s*m^2))] - gets multiplied by light to determine rate
    }


    def __init__(self, pars={}):  # -- Anna changed here for pars to be optional
        mypars = pars.copy()

        # load default parameters for electron transport chain reactions
        for k in NPQmodelParameterSet.default.keys():
            mypars.setdefault(k, NPQmodelParameterSet.default[k])

        for k in mypars.keys():
            setattr(self,k,mypars[k])

        self.setCompositeParameters()
        setattr(self,'KeqPQred',self.Keq_PQred())

    def setCompositeParameters(self):
        setattr(self, 'RT', self.R * self.T)
        setattr(self, 'dG_pH', log(10)*self.RT)
        setattr(self, 'Hstroma', 3.2e4*10**(-self.pHstroma))    # proton concentration in stroma -- Anna
        setattr(self, 'kProtonation', 4.e-3 / self.Hstroma)     # [1/s] converted from 4 * 10^-6 [1/ms] protonation of L
                                                                # depends on pH value in lumen

    def Keq_PQred(self):
        DG1 = -self.E0_QA * self.F
        DG2 = -2. * self.E0_PQ * self.F
        DG = -2. * DG1 + DG2 + 2. * self.pHstroma * self.dG_pH
        K = exp(-DG/self.RT)
        return K
