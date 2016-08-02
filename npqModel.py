#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Universal mathematical model of the non-photochemical quenching

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
from misc import pH, pHinv


class NPQModel(object):

    """ Simple model of non-photochemical quenching (NPQ)
        includes only two lumbed reactions
        NPQ: violaxanthin cycle and LHCII complex
    """

    def __init__(self, par, mechanism=2):
        """
        :param par: object with the parameter space
        :param mechanism: mechanism of quencher operation, 0-additive, 1-cooperative, 2-default, 4state
        """
        self.par = par
        self.mechanism = mechanism

    # --------------------------------------------------------------------------------------------------------------- #
    #                  Quasi-steady-state approximation for reactions within the PSII                                  #
    # --------------------------------------------------------------------------------------------------------------- #

    def ps2states(self, Pred, Q, L):

        """ QSSA, calculates the states of the photosystem II
            accepts values:
            Pred: reduced fraction of the PQ pool (PQH2)
            Q: quencher
            L: light, int or array of the n x 1 dimension, that gives the light intensity

            returns:
            B: array of arrays with the states of PSII; rows: time, columns states: 1 and 3 excited
        """
        Pox = self.par.PQtot - Pred

        k2 = self.par.k2
        kF = self.par.kF

        # 14 April 2015
        # kH = self.par.kH0 + self.par.kH * Q

        kH = self.par.kH * Q

        k3p = self.par.kPQred * Pox
        k3m = self.par.kPQred * Pred / self.par.KeqPQred

        M = np.array([[-L-k3m, kH+kF,       k3p, 0],
                      [L,      -(kH+kF+k2), 0,   0],
                      [0,      0,           L,   -(kH+kF)],
                      [1,      1,           1,   1]])

        A = np.array([0, 0, 0, self.par.PSIItot])

        # B = [B0,B1*,B2,B3*]
        B = np.linalg.solve(M, A)

        return B

    # --------------------------------------------------------------------------------------------------------------- #
    #                                               reaction rates                                                    #
    # --------------------------------------------------------------------------------------------------------------- #

    def vps2(self,P,L,V,PFD):
        """ reaction rate constant for reduction of PQ """
        Q = self.quencher(L,V)
        B = self.ps2states(P, Q, PFD)

        v = self.par.k2 * B[1] / 2
        return v

    def vPQox(self, Pred, Hlf, PFD):
        """
        :param Pred: reduced PQ (PQH2)
        :param ph: lumen pH (free protons)
        :param PFD: light dependency implicit in b6f, because we do not have PSI in model
        :return: rate of oxidation of the PQ pool through cytochrome and PTOX
        """
        ph = pH(Hlf)

        Keq = self.Keq_cytb6f(ph)
        k = self.par.kCytb6f * PFD

        kf = k * Keq/(Keq+1)
        kr = k/(Keq+1)
        kPTOX = self.par.kPTOX * self.par.xO2 # PTOX reaction added

        v = (kf + kPTOX) * Pred - kr * (self.par.PQtot - Pred)
        return v

    def vATPsynthase(self, A, Hlf, Eact):
        ph = pH(Hlf)
        ADP = self.par.APtot - A
        v = Eact * self.par.kATPsynth * (ADP - A / self.Keq_ATP(ph))
        return v

    def vATPactivity(self, PFD, Eact):
        """ replaces rate constant kATPsynth in the raction of ATP synthase
            increases ATP activity in light and decreases in darkness
        """
        ATPasetot = 1
        Ein = ATPasetot - Eact
        switch = 1.0 * (PFD >= 1)
        vf = self.par.kActATPase * switch * Ein
        vr = self.par.kDeactATPase * (1 - switch) * Eact
        return vf - vr

    def vLeak(self, Hlf):
        """ transmembrane proton leak modelled by mass-action kinetics"""
        v = self.par.kLeak * (Hlf - pHinv(self.par.pHstroma))

        return v

    def vATPconsumption(self, A):
        """ ATP consuming reaction modelled by mass-action kinetics"""
        v = self.par.kATPcons * A

        return v

    def vXcyc(self, V, Hlf):
        """
        activity of xantophyll cycle.
        vf: de-epoxidation of violaxanthin, modelled by Hill kinetics
        vf: epoxidation
        """
        nH = self.par.kHillX
        vf = self.par.kDeepoxV * ((Hlf ** nH)/ (Hlf ** nH + pHinv(self.par.kphSat) ** nH)) * V
        #vr = self.par.kEpoxZ *((Hlf ** nH)/ (Hlf ** nH + pHinv(kphSat) ** nH)) *  (self.par.Xtot - V)
        vr = self.par.kEpoxZ * (self.par.Xtot - V)

        return vf - vr

    def vPsbS(self, L, V, Hlf):
        """
        activity of PsbS protein protonation.
        vf: protonation modelled by Hill kinetics
        vf: deprotonation
        """
        nH = self.par.kHillL # shall be changed, for now =5; also half saturation of pH could be change
        vf = self.par.kProtonationL * ((Hlf ** nH)/ (Hlf ** nH + pHinv(self.par.kphSatLHC) ** nH)) * L
        vr = self.par.kDeprotonation * (self.par.PsbStot - L)

        return vf - vr

    def quencher(self, L, V):
        """
        quencher activity depending on the assumed mechanism as parameter
        0: additive by default
        1: cooperative
        """
        Lp = self.par.PsbStot - L
        Z = self.par.Xtot - V
        ZAnt = Z / (Z + self.par.kZSat)
        if self.mechanism == 0:
            # simple additive mechanism
            Q = self.par.alpha * Lp + self.par.beta * Z / (Z + self.par.kZSat)
        elif self.mechanism == 1:
            # no saturation effect of zeaxanthin
            Q = self.par.gamma0 * V * L \
                + self.par.gamma1 * V * Lp \
                + self.par.gamma2 * ZAnt * Lp \
                + self.par.gamma3 * ZAnt * L
        elif self.mechanism == 2:
            Q = self.par.gamma0 * (1-ZAnt) * L \
                + self.par.gamma1 * (1-ZAnt) * Lp \
                + self.par.gamma2 * ZAnt * Lp \
                + self.par.gamma3 * ZAnt * L
        else:
            raise NotImplementedError

        return Q

    # --------------------------------------------------------------------------------------------------------------- #
    #                            model (parameter) dependent helper functions                                         #
    # --------------------------------------------------------------------------------------------------------------- #
    def Keq_ATP(self, pH):
        DG = self.par.DeltaG0_ATP - self.par.dG_pH * self.par.HPR * (self.par.pHstroma - pH)
        Keq = self.par.Pi_mol * np.exp(-DG/self.par.RT)
        return Keq

    def Keq_cytb6f(self, pH):
        DG1 = -2 * self.par.F * self.par.E0_PQ
        DG2 = -self.par.F * self.par.E0_PC

        DG = - (DG1 + 2*self.par.dG_pH * pH) + 2 * DG2 + 2*self.par.dG_pH * (self.par.pHstroma - pH)
        Keq = np.exp(-DG/self.par.RT)
        return Keq

    # ----------------------------------------------------------------------------------------------------------------#

    def rates(self, y, PFD):
        """
        Method calculating reaction rates.
        input parameters:
        :param y:vector of state variables
        :param PFD: photon flux density
        :return: a dict with reaction rates
       """
        P = y[0]        #PQH2
        Hlf = y[1]      #protons
        A = y[2]        #ATP
        L = y[3]        #LHCII Trimer
        V = y[4]        #Violaxanthin
        Eact = y[5]     # ATPsynthase enzyme activity

        v = {
            'ps2': self.vps2(P, L, V, PFD),  # reaction rate constant for reduction of PQ

            'PQox': self.vPQox(P, Hlf, PFD),

            'ATPsynthase': self.vATPsynthase(A, Hlf, Eact),

            'leak': self.vLeak(Hlf),

            'ATPconsumption': self.vATPconsumption(A),

            'Xcyc': self.vXcyc(V, Hlf),

            'PsbS': self.vPsbS(L, V, Hlf),  # protonation depends on state of antennae (faster for violaxanthin)
            
            'ATPactivity': self.vATPactivity(PFD, Eact),

            'quencher': self.quencher(L, V)
        }

        return v

    def npq_model(self, y, PFD):
        """
        Defining the system of equations governing the evolution of PQ
        return: an array [PQH2, lumenal protons, ATP in stroma, LHCII trimeric, Violaxanthin, active ATPase]
        """
        if y[1] <= 0:
            raise ModelError("H went below zero")

        # Reaction rates, v is a dictionary
        v = self.rates(y, PFD)

        # Output from ODE function must be a COLUMN vector, with n rows
        n = len(y)      # 1: implies its a single ODE
        dydt = np.zeros((n, 1))

        # dPQH/dt = -vPSII + vCytb6f - vcyc + vPTOX - vNDH
        dydt[0] = v['ps2'] - v['PQox']
        # dH/dt = water_split + proton_pump - synthase_ATP - leak
        dydt[1] = (2 * v['ps2'] + 4 * v['PQox'] - self.par.HPR * v['ATPsynthase'] - v['leak']) / self.par.bH
        # dATP/dt = vATPsynthase - vATPconsumption
        dydt[2] = v['ATPsynthase'] - v['ATPconsumption']
        # dLhc/dt = -PsbS
        dydt[3] = - v['PsbS']
        # dViola/dt = -vViolaxanthin
        dydt[4] = - v['Xcyc']
        # dE*/dt = vEnzyme
        dydt[5] = v['ATPactivity']

        return dydt


class ModelError(Exception):

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)