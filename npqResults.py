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
import matplotlib.pyplot as plt
import numpy as np
import misc
from simulate import Sim


class NPQResults(Sim):
    """
    routines to extract the results of the simulations and easily plot them
    """

    def __init__(self, s):
        """
        class of methods inheriting from Simulation useful for plotting and analysis
        """
        self.model = s.model
        self.results = s.results
        self.par = s.model.par


    def getT(self, r=None):
        """
        :param r: range of steps of simulation for which results we are interested in
        :return: time of all results as one vector
        """
        if r is None:
            r = range(len(self.results))
            
        T = np.hstack([self.results[i]['t'] for i in r])
        
        return T

    def getVar(self, j, r=None):
        """
        :param j: int of the state variable [0:PQred, ..., 5:ATPsynth]
        :param r: range of steps of simulation for which results we are interested in
        :return: concentrations of variable j as one vector
        """
        if r is None:
            r = range(len(self.results))

        X = np.hstack([self.results[i]['y'][:,j] for i in r])
        
        return X

    def getL(self, t):
        """
        :param t: takes time and returns what was the light intensity at this moment
        :return: int
        """
        r = range(len(self.results))
        for i in range(len(r)+1):
            if self.getT([i])[0] <= t < self.getT([i])[-1]:
                x = i
                break

        L = self.results[x]['lfn'].lightintensity(t)

        return L

    def getRate(self, rate, r=None):
        """
        :param rate: name of the rate
        :param r: range of steps of simulation for which results we are interested in
        :return: rate with name 'rate' of all results as one vector
        """
        if r is None:
            r = range(len(self.results))

        V = np.array([])

        for i in r:
            t = self.results[i]['t']
            y = self.results[i]['y']
            l = self.results[i]['lfn']

            V = np.hstack([V,
                           np.array(
                           [self.model.rates(y[j], l.lightintensity(t[j]))[rate] for j in range(len(t))]
                           )])

        return V

    def fluo(self, r=None):
        """
        :param r: range of steps of simulation for which results we are interested in
        :return: time vector, overall fluorescence from PSII, fluorescence emitted from open reaction centres,
        fluorescence emitted from closed reaction centres and state of the PSII
        """
        if r == None:
            r = range(len(self.results))

        Fm = [] # maximal fluorescence
        Fs = [] # base fluorescence
        F = []  # total fluorescence trace (flashes)
        T = []
        Tm = [] # time for maximal fluorescence (flashes)
        Ts = []
        Fmax = 0
        Bst = np.array([]).reshape(0, 4)

        for i in r:
            t = self.results[i]['t']
            y = self.results[i]['y']
            l = self.results[i]['lfn']
            Pred = y[:,0]
            L = y[:,3]
            V = y[:,4]
            Q = self.model.quencher(L,V)

            B = np.array([self.model.ps2states(Pred[j], Q[j], l.lightintensity(t[j])) for j in range(len(t))])

            T = np.hstack([T, t])
            Bst = np.vstack([Bst, B])

            # The fluorescence signal comes from two states: from state B1 is proportional to kF /(kH + kF + kPchem )
            # and that resulting from B3 is proportional to kF /(kH + kF ).
            # Further, the fluorescence signals are proportional to the occupation of the two respective ground states
            # and both are proportional to the cross-section of PSII.

            Fst0 = self.par.kF / (self.par.kF + self.par.kH * Q + self.par.k2) * B[:, 0]
            Fst2 = self.par.kF / (self.par.kF + self.par.kH * Q) * B[:, 2]

            F = np.hstack([F, Fst0 + Fst2])

            # Base fluorescence F0 is observed when all reaction centres are open
            # Maximal fluorescence FM is observed under saturating flashes

            if l.lightintensity(t[0]) == 5000:
                Tm = np.hstack([Tm, t[0]])
                Fm = np.hstack([Fm, max(Fst0 + Fst2)])
            else:
                Ts = np.hstack([Ts, t[0]])
                Fs = np.hstack([Fs, min(Fst0 + Fst2)])

            if i == 0:
                Fmax = Fst0 + Fst2

        return T, F, max(Fmax), Tm, Fm, Ts, Fs, Bst

    # =============================================================================================================== #
    #                                               Plotting routines                                                 #
    # =============================================================================================================== #

    def plotph(self, r=None):
        """
        :param r: range of steps of simulation for which results we are interested in
        :return: plot of lumenal ph
        """
        if r is None:
            r = range(len(self.results))

        T = self.getT()
        H = self.getVar(1)

        plt.plot(T, misc.pH(H), 'g')
        plt.title('Lumenal ph')


    def plotRel(self, line='--', r=None):
        """
        :return: plot results of integration for PQ, ATP, LHC, Xtot and quencher
        """
        if r is None:
            r = range(len(self.results))

        t = self.getT(r)
        PQH = self.getVar(0, r)
        A = self.getVar(2, r)
        L = self.getVar(3, r)
        V = self.getVar(4, r)

        fig = plt.figure()

        plt.plot(t, PQH/self.par.PQtot, color='black', linewidth=2, label = 'PQred') # PQred
        plt.plot(t, 1-L/self.par.Ltot, 'm', linestyle=line, linewidth=2, label = 'fast component') # PC
        plt.plot(t, 1-V/self.par.Xtot, 'b', linestyle=line, linewidth=2, label = 'slow component') # PC
        plt.plot(t, self.model.quencher(L,V), color='red', linewidth=2, label = 'quencher') # PC
        plt.plot(t, A/self.par.APtot, color='cyan', label='ATP')             # ATP

        plt.xlabel('time')
        plt.title('Temporal evolution of state variables')
        plt.legend(loc='best')

        return fig

    def plotQuencherStates(self, r=None):
        """
        :return: plot with occupation of one of four possible states for antennae
        """
        if r is None:
            r = range(len(self.results))

        t = self.getT(r)
        L = self.getVar(3, r)
        V = self.getVar(4, r)
        Lp = self.par.Ltot - L
        Z = self.par.Xtot - V

        fig = plt.figure()

        plt.plot(t, self.par.gamma0 * V * L, color = 'magenta', ls='--', label ='state 0')
        plt.plot(t, self.par.gamma1 * V * Lp, color = 'blue', ls='-.', label ='state 1')
        plt.plot(t, self.par.gamma2 * Z * L , color = 'red', ls=':', label ='state 2')
        plt.plot(t, self.par.gamma3 * Z * Lp, color = 'black', ls='--', label ='state 3')
        plt.xlabel('time')
        plt.title('State occupation for quencher')
        plt.legend(loc='best')

        return fig

    def plotRate(self, rate, r=None):
        """
        :param rate: name of one of the reactions
        :return: plot reaction rate of a selected rate
        """
        if r is None:
            r = range(len(self.results))

        t = self.getT(r)
        v = self.getRate(rate, r)

        fig = plt.figure()
        plt.plot(t, v)
        plt.xlabel('time')
        plt.title('Rate of ' + rate)

        return fig

    def plotFluo(self, color='b'):
        """
        :return: fluorescence trace, uses function fluo to calculate FM and F0
        """
        T, F, Fmax, _,_,_,_,_ = self.fluo()

        fig = plt.figure()
        plt.plot(T, F/Fmax, c=color)
        plt.xlabel('time')
        plt.title('Normalized fluorescence trace')

        return fig

    def plotNPQ(self, color='g'):
        """
        :return:  Stern-Volmer parameter [Holzwarth et al.2013]
        NPQ trace uses function fluo to calculate FM and F0
        """
        T, F, Fmax, Tm, Fm, _,_,_ = self.fluo()

        fig = plt.figure()
        plt.plot(Tm, (Fmax - Fm)/Fm, c=color, marker='o', linestyle='none')
        plt.xlabel('time')
        plt.title('NPQ')

        return fig


    def plotLight(self, r=None):
        """
        :return: changes in the light used in the simulation
        """
        if r is None:
            r = range(len(self.results))

        L = []

        for i in r:
            t = self.getT([i])
            l = self.results[i]['lfn']
            L = np.hstack([L, [l.lightintensity(t[j]) for j in range(len(t))]])

        fig = plt.figure()
        plt.plot(self.getT(r), L)
        plt.title('Light function')

        return fig

    def plotPS2(self):
        """ plot fluorescence trace, uses function fluo to calculate FM and F0 """
        T, F, Fmax, Tm, Fm, Ts, Fs, Bst = self.fluo()

        fig = plt.figure()
        for i in range(4):
            plt.plot(T, Bst[:,i], label='B'+str(i))
        plt.xlabel('time')
        plt.title('Occupation of states in the photosystem II')
        plt.legend(loc='best')

        return fig

