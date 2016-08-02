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

import dataAnalysis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import misc

plt.rcParams['legend.numpoints'] = 1


class VisualizeExperiment(dataAnalysis.DB):

    """ A class of routines to visualize a single experiment averaging between the replicates.
        Important: each fluorescence measurement is first normalized to the first measurement of the experiment

    Attributes:
        light_intens: An integer representing the intensity of light used for the light phases in the experiment
        dark_dur:
        species:
    """

    def __init__(self,  specie, intensity, duration, dblocation='DataBase/paperdata_new.db'):
        """
        :param specie: name of the specie. Accepts: Arabidopsis, Pothos
        :param intensity: light intensity. Accepts: 100, 300, 900
        :param duration: number of minutes in the darkness. Accepts: 15, 30, 60
        """
        dataAnalysis.DB.__init__(self, dblocation)

        self.specie = specie
        self.intensity = intensity
        self.duration = duration
        self.marker = {'Arabidopsis': ['o', '--'],
                       'Pothos': ['^', '-']}
        self.color = misc.load_colours()

    def getAvgValue(self, value):
        """
        :param value: field name to retrieve. Accepts: 'Fm', 'Ft', 'T', 'qN'
        :return: an array with the avg value for the three replicates in the first array and
                 standard deviation in the second
        """

        if value not in ['Fm', 'Ft', 'T', 'qN', 'NPQ', 'Yield', 'ETR']:
            print 'Not existing value'
        else:
            n = len(self.retrieve_data_sets({'specie':self.specie, 'lightintensity':self.intensity, 'darkduration':self.duration}))
            print n
            x = np.zeros([n, 22])
            x_avg = np.zeros([1, 22])
            x_sd = np.zeros([1, 22])

            ds = self.retrieve_data_sets({'specie': self.specie,
                                          'lightintensity': self.intensity,
                                          'darkduration': self.duration})
            for exp in range(n):
            #for exp in range(len(ds)):
                if value == 'NPQ':
                    # calculate the NPQ value our of the Fm':
                    fluo = getattr(ds[exp], 'Fm')
                    x[exp, :] = (fluo[0] - fluo)/fluo
                else:
                    x[exp, :] = getattr(ds[exp], value)

                # normalization of the fluorescence
                if value in ['Fm', 'Ft']:
                    fm = getattr(ds[exp], 'Fm')[0]
                    x[exp, :] = x[exp, :]/fm

                for i in range(22):
                    # hardcoded number of measured points as the script is crafted for specific data set
                    x_avg[0, i] = np.mean(x[:, i])
                    x_sd[0, i] = np.std(x[:, i])

        return x_avg[0], x_sd[0]

    def setGraphics(self, T, ylim, yheight=0.1):
        """
        :param T: a five element list of time points of the light change
        :param ylim:
        :param yheight: height of the bar, default value set to compose well with the plot with ylim <0,1>
        :return: set up the current axes with the xlabel ticks indicating switch points for turning light on and off
                 and adds color bars on the top visualizing light memory protocol
        """
        axes = plt.gca()
        axes.set_xticks([0, T[1], T[2], T[3], T[4]])
        axes.set_xticklabels([0, '', np.int(T[2]), np.int(T[3]), np.int(T[4])])
        axes.set_xlim(T[0], T[-1])
        axes.set_ylim(0, ylim)
        axes.set_xlabel('time [s]', fontsize=16)

        for p in [
            patches.Rectangle(
                (T[0], ylim-yheight), T[1]-T[0], yheight, fill=True,
                edgecolor="black", color='black'
            ),
            patches.Rectangle(
                (T[1], ylim-yheight), T[2]-T[1], yheight, fill=False,
                edgecolor="black"
            ),
            patches.Rectangle(
                (T[2], ylim-yheight), T[3]-T[2], yheight, fill=True,
                edgecolor="black", color='black'
            ),
            patches.Rectangle(
                (T[3], ylim-yheight), T[4], yheight, fill=False,
                edgecolor="black"
            ),
        ]:
            axes.add_patch(p)

        return axes

    def plotFluorescence(self):
        """
        :return: plot of a single experiment obtained from three replicates
        """
        T_avg, T_sd = self.getAvgValue('T')
        Fm_avg, Fm_sd = self.getAvgValue('Fm')
        Ft_avg, Ft_sd = self.getAvgValue('Ft')

        # Set the y-axis limit depending on the maximal fluorescence measurements + its error
        self.ylim = 1.11 * (max(Fm_avg) + Fm_sd[np.where(Fm_avg == max(Fm_avg))])
        T = [-100, T_avg[1], T_avg[9], T_avg[16], T_avg[-1]]
        axes = self.setGraphics(T, self.ylim)
        axes.set_ylim(0, self.ylim)

        axes.errorbar(T_avg, Fm_avg/Fm_avg[0], Fm_sd, color='k', linestyle='None', marker='^', markersize=8, label='F_${M}\'$')
        axes.errorbar(T_avg, Ft_avg/Fm_avg[0], Ft_sd, color=self.color[5], linestyle='None', markersize=8, marker='^', label='F_${t}$')

        leg = plt.legend(title='Experiment (triplicates)')
        if leg:
            leg.draggable()

        return axes

    def plotNPQ(self, col, graphic='True'):
        """
        :return:
        """
        T_avg, T_sd = self.getAvgValue('T')
        NPQ_avg, NPQ_sd = self.getAvgValue('NPQ')

        ylim = 1.11 * (max(NPQ_avg) + NPQ_sd[np.where(NPQ_avg == max(NPQ_avg))])
        T = [0, T_avg[1], T_avg[9], T_avg[16], T_avg[-1]]

        if graphic == 'True':
            axes = self.setGraphics(T, ylim)
            axes.set_ylim(0, ylim)
            plt.title('NPQ of ' + str(self.specie) +' exposed to light of ' + str(self.intensity) +
                  '$\mu E m^{-2}s^{-1}$ with ' + str(self.duration) + ' m. dark pause', fontsize=16)
            col = self.color[3]
        else:
            pass
        axes = plt.gca()
        axes.errorbar(T_avg, NPQ_avg, NPQ_sd, color=col,
                      marker=self.marker[self.specie][0],
                      linestyle=self.marker[self.specie][1],
                      label=str(self.intensity))

        return axes

if __name__ == '__main__':
    v = VisualizeExperiment('Arabidopsis', 100, 15)
    print('Load the single experiment for Arabidopsis exposed for 15 minutes to low light')

    plt.figure()
    v.plotFluorescence()
    plt.show()



