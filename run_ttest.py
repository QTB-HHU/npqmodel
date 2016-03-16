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
from matplotlib import gridspec
import scipy.stats
import matplotlib
import numpy as np
from matplotlib.colors import ListedColormap


class Ttest():

    def __init__(self, db, intens, dur, spec):
        """
        :param db: name of the data base where the data are stored
        :param intens: one of the parameters required to retrieve the desired data set from the db,
        corresponding to light intensity [minutes]
        :param dur: another parameter, corresponding to relaxation duration [minutes]
        :param spec: another parameter, corresponding to the specie [Arabidopsis or Pothos]
        """
        self.db = db
        self.intens = intens
        self.dur = dur
        self.spec = spec

    def test_f_norm(self):
        """
        :return: an array of result of the Ttest performed for
        """
        intensity = np.zeros([len(self.spec), len(self.intens), len(self.dur)])
        for s in range(len(self.spec)):
            for i in range(len(self.intens)):
                for d in range(len(self.dur)):
                    npq1 = []
                    npq2 = []
                    ds = self.db.retrieveDataSets({'specie': self.spec[s],
                                              'lightintensity': self.intens[i],
                                              'darkduration': self.dur[d]})
                    for exper in range(len(ds)):
                        npq1 = np.hstack([npq1, ds[exper].Fm[2]/ds[exper].Fm[0]])
                        npq2 = np.hstack([npq2, ds[exper].Fm[17]/ds[exper].Fm[0]])
                    intensity[s, i, d] = scipy.stats.ttest_rel(npq1, npq2)[1]
        return intensity

    def plotResults(self):
        fig = plt.figure()
        gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, .1], wspace=.1, hspace=0.1)

        # values x and y give values at z
        xmin = 1; xmax = 4; dx = 1
        ymin = 1; ymax = 4; dy = 1
        x, y = np.meshgrid(np.arange(xmin, xmax, dx), np.arange(ymin, ymax, dy))

        # transform x and y to boundaries of x and y
        x2, y2 = np.meshgrid(np.arange(xmin,xmax+dx,dx)-dx/2.,np.arange(ymin,ymax+dy,dy)-dy/2.)

        # coloring
        cMap = ListedColormap(['red', 'white'])

        for i in range(len(self.test_f_norm())):
            ax = plt.subplot(gs[i])
            pmesh = ax.pcolormesh(x2, y2, self.test_f_norm()[i], vmin=0.00,vmax=0.15, cmap=cMap)
            plt.title(str(spec[i]))

            # set the labels
            plt.axis([x2.min(),x2.max(),y2.min(),y2.max()])
            plt.xticks(np.arange(xmin,xmax,dx), [15,30,60])
            plt.yticks(np.arange(ymin,ymax,dy), [100, 300,900])
            # hide the yaxis labels on the second plot
            if i != 0:
                plt.setp(ax.get_yticklabels(), visible=False)

            # annotate the heat map with the p-values
            for y in range(self.test_f_norm()[0].shape[0]):
                for x in range(self.test_f_norm()[0].shape[1]):
                    plt.text(x + 1, y + 1, '%.4f' % self.test_f_norm()[i, y, x],
                         horizontalalignment='center',
                         verticalalignment='center',
                         )

            ax3 = plt.subplot(gs[2])
            cbar = matplotlib.colorbar.ColorbarBase(ax3, cmap=cMap, orientation='vertical')
            plt.setp(ax3.get_xticklabels())
            plt.tick_params(
                    axis='y', left='on', top='off', bottom='off', labelleft='off')

            # shift the axis ticks properly
            cbar.ax.get_yaxis().set_ticks([])

            for j, lab in enumerate(['$<0.05$', '$>0.05$']):
                cbar.ax.text(.5, j/2.+0.15, lab, ha='center', va='center')
            cbar.ax.get_yaxis().labelpad = 20
            cbar.ax.set_ylabel('p-value', rotation=270)
        return fig

if __name__ == '__main__':
    import dataAnalysisPaperData
    import matplotlib.pyplot as plt

    db = dataAnalysisPaperData.DB()
    intens = [100, 300, 900]
    dur = [15, 30, 60]
    spec = ['Arabidopsis', 'Pothos']

    ttest = Ttest(db, intens, dur, spec)

    ttest.plotResults()
    plt.show()


