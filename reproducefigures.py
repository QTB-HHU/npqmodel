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
import dataAnalysis
import visualizeexperiment
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import gridspec

import parameters
import npqModel
import simulate
import misc
import npqResults
import cpfd
import loadspecie


db = dataAnalysis.DB()

class ReproduceFigures():

    def __init__(self, PlotSave='Display'):
        self.light_intens = ['100', '300', '900']
        self.dark_dur = ['15', '30', '60']

        # common attributes for plotting consistency
        self.col = misc.load_colours()
        self.species = ['Arabidopsis', 'Pothos']
        self.marker = {'Arabidopsis': ['o', '--'],
                       'Pothos': ['^', '-']}
        self.PlotSave = PlotSave

        self.bbox_props = dict(boxstyle="square,pad=0.3", fc="white", ec="k", lw=2)

    def save_plot(self, figure, name):
        """ Plot the figure and display or save """

        if self.PlotSave == "Display":
            figure.show()
        elif self.PlotSave == "Save":
            FileName = 'Matuszynska_' + name
            figure.savefig(FileName)
            plt.close(figure)
        else:
            print("Error in plotting")


    def flashes(self, v, PFD):
        """
        :param v: object with data of specific experiment
        :param PFD: light intensity of the light phase
        :return: experiment specific time vector when saturating pulses of light were applied
                and corresponding vector of light intensities
        """

        flashes = v.getAvgValue('T')[0]
        flashes[0] = 1

        light = [0, 0, PFD, PFD, PFD, PFD, PFD, PFD, PFD, PFD, 0, 0, 0, 0, 0, 0, 0, PFD, PFD, PFD, PFD, PFD]
        Tfls = []
        PFDs = []

        for t in flashes:
            Tfls = np.hstack([Tfls, t, t+0.8])
        for l in range(len(light)):
            PFDs = np.hstack([PFDs, light[l], 5000])

        return Tfls, PFDs


    def setMultiFigure(self, n_rows, n_cols, xlabels, ylabels, figsize=(30,15)):
        """ set up multiple graphs figure with shared y-axis
            each row corresponds to specific experiment """
        fig, ax = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=figsize)

        for i in range(n_cols):
            ax[0][i].text(x=190, y=1.12, s=str(xlabels[i])+' $\mu \mathrm{E m^{-2}s^{-1}}$', fontsize=20, weight='bold',
                          ha="center", bbox=self.bbox_props)
        for i in range(n_rows):
            ax[i][0].text(x=-80, y=0.5, s=ylabels[i], weight='bold', fontsize=18, ha="center", bbox=self.bbox_props)

        return fig, ax

    def fig1b(self, specie='Arabidopsis', cut='Yes'):
        """
        :param specie: default Arabidopsis.
        :param cut: default Yes; cuts the period of light to only 5 minutes as the first period of light is much longer
        :return: Comparison of Fm' measured in the first 5 minutes of first and second phase of light.
        """
        row = 0
        fig, axes = self.setMultiFigure(3, 3, [100, 300, 900], [15, 30, 60], (15,15))
        fig.subplots_adjust(hspace=0, wspace=0)

        for d_d in self.dark_dur:
            column = 0
            for l_i in self.light_intens:

                # calculate the averages and sd for each experiment
                v = visualizeexperiment.VisualizeExperiment(specie, l_i, d_d)
                T, T_sd = v.getAvgValue('T')
                Fm, Fm_sd = v.getAvgValue('Fm')

                # For the visual reasons the figure for the manuscript presents only first 5 min of light phases
                if cut == 'Yes':
                    cutT = 6
                else:
                    cutT = 9

                s1 = axes[row][column].errorbar(T[0:cutT], Fm[0:cutT], Fm_sd[0:cutT], color='b',
                                                marker=self.marker[specie][0], markersize=10, linestyle='None')
                s2 = axes[row][column].errorbar(T[15:] - T[15], Fm[15:], Fm_sd[15:], color=self.col[6],
                                                marker=self.marker[specie][0], markersize=10, linestyle='None')

                # set the labels and ticks
                axes[row][column].set_ylim(0, 1.05)
                axes[row][column].set_yticks([0.25, 0.5, 0.75, 1.])
                axes[row][column].set_yticklabels([0.25, 0.5, 0.75, 1.], fontsize=14)
                axes[row][column].set_xticks(T[0:7])
                axes[row][column].set_xlim(-10, 390)
                axes[row][column].set_xticklabels([0, '', 61, '', '', '', ''], fontsize=18)

                for p in [patches.Rectangle((0, 0), 30, 1.05, fill=True, edgecolor="black", color='grey', alpha=0.3)]:
                    axes[row][column].add_patch(p)
                column += 1
            row += 1

        # mark the end time only for teh last box
        #axes[2][2].set_xticklabels([0, '', 61, '', '', '', 384], fontsize=18)

        leg = plt.legend((s1, s2), ('$\mathrm{F\'_{M}}$ training phase', '$\mathrm{F\'_{M}}$ memory phase'),
                         loc='upper right', fontsize=20)
        if leg:
            leg.draggable()

        fig.text(0.5, 0.04, 'time [s]', ha='center', fontsize=22)
        fig.text(0.045, 0.5, 'dark duration', va='center', fontsize=22, rotation='vertical')

        self.save_plot(fig, 'Fig1')

    def fig4(self):
        """
        :return: three plots with PAM traces simulated and measured experimentally
        """
        timesft = 5    # shift the simulation by an interval to see the experimental points

        fig = plt.figure(figsize=(17, 5))
        gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1])

        for i in range(len(self.light_intens)):
            ax = plt.subplot(gs[i])
            v = visualizeexperiment.VisualizeExperiment('Arabidopsis', self.light_intens[i], self.dark_dur[i])
            v.plotFluorescence()

            # ===================================================== #
            # read out time of flashes as applied in the experiment #
            # ===================================================== #
            light = cpfd.cpfd('Arabidopsis', int(self.light_intens[i]))
            Tfls, PFDs = self.flashes(v, light)

            # ===================================================== #
            #                   integrate over the time             #
            # ===================================================== #
            model, y0 = loadspecie.loadspecie('Arabidopsis', 'wt')
            s = simulate.Sim(model)
            s.piecewiseConstant(PFDs[0:], Tfls[0:], y0)

            res = npqResults.NPQResults(s)

            T, F, Fmax, Tm, Fm, Ts, Fs, Bst = res.fluo()

            # plot the simulation
            plt.plot(T + timesft, F/max(F), 'orchid', label='simulation')

            # add the legend for the experiment
            ax = plt.gca()
            T_avg, T_sd = v.getAvgValue('T')
            Fm_avg, Fm_sd = v.getAvgValue('Fm')

            l1 = ax.errorbar(T_avg, Fm_avg/Fm_avg[0], Fm_sd, color='k', linestyle='None', marker='^', label='Fm\'')
            l3 = plt.Line2D([], [], linewidth=3, color="orchid")

            leg = plt.legend([l1, l3], ["Experiment", "Simulation"])
            if leg:
                leg.draggable()

            if i == 0:
                ax.annotate(str(self.light_intens[i]) + ' $\mu$Em$^{-2}$s$^{-1}$', (T_avg[9]/2, 1.05), color='k', weight='bold',
                            fontsize=16, ha='center', va='center')
                ax.set_ylabel('Fluorescence ($\mathrm{F_M\'}/\mathrm{F_M}$)',  fontsize=18)
            else:
                ax.annotate(str(self.light_intens[i]), (T_avg[9]/2, 1.05), color='k', weight='bold',
                                fontsize=16, ha='center', va='center')
            ax.annotate(str(self.dark_dur[i]) +' min', (T_avg[9]+(T_avg[16]-T_avg[9])/2, 1.05), color='w', weight='bold',
                            fontsize=16, ha='center', va='center')

            if i == 0:
                ax.set_xticks([0, 840, 1800, 2100] )
                ax.set_xticklabels([0, '14', '30', '35'], fontsize=16)
            elif i ==1:
                ax.set_xticks([0, 840, 2700, 3000] )
                ax.set_xticklabels([0, '14', '45', '50'], fontsize=16)
            else:
                ax.set_xticks([0, 840, 4500, 4800] )
                ax.set_xticklabels([0, '14', '75', '80'], fontsize=16)
            ax.set_xlabel('time [min]', fontsize=18)

        plt.tight_layout()

        self.save_plot(fig, 'Fig4')

    def fig5(self):
        """
        :return: figure visualizing pH dependance of quencher and its components.
                 Upper plot: phase plane trajectories
                 Bottom plots: pH and quencher components for 3 light intensities
        """
        color = list(self.col[i] for i in [1,3,5])  # select 3 distinct colours
        mark = ['x', 'o', 's']
        v = visualizeexperiment.VisualizeExperiment('Arabidopsis', 100, 15)

        # set up the figure
        fig = plt.figure(figsize=(7, 10))

        #make outer gridspec
        outer_grid = gridspec.GridSpec(2, 1, hspace=0.2) # gridspec with two adjacent horizontal cellstwo rows
        ax = plt.subplot(outer_grid[0])
        ax.set_ylabel('lumen pH', fontsize=18)
        ax.set_xlabel('quencher activity [Q]', fontsize=18)
        ax.set_ylim(3, 8)
        ax.set_xlim(0.05, 0.4)
        ax.set_yticks([3,4,5,6,7,8])
        ax.set_yticklabels(['',4,5,6,7,8], fontsize=14)
        ax.set_xticks([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4])
        ax.set_xticklabels(['','0.1','','0.2','','0.3','','0.4'], fontsize=14)

        # add the steady state
        m, y0 = loadspecie.loadspecie('Arabidopsis', 'wt')
        ss = simulate.Sim(m)

        scan = [10, 50, 100, 120, 150, 200, 220, 300, 320, 330, 340, 350, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500]
        Y = ss.steadyStateLightScan(scan, y0)

        ax.plot(m.quencher(Y[:,3], Y[:,4]), misc.pH(Y[:,1]), c='r', linestyle='-')

        for i in[6, 8, 17]:
            ax.plot(m.quencher(Y[i,3], Y[i,4]), misc.pH(Y[i,1]), 'o', linestyle='None',
                markersize=15, mfc='r')

        upper_cell = outer_grid[1]
        inner_grid = gridspec.GridSpecFromSubplotSpec(2, 1, upper_cell, hspace=0.0)

        # From here we can plot using inner_grid's SubplotSpecs
        ax2 = plt.subplot(inner_grid[0, 0])
        # add the three points for the calibrated light intensity
        PFDs = [cpfd.cpfd('Arabidopsis', 100), cpfd.cpfd('Arabidopsis', 300), cpfd.cpfd('Arabidopsis', 900)]
        time = [0, 30, 930, 1830, 2130]
        nsteps = [31, 901, 901, 301]    # nr of steps for each integration phase
        v.setGraphics(time, 9, 1)
        ax2.set_ylabel('lumen pH', fontsize=14)
        ax2.set_ylim(2, 9)
        ax2.get_xaxis().set_visible(False)
        ax2.get_yaxis().set_visible(True)
        ax2.set_yticks([2,3,4,5,6,7,8])
        ax2.set_yticklabels(['',3,4,5,6,7,8], fontsize=14)

        ax3 = plt.subplot(inner_grid[1, 0])
        ax3.set_ylabel('Q components', fontsize=16)
        ax3.set_xlabel('time [min]', fontsize=18)
        ax3.set_xlim([0,time[-1]])
        ax3.set_xticks(time[1:])
        ax3.set_xticklabels(['30', '14', '30', '35'], fontsize=14)
        ax3.set_yticks([0, 0.2, 0.5, 0.8, 1.])
        ax3.set_yticklabels([0, 0.2, 0.5, 0.8, ''], fontsize=14)

        for i in range(len(PFDs)):
            s = simulate.Sim(m)
            s.clearResults()
            PFD = PFDs[i]
            light = [0, PFD, 0, PFD]

            s.piecewiseConstant(light, time[1:], y0, nsteps)
            res = npqResults.NPQResults(s)

            ax.plot(m.quencher(res.getVar(3), res.getVar(4)), misc.pH(res.getVar(1)), c=color[i],
                    marker=mark[i], label=str(self.light_intens[i]))

            ax2.plot(res.getT(), misc.pH(res.getVar(1)), c=color[i], label=str(self.light_intens[i]))

            ax3.plot(res.getT(), 1-res.getVar(3), c=color[i], linestyle='--', label=str(self.light_intens[i]))
            ax3.plot(res.getT(), 1-res.getVar(4), c=color[i], linestyle='-', label=str(self.light_intens[i]))

        handles, labels = ax.get_legend_handles_labels()
        leg1 = ax.legend(handles, labels,
                                      bbox_to_anchor=(0., 1.02, 1., .102), loc=10,
                                      ncol=1,
                                      borderaxespad=0.1)
        if leg1:
            leg1.draggable()

        handles, labels = ax3.get_legend_handles_labels()
        leg = ax3.legend([handles[0], handles[2], handles[4], handles[1], handles[3], handles[5]],
                         [labels[0], labels[2], labels[4], labels[1], labels[3], labels[5]],
                                      bbox_to_anchor=(0., 1.02, 1., .102), loc=10,
                                      title ='PsbS$^P$        Zx',
                                      ncol=2,
                                      borderaxespad=0.1)
        if leg:
            leg.draggable()

        self.save_plot(fig, 'Fig5')

    def fig6(self):
        """
        :return: figure with fluorescence dynamics for Pothos simulated and measured
                 and simulated and measured photosynthetic yield
        """
        light = cpfd.cpfd('Pothos', 100)

        fig = plt.figure(figsize=(7, 3))

        ax = plt.subplot(211)
        v = visualizeexperiment.VisualizeExperiment('Pothos', 100, 15)
        v.plotFluorescence()

        Tfls, PFDs = self.flashes(v, light)

        # load the Pothos model, automatically sets the gamma2 parameter to higher value
        model, y0 = loadspecie.loadspecie('Pothos', 'wt')
        s = simulate.Sim(model)

        s.piecewiseConstant(PFDs[0:], Tfls[0:], y0)

        # load the results into the res object that allow for easier data ansalysis
        res = npqResults.NPQResults(s)

        T, F, Fmax, Tm, Fm, Ts, Fs, Bst = res.fluo()
        timesft = 20    # shift the simulation by an interval to see the experimental points
        ax.plot(T + timesft, F/max(F), 'r')

        # collect handles to set the legend with experimental results
        T_avg, T_sd = v.getAvgValue('T')
        Fm_avg, Fm_sd = v.getAvgValue('Fm')

        l1 = plt.errorbar(T_avg, Fm_avg/Fm_avg[0], Fm_sd, color=self.col[5],
                          linestyle='None', marker='^', label='Experiment')
        l3 = plt.Line2D([], [], linewidth=3, color='r')

        ax.set_ylabel('Fluorescence (a.u.)', fontsize=18)
        ax.set_xlabel('')
        ax.set_xticks([0, 840, 1800, 2100] )
        ax.set_xticklabels([0, '14', '30', '35'])

        leg = plt.legend([l1, l3], ["Experiment", "Simulation"])
        if leg:
            leg.draggable()

        ax2 = plt.subplot(212)
        v.setGraphics([-100, 0, 841,1792,2147], 1, yheight=0.1)

        ax2.plot(v.getAvgValue('T')[0], (v.getAvgValue('Fm')[0] - v.getAvgValue('Ft')[0]) /v.getAvgValue('Fm')[0],
                 color=self.col[5], marker='^', label='experiment')
        ax2.plot(Tm, (Fm - Fs)/Fm, color='r', marker='o', label='simulation')
        ax2.set_xticks([0, 840, 1800, 2100] )
        ax2.set_xticklabels([0, '14', '30', '35'])
        ax2.set_xlabel('time [min]', fontsize=16)
        ax2.set_ylabel('$\Phi$ PSII', fontsize=18)

        leg = plt.legend()
        if leg:
            leg.draggable()

        self.save_plot(fig, 'Fig6')

    def figs3(self, dark=[15]):
        """
        :param dark: list of dark duration
        :return: figure with NPQ traces for Arabidopsis and Pothos.
                 Call with dark == self.dark_dur to reproduce for all three light intensities
        """
        fig = plt.figure(figsize=(12, 15))
        fig.subplots_adjust(hspace=0.5)
        index = len(dark)
        cnt=1
        for d_d in dark:
            # sets up axes depending on the length of the dark list
            axes = fig.add_subplot(index, 1, cnt)
            for s in self.species:
                i = 0
                for l_i in self.light_intens:
                    v = visualizeexperiment.VisualizeExperiment(s, l_i, d_d)
                    v.plotNPQ(self.col[i], 'False')
                    i += 1
                handles, labels = axes.get_legend_handles_labels()
                leg = axes.legend(handles, labels,
                                  bbox_to_anchor=(0., 1.02, 1., .102), loc=10,
                                  title ='Arabidopsis   Pothos  ',
                                  ncol=2,
                                  borderaxespad=0.1)
                if leg:
                    leg.draggable()
                plt.title('NPQ induction for experiment with dark period=' + str(d_d), fontsize=16)

            T_avg, Tsd = v.getAvgValue('T')
            axes.set_ylim(0, 2.7)
            axes.set_xticks([T_avg[0], T_avg[1], T_avg[9], T_avg[16], T_avg[-1]])
            axes.set_xlim(0, 1.01 * T_avg[-1])
            axes.set_xlabel('time [s]', fontsize=16)
            axes.set_ylabel(r'NPQ = $\frac{\mathrm{F_M} - \mathrm{F_M}}{\mathrm{F_M}}$', fontsize=16)

            # Add the light banner
            v.setGraphics([T_avg[0], T_avg[1], T_avg[9], T_avg[16], T_avg[-1]], 2.7, 0.3)
            cnt+=1

        self.save_plot(fig, 'FigS3')

    def figs4(self):
        """ return figure with the experimental results of the extent of relaxation in the dark phase
            for both Arabidopsis and Pothos
            depending on the time spent in the darkness
        """
        textt = ['A 15 minutes', 'B 30 minutes', 'C 60 minutes']


        fig = plt.figure(figsize=(10, 5))
        gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1], wspace=0.0, hspace=0.0)

        for l in range(len(self.dark_dur)):
            # set up the figure
            ax = plt.subplot(gs[l])
            sol = np.zeros([2, 3, 2])   # to collect the experimental points

            cnt = 0

            for i in self.species:
                for d in range(len(self.light_intens)):
                    relax = []
                    ds = db.retrieve_data_sets({'specie': i,
                                                'lightintensity': self.light_intens[d],
                                                'darkduration': self.dark_dur[l]})
                    for exper in range(len(ds)):
                        relax = np.hstack([relax, ds[exper].Fm[15]/ds[exper].Fm[0]])
                    sol[cnt, d, 0] = np.mean(relax)
                    sol[cnt, d, 1] = np.std(relax)

                ax.errorbar([1, 2, 3], sol[cnt, :, 0], sol[cnt, :, 1], marker=self.marker[i][0],
                            linestyle='--', color=self.col[cnt], label=i)
                cnt += 1

            ax.set_xticks([1, 2, 3])
            ax.set_xticklabels([100, 300, 900], fontsize=14)
            ax.set_xlim(0, 4)
            ax.set_ylim(0.7, 1.05)
            ax.set_yticks([0.8, 0.85, 0.9, 0.95, 1.])
            ax.set_yticklabels([0.8,0.85,  0.9, 0.95, 1.], fontsize=14)
            ax.text(0.25, 1.025, textt[l][0], ha="left", va="center", size=15, bbox=self.bbox_props)
            ax.text(2.25, 1.025, textt[l][1:], ha="center", va="center", size=15)

        leg = ax.legend(title='Measurement')
        if leg:
            leg.draggable()

        fig.text(-1, 1, 'Relaxation time (min)', ha='center', rotation='vertical', fontsize=18)
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)

        fig.text(0.5, 0.04, 'Illumination intensity($\mu \mathrm{Em}^{-2}\mathrm{s}^{-1}$)', ha='center', fontsize=18)
        plt.show()

        self.save_plot(fig, 'FigS4')

    def figs6(self, PFDs=[0.1, 10, 50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000]):
        """ returns the steady state fraction of the reduced PQ pool for a given list of different light intensities
            from dim light to high light """
        model, y0 = loadspecie.loadspecie('Arabidopsis', 'wt')
        s = simulate.Sim(model)
        # integrate over time
        Y = s.steadyStateLightScan(PFDs, y0)
        # visualize results
        fig = plt.figure()
        ax = plt.subplot(111)

        plt.plot(PFDs, Y[:, 0]/model.par.PQtot, 'sr', linestyle='-', markersize=10)
        plt.xlabel('light intensity', fontsize=16)
        plt.ylabel('PQH$_2$/(PQ+PQH$_2$)', fontsize=16)

        ax.set_yticks([0, 0.5,1])
        ax.set_yticklabels(['0%','50%','100%'])
        ax.set_xticks([0, 220, 320, 900])
        ax.set_xticklabels(['dim', 100, 300, 900])

        self.save_plot(fig, 'FigS6')

    def figs7(self):
        """ returns figure with NPQ trace for Arabidopsis WT simulated, measured and simulated npq4 mutant"""
        fig = plt.figure()
        ax = plt.subplot(111)
        v = visualizeexperiment.VisualizeExperiment('Arabidopsis', 900, 15)
        v.plotNPQ('r')
        light = cpfd.cpfd('Arabidopsis', 900)
        Tfls, PFDs = self.flashes(v, light)
        # ===================================================== #
        #                   integrate over the time             #
        # ===================================================== #
        model, y0 = loadspecie.loadspecie('Arabidopsis', 'wt')
        s = simulate.Sim(model)
        s.piecewiseConstant(PFDs[0:], Tfls[0:], y0)
        res = npqResults.NPQResults(s)

        T, F, Fmax, Tm, Fm, Ts, Fs, Bst = res.fluo()
        timesft = 20    # shift the simulation by an interval to see the experimental points
        ax.plot(Tm + timesft, (Fm[0] - Fm)/Fm, c=self.col[6], marker='s', linestyle='-', label='wt')

        model2, y0 = loadspecie.loadspecie('Arabidopsis', 'npq4')
        s2 = simulate.Sim(model2)
        s2.piecewiseConstant(PFDs[0:], Tfls[0:], y0)
        res2 = npqResults.NPQResults(s2)

        T, F, Fmax, Tm, Fm, Ts, Fs, Bst = res2.fluo()
        timesft = 20    # shift the simulation by an interval to see the experimental points
        ax.plot(Tm + timesft, (Fm[0] - Fm)/Fm, c=self.col[1], marker='s', linestyle='-', label='npq4')

        # add light banner
        v.setGraphics([-30, 0, 841,1792,2147], 1.7, yheight=0.2)

        handles, labels = ax.get_legend_handles_labels()
        leg = ax.legend(handles, ['simulated wt', 'simulated npq4', 'measured wt'],
                                  bbox_to_anchor=(0., 1.02, 1., .102), loc=10,
                                  borderaxespad=0.1)
        if leg:
            leg.draggable()

        plt.title('')

        self.save_plot(fig, 'FigS7')


if __name__ == '__main__':
    r = ReproduceFigures()

    print('Reproduce the first figure of the manuscript. You can move the legend around')
    r.fig1b()

    plt.show()