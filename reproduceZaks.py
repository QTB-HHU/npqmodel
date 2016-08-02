__author__ = 'anna'


import dataAnalysis
import visualizeexperiment
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import gridspec
import pickle

import parameters
import npqModel
import simulate
import misc
import npqResults
import cpfd
import loadspecie

w, y0w = loadspecie.loadspecie('Arabidopsis', 'wt')
m, y0m = loadspecie.loadspecie('Arabidopsis', 'npq4')

yheight=1.
T = [0, 200, 800, 1200]
ylim = 7

# ================================================== #
# Fig S1: npq4 mutant, no qE
# ================================================== #
PFD = [100, 1000]

for i in range(len(PFD)):
    plt.subplot(2,1,i+1)
    plt.ylabel('FI. Yield (Relative)', fontsize=18)

    axes = plt.gca()

    axes.set_xticks([0, 200, 800, 1200])
    axes.set_xticklabels([0, 200, 800, 1200])
    axes.set_xlim(T[0], T[-1])
    axes.set_ylim(0, 7)
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
        )
    ]:
        axes.add_patch(p)

    axes.annotate(str(PFD[i]) + '$ \mu \mathrm{mol photons/(m^{-2}s)}$', (450, 6.5), color='k',
                            fontsize=16, ha='center', va='center')


    light = cpfd.cpfd('Arabidopsis', PFD[i])

    sw = simulate.Sim(w)
    sm = simulate.Sim(m)

    lT = [0, 5000, 0, 5000,light,
          5000, light, 5000, light, 5000, light, 5000, light, 5000, light, 5000, light, 5000, light,
          5000, light, 5000, light, 5000, light, 5000, light, 5000, light, 5000, light, 5000, light, 5000, light, 5000, 0,
          5000, 0, 5000, 0, 5000, 0, 5000, 0, 5000, 0, 5000, 0, 5000]
    tT = [60, 60.8,120,120.8, 130,
          130.8, 150, 150.8, 180, 180.8, 200, 200.8, 220, 220.8, 240, 240.8, 260, 260.8, 300,
          300.8, 360, 360.8, 420, 420.8, 480, 480.8, 540, 540.8, 600, 600.8, 660, 660.8, 720, 720.8, 780, 780.8, 840,
          840.8, 900, 900.8, 960, 960.8, 1020, 1020.8, 1080, 1080.8, 1140, 1140.8, 1200, 1200.8]

    sw.piecewiseConstant(lT, tT,y0w)
    resw = npqResults.NPQResults(sw)
    Tw, Fw, Fmaxw, _,_,_,_,_ = resw.fluo()
    plt.plot(Tw+10, Fw/Fmaxw, 'k', linewidth=2, label='wt')

    sm.piecewiseConstant(lT, tT,y0m)
    resm = npqResults.NPQResults(sm)
    Tm, Fm, Fmaxm, _,_,_,_,_ = resm.fluo()
    plt.plot(Tm, Fm/Fmaxm, 'r', linewidth=2, label='npq4 mutant')

    leg = plt.legend()
    if leg:
        leg.draggable

plt.suptitle('Simulation', fontsize=18)

# ================================================== #
# Fig S3: simulated lumen pH values
# ================================================== #