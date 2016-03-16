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
import parameters
import npqModel
import simulate
import npqResults
import misc as misc
import lightProtocol

p = parameters.NPQmodelParameterSet()
col = misc.load_colours()

PFDs = [20, 50, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 400, 600, 900]

T = np.linspace(0,500,1000)
y0 = np.array([0, p.bH*misc.pHinv(p.pHstroma), p.APtot/2, 1, 1., 0])

# Y = [Q, pH, PQ, Fs]

Pref = [p.gamma0, p.gamma1, p.gamma2, p.gamma3]
Xref = np.zeros([len(PFDs), 4])

for i in range(len(PFDs)):

    m = npqModel.NPQModel(p, 2)
    s = simulate.Sim(m)
    s.clearResults()

    l = lightProtocol.LightProtocol({'prot':'const', 'PFD':PFDs[i]})

    s.timeCourse(l, T, y0)

    res = npqResults.NPQResults(s)

    T, F, Fmax, Tm, Fm, Ts, Fs, Bst = res.fluo()

    Xref[i] = [m.quencher(res.getVar(3), res.getVar(4))[-1], misc.pH(res.getVar(1))[-1], res.getVar(0)[-1], F[-1]]

# change the value of gamma parameter by 1%
perturb = [1.01, 0.99]

# collect new calculated values for perturbated parameters for different light intensities
Xgl = np.zeros([len(Pref), 2, len(PFDs), 4])
Rgl = np.zeros([len(Pref), len(PFDs), 4])

for g in range(len(Pref)):
    p = parameters.NPQmodelParameterSet()
    m = npqModel.NPQModel(p, 2)
    for pert in range(len(perturb)):
        if g==0:
            m.par.gamma0 = Pref[g] * perturb[pert]
        elif g==1:
            m.par.gamma1 = Pref[g] * perturb[pert]
        elif g==2:
            m.par.gamma2 = Pref[g] * perturb[pert]
        else:
            m.par.gamma3 = Pref[g] * perturb[pert]

        for l in range(len(PFDs)):

            s = simulate.Sim(m)

            light = lightProtocol.LightProtocol({'prot':'const', 'PFD':PFDs[l]})

            s.timeCourse(light, T, y0)

            res = npqResults.NPQResults(s)

            T, F, Fmax, Tm, Fm, Ts, Fs, Bst = res.fluo()

            Xgl[g, pert, l, :] = [m.quencher(res.getVar(3), res.getVar(4))[-1],
                                  misc.pH(res.getVar(1))[-1],
                                  res.getVar(0)[-1],
                                  F[-1]]
            print Xgl[g, pert, l, :]

for g in range(len(Pref)):
    for l in range(len(PFDs)):
        Rgl[g,l,:] =  ((Xgl[g,0,l] - Xgl[g,1,l]) / Xref[l]) / ((Pref[g]  * perturb[0] - Pref[g]  * perturb[1]) / Pref[g])


fig, ax = plt.subplots(2, 2, sharex=True, sharey=True)
plt.suptitle('Effect of $\gamma$ parameters on key internal variables', fontsize=18)
cnt = 0

labels=['Q', 'pH', 'PQ', 'Fs']

for l in range(len(PFDs)):
    for j in range(len(labels)):
        ax[0][0].plot(PFDs, Rgl[0,:,j], color = col[j], label=labels[j] if l==0  else None)
ax[0][0].set_xticks([220, 320, 900])
ax[0][0].set_xticklabels([100,300,900])
ax[0][0].set_title('$\gamma_0$', fontsize=18)
ax[0][0].set_ylim(-1,1)
leg = ax[0][0].legend()
if leg:
    leg.draggable()

for l in range(len(PFDs)):
    for j in range(len(labels)):
        ax[0][1].plot(PFDs, Rgl[1,:,j], color = col[j], label=labels[j] if l==0 else None)
ax[0][1].set_xticks([220, 320, 900])
ax[0][1].set_xticklabels([100,300,900])
ax[0][1].set_title('$\gamma_1$', fontsize=18)

for l in range(len(PFDs)):
    for j in range(len(labels)):
        ax[1][0].plot(PFDs, Rgl[2,:,j], color = col[j], label=labels[j] if l==0 else None)
ax[1][0].set_xticks([220, 320, 900])
ax[1][0].set_xticklabels([100,300,900])
ax[1][0].set_title('$\gamma_2$', fontsize=18)

for l in range(len(PFDs)):
    for j in range(len(labels)):
        ax[1][1].plot(PFDs, Rgl[3,:,j], color = col[j], label=labels[j] if l==0 else None)
ax[1][1].set_xticks([220, 320, 900])
ax[1][1].set_xticklabels([100,300,900])
ax[1][1].set_title('$\gamma_3$', fontsize=18)
