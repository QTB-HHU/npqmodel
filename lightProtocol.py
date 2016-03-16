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


class LightProtocol:
    """ class containing protocols used for investigation:
        ldl: Light-Dark_Light
        lightmemory: Soma's protocol, without 30 s of darkness at the beginning
        doublelightmemory: Ioannis's protocol, with double light-dark cycle
    """


    def __init__(self, options={'PFD': 0, 'protocol': 'constant'}):

        PFD = options.setdefault('PFD', 0)
        Tmax = options.setdefault('Tmax', 1800)

        lightFnKey = 'lightintensity'
        timeFnKey = 'tfls'
        lightchangeFnKey = 'pfds'

        proto = options.setdefault('protocol', 'const')

        if proto.startswith('user'):
            # user defined light function
            # keys in options:
            # 'LightFn': <user_defined_function>
            _userDefinedLightFn = options.setdefault('lightFn', lambda t:0)
            setattr(self,lightFnKey,_userDefinedLightFn)

        elif proto.startswith('const'):
            """ constant light for the whole duration of the experiment """
            setattr(self, lightFnKey, lambda t: PFD)

        elif proto == 'ldl':
            """ simple light-dark-light protocol with no flashes
                requires Tof to switch off the light and Ton to switch it on again
            """
            self.Toff = options.setdefault('Toff', 600)
            self.Ton = options.setdefault('Ton', 1200)

            space = np.hstack([self.Toff, self.Ton, Tmax])

            light = np.hstack([PFD, 0.000001, PFD])

            setattr(self, timeFnKey, space)
            setattr(self, lightchangeFnKey, light)

            setattr(self,lightFnKey,
                     lambda t: PFD * ((t < self.Toff) or (t >= self.Ton)) + 0.0001 * ((t >= self.Toff) and (t < self.Ton))
                     )

        else:
            raise NameError('Protocol \''+options['protocol']+'\' not yet implemented.') 
