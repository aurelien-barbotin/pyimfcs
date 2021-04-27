# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 14:31:37 2021

@author: abarbotin
"""

import numpy as np

# amplitude 100 percent curve 1.16*10**-4

c = 20*10**-9   # mol/L
nav = 6.022*10**23 #mol-1
dx= 64*10**-9 #m, pixel ize
dz = 180*10**-9

Nbin = 4
dV = (Nbin*dx)**2*dz*10**3 #L

Nexpected = c*dV*nav
amplitude_expected = 1/Nexpected
print("N expected {}, amplitude expected {}".format(Nexpected, amplitude_expected))


D = 4 #um2/s 100 nm bead in water

D_expected = 4/10 #um2/s, factor 6 bc of sucrose
w = 0.2
tau = (dx*Nbin*10**6)**2/(8*np.log(2)*D_expected)
tauz = (dz*10**6)**2/(8*np.log(2)*D_expected)*10**3

