# licensed to the apache software foundation (asf) under one
# or more contributor license agreements.  see the notice file
# distributed with this work for additional information
# regarding copyright ownership.  the asf licenses this file
# to you under the apache license, version 2.0 (the
# "license"); you may not use this file except in compliance
# with the license.  you may obtain a copy of the license at
# 
#   http://www.apache.org/licenses/license-2.0
# 
# unless required by applicable law or agreed to in writing,
# software distributed under the license is distributed on an
# "as is" basis, without warranties or conditions of any
# kind, either express or implied.  see the license for the
# specific language governing permissions and limitations
# under the license.

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import potential
import xstabpy


# number of points in tab files (for EAM tabulation, the number of
# points in r and s arrays does not need to be the same)
N = 1000



def filename(pot, n):
    '''Utility function to generate a filename. Arg pot must be a string,
       and n should be a positive interger.
    '''
    return 'tab_' + pot + '_' + str(n) + '.yaml'


def ljtab():
    ''' Lennard-Jones -- parameters values are for copper -- IS units
    -- bounds for array r are arbitrary
    '''

    rcut = 5.68E-10
    sigma = 2.27E-10
    epsilon = 9.34E-20

    p = potential.lennardjones(rcut, sigma, epsilon)
    r = np.linspace(1.0E-10, 1.1 * rcut, N)
    f = filename('lj', N)

    xstabpy.tab_and_write_pair(p, r, f)


def sctab():
    ''' Sutton-Chen -- parameters values are for copper -- IS units --
    bounds for array r and s are arbitrary
    '''

    rcut = 0.729E-09
    c = 3.317E+01
    epsilon = 3.605E-21
    a0 = 0.327E-09
    n = 9.050E+00
    m = 5.005E+00

    p = potential.suttonchen(rcut, c, epsilon, a0, n, m)
    r = np.linspace(1.0E-10, 1.1 * rcut, N)
    s = np.logspace(-3., +3., N)
    f = filename('sc', N)

    xstabpy.tab_and_write_eam(p, r, s, f)


def vniitftab():
    ''' Eam Vniitf -- parameters values are for tin -- IS units --
    bounds for array r and s are arbitrary
    '''

    rcut = 5.599E-10
    rmax = 5.599E-10
    rmin = 1.000E-10
    rt0 = 3.437E-10
    Ecoh = 2.956031e-19
    E0 = 5.15003855e-20
    beta = 6.000E+00
    A = 1.401E+00
    Z = 7.618E+00
    n = 0.724E+00
    alpha = 3.072E+00
    D = 1.450E-01
    eta = 2.720E+00
    mu = -1.87E+00

    p = potential.vniitf(rcut, rmax, rmin, rt0, Ecoh, E0, beta, A, Z, n, alpha, D, eta, mu)
    r = np.linspace(0.9 * rmin, 1.1 * rmax, N)
    s = np.logspace(-3., +3., N)
    f = filename('ev', N)

    xstabpy.tab_and_write_eam(p, r, s, f)


''' Run above functions
'''

ljtab()
sctab()
vniitftab()
