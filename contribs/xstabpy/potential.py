#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


class lennardjones(object):
    ''' Lennard-Jones potential class
    '''

    def __init__(self, rcut=None, sigma=None, epsilon=None, zerocut=False):
        '''Constructor with args (rcut, sigma, epsilon). If zerocut is set to
        True, all energy values are shifted so that energy at rcut is
        equal to zero.
        '''
        self.rcut = rcut
        self.ecut = 0.
        self.sigma = sigma
        self.epsilon = epsilon
        if zerocut:
            self.ecut, _ = self.__call__(rcut)

    def __call__(self, r):
        ''' Energy function. Returns energy and its derivative as a tuple. Both 
        elements should have tyhe same size as r.
        '''
        ratio = self.sigma / r
        e = 4. * self.epsilon * (ratio**12 - ratio**6) - self.ecut
        de = -24. * self.epsilon * (2. * ratio**12 - ratio**6) / r
        return e, de


class suttonchen(object):
    ''' Sutton-Chen potential
    '''

    def __init__(self, rcut=None, c=None, epsilon=None, a0=None, n=None, m=None, zerocut=False):
        ''' Constructor with args (...)
        '''
        self.rcut = rcut
        self.phicut = 0.
        self.rhocut = 0.
        self.c = c
        self.epsilon = epsilon
        self.a0 = a0
        self.n = n
        self.m = m
        if zerocut:
            self.phicut, _ = self.phi(rcut)
            self.rhocut, _ = self.rho(rcut)

    def phi(self, r):
        ''' Phi function
        '''
        q = self.epsilon * (self.a0 / r) ** self.n
        dq = -1. * self.n * q / r
        q -= self.phicut
        return q, dq

    def rho(self, r):
        ''' Rho function
        '''
        q = (self.a0 / r) ** self.m
        dq = -1. * self.m * q / r
        q -= self.rhocut
        return q, dq

    def f(self, rho):
        ''' Fembed function
        '''
        q = -1. * self.c * self.epsilon * np.sqrt(rho)
        dq = 0.5 * q / rho
        return q, dq

    
def spline3(x):
    ''' Spline 3 function used in vniitf potential
    '''
    return x**4 * (-20. * x**3 + 70. * x**2 - 84. * x + 35.)


def dspline3(x):
    ''' Spline 3 function (derivative) used in vniitf potential
    '''
    return 140. * x**3 * (-1. * x**3 + 3. * x**2 - 3. * x + 1.)


class vniitf(object):
    ''' VNIITF potential
    '''

    def __init__(self, rcut=None, rmax=None, rmin=None, rt0=None, Ecoh=None,
                 E0=None, beta=None, A=None, Z=None, n=None, alpha=None,
                 D=None, eta=None, mu=None, zerocut=False):
        ''' Constructor with args (...)
        '''
        self.rcut = rcut
        self.phicut = 0.
        self.rhocut = 0.
        self.rmax = rmax
        self.rmin = rmin
        self.rt0 = rt0
        self.Ecoh = Ecoh
        self.E0 = E0
        self.beta = beta
        self.A = A
        self.Z = Z
        self.n = n
        self.alpha = alpha
        self.D = D
        self.eta = eta
        self.mu = mu
        if zerocut:
            self.phicut, _ = self.phi(rcut)
            self.rhocut, _ = self.rho(rcut)

    def phi(self, r):
        ''' Phi function
        '''
        ir = 1 / r
        irt0 = 1 / self.rt0
        dr = r * irt0 - 1.
        dr2 = dr * dr
        a = -2. * self.Ecoh / self.Z
        b = self.alpha**3 * self.D * self.rt0
        f1 = a * (1. + self.alpha * dr + self.eta * dr2 + (self.mu + b * ir) *
                  dr2 * dr)
        df1 = a * (self.alpha * irt0 + 2. * self.eta * irt0 * dr + 3. * self.mu *
                   irt0 * dr2 + b * (3. * irt0 - dr * ir) * dr2 * ir)
        f2 = np.exp(-self.alpha * dr)
        df2 = -1. * self.alpha * irt0 * f2
        drS = (self.rmax - r) / (self.rmax - self.rmin)
        S = spline3(drS)
        dS = dspline3(drS) / (self.rmin - self.rmax)
        if np.isscalar(r):
            if drS < 0.:
                S = 0.
                dS = 0.
            elif drS > 1.:
                S = 1.
                dS = 0.
        else:
            S[drS < 0.] = 0.
            dS[drS < 0.] = 0.
            S[drS > 1.] = 1.
            dS[drS > 1.] = 0.
        q = (self.E0 + f1*f2) *  S
        dq = (self.E0 + f1*f2) * dS + (f1*df2 + f2*df1) * S
        q -= self.phicut
        return q, dq

    def rho(self, r):
        ''' Rho function
        '''
        irt0 = 1. / self.rt0
        F = np.exp(-self.beta * (r * irt0 - 1.)) / self.Z
        dF = -self.beta * F * irt0
        drS = (self.rmax - r) / (self.rmax - self.rmin)
        S = spline3(drS)
        dS = dspline3(drS) / (self.rmin - self.rmax)
        if np.isscalar(r):
            if drS < 0.:
                S = 0.
                dS = 0.
            elif drS > 1.:
                S = 1.
                dS = 0.
        else:
            S[drS < 0.] = 0.
            dS[drS < 0.] = 0.
            S[drS > 1.] = 1.
            dS[drS > 1.] = 0.
        q  = F * S;
        dq = F * dS + S * dF;
        q -= self.rhocut
        return q, dq

    def f(self, rho):
        ''' Fembed function
        '''
        a = rho**self.n
        b = self.A * self.Ecoh * a
        c = np.log(a)
        q = b * (c - 1.)
        dq = self.n * b * c / rho
        if np.isscalar(rho) and rho < 0.:
            q = 0.
            dq = 0.
        else:
            q[rho < 0.] = 0.
            dq[rho < 0.] = 0.
        return q, dq
