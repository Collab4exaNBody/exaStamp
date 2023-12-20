#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import yaml


def tabulate_pair(pot, r, addrcut=True):
    ''' Pair potentials tabulation
    -- 
    This function takes a pair potential (pot) and an array of
    positions (r), then returns a dictionnary containing data ready to
    be written as a yaml file
    '''
    if len(r.shape) > 1:
        NameError('Tabulate method only works with 1D input!')
    s = np.copy(r)
    if addrcut:
        s = np.append(s, pot.rcut)
    s = np.unique(s)
    e, de = pot(s)
    return {'r' : s.tolist(), 'e' : e.tolist(), 'de' : de.tolist()}


def tab_and_write_pair(pot, r, filename):
    ''' Tabulate potential and write file
    --
    This function takes a pair potential (pot), an array of positions
    (r) and a string (filename), then write a yaml file containing
    tabulated data meant to be used in exastamp
    '''
    data = tabulate_pair(pot, r)
    strm = open(filename, 'w')
    strm.write('# ' + 'Tabulated pair potential -- r (m), e (J), de (N)' + '\n')
    yaml.dump(data, strm, encoding=None)
    strm.close()

    
def tabulate_eam(pot, r, rhof, addrcut=True):
    ''' Eam potentials tabulation
    --
    This function takes an eam potential (pot) and an array of
    positions (r), then returns a dictionnary containing data ready to
    be written as a yaml file 
    '''
    if len(r.shape) > 1 or len(rhof.shape) > 1:
        NameError('Tabulate method only works with 1D input!')
    s = np.copy(r)
    if addrcut:
        s = np.append(s, pot.rcut)
    s = np.unique(s)
    t = np.copy(rhof)
    t = np.unique(t)
    qp, dqp = pot.phi(s)
    qr, dqr = pot.rho(s)
    qf, dqf = pot.f(t)
    return {'r' : s.tolist(),
            'phi' : qp.tolist(),
            'dphi' : dqp.tolist(),
            'rho' : qr.tolist(),
            'drho' : dqr.tolist(),
            'rhof' : t.tolist(),
            'f' : qf.tolist(),
            'df' : dqf.tolist()}


def tab_and_write_eam(pot, r, s, filename):
    ''' Tabulate potential and write file
    --
    This function takes an eam potential (pot), an array of positions
    (r) and a string (filename), then write a yaml file containing
    tabulated data meant to be used in exastamp
    '''
    data = tabulate_eam(pot, r, s)
    strm = open(filename, 'w')
    strm.write('# ' + 'Tabulated EAM potential -- r (m), rho (none), ' + 
               'drho (m^-1), phi (J), dphi (N), rhof (none), f (J), df (N)' +
               '\n')
    yaml.dump(data, strm, encoding=None)
    strm.close()

