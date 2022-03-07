# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 09:22:57 2022

@author: Lara Bohnenblust
"""

import numpy as np

def imstep(deriv,r,p,dt,m1,m2):
    """
    Does one step with the implicit midpoint method up to machine precision

    Parameters
    ----------
    deriv : Equation of motion: dH(r,p,m1,m2): return dr/dt, dp/dt
    z : array(r,p)
    dt : time-step
    m1,m2 : Masses of the two CO or mu(reduced mass) and m(total mass)

    Returns
    -------
    zprime : new z = (r,p)
    len(hist) : Number of steps it took to reach machine precision

    """
    rprime,pprime = r,p
    hist = []
    while True:
        dr,dp = deriv((r+rprime)/2,(p+pprime)/2,m1,m2)
        dr,dp = dr*dt, dp*dt
        rem = np.sum(abs(dr)) +np.sum(abs(dp))  # or any kind of hash
        if rem in hist:  #checks if we have reached machine precision
            break
        hist.append(rem)
        rprime = r + dr
        pprime = p + dp
    return rprime,pprime,len(hist)

def imstepadaptive(deriv,r,p,ds,m1,m2,phi):
    """
    Does one step with the implicit midpoint method up to machine precision,
    but does an adaptive timestep!
    See arXiv: 0906.2226

    Parameters
    ----------
    deriv : Equation of motion: dH(r,p,m1,m2): return dr/dt, dp/dt
    z : array(r,p)
    ds : time-step
    m1,m2 : Masses of the two CO or mu(reduced mass) and m(total mass)

    Returns
    -------
    zprime : new z = (r,p)
    len(hist) : Number of steps it took to reach machine precision

    """
    rprime,pprime = r,p
    hist = []
    while True:
        dr,dp = deriv((r+rprime)/2,(p+pprime)/2,m1,m2)
        dr,dp = dr*ds/phi, dp*ds/phi
        rem = np.sum(abs(dr)) +np.sum(abs(dp))  # or any kind of hash
        if rem in hist:  #checks if we have reached machine precision
            break
        hist.append(rem)
        rprime = r + dr
        pprime = p + dp
    return rprime,pprime,len(hist)


