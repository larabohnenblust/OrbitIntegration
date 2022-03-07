# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 17:22:46 2022

@author: Lara Bohnenblust
"""

import numpy as np
from numpy.linalg import norm
from numpy import dot, sqrt,pi
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.constants as const
import astropy.constants as astroconst
from math import atan2
import time as time

from integrate import orbit,phase

c = const.c

def countorbits(phi):
    occ1 = phi[:-1]*phi[1:] < 0
    occ2 = phi[:-1]*phi[1:] == 0
    tot1 = occ1.sum()
    tot2 = occ2.sum()/2
    tot = (tot1 + tot2 -1)/2
    return tot

def evaluate(r01,r02,v01,v02,m1,m2,t0,tend,N):
    """
    m1 is the primary, m2 is the secondary!
    """
    start = time.time()
    r0PN1,r0PN2,v0PN1,v0PN2,r0PN,_,E0PN,T0PN = orbit(r01,r02,v01,v02,m1,m2,t0,tend,N, "PN",0)
    r1PN1,r1PN2,v1PN1,v1PN2,r1PN,_,E1PN,T1PN = orbit(r01,r02,v01,v02,m1,m2,t0,tend,N, "PN",1)
    r2PN1,r2PN2,v2PN1,v2PN2,r2PN,_,E2PN,T2PN = orbit(r01,r02,v01,v02,m1,m2,t0,tend,N, "PN",2)
    r3PN1,r3PN2,v3PN1,v3PN2,r3PN,_,E3PN,T3PN = orbit(r01,r02,v01,v02,m1,m2,t0,tend,N, "PN",3)
    r1PM1,r1PM2,v1PM1,v1PM2,r1PM,_,E1PM,T1PM = orbit(r01,r02,v01,v02,m1,m2,t0,tend,N, "PM",1)
    r2PM1,r2PM2,v2PM1,v2PM2,r2PM,_,E2PM,T2PM = orbit(r01,r02,v01,v02,m1,m2,t0,tend,N, "PM",2)
    r3PM1,r3PM2,v3PM1,v3PM2,r3PM,_,E3PM,T3PM = orbit(r01,r02,v01,v02,m1,m2,t0,tend,N, "PM",3)
    stop = time.time()
    print("Time to simulate orbits in [s]:", stop-start)
    
    v0PN1norm = np.linalg.norm(v0PN1,axis = 0)/c
    v1PN1norm = np.linalg.norm(v1PN1,axis = 0)/c
    v2PN1norm = np.linalg.norm(v2PN1,axis = 0)/c
    v3PN1norm = np.linalg.norm(v3PN1,axis = 0)/c
    v1PM1norm = np.linalg.norm(v1PM1,axis = 0)/c
    v2PM1norm = np.linalg.norm(v2PM1,axis = 0)/c
    v3PM1norm = np.linalg.norm(v3PM1,axis = 0)/c
    
    v0PN2norm = np.linalg.norm(v0PN2,axis = 0)/c
    v1PN2norm = np.linalg.norm(v1PN2,axis = 0)/c
    v2PN2norm = np.linalg.norm(v2PN2,axis = 0)/c
    v3PN2norm = np.linalg.norm(v3PN2,axis = 0)/c
    v1PM2norm = np.linalg.norm(v1PM2,axis = 0)/c
    v2PM2norm = np.linalg.norm(v2PM2,axis = 0)/c
    v3PM2norm = np.linalg.norm(v3PM2,axis = 0)/c
    
    phi0PN = phase(r0PN[0,:],r0PN[1,:])
    phi1PN = phase(r1PN[0,:],r1PN[1,:])
    phi2PN = phase(r2PN[0,:],r2PN[1,:])
    phi3PN = phase(r3PN[0,:],r3PN[1,:])
    phi1PM = phase(r1PM[0,:],r1PM[1,:])
    phi2PM = phase(r2PM[0,:],r2PM[1,:])
    phi3PM = phase(r3PM[0,:],r3PM[1,:])
    
    #Number of orbits
    no0PN = countorbits(phi0PN)
    no1PN = countorbits(phi1PN)
    no2PN = countorbits(phi2PN)
    no3PN = countorbits(phi3PN)
    no1PM = countorbits(phi1PM)
    no2PM = countorbits(phi2PM)
    no3PM = countorbits(phi3PM)
       
    rcart = [[r0PN1, r0PN2], [r1PN1, r1PN2],[r2PN1, r2PN2],[r3PN1, r3PN2],[r1PM1, r1PM2],[r2PM1, r2PM2],[r3PM1, r3PM2]]
    r =[r0PN,r1PN,r2PN,r3PN,r1PM,r2PM,r3PM]
    E = [E0PN,E1PN,E2PN,E3PN,E1PM,E2PM,E3PM]
    T = [T0PN,T1PN,T2PN,T3PN,T1PM,T2PM,T3PM]
    v1 = [v0PN1norm,v1PN1norm,v2PN1norm,v3PN1norm,v1PM1norm,v2PM1norm,v3PM1norm]
    v2 = [v0PN2norm,v1PN2norm,v2PN2norm,v3PN2norm,v1PM2norm,v2PM2norm,v3PM2norm]
    phi = [phi0PN,phi1PN,phi2PN,phi3PN,phi1PM,phi2PM,phi3PM]
    no = [no0PN,no1PN,no2PN,no3PN,no1PM,no2PM,no3PM]
    
    return rcart,r,E,T,v1,v2,phi,no