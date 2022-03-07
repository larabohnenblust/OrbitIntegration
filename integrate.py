# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 08:51:49 2022

@author: Lara Bohnenblust
"""
import numpy as np
from numpy.linalg import norm
from numpy import dot, sqrt, pi,arcsinh
import scipy.constants as const
from sympstep import imstep,imstepadaptive
from math import atan2

from MakeFunction import build_functionHPN0, build_functionHPN,build_functiondHPN,\
    build_functionHPM,build_functiondHdrPM,build_functiondHdpPM

#Constants:
Msolar = 1.9885E30

#Gravitational Constant in solar Masses
G = const.G*Msolar
g=G
c = const.c
pc = 3.085677581E16

def transform(x1,x2,v1,v2,m1,m2):
    """(For the PN integration)
    The function turns the Cartesian coordinates into a form suitable for the
    PN-Hamiltonian
    INPUT:
        x1,x2: The two positions of the objects as d-vector in Meters
        v1,v2: The two velocities of the objects as d-vectors in m/s
        m1,m2: The masses of the objects in solar masses.
    OUTPUT:
        r: relative reduced distance
        p: relative momentum
        mu: reduced mass
        M: total mass
    See 1805.07240, Equation (6.4)
    """
    m = m1 + m2
    mu = m1*m2/m
    
    r = (x1-x2)/(G*m)
    vcom = 1/m*(m1*v1 + m2*v2)
    v1new = v1-vcom
    #v2new = v2-vcom
    p = m1*v1new/mu
    return r,p,mu,m

def CartesianCoord(r,p,m1,m2):
    """
    Turns Relative Coordinates back into Cartesian coordinates

    Parameters
    ----------
    r : Distance between objects
    p : c.o.m. momentum
    m1,m2 : masses 

    Returns
    -------
    x1,x2 : Cartesian Coordinates of two objects
    v1,v2 : velocities of the two objects
    """
    mu = m1*m2/(m1+m2)
    x1 = mu/m1*r
    x2 = -mu/m2*r
    v1 = p/m1
    v2 = -p/m2
    return x1,x2,v1,v2

def com(x1,x2,v1,v2,m1,m2):
    """
    Turns Cartesian coordinates into com coordinates

    Parameters
    ----------
    x1, x2, v1, v2 : Cartesian coordinates (array) and velocities.
    m1, m2 : Masses of the two CO.

    Returns
    -------
    r, p : relative coordinate and c.o.m. momentum.
    mu, m : reduced mass and total mass.

    """
    r = x1-x2
    m = m1+m2
    vcom = 1/m*(m1*v1+m2*v2)
    v1new = v1-vcom
    p = m1*v1new
    mu = m1*m2/m
    return r,p,mu,m

def expeuler(deriv,phi,r,p,ds,m1,m2):
    """
    Does one step with explicit Euler step for phi

    Parameters
    ----------
    deriv : right-hand side of ode (returns dr/dt, dp/dt)
    phi : function for stretching of timestep
    r : initial position
    p : initial momentum
    ds : non-stretched timestep
    m1,m2 : masses of the two compact objects

    Returns
    -------
    
        DESCRIPTION.

    """
    dr,_ = deriv(r,p,m1,m2)
    dphi = - dot(dr,r)/dot(r,r)
    return phi + dphi*ds

#Let psi = r0/r
def dxlnpsi(r):
    return -r/dot(r,r)

#r0,p0,t0,tend,N,deriv,Ham,m1,m2
def adaptiveimplicitmidpoint(r0,p0,t0,tend,N,deriv,Ham,m1,m2):
    if(tend <= t0 or N<1):
        print("Enter valid time interval and number of timesteps!")
        return 0
    r = [r0]
    p = [p0]
    ds = (tend-t0)/N
    phi_init = 1
    phi = expeuler(deriv,phi_init,r0,p0,ds/2,m1,m2)
    #phi_array = [phi]
    T = [t0]
    H= [Ham(r0,p0,m1,m2)]
    dt = 1/phi*ds
    T.append(T[-1]+dt)
    while (T[-1]<tend):
        rk,pk,n = imstepadaptive(deriv,r[-1],p[-1],ds,m1,m2,phi)
        H.append(Ham(rk,pk,m1,m2))
        r.append(rk)
        p.append(pk)
        phi = expeuler(deriv,phi,rk,pk,ds,m1,m2)
        dt = 1/phi*ds
        T.append(T[-1]+dt)
        #phi_array.append(phi)
    T[-1] = tend
    dt=tend-T[-2]
    rk,pk,n = imstep(deriv,r[-2],p[-2],dt,m1,m2)
    H.append(Ham(rk,pk,m1,m2))
    r.append(rk)
    p.append(pk)
    return np.array(r).T,np.array(p).T,np.array(H),np.array(T)


"""EoM from Hamiltonians:
    dr/dt = dH/dp
    dp/dt = -dH/dr
    
    dH = (dH/dp, -dH/dr)
"""

"""
The following Hamiltoniasn are imported from a Mathematica Output.
At different orders, they have differnt dependencies on the variables ...
"""
#Post-Newtonian: From 1805.07240, page 48

#H0PN(r,p,mu,m)
H0PN = build_functionHPN0("E0PN.txt")

#H1PN_temp(r,p,mu,m,nu)
H1PN_temp = build_functionHPN("E1PN.txt")
def H1PN(r,p,mu,m):
    nu = mu/m
    return H1PN_temp(r,p,mu,m,nu)

#H2PN_temp(r,p,mu,m,nu)
H2PN_temp = build_functionHPN("E2PN.txt")
def H2PN(r,p,mu,m):
    nu = mu/m
    return H2PN_temp(r,p,mu,m,nu)

#H3PN_temp(r,p,mu,m,nu)
H3PN_temp = build_functionHPN("E3PN.txt")
def H3PN(r,p,mu,m):
    nu = mu/m
    return H3PN_temp(r,p,mu,m,nu)

"""
Derivatives of the Hamiltonians:
    Imported from txt-files generated in Mathematica
"""


#0PN
dHdp0PN = build_functiondHPN("dHdp0PN.txt")
dHdr0PN = build_functiondHPN("dHdr0PN.txt")
def dH0PN(r,p,mu,m):
    nu = mu/m
    dr = dHdp0PN(r,p,nu)
    dp = -dHdr0PN(r,p,nu)
    return dr,dp

#1PN
dHdp1PN = build_functiondHPN("dHdp1PN.txt")
dHdr1PN = build_functiondHPN("dHdr1PN.txt")
def dH1PN(r,p,mu,m):
    nu = mu/m
    dr = dHdp1PN(r,p,nu)
    dp = -dHdr1PN(r,p,nu)
    return dr,dp

#2PN
dHdp2PN = build_functiondHPN("dHdp2PN.txt")
dHdr2PN = build_functiondHPN("dHdr2PN.txt")
def dH2PN(r,p,mu,m):
    nu = mu/m
    dr = dHdp2PN(r,p,nu)
    dp = -dHdr2PN(r,p,nu)
    return dr,dp

#3PN
dHdp3PN = build_functiondHPN("dHdp3PN.txt")
dHdr3PN = build_functiondHPN("dHdr3PN.txt")
def dH3PN(r,p,mu,m):
    nu = mu/m
    dr = dHdp3PN(r,p,nu)
    dp = -dHdr3PN(r,p,nu)
    return dr,dp

"""
The same for the PM-Hamiltonians
"""
    
#Post-Minkowskian, Bern et. al
#PM Hamiltonians

#0PM(This is actually 0PN ... as it is expanded in v)
def H0PM(r,p,m1,m2):
    return (m1+m2)*c**2+norm(p)**2/(2*m1)+norm(p)**2/(2*m2)-m1*m2*G/norm(r)
    

#1PM
#H1PM_temp(r,p,E1,E2,m,nu,sig,gamma,xi)
H1PM_temp = build_functionHPM("E1PM.txt")
def H1PM(r,p,m1,m2):
    E1 = np.sqrt(norm(p)**2*c**2+m1**2*c**4)
    E2 = np.sqrt(norm(p)**2*c**2+m2**2*c**4)
    E = E1+E2
    m = m1+m2
    nu = m1*m2/m**2
    sig = 1/(m1*m2*c**2)*(E1*E2/c**2+norm(p)**2)
    gamma = E/(m*c**2)
    xi = E1*E2/E**2
    return H1PM_temp(r,p,E1,E2,m,nu,sig,gamma,xi)




#2PM
#H2PM_temp(r,p,E1,E2,m,nu,sig,gamma,xi)
H2PM_temp = build_functionHPM("E2PM.txt")
def H2PM(r,p,m1,m2):
    E1 = np.sqrt(norm(p)**2*c**2+m1**2*c**4)
    E2 = np.sqrt(norm(p)**2*c**2+m2**2*c**4)
    E = E1+E2
    m = m1+m2
    nu = m1*m2/m**2
    sig = 1/(m1*m2*c**2)*(E1*E2/c**2+norm(p)**2)
    gamma = E/(m*c**2)
    xi = E1*E2/E**2
    return H2PM_temp(r,p,E1,E2,m,nu,sig,gamma,xi)
    
#3PM
#H3PM_temp(r,p,E1,E2,m,nu,sig,gamma,xi)
H3PM_temp = build_functionHPM("E3PM.txt")
def H3PM(r,p,m1,m2):
    E1 = np.sqrt(norm(p)**2*c**2+m1**2*c**4)
    E2 = np.sqrt(norm(p)**2*c**2+m2**2*c**4)
    E = E1+E2
    m = m1+m2
    nu = m1*m2/m**2
    sig = 1/(m1*m2*c**2)*(E1*E2/c**2+norm(p)**2)
    gamma = E/(m*c**2)
    xi = E1*E2/E**2
    return H3PM_temp(r,p,E1,E2,m,nu,sig,gamma,xi)



#0PM
def dH0PM(r,p,m1,m2):
    dr = p/(m1+m2)
    dp = -m1*m2*G/norm(r)**3*r
    return dr,dp

#1PM
dHdp1PM = build_functiondHdpPM("dHdp1PM.txt")
dHdr1PM = build_functiondHdrPM("dHdr1PM.txt")
#1PM
def dH1PM(r,p,m1,m2):
    E1 = np.sqrt(norm(p)**2*c**2+m1**2*c**4)
    E2 = np.sqrt(norm(p)**2*c**2+m2**2*c**4)
    dE1 = c**2/E1*p
    dE2 = c**2/E2*p
    E = E1+E2
    m = m1+m2
    nu = m1*m2/m**2
    sig = 1/(m1*m2*c**2)*(E1*E2/c**2+norm(p)**2)
    gamma = E/(m*c**2)
    xi = E1*E2/E**2
    dr = dHdp1PM(r,p,m1,m2,m,nu,E1,E2,E,dE1,dE2,gamma,xi,sig)
    dp = -dHdr1PM(r,p,m,nu,sig,gamma,xi)
    return dr,dp


#2PM
dHdp2PM = build_functiondHdpPM("dHdp2PM.txt")
dHdr2PM = build_functiondHdrPM("dHdr2PM.txt")
#2PM
def dH2PM(r,p,m1,m2):
    E1 = np.sqrt(norm(p)**2*c**2+m1**2*c**4)
    E2 = np.sqrt(norm(p)**2*c**2+m2**2*c**4)
    dE1 = c**2/E1*p
    dE2 = c**2/E2*p
    E = E1+E2
    m = m1+m2
    nu = m1*m2/m**2
    sig = 1/(m1*m2*c**2)*(E1*E2/c**2+norm(p)**2)
    gamma = E/(m*c**2)
    xi = E1*E2/E**2
    dr = dHdp2PM(r,p,m1,m2,m,nu,E1,E2,E,dE1,dE2,gamma,xi,sig)
    dp = -dHdr2PM(r,p,m,nu,sig,gamma,xi)
    return dr,dp

#3PM
dHdp3PM = build_functiondHdpPM("dHdp3PM.txt")
dHdr3PM = build_functiondHdrPM("dHdr3PM.txt")
#2PM
def dH3PM(r,p,m1,m2):
    E1 = np.sqrt(norm(p)**2*c**2+m1**2*c**4)
    E2 = np.sqrt(norm(p)**2*c**2+m2**2*c**4)
    dE1 = c**2/E1*p
    dE2 = c**2/E2*p
    E = E1+E2
    m = m1+m2
    nu = m1*m2/m**2
    sig = 1/(m1*m2*c**2)*(E1*E2/c**2+norm(p)**2)
    gamma = E/(m*c**2)
    xi = E1*E2/E**2
    dr = dHdp3PM(r,p,m1,m2,m,nu,E1,E2,E,dE1,dE2,gamma,xi,sig)
    dp = -dHdr3PM(r,p,m,nu,sig,gamma,xi)
    return dr,dp
    
    

def integratePN(r1,r2,v1,v2,m1,m2,t0,tend,N, order = 0):
    #Go to com frame and scale r and p with G*M/mu
    r0,p0,mu,M = transform(r1,r2,v1,v2,m1,m2)
    #rescale time according to 1805.07240 (6.4)
    t0 = t0/(G*M)
    tend = tend/(G*M)
    if (order == 0):
        r,p,H,T = adaptiveimplicitmidpoint(r0,p0,t0,tend,N,dH0PN,H0PN,mu,M)
    elif (order == 1):
        r,p,H,T = adaptiveimplicitmidpoint(r0,p0,t0,tend,N,dH1PN,H1PN,mu,M)
    elif (order == 2):
        r,p,H,T = adaptiveimplicitmidpoint(r0,p0,t0,tend,N,dH2PN,H2PN,mu,M)
    elif (order == 3):
        r,p,H,T = adaptiveimplicitmidpoint(r0,p0,t0,tend,N,dH3PN,H3PN,mu,M)
    else:
        print("NO VALID PN ORDER")
        return 0
    #Scale Variables back
    r = r*G*M
    p = mu*p
    T = T*G*M
    r1_res, r2_res, v1_res,v2_res = CartesianCoord(r,p,m1,m2)
    return r1_res,r2_res,v1_res,v2_res,r,p,H,T


def integratePM(r1,r2,v1,v2,m1,m2,t0,tend,N, order=0):
    #Transform coordinates to the COM frame
    r0,p0,mu,M = com(r1,r2,v1,v2,m1,m2)
    if(order ==0):
        r,p,H,T = adaptiveimplicitmidpoint(r0,p0,t0,tend,N,dH0PM, H0PM,m1,m2)
    elif(order==1):
        r,p,H,T = adaptiveimplicitmidpoint(r0,p0,t0,tend,N,dH1PM, H1PM,m1,m2)
    elif(order==2):
        r,p,H,T = adaptiveimplicitmidpoint(r0,p0,t0,tend,N,dH2PM, H2PM,m1,m2)
    elif(order==3):
        r,p,H,T = adaptiveimplicitmidpoint(r0,p0,t0,tend,N,dH3PM, H3PM,m1,m2)

    else:
        print("NO VALID PM ORDER")
        return 0
    #Transform coordinates back to Cartesian Coordinates
    r1_res,r2_res,v1_res,v2_res = CartesianCoord(r,p,m1,m2)
    return r1_res,r2_res,v1_res,v2_res,r,p,H,T

def orbit(x1,x2,v1,v2,m1,m2,t0,tend,N,method = "PN", order = 1):
    """
    INPUT:
        x1,x2: Initial Positions as arrays in Cartesian Coordinates
        v1,v2: Initial Velocities as arrays
        m1,m2: Masses of the COs
        t0,tend: Initial and final time
        N: Number of steps
        method: PN for Post-Newtonian Expansion
                PM for Post-Minkowskian Expansion
        order: Order of PN/PM Expansion that one wants
    OUTPUT:
        r1,r2: Orbits in Cartesian Coordinates
        v1,v2: Velocity during Orbits in Cartesian Coordinates
        r: Relative Coordinate of Orbit
        p: c.o.m. momentum during orbit
        H: Array with energies at each point of T
        T: Time-array
    """
    if(method == "PN"):
        r1,r2,v1,v2,r,p,H,T = integratePN(x1,x2,v1,v2,m1,m2,t0,tend,N,order)
    elif(method == "PM"):   
        r1,r2,v1,v2,r,p,H,T = integratePM(x1,x2,v1,v2,m1,m2,t0,tend,N,order)
    else:
        print("NO VALID EXPANSION GIVEN!")
        return 0
    return r1,r2,v1,v2,r,p,H,T

def phase(x,y):
    """
    Parameters
    ----------
    r : d x N array. r[0,:] = x, r[1,:] = y

    Returns
    -------
    phase of points r[i,:]

    """
    
    temp = [atan2(a,b) for b,a in zip(x,y)]
    return np.array(temp)
    
    
    

