  # -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 08:28:36 2022

@author: Lara Bohnenblust
"""

#Define functions for different numbers of Arguments

from numpy import dot, sqrt, pi, arcsinh
import scipy.constants as const

#Constants:
Msolar = 1.9885E30

#Gravitational Constant in solar Masses
G = const.G*Msolar
g=G
c = const.c
pc = 3.085677581E16
    
def build_functionHPN0(filename):
    with open(filename, 'r') as f:
        eqn = f.read().lower().strip()
        exec("def fcn(r,p,mu,m):\n return ({})".format(eqn))
        return locals()['fcn'] 
    
def build_functionHPN(filename):
    with open(filename, 'r') as f:
        eqn = f.read().lower().strip()
        exec("def fcn(r,p,mu,m,nu):\n return ({})".format(eqn))
        return locals()['fcn'] 
    
def build_functiondHPN(filename):
    with open(filename, 'r') as f:
        eqn = f.read().lower().strip()
        exec("def fcn(r,p,nu):\n return ({})".format(eqn))
        return locals()['fcn'] 
    
def build_functionHPM(filename):
    with open(filename, 'r') as f:
        eqn = f.read().lower().strip()
        exec("def fcn(r,p,e1,e2,m,nu,sig,gamma,xi):\n return ({})".format(eqn))
        return locals()['fcn']  
    
    
    

def build_functiondHdrPM(filename):
    with open(filename, 'r') as f:
        eqn = f.read().lower().strip()
        exec("def fcn(r,p,m,nu,sig,gamma,xi):\n return ({})".format(eqn))
        return locals()['fcn']  
    
def build_functiondHdpPM(filename):
    with open(filename, 'r') as f:
        eqn = f.read().lower().strip()
        exec("def fcn(r,p,m1,m2,m,nu,e1,e2,ee,de1,de2,gamma,xi,sig):\n return ({})".format(eqn))
        return locals()['fcn']
    