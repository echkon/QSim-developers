# This file is a part of QSim: Quantum computer Simulation code.
# Copyright (c) 2019 and later
# Authors: Binho Le 
# All rights reserved.

__all__ = ['ghz', 'w', 'dicke', 'qubit', 'randoms']

import numpy as np
from math import *
from qsim.utilities import ide,ran,dag

def ghz(n, fid=1.0):
    """ to generate GHZ state
    Parameters
    ----------
    n: number of qubits
    fid: default fidelity

    Return: GHZ state, 
    ie. (|00...0> + |11...1>)/sqrt(2)

    """
    dim = 2**n
    err = dim*(1.0-fid)/(dim-1.0)
    up,down = _up_down_()
    ups,upd = up,down
    for i in range(n-1):
        ups = np.kron(ups,up)
        upd = np.kron(upd,down)
    GHZ = (ups+upd)/sqrt(2.0)
    rho = np.dot(GHZ,GHZ.transpose())
    rho = (1.0-err)*rho + err*ide(dim)/dim
    return rho

def w(n,fid=1.0):
    """ to generate W state
    Parameters
    ----------
    n: number of qubits
    fid: default fidelity

    Return: W state, 
    ie. |00...0> +  

    """
    dim = 2**n
    err = dim*(1.0-fid)/(dim-1.0)
    W = np.zeros((dim,1))
    W[0] = 0.0
    W[1] = 1.0
    W[2] = 1.0
    W[3] = 0.0
    for i in range(3,n+1):
        W[2**(i-1)]= 1.0
    W = W/sqrt(n)
    rho = np.dot(W,W.transpose())
    rho = (1.0-err)*rho + err*ide(dim)/dim
    return rho

def dicke(n,e,fid=1.0):
    """ to generate Dicke state
    Parameters
    ----------
    n: number of qubits
    e: excited qubits
    fid: default fidelity

    Return: dicke state, 
    ie. |00...0> +  

    """
    dim = 2**n
    err = dim*(1.0-fid)/(dim-1.0)
    up,down = _up_down_()
    if e==0:
       temp_vec = np.zeros((dim,1))
       temp_vec[dim-1] = 1.0
    if n==e:
       temp_vec = np.zeros((dim,1))
       temp_vec[0] = 1.0
    if (e!=0) and (n!=e):
       temp_vec = up
       for i in range(1,e):
           temp_vec0 = np.kron(temp_vec,up)
           temp_vec = temp_vec0
       for i in range(1,n-e+1):
           temp_vec0 = np.kron(temp_vec,down)
           temp_vec = temp_vec0
    rho = temp_vec
    rho = (1.0-err)*rho + err*ide(dim)/dim
    return rho

def qubit():
    """ to generate a qubit |0>
    Return: |0>
    
    """
    up, down = _up_down_()
    return up

def _up_down_():
    # to generate up and down state
    up = np.zeros((2,1))
    up[0] = 1.0
    down = np.zeros((2,1))
    down[1] = 1.0
    return up,down

def randoms(d,fid=1.0):
    """ to generate a random state
    Parameters:
    ----------
    d: dimension
    fid: defaut fidelity

    Return: random state

    """
    err = d*(1.0-fid)/(d-1.0)
    rpsi = np.zeros((d,1))
    ipsi = np.zeros((d,1))
    for i in range (d):
        rpsi[i] = ran()
        ipsi[i] = ran()
    ppsi = rpsi + 1j*ipsi
    ppsi /= np.linalg.norm(ppsi)
    if fid == 1:
       return ppsi
    else:
       rho = np.dot(ppsi,dag(ppsi))
       rho = (1.0-err)*rho + err*ide(d)/d
       return rho
