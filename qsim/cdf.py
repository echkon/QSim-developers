# This file is a part of QSim: Quantum computor Simulation code.
# Copyright (c) 2019 and later
# Authors: Binho Le  
# All rights reserved.

__all__ = ['cdf']

import numpy as np
from math import *
from qsim.utilities import ide,ran
import matplotlib.pyplot as plt

def cdf(f):
    """ comulative distribution function
    input: a function f
    return: average of original function

    """
    dim = np.int(len(f))
    nbin = 100
    binh,binb = np.histogram(f,nbin)
    
    s = 0.0 
    for i in range(np.int(nbin/2)):
         s += binh[i] 
    return s/(f.max()/nbin*dim)/(nbin/2)
