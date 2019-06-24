# This file is a part of QSim: Quantum computer Simulation code.
# Copyright (c) 2019 and later
# Authors: Binho Le 
# All rights reserved.

import numpy as np
from math import *
from random import random

def ide(n):
    # to gererate an n x n identity matrix
    return np.identity(n)

def ran():
    # to generate a random number [0,1]
    return random()

def dag(s):
    # to generate conjugate transpost (dagger)
    return s.conjugate().T

def kdel(i,j):
    # to generate kronecker delta
    return 1 if i==j else 0

def gtrace(a,b):
    # to get trace of a-b
    return np.trace(np.abs(a-b))/2.0
