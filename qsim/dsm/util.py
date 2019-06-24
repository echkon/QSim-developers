# This file is a part of QSim: Quantum computer Simulation code.
# Copyright (c) 2019 and later
# Authors: Binho Le 
# All rights reserved.

import numpy as np
from math import *
from random import random
from qsim.utilities import dag

def gtrace(a,b):
    # to get trace of a-b
    return np.trace(np.abs(a-b))/2.0

def gfide(a,b):
    # to get fidelity of a and original b
    return np.real(np.trace(np.dot(a,b)))
