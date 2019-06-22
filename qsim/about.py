# This file is part of QSim: Quantum computor Simulation code.
# Copyright (c) 2019 and later
# Authors: Binho Le 
# All rights reserved.

"""
Command line output of information on QSim and dependencies.
"""

__all__ = ['about']

import sys
import os
import platform
import numpy
import scipy
import inspect
import qsim.version
from qutip.hardware_info import hardware_info

def about():
    """
    About box for QSim. Gives version numbers for
    QSim, NumPy, SciPy, Cython, and MatPlotLib.
    """
    print("")
    print("QSim: Quantum computor Simulation code")
    print("Copyright (c) 2019 and later.")
    print("Authors: Binho Le ")
    print("")
    print("QSim Version:       %s" % qsim.__version__)
    print("Numpy Version:      %s" % numpy.__version__)
    print("Scipy Version:      %s" % scipy.__version__)
    try:
        import Cython
        cython_ver = Cython.__version__
    except:
        cython_ver = 'None'
    print("Cython Version:     %s" % cython_ver)
    try:
        import matplotlib
        matplotlib_ver = matplotlib.__version__
    except:
        matplotlib_ver = 'None'
    print("Matplotlib Version: %s" % matplotlib_ver)
    print("Python Version:     %d.%d.%d" % sys.version_info[0:3])
    print("Number of CPUs:     %s" % hardware_info()['cpus'])
    #print("BLAS Info:          %s" % _blas_info())
    #print("OPENMP Installed:   %s" % str(qutip.settings.has_openmp))
    #print("INTEL MKL Ext:      %s" % str(qutip.settings.has_mkl))
    print("Platform Info:      %s (%s)" % (platform.system(),
                                           platform.machine()))
    #qsim_install_path = os.path.dirname(inspect.getsourcefile(qsim))
    #print("Installation path:  %s" % qsim_install_path)
    print("")

if __name__ == "__main__":
    about()
