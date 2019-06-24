#!/usr/bin/env python
"""
QSim: A Quantum computor Simulation code

QSim is open-source software for simulating a quantum computer. 
The QSim library depends on the  Numpy and Scipy numerical packages. 
QSim aims to provide user-friendly and efficient numerical
simulations of a wide variety of quantum computer circuits, including those
with quantum state tomography.
"""
import os
import sys

from setuptools import setup

# all information about QSim is here
MAN = 1
SUB = 0
SUBS= 1
VERSION = '%d.%d.%d' % (MAN,SUB,SUBS)
REQUIRES = ['numpy (>=1.8)', 'scipy (>=0.15)', 'cython (>=0.21)']
PACKAGES = ['qsim', 'qsim/dsm']

NAME = "qsim"
AUTHOR = ("Binho Le")
AUTHOR_EMAIL = ("binho@kindai.ac.jp")
LICENSE = "GNU"
DESCRIPTION = "A Quantum computor Simulation code"
KEYWORDS = "quantum computor simulation"
URL = "http://"
PLATFORMS = ["Linux", "Mac OSX", "Unix", "Windows"]

def write_version_py(filename='qsim/version.py'):
    cnt = """\
# this file is generated from QSim setup.py
version = '%(version)s'
"""
    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION})
    finally:
        a.close()

# if exists then remove
if os.path.exists('qsim/version.py'):
    os.remove('qsim/version.py')

write_version_py()

setup(name = NAME,
      version = VERSION,
      description=DESCRIPTION,
      url=URL,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      packages=PACKAGES,
      zip_safe=False)

