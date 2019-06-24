# This file is a part of QSim: Quantum computer Simulation code.
# Copyright (c) 2019 and later
# Authors: Binho Le 
# All rights reserved.

import numpy as np
from qsim.dsm.dsmStrong import execu as exeStr
from qsim.dsm.dsmWeak import execu as exeWeak
from qsim.dsm.dsmProb import execu as exeProb

class execute:
   # to run the program
   
   def __init__(self,state,pstate,fid,niter,nsamp,theta,pHist):
      self.state = state
      self.pstate = pstate
      self.fid = fid
      self.niter = niter
      self.nsamp = nsamp
      self.theta = theta
      self.pHist = pHist

   def job(self,name):
      if name == 'strong':
         model = exeStr(self.state,self.pstate,\
                 self.fid,self.niter,self.nsamp,\
                 self.theta,self.pHist)
      if name == 'weak':
         model = exeWeak(self.state,self.pstate,\
                 self.fid,self.niter,self.nsamp,\
                 self.theta,self.pHist)
      if name == 'prob':
         model = exeProb(self.state,self.pstate,\
                 self.fid,self.niter,self.nsamp,\
                 self.theta,self.pHist)
      return model

   def info(self):
      print("This is to run a DSM code")
