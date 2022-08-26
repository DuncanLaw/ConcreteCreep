# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 13:33:19 2022

@author: Duncan.Law
"""

import numpy as np 

### Initialise paramters
# environmental parameters as per CreepShrinkageEC2
h = 300 # effective thickness in mm
RH = 50 # rel. humidity
t_s = 7 # days after casting that drying commences, normally at the end of curing
fck = 40 # cylinder strength of the concrete at 28 days in MPa

f_cm = 48
t_0adj = 50

m = 9 # number of retardation times
Gamma_m = [0.01, 0.1, 1, 5, 10, 50, 100, 500, 1000] # retardation times (9)
#tau_0 = np.linspace(0.025, 9125, 50) # loading ages (50)
tau = 40 # trial loading age - if we vary tau and t would this not be a 3D graph?
n = 60 # number of observation times
t = np.linspace(0.05, 9125, n) # observation times
###

### Initialise empty zero matrices to be filled
A = np.zeros((n, m))
a = np.zeros((m, 1))
b = np.zeros((n, 1))
C = np.zeros((n, 1)) # numerical creep coefficient vector
phi = np.zeros((n, 1)) # FIB code (analytical) creep coefficeint vector
asympTime = 10e100
###

tau_0 = np.linspace(0.025, 9125, 50) # loading ages (50)
output = np.zeros((n, 1))
for tau_i in tau_0:
    for i, t_i in enumerate(t):
        for j in range(m):
            for idx, Gamma in enumerate(Gamma_m):
                if idx == j:
                    if t_i <= tau_i:
                        A[i][j] = 0
                    else:
                        A[i][j] = 1 - np.exp(-(t_i - tau_i)/Gamma)
    a = np.linalg.lstsq(A, phi, rcond=None)[0]
    C_ans = A@a
    output = np.c_[output, C_ans]
np.savetxt('creepOutput.csv', output[:,1:], delimiter=',')