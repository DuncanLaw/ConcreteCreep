# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 09:44:15 2022

@author: Duncan.Law
"""

import numpy as np 
from scipy.optimize import nnls
import matplotlib.pyplot as plt

t1 = 1
t2 = 2
t3 = 3
phi1 = 1
phi2 = 3
phi3 = 7
tau = 0.5

### Test 1
# t = [t1, t2]
# A = np.array([[t1, 1], [t2, 1]])
# b = np.array([phi1, phi2])
# a = np.linalg.lstsq(A, b, rcond=None)[0]
# solution = A@a
# plt.plot(t, [phi1, phi2], 'bx')
# plt.plot(t, solution, 'rx')
# plt.show()

### Test 2
# t = [t1, t2, t3]
# A = np.array([[t1, 1], [t2, 1], [t3, 1]])
# b = np.array([phi1, phi2, phi3])
# a = np.linalg.lstsq(A, b, rcond=None)[0]
# solution = A@a
# plt.plot(t, [phi1, phi2, phi3], 'bx')
# plt.plot(t, solution, 'r')
# plt.show()

### Test 3
def findX(index, time):
    t = []