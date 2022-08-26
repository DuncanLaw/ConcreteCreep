# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 13:28:40 2022

@author: Duncan.Law
"""

import numpy as np 
import matplotlib.pyplot as plt

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

### Calculate phi
def findPhi(time):
    if time <= tau:
        return 0
    else:
        beta_bcf = 1.8/(f_cm)**0.7 # FIB code eq 5.1-65
        beta_bct = np.log((30/t_0adj + 0.035)**2 * (time - tau) + 1) # FIB code eq 5.1-66
        phi_bc = beta_bcf * beta_bct  #**basic creep coefficient
        
        beta_dcf = 412/(f_cm)**1.4 # FIB code eq 5.1-68
        beta_RH = (1 - RH/100)/(0.1 * h/100)**(1/3) # FIB code eq 5.1-69
        beta_dct0 = 1/(0.1 + t_0adj**0.2) # FIB code eq 5.1-70
        gamma_t0 = 1/(2.3 + 3.5 / (t_0adj)**0.5)
        alpha_fcm = (35/f_cm)**0.5
        beta_h = 1.5*h + 250*alpha_fcm
        beta_dct = ((t_i - tau)/(beta_h + (time - tau)))**gamma_t0 # FIB code eq 5.1-71a
        phi_dc = beta_dcf * beta_RH * beta_dct0 * beta_dct # dry creep coefficient
        return phi_bc + phi_dc

for i, t_i in enumerate(t):
    phi[i] = findPhi(t_i)
### 

# ### Calculate C
def findC(time):
    for i, t_i in enumerate(time):
        for j in range(m):
            for idx, Gamma in enumerate(Gamma_m):
                if idx == j:
                    if t_i <= tau:
                        A[i][j] = 0
                    else:
                        A[i][j] = 1 - np.exp(-(t_i - tau)/Gamma)
    a = np.linalg.lstsq(A, phi, rcond=None)[0]
    return A@a

C = findC(t)
##     


## Plot the analytical solution against the FIB code
plt.xlabel('time (days)')
plt.ylabel('creep coefficient')
plt.plot(t, C, 'r', label='numerical')
plt.plot(t, phi, 'bx', label='FIB')
plt.legend(loc='lower right')
plt.show()
##

## Compare asymptotes
print(f'The asymptotic creep coefficient of the numerical solution is {findC([asympTime])[-1][0]}')
print(f'The asymptotic creep coefficient of the FIB solution is {findPhi(asympTime)}')
###