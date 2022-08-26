# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 14:45:58 2022

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
n = 1000 # number of observation times
t = np.linspace(0.05, 400, n) # observation times
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


### Applying single constant step pressure to get  strain
def find_delEsp_i(delT, delSig, a_i, g_i, Gamma_i):
    expExpression = 1 - np.exp(- delT / Gamma_i)
    delEps = a_i * delSig - a_i * Gamma_i * expExpression * delSig/delT + g_i * expExpression # Manual eq (4.10-27)
    return delEps

# Initialise constant step app'd stress(es)
sigma1 = []
t_sigma1 = tau

for t_i in t:
    if t_i < t_sigma1:
        sigma1.append(0)
    else:
        sigma1.append(1)
### Test 1 
eps1 = np.zeros((n, ))
g = []
eps1_tot = np.zeros((n, ))

for i, t_i in enumerate(t):
    for j in range(m):
        for idx, Gamma in enumerate(Gamma_m):
            if idx == j:
                if t_i <= tau:
                    A[i][j] = 0
                else:
                    A[i][j] = 1 - np.exp(-(t_i - tau)/Gamma)
a = np.linalg.lstsq(A, phi, rcond=None)[0]

for j in range(m):
    Gamma = Gamma_m[j]
    for i, t_i in enumerate(t):
        delT = t[i] - t[i-1]
        if t_i < t_sigma1:
            eps1[i] = 0
            g.append(0)
        else:
            sigma_1 = sigma1[i - 1]
            sigma_2 = sigma1[i]
            delSigma = sigma_2 - sigma_1
            g_i_delT = np.exp(-delT/Gamma)*g[i - 1] + a[j][0]*(delSigma/delT)*Gamma*(1 - np.exp(-delT/Gamma))
            g.append(g_i_delT)
            delEps = find_delEsp_i(delT, delSigma, a[j][0], g[i - 1], Gamma)
            eps1[i] = (delEps + eps1[i - 1])
            eps1_tot[i] = eps1_tot[i] + eps1[i]
    g = []

plt.plot(t, eps1_tot, 'r')
plt.show()

np.savetxt('strain_test1.csv', eps1_tot, delimiter=',')
###