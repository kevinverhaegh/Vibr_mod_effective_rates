#%%
threads = 16
import os
os.environ["MKL_NUM_THREADS"] = str(threads)
os.environ["NUMEXPR_NUM_THREADS"] = str(threads)
os.environ["OMP_NUM_THREADS"] = str(threads)

import Functions as fn
import numpy as np
import CRUMPET
import matplotlib.pyplot as plt
import greenland_criterion as gc
import fc_mapping as fc
import scipy.interpolate as interpolate

Te = 1
ne = 1e14
iso_mass = 2
tau = 1e-4

crm = CRUMPET.Crumpet('input_fin.dat')
source = np.zeros(15)
# source[0] = 1e16



# Initializing the indices for different species
indx_X1 = np.append(0,np.arange(1,14))
indx_B1 = np.append(14,np.arange(15,29))
indx_C1 = np.append(29,np.arange(30,44))
indx_EF1 = np.append(44,np.arange(45,59))
indx_a3 = np.append(59,np.arange(60,74))
indx_c3 = np.append(74,np.arange(75,89))
indx_d3 = np.append(89,np.arange(90,104))
indx_H2Plus = len(crm.species)-2
indx_Hmin = len(crm.species)-1
# T_vibr = 1e0
# E_X1 = fn.get_energies('X1Sg')
# for i,p in enumerate(indx_X1):
#     source[i] = 1e18*np.exp(-E_X1[i]/T_vibr)

# crm.source[indx_X1] = source

# f_H2 = crm.steady_state(Te,ne,Ti=Te/iso_mass,plot=True,dt=True, n0=np.ones(len(crm.species)))

# n = np.zeros(len(crm.species))
# f_H2 = crm.solve_crm(tau,Te,ne,Ti=Te/iso_mass, gl=False, n=n)

#%%
Tg = 10000
k0 = 4.23e-15*(300/Tg)**(1/3)
d = 0.21*(Tg/300)**(0.5)
D1 = 0.236*(Tg/300)**0.25
D2 = 0.0572*(300/Tg)**(1/3)

test = np.zeros(5000)


n0 = 1e14

n=np.zeros(len(crm.species))
for i in range(100):

    # Calculate time step with n0=n
    f_H2 = crm.solve_crm(tau,Te,ne,Ti=Te/iso_mass, gl=False, n=n, densonly=True).y[:,-1]
    f = np.append(1,f_H2[indx_X1])

    # Calculate corresponding source vector
    S = np.zeros(15)
    for w in range(14):
        for v in range(w+1,14):
            rate = n0*f[v]*f[w+1]*(v+1)*(w+1)*k0*(1.5-0.5*np.exp(-d*(v-w)))*np.exp(D1*(v-w)-D2*(v-w)**2)
            S[v] += -rate
            S[w+1] += -rate
            S[v+1] += rate
            S[w] += rate

    # Update source in CRM
    # S = np.zeros(len(indx_X1)+1)
    crm.source[indx_X1] = S[1:]

    # Update density
    n = f_H2

    # test[i] = S[0]
    # plt.plot(f_H2[indx_d3])
    print(i)

#%%
crm.source[:]=0
f_H2_old = crm.steady_state(Te,ne,Ti=Te/iso_mass,plot=False,dt=True)

plt.figure()
plt.plot(f_H2_old[indx_d3],'o', label = 'Without VVE')
plt.plot(f_H2[indx_d3],'o', label = 'With VVE')
plt.legend()

plt.figure()
plt.plot(f_H2_old[indx_X1],'o', label = 'Without VVE')
plt.plot(f_H2[indx_X1],'o', label = 'With VVE')
plt.legend()



plt.show()
# %%
