import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.optimize import curve_fit


def dissipationRateFunction(t, E):
    keep = int(len(t)/5)
    a, b = np.shape(t)[:-1]

    dissipationRate = np.zeros((a, b))
    for i in range(a):
        for j in range(b):
            dissipationRate[i, j] = curve_fit(linear, t[i, j, keep:], energyForDissipation[i, j, keep:])[0][0]
    return dissipationRate


def linear(x, a, b):
    return a*x + b

with open(r"lammps data/data.pkl", "rb") as f:
    _, syncVariance, Ez, sync, Et = pickle.load(f)
with open(r"lammps data/dataTimed.pkl", "rb") as f:
    t, syncKuramotoWithTime = pickle.load(f)
    
with open(r"lammps data/dataEnergyTimed.pkl", "rb") as f:
    tEnergy, energyForDissipation = pickle.load(f)

E = np.mean(Et[:, :, :, -10:], axis = 3)
syncKuramoto = np.mean(syncKuramotoWithTime[:, :, -100:], axis = 2)
dissipationRate = dissipationRateFunction(tEnergy, energyForDissipation)



"""
PURELY DEAD STATE (NO XY VELOCITY):
    TAB[heigh, amp]
    
    syncKuramotoWithTime = synchronization according to theta for a given (h, a) according to time
    
    syncKuramoto = synchronization according to theta for a given (h, a)
    syncVariance = synchronizeation accorgind to variance for a given (h, a)
    Ez = Energy in the z direction for a given (h, a)
    ________________
    
DISSIPATION RATE:
    tEnergy = time of simulation
    energyForDissipation = evolution of energy of non interacting particle according to tEnergy
    dissipationRate = rate of dissipation for a given (h, a)
    
ACTIVE STATE
    TAB[phi, heigh, amp]

    Et = Energy over time at a given (phi, h, a)
    E = Energy of the steady state for a given (phi, h, a)
    sync = dynamics synchronization according to theta for a given (phi, h, a)
"""






phi = np.array([0.05, 0.06, 0.07, 0.08, 0.09, 0.1 , 0.11, 0.12, 0.13, 0.14, 0.15,
       0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28,
       0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4 , 0.41,
       0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5 , 0.51, 0.52,
       0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6 , 0.61, 0.62, 0.63,
       0.64])

amp = np.array([0.045, 0.047, 0.049, 0.051, 0.053, 0.055, 0.057, 0.059, 0.061,
       0.063, 0.065, 0.067, 0.069, 0.071, 0.073, 0.075, 0.077, 0.079,
       0.081, 0.083, 0.085, 0.087, 0.089, 0.091, 0.093, 0.095, 0.097,
       0.099, 0.101, 0.103, 0.105, 0.107, 0.109, 0.111, 0.113, 0.115,
       0.117, 0.119, 0.121])

h = np.array([1.2   , 1.2615, 1.3231, 1.3846, 1.4462, 1.5077, 1.5692, 1.6308,
       1.6923, 1.7538, 1.8154, 1.8769, 1.9385, 2. ])



PHI = {k: v for v, k in enumerate(phi)}
AMP = {k: v for v, k in enumerate(amp)}
H = {k: v for v, k in enumerate(h)}

plt.figure()
plt.imshow(E[18], extent = [min(amp), max(amp), min(h), max(h)], origin = 'lower', cmap = "magma", aspect = (np.max(amp) - np.min(amp))/(np.max(h) - np.min(h)))
plt.xlabel(r"amp")
plt.ylabel(r"height")
plt.title(rf"$E$ at $\phi = {phi[18]}$")
clb=plt.colorbar()
clb.ax.set_title(r'$E$',fontsize=15)

plt.figure()
plt.scatter(phi, E[:, H[1.5077], AMP[0.085]])
plt.xlabel(r"$\phi$")
plt.ylabel(r"E")

plt.figure()
plt.imshow(Ez, extent = [min(amp), max(amp), min(h), max(h)], origin = 'lower', cmap = "magma", aspect = (np.max(amp) - np.min(amp))/(np.max(h) - np.min(h)))
plt.xlabel(r"amp")
plt.ylabel(r"height")
plt.title(r'$E_z$ in the dead state')
clb=plt.colorbar()
clb.ax.set_title(r'$E$',fontsize=15)

plt.figure()
phi_critical = 5e-6
phi_c = phi[(E > phi_critical).argmax(axis = 0)]
phi_c[phi_c == np.min(phi)] = np.nan
plt.imshow(phi_c, extent = [min(amp), max(amp), min(h), max(h)], origin = 'lower', cmap = "magma", aspect = (np.max(amp) - np.min(amp))/(np.max(h) - np.min(h)))
plt.xlabel(r"amp")
plt.ylabel(r"height")
clb=plt.colorbar()
clb.ax.set_title(r'$\phi_c$',fontsize = 15)

plt.figure()
plt.imshow(syncKuramoto, extent = [min(amp), max(amp), min(h), max(h)], origin = 'lower', cmap = "magma", aspect = (np.max(amp) - np.min(amp))/(np.max(h) - np.min(h)))
plt.xlabel(r"amp")
plt.ylabel(r"height")
clb=plt.colorbar()
clb.ax.set_title(r'$sync$',fontsize = 15)

plt.figure()
plt.plot(amp, syncKuramoto[H[1.5077]], label = "sync param")
plt.plot(amp, phi_c[H[1.5077]], label = "critical_phi" )
plt.legend()
plt.xlabel("amp")
plt.title("h = 1.5077")
plt.ylabel(r"sync or $\phi_c$")


plt.figure()
plt.plot(t, syncKuramotoWithTime[H[1.5077], AMP[0.085]], label = "async")
plt.plot(t, syncKuramotoWithTime[H[2], AMP[0.085]], label = "sync")
plt.legend()
plt.xlabel("time")
plt.ylabel("Synchronization")

plt.figure()
plt.imshow(dissipationRate, extent = [min(amp), max(amp), min(h), max(h)], origin = 'lower', cmap = "magma", aspect = (np.max(amp) - np.min(amp))/(np.max(h) - np.min(h)))
plt.xlabel(r"amp")
plt.ylabel(r"height")
plt.title(r'dissipation rate')
clb=plt.colorbar()
clb.ax.set_title(r'$E$',fontsize=15)
