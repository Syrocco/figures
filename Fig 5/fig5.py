import numpy as np
import matplotlib.pyplot as plt
from functions2 import dataArray
import glob
import os
from scipy.optimize import fsolve

def wo(phi):
    g = (((1 + phi**2/8 - phi**4/10)/(1 - phi)**2) - 1)/(2*phi)
    return g*4*phi*np.sqrt(1/np.pi)

def meanFreePath(phi):
    return np.sqrt(np.pi/2)/wo(phi)

def whereDoesItBreak(delta, gamma):
    return fsolve(lambda x: meanFreePath(x) - delta/gamma, 0.1)[0]

def getEnergy(e, l):
    if len(e) == l:
        return np.mean(e[int(l/10):])
    return 0


def equation(v0, T, gamma, l):
    return 4*T - v0**2/np.log(1-gamma*l/v0)*((1-gamma*l/v0)**2 - 1)

def getV0asT(T, gamma, l):
    return fsolve(equation, 0.031, args = (T, gamma, l))[0]

def simpleFreq(T, gamma, phi):
    l = mfp(phi)
    v0 = getV0asT(T, gamma, l)/(2/np.sqrt(np.pi))
    return -gamma/np.log(1 - gamma*l/v0)

allData = dataArray("/home/syrocco/Documents/Data/Fig 5/data/")

ts = np.sort(np.array(list(set([data.T for data in allData]))))
gamma = np.sort(np.array(list(set([data.gamma for data in allData]))))
phi =  np.sort(np.array(list(set([data.phi for data in allData]))))
maxL = max([len(data.E) for data in allData])

lts = len(ts)
lg = len(gamma)
lp = len(phi)

E = np.zeros((lp, lts, lg))
Phi_c = np.zeros((lts, lg))

for data in allData:
    TS = np.where(data.T == ts)[0][0]
    GAMMA = np.where(data.gamma == gamma)[0][0]
    PHI = np.where(data.phi == phi)[0][0]
    E[PHI, TS, GAMMA] = getEnergy(data.E, maxL)




phi_critical = 5e-6
phi_c = phi[(E > phi_critical).argmax(axis = 0)]
plt.figure(figsize = (15, 10))
plt.subplot(121)
color = ["C0", "C1", "C2", "C3", "C4"]
marker = ["o", "x", "<", "s", "p"]
for i in range(5):
    for j in range(5):
        plt.scatter(phi, E[:, i, j], c = color[i], marker=  marker[j], s = 15)

plt.title(r"symbol $\propto\gamma$ AND color $\propto \tau_s$")
plt.xlabel(r"$\phi$")
plt.ylabel("E/?")
plt.subplot(122)
plt.imshow(phi_c, extent = [min(gamma), max(gamma), min(ts), max(ts)], origin = 'lower', cmap = "magma", aspect = (np.max(gamma) - np.min(gamma))/(np.max(ts) - np.min(ts)))
plt.xlabel(r"$\gamma$")
plt.ylabel(r"$\tau_s$")
clb=plt.colorbar()
clb.ax.set_title(r'$\phi_c$',fontsize = 15)

plt.figure(figsize = (15, 10))
plt.subplot(121)
for i in range(5):
    for j in range(5):
        phic = whereDoesItBreak(0.03, gamma[j])
        phi_c[i, j] -= phic
        plt.scatter(phi - phic, E[:, i, j], c = color[i], marker=  marker[j], s = 15)

plt.title(r"symbol $\propto\gamma$ AND color $\propto \tau_s$")
plt.xlabel(r"$\phi$")
plt.ylabel("E/?")
plt.subplot(122)
plt.imshow(phi_c, extent = [min(gamma), max(gamma), min(ts), max(ts)], origin = 'lower', cmap = "magma", aspect = (np.max(gamma) - np.min(gamma))/(np.max(ts) - np.min(ts)))
plt.xlabel(r"$\gamma$")
plt.ylabel(r"$\tau_s$")
clb=plt.colorbar()
clb.ax.set_title(r'$\phi_c$',fontsize = 15)

        





if 0:
    
    a = glob.glob("/home/syrocco/Documents/Data/Fig 7/data/*v_2*")
    for string in a:
        temp = string.split("v_2")
        string2 = temp[0] + "v_1" + temp[1]
        os.remove(string2)
    