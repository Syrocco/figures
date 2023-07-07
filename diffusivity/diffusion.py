import numpy as np
from functions2 import dataArray, Data
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def f(x, a):
    return a*x

"""
def computeMSD(data):
    dump = data.dump
    n = dump.nframes
    t = data.ticks - data.ticks[0]
    dump.jump_to_frame(0)
    xi = dump.get_atompropf("xu")
    yi = dump.get_atompropf("yu")
    MSD = []
    for i in range(n):
        dump.jump_to_frame(i)
        val = (dump.get_atompropf("xu") - xi)**2 + (dump.get_atompropf("yu") - yi)**2
        MSD.append(np.mean(val/(2*t[i])))
    return np.array(MSD)




allData = np.array(dataArray("/home/syrocco/Documents/Data/diffusivity/temp/"))

phi = np.array([data.phi for data in allData])
allData = allData[np.argsort(phi)]

t = allData[0].ticks
t = t - t[0]

D = np.zeros((len(allData), len(t)))
E = np.zeros((len(allData), len(t)))
phi = np.zeros((len(allData)))
for i, data in enumerate(allData):
    D[i] = computeMSD(data)
    E[i] = data.E/data.N
    phi[i] = data.phi

np.savez('MSD.npz', t = t, D = D, phi = phi, E = E)

"""
arr = np.load("MSD.npz")
t = arr['t']
D = arr['D']
phi = arr['phi']
Earr = arr["E"]




N = len(D)

Dav = np.zeros(N)
D_err = np.zeros(N)
E = np.zeros(N)
E_err = np.zeros(N)


w = t > 5

for i in range(N):
    Dav[i] = np.mean(D[i][w])
    D_err[i] = np.std(D[i][w])
    E[i] = np.mean(Earr[i])
    E_err[i] = np.std(Earr[i])


fig, ax1 = plt.subplots()

c1 = "black"
ax1.set_xlabel(r'$\phi$')
ax1.set_ylabel('E', color = c1)
ax1.errorbar(phi, E, yerr = E_err, color = c1, fmt = "o")
ax1.tick_params(axis='y', labelcolor = c1)

ax2 = ax1.twinx()  

c2 = "grey"
ax2.set_ylabel('D', color = c2)
ax2.errorbar(phi, Dav, yerr = D_err, fmt = "o", color = c2)
ax2.tick_params(axis='y', labelcolor = c2)

fig.tight_layout()  
plt.show()