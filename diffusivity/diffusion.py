import numpy as np
from functions2 import dataArray, Data
import matplotlib.pyplot as plt
from scipy.stats import linregress



"""
def computeMSD(data):
    dump = data.dump
    n = dump.nframes
    
    dump.jump_to_frame(0)
    xi = dump.get_atompropf("xu")
    yi = dump.get_atompropf("yu")
    MSD = []
    for i in range(n):
        dump.jump_to_frame(i)
        MSD.append(np.mean((dump.get_atompropf("xu") - xi)**2 + (dump.get_atompropf("yu") - yi)**2))
    return np.array(MSD)
    

allData = np.array(dataArray("/home/syrocco/Documents/Data/diffusivity/temp/"))

phi = np.array([data.phi for data in allData])
allData = allData[np.argsort(phi)]

t = allData[0].ticks
t = t - t[0]

MSD = np.zeros((len(allData), len(t)))
E = np.zeros((len(allData), len(t)))
phi = np.zeros((len(allData)))
for i, data in enumerate(allData):
    MSD[i] = computeMSD(data)
    E[i] = data.E/data.N
    phi[i] = data.phi

np.savez('MSD.npz', t = t, MSD = MSD, phi = phi, E = E)
"""

arr = np.load("MSD.npz")
t = arr['t']
MSD = arr['MSD']
phi = arr['phi']
Earr = arr["E"]




N = len(MSD)

D = np.zeros(N)
D_err = np.zeros(N)
E = np.zeros(N)
E_err = np.zeros(N)

for i in range(N):
    res = linregress(t, MSD[i])
    D[i] = res[0]/2
    D_err[i] = res[-1]/2
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
ax2.errorbar(phi, D, yerr = D_err, fmt = "o", color = c2)
ax2.tick_params(axis='y', labelcolor = c2)

fig.tight_layout()  
plt.show()