import numpy as np
from functions2 import dataArray, Data
import matplotlib.pyplot as plt
from scipy.stats import linregress



data = Data("/home/syrocco/Documents/Data/diffusivity/phi_0.1367839freq_53T_0.085h_1.5077.dumpL")
allData = np.array(dataArray(loc))

phi = np.array([data.phi for data in allData])
allData = allData[np.argsort(phi)]

t = allData[0].ticks

MSD = np.zeros((len(allData), len(t)))
E = np.zeros((len(allData), len(t)))
phi = np.zeros((len(allData)))
for i, data in enumerate(allData):
    MSD[i] = data.D
    E[i] = data.E/data.N
    phi[i] = data.phi

np.save('MSD.npy', [t, MSD, phi, E])

arr = np.load("MSD.txt")
t = arr[0]
MSD = arr[1]
phi = arr[2]
Earr = arr[3]


# time after which we start to compute the observables
where = t > 4

N = len(MSD)

D = np.zeros(N)
D_err = np.zeros(N)
E = np.zeros(N)
E_err = np.zeros(N)

for i in range(N):
    res = linregress(t[where], MSD[i][where])
    D[i] = res[0]/2
    D_err[i] = res[-2]/2
    E[i] = np.mean(Earr[i][where])
    E_err[i] = np.std(Earr[i][where])


fig, ax1 = plt.subplots()

ax1.set_xlabel(r'$\phi$')
ax1.set_ylabel('E')
ax1.plot(phi, E)

ax2 = ax1.twinx()  

ax2.set_ylabel('D')
ax2.plot(phi, D)

fig.tight_layout()  
plt.show()