import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from functools import partial

def exp(x, a, maxim, minim):
    return maxim*np.exp(-x/a) + minim


def tempMax(a, delta):
    return np.pi*a**2/(4*(1-a**2)**2)*(1 + np.sqrt(1 + (4*(1-a**2))/(np.pi*a**2)))**2*delta**2
"""

alpha = 0.95 gamma = 0.01 

delta1 = 0.03
delta2 = 0.035
delta3 = 0.025 

"""

t = np.linspace(0, 20000, 2000)
phi = np.loadtxt(r"delta model example/phi.txt")
E1 = np.loadtxt(r"delta model example/E1.txt")
E3 = np.loadtxt(r"delta model example/E3.txt")
E2 = np.loadtxt(r"delta model example/E2.txt")


ME1 = np.array([np.mean(E[-100:]) for E in E1])
ME2 = np.array([np.mean(E[-100:]) for E in E2])
ME3 = np.array([np.mean(E[-100:]) for E in E3])


plt.scatter(phi, ME3/tempMax(0.95, 0.025), label = r"$\Delta/\beta = 2.5$", s = 7)
plt.scatter(phi, ME1/tempMax(0.95, 0.03), label = r"$\Delta/\beta = 3.0$", marker = "x", s = 7)
plt.scatter(phi, ME2/tempMax(0.95, 0.035), label = r"$\Delta/\beta = 3.5$", marker = "s", s = 7)
plt.xlabel(r"$\phi$")
plt.ylabel(r"$T_{ss}/T_{ss}(\gamma = 0)$")
plt.legend()
plt.grid()

plt.figure()
time1 = np.array([curve_fit(partial(exp, maxim = E1[i][0], minim = np.mean(E1[i][-100:])), t, E1[i])[0] for i in range(300)])
time2 = np.array([curve_fit(partial(exp, maxim = E2[i][0], minim = np.mean(E2[i][-100:])), t, E2[i])[0] for i in range(300)])
time3 = np.array([curve_fit(partial(exp, maxim = E3[i][0], minim = np.mean(E3[i][-100:])), t, E3[i])[0] for i in range(300)])

plt.scatter(phi, time3*0.01, label = r"$\Delta/\beta = 2.5$", s = 7)
plt.scatter(phi, time1*0.01, label = r"$\Delta/\beta = 3.0$", marker = "x", s = 7)
plt.scatter(phi, time2*0.01, label = r"$\Delta/\beta = 3.5$", marker = "s", s = 7)
plt.xlabel(r"$\phi$")
plt.ylabel(r"$\chi\gamma$")
plt.legend()
plt.grid()


"""
def sortAccordingToA(A, *args):
    arg = np.argsort(A)
    # [:] allows to change by pointer the arrays (and not overwrite them), so no need to return them
    A[:] = A[arg]
    for arr in args:
        arr[:] = arr[arg]



allData1 = dataArray("/home/syrocco/Documents/absorbing/Data/simple delta/0.03/")
allData2 = dataArray("/home/syrocco/Documents/absorbing/Data/simple delta/0.035/")
allData3 = dataArray("/home/syrocco/Documents/absorbing/Data/simple delta/0.025/")

phi1 = np.array([data.phi for data in allData1])
E1 = np.array([data.E for data in allData1])
phi2 = np.array([data.phi for data in allData2])
E2 = np.array([data.E for data in allData2])
phi3 = np.array([data.phi for data in allData3])
E3 = np.array([data.E for data in allData3])

sortAccordingToA(phi1, E1)
sortAccordingToA(phi2, E2)
sortAccordingToA(phi3, E3)"""