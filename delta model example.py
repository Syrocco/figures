import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from functools import partial

def exp(x, a, maxim, minim):
    return maxim*np.exp(-x/a) + minim

def wo(phi):
    g = (((1 + phi**2/8 - phi**4/10)/(1 - phi)**2) - 1)/(2*phi)
    return g*4*phi*np.sqrt(1/np.pi)

def meanFreePath(phi):
    return np.sqrt(np.pi/2)/wo(phi)


def tempMax(a, delta): # == theoDelta(gamma = 0)
    return np.pi*a**2/(4*(1-a**2)**2)*(1 + np.sqrt(1 + (4*(1-a**2))/(np.pi*a**2)))**2*delta**2

def theoDelta(phi, delta = 0.03, gamma = 0.01, a = 0.95):
    w = wo(phi)
    temp = a*delta*np.sqrt(np.pi) - 4*gamma/w
    return ((+temp + np.sqrt(temp**2 + (4*delta**2*(1-a**2))))/(2*(1-a**2)))**2

def whereDoesItBreak(delta, gamma):
    return fsolve(lambda x: meanFreePath(x) - delta/gamma, 0.1)[0]
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

ME1 = np.array([np.mean(E[-100:]) for E in E1])/tempMax(0.95, 0.03)
ME2 = np.array([np.mean(E[-100:]) for E in E2])/tempMax(0.95, 0.035)
ME3 = np.array([np.mean(E[-100:]) for E in E3])/tempMax(0.95, 0.025)

plt.scatter(phi, ME3, label = r"$\Delta/\beta = 2.5$", s = 7)
plt.scatter(phi, ME1, label = r"$\Delta/\beta = 3.0$", marker = "x", s = 7)
plt.scatter(phi, ME2, label = r"$\Delta/\beta = 3.5$", marker = "s", s = 7)
plt.xlabel(r"$\phi$")
plt.ylabel(r"$T_{ss}/T_{ss}(\gamma = 0)$")
plt.legend()
plt.grid()

plt.figure()

ME1t = theoDelta(phi, delta = 0.03)/tempMax(0.95, 0.03)
ME2t = theoDelta(phi, delta = 0.035)/tempMax(0.95, 0.035)
ME3t = theoDelta(phi, delta = 0.025)/tempMax(0.95, 0.025)

plt.scatter(phi, ME3, label = r"sim $\Delta/\beta = 2.5$", s = 7, c = "C0")
plt.scatter(phi, ME1, label = r"sim $\Delta/\beta = 3.0$", marker = "x", s = 7, c = "C1")
plt.scatter(phi, ME2, label = r"sim $\Delta/\beta = 3.5$", marker = "s", s = 7, c = "C2")

plt.plot(phi, ME3t, label = r"theory $\Delta/\beta = 2.5$", c = "C0")
plt.plot(phi, ME1t, label = r"theory $\Delta/\beta = 3.0$", c = "C1")
plt.plot(phi, ME2t, label = r"theory $\Delta/\beta = 3.5$", c = "C2")


plt.vlines(whereDoesItBreak(0.025, 0.01), 0, 1, colors = "C0", linestyles="--")
plt.vlines(whereDoesItBreak(0.03, 0.01), 0, 1, colors = "C1", linestyles="--")
plt.vlines(whereDoesItBreak(0.035, 0.01), 0, 1, colors = "C2", linestyles="--")


plt.xlabel(r"$\phi$")
plt.ylabel(r"$T_{ss}/T_{ss}(\gamma = 0)$")
plt.yscale("log")
plt.legend()

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
plt.figure()
time4 = np.array([curve_fit(partial(exp, maxim = E4[i][0], minim = np.mean(E1[i][-100:])), t[:50], E4[i][:50])[0] for i in range(len(E4))])
time5 = np.array([curve_fit(partial(exp, maxim = E5[i][0], minim = np.mean(E2[i][-100:])), t[:50], E5[i][:50])[0] for i in range(len(E5))])

plt.scatter(phi[:-1], time4*0.005, label = r"$\gamma = 0.005$", s = 7)
plt.scatter(phi, time1*0.01, label = r"$\gamma = 0.01$", marker = "x", s = 7)
plt.scatter(phi[:-2], time5*0.02, label = r"$\gamma = 0.02$", marker = "s", s = 7)
plt.xlabel(r"$\phi$")
plt.ylabel(r"$\chi\gamma$")
plt.legend()
plt.grid()
"""

plt.figure()
a1 = 0.1509
a3 = 0.1700
a2 = 0.1339
plt.scatter((phi - a3)/a3, ME3, label = r"$\Delta/\beta = 2.5$", s = 7)
plt.scatter((phi - a1)/a1, ME1, label = r"$\Delta/\beta = 3.0$", marker = "x", s = 7)
plt.scatter((phi - a2)/a2, ME2, label = r"$\Delta/\beta = 3.5$", marker = "s", s = 7)
plt.xlabel(r"$(\phi - \phi_c)/\phi_c$")
plt.ylabel(r"$T_{ss}/T_{ss}(\gamma = 0)$")
plt.legend()
plt.grid()
plt.title("$\phi_c$ = critical packing fraction obtained with mean free path argument")

plt.figure()
a1 = phi[np.argmax(time1)]
a3 = phi[np.argmax(time3)]
a2 = phi[np.argmax(time2)]
plt.scatter((phi - a3)/a3, ME3, label = r"$\Delta/\beta = 2.5$", s = 7)
plt.scatter((phi - a1)/a1, ME1, label = r"$\Delta/\beta = 3.0$", marker = "x", s = 7)
plt.scatter((phi - a2)/a2, ME2, label = r"$\Delta/\beta = 3.5$", marker = "s", s = 7)
plt.xlabel(r"$(\phi - \phi_c)/\phi_c$")
plt.ylabel(r"$T_{ss}/T_{ss}(\gamma = 0)$")
plt.legend()
plt.grid()
plt.title(r"$\phi_c$ = critical packing fraction obtained from max of time scale")




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
sortAccordingToA(phi3, E3)

allData2 = dataArray("/home/syrocco/Documents/absorbing/Data/simple delta/0.03 g=0.005/")
allData3 = dataArray("/home/syrocco/Documents/absorbing/Data/simple delta/0.03 g=0.02/")

phi2 = np.array([data.phi for data in allData2])
E2 = np.array([data.E for data in allData2])
phi3 = np.array([data.phi for data in allData3])
E3 = np.array([data.E for data in allData3])

sortAccordingToA(phi2, E2)
sortAccordingToA(phi3, E3)

np.savetxt("delta model example/E4.txt", E2)
np.savetxt("delta model example/E5.txt", E3)"""