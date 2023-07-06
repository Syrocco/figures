import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from functools import partial
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def exp(x, a, maxim, minim):
    return maxim*np.exp(-x/a) + minim

def wo(phi):
    g = (((1 + phi**2/8 - phi**4/10)/(1 - phi)**2) - 1)/(2*phi)
    return g*4*phi*np.sqrt(1/np.pi)

def meanFreePath(phi):
    return np.sqrt(np.pi/2)/wo(phi)

def mfp(phi):
    return (np.sqrt(np.pi/2))/freq(1, phi)

def equation(v0, T, gamma, l):
    return 4*T - v0**2/np.log(1-gamma*l/v0)*((1-gamma*l/v0)**2 - 1)

def getV0asT(T, gamma, l):
    return fsolve(equation, 0.031, args = (T, gamma, l))[0]

def simpleFreq(T, gamma, phi):
    l = mfp(phi)
    v0 = getV0asT(T, gamma, l)/(2/np.sqrt(np.pi))
    return -gamma/np.log(1 - gamma*l/v0)


def freq(T, phi):
    g = (((1 + phi**2/8 - phi**4/10)/(1 - phi)**2) - 1)/(2*phi)
    return g*4*phi*np.sqrt(1/np.pi)*np.sqrt(T)

def tempMax(a, delta): # == theoDelta(gamma = 0)
    return np.pi*a**2/(4*(1-a**2)**2)*(1 + np.sqrt(1 + (4*(1-a**2))/(np.pi*a**2)))**2*delta**2

def theoDelta(phi, delta = 0.03, gamma = 0.01, a = 0.95):
    w = wo(phi)
    temp = a*delta*np.sqrt(np.pi) - 4*gamma/w
    return ((+temp + np.sqrt(temp**2 + (4*delta**2*(1-a**2))))/(2*(1-a**2)))**2

def whereDoesItBreak(delta, gamma):
    return fsolve(lambda x: meanFreePath(x) - delta/gamma, 0.1)[0]

def cutPhi(phi, E):
    loc = np.where(phi>0.6)[0][0]
    return E[:loc]


def Gdelta(T, phi, delta, gamma, a):
    w = simpleFreq(T, gamma, phi)
    return w/2*(delta**2 + a*delta*np.sqrt(np.pi*T) - T*(1-a**2)) - 2*gamma*T


def theo(phi, delta = 0.03, gamma = 0.01, a = 0.95):
    def find0Delta(phi, delta = 0.03, gamma = 0.01, a = 0.95):
        res, _, b, _ = fsolve(Gdelta, 0.2, args = (phi, delta, gamma, a), full_output = True)
        if b == 1:
            return res[0]
        return 0
    return np.array([find0Delta(p, delta, gamma, a) for p in phi])

def getE(data, maxl):
    
    E = data.E
    if len(E) < maxl:
        return 0
    b = np.mean(E[int(maxl/10):])
    return b

fig, ax = plt.subplots(1, 1, figsize = (15, 10))

"""
alpha = 0.95 gamma = 0.01 

delta1 = 0.03
delta2 = 0.035
delta3 = 0.025 

"""

t = np.linspace(0, 20000, 2000)
phi = np.loadtxt(r"phi.txt")
E1 = np.loadtxt(r"E1.txt")
E3 = np.loadtxt(r"E3.txt")
E2 = np.loadtxt(r"E2.txt")
E1 = cutPhi(phi, E1)
E2 = cutPhi(phi, E2)
E3 = cutPhi(phi, E3)
phi = cutPhi(phi, phi)
ME1 = np.array([np.mean(E[-100:]) for E in E1])/tempMax(0.95, 0.03)
ME2 = np.array([np.mean(E[-100:]) for E in E2])/tempMax(0.95, 0.035)
ME3 = np.array([np.mean(E[-100:]) for E in E3])/tempMax(0.95, 0.025)
"""
plt.scatter(phi, ME3, label = r"$\Delta/\beta = 2.5$", s = 7)
plt.scatter(phi, ME1, label = r"$\Delta/\beta = 3.0$", marker = "x", s = 7)
plt.scatter(phi, ME2, label = r"$\Delta/\beta = 3.5$", marker = "s", s = 7)
plt.xlabel(r"$\phi$")
plt.ylabel(r"$T_{ss}/T_{ss}(\gamma = 0)$")
plt.legend()
plt.grid()

plt.figure()"""


plt.xlabel(r"$\phi$")
plt.ylabel(r"$T_{ss}/T_{ss}(\gamma = 0)$")
plt.yscale("log")




ME1t = theoDelta(phi, delta = 0.03)/tempMax(0.95, 0.03)
ME2t = theoDelta(phi, delta = 0.035)/tempMax(0.95, 0.035)
ME3t = theoDelta(phi, delta = 0.025)/tempMax(0.95, 0.025)

ME1t2 = theo(phi, delta = 0.03)/tempMax(0.95, 0.03)
ME2t2 = theo(phi, delta = 0.035)/tempMax(0.95, 0.035)
ME3t2 = theo(phi, delta = 0.025)/tempMax(0.95, 0.025)


ax.scatter(phi, ME3, label = r"sim $\Delta/\beta = 2.5$", s = 7, c = "C0")
ax.scatter(phi, ME1, label = r"sim $\Delta/\beta = 3.0$", marker = "x", s = 7, c = "C1")
ax.scatter(phi, ME2, label = r"sim $\Delta/\beta = 3.5$", marker = "s", s = 7, c = "C2")

ax.plot(phi, ME3t, label = r"simple Theory", linestyle = "-", c = "C0")
ax.plot(phi, ME1t, linestyle = "-", c = "C1")
ax.plot(phi, ME2t, linestyle = "-", c = "C2")

ax.plot(phi, ME3t2, label = "corrected theory", linestyle = "--",  c = "C0")#, label = r"corrected theory $\Delta/\beta = 2.5$", c = "black")
ax.plot(phi, ME1t2, linestyle = "--", c = "C1")#, label = r"corrected theory $\Delta/\beta = 3.0$", c = "black")
ax.plot(phi, ME2t2, linestyle = "--", c = "C2")#, label = r"corrected theory $\Delta/\beta = 3.5$", c = "black"

ax.vlines(whereDoesItBreak(0.025, 0.01), 0, 1, colors = "black", linestyles="--", label = "Mean free path argument")
ax.vlines(whereDoesItBreak(0.03, 0.01), 0, 1, colors = "black", linestyles="--")
ax.vlines(whereDoesItBreak(0.035, 0.01), 0, 1, colors = "black", linestyles="--")

plt.legend()

axins3 = inset_axes(ax, width="40%", height="65%", loc=4, borderpad=3)

plt.xlabel(r"$\phi$")
plt.ylabel(r"$T_{ss}/T_{ss}(\gamma = 0)$")

axins3.scatter(phi[::3], ME3[::3], label = r"sim $\Delta/\beta = 2.5$", s = 7, c = "C0")
axins3.scatter(phi[::3], ME1[::3], label = r"sim $\Delta/\beta = 3.0$", marker = "x", s = 7, c = "C1")
axins3.scatter(phi[::3], ME2[::3], label = r"sim $\Delta/\beta = 3.5$", marker = "s", s = 7, c = "C2")

axins3.plot(phi, ME3t, label = r"theory $\Delta/\beta = 2.5$", c = "C0")
axins3.plot(phi, ME1t, label = r"theory $\Delta/\beta = 3.0$", c = "C1")
axins3.plot(phi, ME2t, label = r"theory $\Delta/\beta = 3.5$", c = "C2")


"""
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