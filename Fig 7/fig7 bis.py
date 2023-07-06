import numpy as np
import matplotlib.pyplot as plt
from functions2 import dataArray
from scipy.optimize import fsolve
import pickle
import glob
import matplotlib


def wo(phi, T):
    g = (((1 + phi**2/8 - phi**4/10)/(1 - phi)**2) - 1)/(2*phi)
    return g*4*phi*np.sqrt(1/np.pi)*np.sqrt(T)

def G(T, Δ, γ, a, w):
    return w/2*(Δ**2 + a*Δ*np.sqrt(np.pi*T) - T*(1-a**2)) - 2*γ*T

def system0(Y, ΔS, ΔA, τ, Φ, γ, a):
    P, T = Y[0], Y[1]
    w = wo(Φ, T)
    dTdt = (2*P - P**2)*G(T, ΔA, γ, a, w) + (1 - P)**2*G(T, ΔS, γ, a, w)
    dPdt = w*P*(1 - P) - np.exp(-w*τ)/τ*P
    
    return np.array([dPdt, dTdt])

def solve(p, t):
    ΔA = 0.1
    ΔS = 0
    a = 0.95
    g = 0.01
    Y0 = np.array([1, 1.3])
    res = fsolve(system0, Y0,  args = (ΔS, ΔA, t, p, g, a), full_output = 1)
    if res[2] == 1:
        return res[0][1]
    else:
        return 0

def getEnergy(e, l):
    if len(e) == l:
        return np.mean(e[int(l/10):])
    return 0


def theoDelta(phi, delta = 0.03, gamma = 0.01, a = 0.95):
    w = wo(phi, 1)
    temp = a*delta*np.sqrt(np.pi) - 4*gamma/w
    return ((+temp + np.sqrt(temp**2 + (4*delta**2*(1-a**2))))/(2*(1-a**2)))**2

allData = dataArray("/home/syrocco/Documents/Data/Fig 7/Data higher delta/")

ts = np.sort(np.array(list(set([data.T for data in allData]))))
phi =  np.sort(np.array(list(set([data.phi for data in allData]))))
maxL = max([len(data.E) for data in allData])

lts = len(ts)
lp = len(phi)

E = np.zeros((lp, lts))
Etheo = np.zeros((lp, lts))

        
        
for data in allData:
    TS = np.where(data.T == ts)[0][0]
    PHI = np.where(data.phi == phi)[0][0]
    E[PHI, TS] = getEnergy(data.E, maxL)
    Etheo[PHI, TS] = solve(data.phi, data.T)



for i in range(len(ts)):
    plt.scatter(phi, E[:, i])
    plt.plot(phi, Etheo[:, i])

plt.plot(phi, theoDelta(phi, 0.1))