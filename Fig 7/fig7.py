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
    dPdt = w*P*(1 - P) - np.exp(-w*τ)*w*P
    
    return np.array([dPdt, dTdt])

def makeData(loc, fr):
    def getE(E, maxLenght):
        if len(E) < maxLenght:
            return 0
        else:
            return np.mean(E[-int(maxLenght/5):])
        
    allData = dataArray(loc)
    E = [data.E for data in allData]
    lenght = [len(e) for e in E]
    maxLenght = max(lenght)
    print(E)
    Em = np.array([getE(e, maxLenght) for e in E])
   
    phi = np.array([data.phi for data in allData])
    arg = np.argsort(phi)
    phi = phi[arg]
    Em = Em[arg]
    
    with open(f"data{fr}.pkl", "wb") as f:
        pickle.dump((fr, phi, Em), f)
    return phi, Em

"""
makeData("/home/syrocco/Documents/Data/Fig 5/6/", 1/6)
makeData("/home/syrocco/Documents/Data/Fig 5/inf/", 0)
makeData("/home/syrocco/Documents/Data/Fig 5/7.7/", 1/7.7)
"""
#makeData("/home/syrocco/Documents/Data/Fig 5/100/", 1/100)
if 1:
    F = []
    PHI = []
    E = []
    for file in glob.glob("*.pkl"):
        with open(file, "rb") as f:
            fr, phi, e = pickle.load(f)
        F.append(fr)
        PHI.append(phi)
        E.append(e)
    
    F = np.array(F)
    PHI = np.array(PHI)
    E = np.array(E)
    arg = np.argsort(F)
    arg = arg[::-1]
    F = F[arg]
    PHI = PHI[arg]
    E = E[arg]
    
    Etheo = np.zeros(np.shape(E))
    
    ΔA = 0.05
    ΔS = 0
    a = 0.95
    g = 0.01
    Y0 = np.array([1, 0.23])
    
    size = len(F)
    Fcolor = F - np.min(F)
    Fcolor = Fcolor/np.max(Fcolor)
    
    for i in range(size):
        for j in range(len(phi)):
            p = PHI[0, j]
            t = 1/F[i]
            res = fsolve(system0, Y0,  args = (ΔS, ΔA, t, p, g, a), full_output = 1)
            if res[2] == 1:
                res = res[0][1]
            else:
                res = 0
            Etheo[i, j] = res
    cmap = matplotlib.colormaps['cool']
    
    alpha = np.ones(size)*0.25
    alpha[0] = 1
    alpha[-1] = 1
    
    fig = plt.figure(figsize=(25, 17))
    for i in range(size - 1):
        scatter = plt.scatter(PHI[i], E[i], color = cmap(Fcolor[i]), alpha = alpha[i])#, label = fr"$f_s = {F[i]:.2f}$")
        plt.plot(PHI[i], Etheo[i], color = cmap(Fcolor[i]))
    norm = plt.Normalize(F[-2], 0.25)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=None, orientation='vertical') 
    cbar.ax.set_title(r'$1/\tau$',fontsize=15)
    i = size - 1
    plt.plot(PHI[i], E[i], color = "black", alpha = alpha[i], label = r"Simple $\Delta$ model")#, label = fr"$f_s = {F[i]:.2f}$")
    plt.xlabel(r"$\phi$")
    plt.ylabel("E")
    plt.legend()
    