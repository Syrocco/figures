import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.optimize import curve_fit

with open(r"continuous delta model/data.pkl", "rb") as f:
    t, evolution_of_E, evolution_of_nactive, fluctuation_E, fluctuation_nactive, time_to_reach_AS = pickle.load(f)
    
phi = np.array([0.1495, 0.1496, 0.1497, 0.1498, 0.1499, 0.15  , 0.1501, 0.1502,
       0.1503, 0.1504, 0.1505, 0.1506, 0.1507, 0.1508, 0.1509, 0.151 ,
       0.1511, 0.1512, 0.1513, 0.1514, 0.1515, 0.1516, 0.1517, 0.1518,
       0.1519, 0.152 , 0.1521, 0.1522, 0.1523, 0.1524, 0.1525, 0.1526,
       0.1527, 0.1528, 0.1529, 0.153 ])

N = np.array([10000, 10769, 11538, 12307, 13076, 13846, 14615, 15384, 16153,
       16923, 17692, 18461, 19230, 20000, 20769, 21538, 22307, 23076,
       23846, 24615, 25384, 26153, 26923, 27692, 28461, 29230, 30000])

plt.figure()
num = -4
for i in range(0, len(phi) - 10, 2):
    plt.plot(t, evolution_of_nactive[i, num], label = phi[i])
plt.xlabel("time, did I have once good units for time? No..")
plt.ylabel("n. of active particles")
plt.title(f"N = {N[num]}")
plt.xscale("log")
plt.yscale("log")
plt.legend()

plt.figure()
for i in range(0, len(N), 2):
    plt.plot(phi, time_to_reach_AS[:, i], label = f"N = {N[i]}")
plt.xlabel("phi")
plt.ylabel("time to reach AS")
plt.xscale("log")
plt.yscale("log")
plt.legend()



plt.figure()
for i in range(0, len(N) -2, 3):
    plt.subplot(121)
    plt.plot(t, evolution_of_nactive[3, i], label = N[i])
    plt.xlabel("time, did I have once good units for time? No..")
    plt.ylabel("n. of active particles")
    plt.title(f"phi = {phi[18]}")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()

    plt.subplot(122)
    plt.plot(t, evolution_of_nactive[25, i], label = N[i])
    plt.xlabel("time, no good units? have i ever?")
    plt.ylabel("n. of active particles")
    plt.title(f"phi = {phi[25]}")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()

plt.figure()
plt.imshow(time_to_reach_AS, extent = [min(N), max(N), min(phi), max(phi)], origin = 'lower', cmap = "magma", aspect = (np.max(N) - np.min(N))/(np.max(phi) - np.min(phi)))
plt.xlabel("N")
plt.ylabel(r"$\phi$")
plt.title("time to reach ss")
plt.colorbar()

plt.figure()
plt.imshow(fluctuation_nactive, extent = [min(N), max(N), min(phi), max(phi)], origin = 'lower', cmap = "magma", aspect = (np.max(N) - np.min(N))/(np.max(phi) - np.min(phi)))
plt.xlabel(r"N")
plt.ylabel(r"$\phi$")
plt.title(r"color = $N(\langle \rho_a^2\rangle - \langle \rho_a\rangle^2)$")
plt.colorbar()

plt.figure()
plt.imshow(fluctuation_E, extent = [min(N), max(N), min(phi), max(phi)], origin = 'lower', cmap = "magma", aspect = (np.max(N) - np.min(N))/(np.max(phi) - np.min(phi)))
plt.xlabel(r"N")
plt.ylabel(r"$\phi$")
plt.title(r"color = $N(\langle E^2\rangle - \langle E\rangle^2)$")
plt.colorbar()


def linearL(x, a, b, c):
    return a - b*np.log(x - c)

def power(x,a, b, c):
    return a*(x - c)**(-b)



def critical(x, y):
    where = np.logical_not(np.isnan(y))
    a, b = curve_fit(linearL, x[where], np.log(y[where]), p0 = (-6, 1.5, 0.149))
    return a

plt.figure()
plt.plot(phi, fluctuation_nactive[:, -7])
plt.xscale('log')
plt.yscale("log")
plt.xlabel(r"$\phi$")
plt.ylabel(r"variance of active particles")

plt.figure()
plt.plot(phi, fluctuation_E[:, -7])
plt.xscale('log')
plt.yscale("log")
plt.xlabel(r"$\phi$")
plt.ylabel(r"variance of energy BY particle")

critical_arr = []
phic = []
plt.figure()
for i in range(6, 10):
    plt.plot(phi, fluctuation_nactive[:, -i], label = N[-i])
    setdata = critical(phi, fluctuation_nactive[:, -i])
    critical_arr.append(setdata[1])
    phic.append(setdata[2])
plt.xscale('log')
plt.yscale("log")
plt.xlabel(r"$\phi$")
plt.ylabel(r"variance of active particles")
phi_c = np.mean(np.array(phic))
critical_val = np.mean(np.array(critical_arr))
print(f"phi_c = {phi_c} and exponent = {critical_val}")