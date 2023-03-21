import numpy as np
import matplotlib.pyplot as plt
from functions2 import Data, energy

"""
data = Data("/home/syrocco/Documents/LAMMPS/Script/phi_0.5000freq_53T_0.1h_1.8.dumpL")
E = energy(data)[1][120:2000]
data2 = Data("/home/syrocco/Documents/LAMMPS/Script/nice2/phi_0.5000freq_53T_0.1h_1.8.dumpL")
E2 = energy(data2)[1][:50000]

v = []
dump = data2.dump
n = dump.nframes
for j in range(0, 50000, 10):
    dump.jump_to_frame(j)
    v.append(dump.get_atompropf("vx")**2 +  dump.get_atompropf("vy")**2)
np.array(v)  
"""



"""
h = 1.8
amp = 0.1
"""

E = np.loadtxt(r"collifing particles k and mu/E first run.txt")
E2 = np.loadtxt(r"collifing particles k and mu/E second run.txt")
v = np.loadtxt(r"collifing particles k and mu/v second run.txt")
plt.xlabel("timestep")
plt.ylabel("energy")
plt.subplot(221)
plt.plot(E)
plt.title("First run")
plt.ylabel("energy")
plt.subplot(222)
plt.plot(E2)
plt.title("Second run")
plt.subplot(223)
plt.plot(E)
plt.title("First run log scale")
plt.xlabel("timestep")
plt.ylabel("energy")
plt.yscale("log")
plt.subplot(224)
plt.plot(E2)
plt.title("Second run log scale")
plt.yscale("log")
plt.xlabel("timestep")
plt.ylabel("energy")

plt.figure()
plt.plot(v[:2500])
plt.xlabel("timestep/10")
plt.ylabel(r"$v^2$ for each particles")