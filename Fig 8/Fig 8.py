import numpy as np
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib


cmap = matplotlib.colormaps["cool"]


with open("EDMD continuous.pkl", "rb") as f:
    ec_N, ec_phi, ec_E, ec_std, ec_jump = pickle.load(f)
with open("EDMD discontinuous.pkl", "rb") as f:
    ed_N, ed_phi, ed_E, ed_std, ed_jump = pickle.load(f)

with open("continuous lammps.pkl", "rb") as f:
    lc_N, lc_phi, lc_E, lc_std, lc_jump = pickle.load(f)
with open("discontinuous lammps.pkl", "rb") as f:
    ld_N, ld_phi, ld_E, ld_std, ld_jump = pickle.load(f)

fig, axes = plt.subplots(2, 2, figsize = (15, 10))
fig.set_tight_layout(True)

#EDMD continuous
ax = axes[0][0]
ec_E[-1, np.logical_and(ec_phi > 0.1512, ec_E[-1] == 0)] = None
"""
ec_c = ec_N - ec_N[0]
ec_c = (ec_c/ec_c[-1])**0.4

for i in range(len(ec_N) - 1, 3, -1):
    ec_E[i][ec_phi > 0.151367] = None
    ax.errorbar(ec_phi, ec_E[i], ec_std[i], fmt="o", label = ec_N[i], alpha = 0.7, c = cmap(ec_c[i]))
"""

for i in [len(ec_N) - 1, int(len(ec_N)/2), 3]:
    ec_E[i][ec_phi > 0.15183] = None
    ax.errorbar(ec_phi, ec_E[i], ec_std[i], fmt="o", label = ec_N[i], alpha = 0.7)

#ax.set_xlabel(r"$\phi$")
ax.set_ylabel("Energy")
ax.set_xticks([])
ax.set_yticks([])


axin = inset_axes(ax, width="40%", height="50%", loc=2, borderpad=0.25)
#axin.scatter(ec_N, ec_jump, c = cmap(ec_c))
ec_c = ["black"]*len(ec_N)
ec_c[len(ec_N) - 1] = "C0" 
ec_c[int(len(ec_N)/2)] = "C1"
ec_c[3] = "C2"
axin.scatter(ec_N, ec_jump, c = ec_c)
axin.yaxis.tick_right()
axin.yaxis.set_label_position("right")
axin.set_xlabel(r"$N$",fontsize = 10)
axin.set_ylabel("Gap", fontsize = 10)
axin.tick_params(axis = 'y', which = 'both', labelsize = 10)

for item in ([axin.title, axin.xaxis.label, axin.yaxis.label]):
    item.set_fontsize(10)

for ti in (axin.get_xticklabels()+axin.get_yticklabels()): 
	ti.set_fontsize(10)  
axin.set_yscale("log")
axin.set_xscale("log")


#LAMMPS continuous
ax = axes[0][1]
lc_c = lc_N - lc_N[0]
lc_c = lc_c/lc_c[-1]
for i in range(len(lc_N) - 1, -1, -1):
    ax.errorbar(lc_phi, lc_E[i], lc_std[i], fmt="o", label = lc_N[i], alpha = 0.7, c = cmap(lc_c[i]))

#ax.set_xlabel(r"$\phi$")
#ax.set_ylabel("Energy")



axin = inset_axes(ax, width="40%", height="50%", loc=2, borderpad=0.25)
axin.scatter(lc_N, lc_jump, c =  cmap(lc_c))
axin.yaxis.tick_right()
axin.yaxis.set_label_position("right")
axin.set_xlabel(r"$N$",fontsize = 10)
axin.set_ylabel("Gap", fontsize = 10)
axin.tick_params(axis = 'y', which = 'both', labelsize = 10)

for item in ([axin.title, axin.xaxis.label, axin.yaxis.label]):
    item.set_fontsize(10)

for ti in (axin.get_xticklabels()+axin.get_yticklabels()): 
	ti.set_fontsize(10)  
axin.set_yscale("log")
axin.set_xscale("log")

#EDMD dicontinuous
ax = axes[1][0]
"""
ed_c = ed_N - ed_N[0]
ed_c = (ed_c/ed_c[-1])**0.6
for i in range(0, len(ed_N)):
    ed_E[i][ed_phi > 0.2839] = None
    ax.errorbar(ed_phi, ed_E[i], ed_std[i], fmt="o", label = ed_N[i], alpha = 0.7, c = cmap(ed_c[i]))


ed_c = ed_N - ed_N[0]
ed_c = (ed_c/ed_c[-1])**0.6"""

ed_c = ["C1", "C2", "C0"]
s = 0
for i in [int(len(ed_N)/2), len(ed_N) - 1, 0]:
    ed_E[i][ed_phi > 0.284] = None
    
    ax.errorbar(ed_phi, ed_E[i], ed_std[i], fmt="o", label = ed_N[i], alpha = 0.7, c = ed_c[s])
    s+= 1
    
ax.set_xlabel(r"$\phi$")
ax.set_ylabel("Energy")
ax.set_xticks([])
ax.set_yticks([])


axin = inset_axes(ax, width="40%", height="50%", loc=2, borderpad=0.25)
ed_c = ["black"]*len(ed_N)
ed_c[len(ed_N) - 1] = "C2" 
ed_c[int(len(ed_N)/2)] = "C1"
ed_c[0] = "C0"
axin.scatter(ed_N, ed_jump, c = ed_c)
#axin.scatter(ed_N, ed_jump, c = cmap( ed_c))
axin.set_yscale("log")
axin.set_xscale("log")
axin.yaxis.tick_right()
axin.yaxis.set_label_position("right")
axin.set_xlabel(r"$N$",fontsize = 10)
axin.set_ylabel("Gap", fontsize = 10)
axin.tick_params(axis = 'y', which = 'both', labelsize = 10)

for item in ([axin.title, axin.xaxis.label, axin.yaxis.label]):
    item.set_fontsize(10)

for ti in (axin.get_xticklabels()+axin.get_yticklabels()): 
	ti.set_fontsize(10)  

#LAMMPS discontinuous
ax = axes[1][1]
ld_c = ld_N - ld_N[0]
ld_c = ld_c/ld_c[-1]
for i in range(len(ld_N)):
    
    ld_E[i][ld_phi > 0.218] = None
    ax.errorbar(ld_phi, ld_E[i], ld_std[i], fmt="o", label = ld_N[i], alpha = 0.7, c = cmap(ld_c[i]))

ax.set_xlabel(r"$\phi$")
#ax.set_ylabel("Energy")
ax.set_xticks([])
ax.set_yticks([])


axin = inset_axes(ax, width="40%", height="50%", loc=2, borderpad=0.25)
axin.scatter(ld_N, ld_jump, c =  cmap(ld_c))
axin.yaxis.tick_right()
axin.yaxis.set_label_position("right")
axin.set_xlabel(r"$N$",fontsize = 10)
axin.set_ylabel("Gap", fontsize = 10)
axin.tick_params(axis = 'y', which = 'both', labelsize = 10)

for item in ([axin.title, axin.xaxis.label, axin.yaxis.label]):
    item.set_fontsize(10)

for ti in (axin.get_xticklabels()+axin.get_yticklabels()): 
	ti.set_fontsize(10)  
axin.set_yscale("log")
axin.set_xscale("log")

fig.savefig('first vs second.pdf', dpi=fig.dpi)