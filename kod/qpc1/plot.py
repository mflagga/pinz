# pyright: reportUndefinedVariable=false
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({"text.usetex":True, "font.size": 16})

def wczytaj(*nazwy): # funkcja do wczytywania plikow od clauda
    for nazwa in nazwy:
        globals()[nazwa] = np.loadtxt(nazwa + ".csv", delimiter=',')

# wczytanie plików
wczytaj("dyspfile","kxfile")

size = 8

plt.figure(figsize=(1.2*size,size))
for j in range(len(dyspfile[0])):
    plt.plot(kxfile,dyspfile[:,j],c="tab:blue")
plt.grid(ls=":")
plt.xlabel(rf"$k_x\ \left[\frac{{1}}{{nm}}\right]$")
plt.ylabel(rf"$E\ \left[eV\right]$")
plt.tight_layout()
plt.savefig("dyspfile.png",dpi=150)
plt.close()
