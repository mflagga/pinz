# pyright: reportUndefinedVariable=false
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({"text.usetex":True, "font.size": 16})

def wczytaj(*nazwy): # funkcja do wczytywania plikow od clauda
    for nazwa in nazwy:
        globals()[nazwa] = np.loadtxt(nazwa + ".csv", delimiter=',')

# wczytanie plików
wczytaj("dyspfile","kxfile","lfile","ufile")

Ny = len(lfile)//2

uRe = ufile[:,0].reshape(2*Ny,2*Ny)
uIm = ufile[:,1].reshape(2*Ny,2*Ny)

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

plt.figure(figsize=(size,size))
plt.scatter(range(len(lfile)),lfile)
plt.savefig("lfile.png")
plt.close()

plt.figure(figsize=(1.25*size,size))
plt.imshow(uRe,origin='lower',cmap='bwr')
plt.colorbar()
plt.tight_layout()
plt.savefig("uRe.png")
plt.close()

plt.figure(figsize=(1.25*size,size))
plt.imshow(uIm,origin='lower',cmap='bwr')
plt.colorbar()
plt.tight_layout()
plt.savefig("uIm.png")
plt.close()
